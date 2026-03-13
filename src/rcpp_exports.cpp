// rcpp_exports.cpp  –  Rcpp bindings for copc4R
//
// Exposes two C++ functions to R:
//   cpp_read_copc_header(path_or_url)
//   cpp_read_copc(path_or_url, bbox, zrange, select, max_points)
//
// The R wrappers live in R/read_copc.R.

#include <Rcpp.h>
#include "range_reader.h"
#include "las14_header.h"
#include "copc_hierarchy.h"
#include "laszip_decompress.h"
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <future>
#include <thread>

using namespace Rcpp;

// ═══════════════════════════════════════════════════════════════════════
// Helper: convert LASHeader to an rlas-style named R list
// ═══════════════════════════════════════════════════════════════════════
static Rcpp::List header_to_rlas_list(const copc4r::LASHeader& h) {
    const auto& p = h.pub;

    Rcpp::List out;

    out["File Signature"]           = std::string(p.file_signature, 4);
    out["File Source ID"]           = static_cast<int>(p.file_source_id);

    // Global Encoding: decompose bitfield into named boolean list (rlas convention)
    {
        uint16_t ge = p.global_encoding;
        Rcpp::List ge_list;
        ge_list["GPS Time Type"]                  = (bool)(ge & 0x01);
        ge_list["Waveform Data Packets Internal"] = (bool)(ge & 0x02);
        ge_list["Waveform Data Packets External"] = (bool)(ge & 0x04);
        ge_list["Synthetic Return Numbers"]       = (bool)(ge & 0x08);
        ge_list["WKT"]                            = (bool)(ge & 0x10);
        ge_list["Aggregate Model"]                = (bool)(ge & 0x20);
        out["Global Encoding"] = ge_list;
    }

    // GUID: format as "xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx" string (rlas convention)
    {
        char guid_str[37];
        snprintf(guid_str, sizeof(guid_str),
                 "%08x-%04x-%04x-%02x%02x-%02x%02x%02x%02x%02x%02x",
                 p.guid_data1, p.guid_data2, p.guid_data3,
                 p.guid_data4[0], p.guid_data4[1],
                 p.guid_data4[2], p.guid_data4[3],
                 p.guid_data4[4], p.guid_data4[5],
                 p.guid_data4[6], p.guid_data4[7]);
        out["Project ID - GUID"] = std::string(guid_str);
    }

    out["Version Major"]            = static_cast<int>(p.version_major);
    out["Version Minor"]            = static_cast<int>(p.version_minor);
    out["System Identifier"]        = std::string(p.system_identifier,
        strnlen(p.system_identifier, 32));
    out["Generating Software"]      = std::string(p.generating_software,
        strnlen(p.generating_software, 32));
    out["File Creation Day of Year"]= static_cast<int>(p.creation_day);
    out["File Creation Year"]       = static_cast<int>(p.creation_year);
    out["Header Size"]              = static_cast<int>(p.header_size);
    out["Offset to point data"]     = static_cast<double>(p.offset_to_point_data);
    out["Number of variable length records"] = static_cast<int>(p.num_vlrs);
    // Strip LASzip compression bit (bit 7) so PDRF reads 6/7/8 not 134/135/136
    out["Point Data Format ID"]     = static_cast<int>(p.point_data_format & 0x7F);
    out["Point Data Record Length"]  = static_cast<int>(p.point_data_record_length);

    // prefer LAS 1.4 64-bit counts
    double npts = (p.point_count_14 > 0)
        ? static_cast<double>(p.point_count_14)
        : static_cast<double>(p.legacy_point_count);
    out["Number of point records"] = npts;

    // Number of points by return (combine legacy + 1.4)
    Rcpp::NumericVector pbr(15);
    for (int i = 0; i < 15; i++)
        pbr[i] = static_cast<double>(p.points_by_return_14[i]);
    out["Number of points by return"] = pbr;

    out["X scale factor"] = p.x_scale;
    out["Y scale factor"] = p.y_scale;
    out["Z scale factor"] = p.z_scale;
    out["X offset"] = p.x_offset;
    out["Y offset"] = p.y_offset;
    out["Z offset"] = p.z_offset;
    out["Max X"] = p.max_x;
    out["Min X"] = p.min_x;
    out["Max Y"] = p.max_y;
    out["Min Y"] = p.min_y;
    out["Max Z"] = p.max_z;
    out["Min Z"] = p.min_z;

    // ── VLRs in rlas format ───────────────────────────────────────────
    // rlas VLR format: reserved, user ID, record ID, length after header,
    // description, and potentially parsed content (e.g. WKT string).
    auto format_vlr = [](const copc4r::ParsedVLR& vlr) -> Rcpp::List {
        Rcpp::List v;
        v["reserved"]            = static_cast<int>(vlr.reserved);
        v["user ID"]             = vlr.user_id;
        v["record ID"]           = static_cast<int>(vlr.record_id);
        v["length after header"] = static_cast<int>(vlr.data.size());
        v["description"]         = vlr.description;

        // Populate well-known VLR data content
        if (vlr.user_id == "LASF_Projection" && vlr.record_id == 2112 &&
            !vlr.data.empty()) {
            // WKT CRS string
            std::string wkt(vlr.data.begin(), vlr.data.end());
            // Strip trailing null if present
            if (!wkt.empty() && wkt.back() == '\0') wkt.pop_back();
            v["WKT OGC COORDINATE SYSTEM"] = wkt;
        } else if (vlr.user_id == "LASF_Projection" && vlr.record_id == 34735 &&
                   vlr.data.size() >= 8) {
            // GeoTIFF GeoKeyDirectoryTag
            Rcpp::IntegerVector tags(vlr.data.size() / 2);
            for (size_t j = 0; j < vlr.data.size() / 2; ++j) {
                uint16_t val;
                std::memcpy(&val, vlr.data.data() + j * 2, 2);
                tags[j] = static_cast<int>(val);
            }
            v["GeoKeyDirectoryTag"] = tags;
        }
        // Other VLR data stored as raw bytes
        return v;
    };

    // Assign names to VLR list elements matching rlas conventions
    auto vlr_name = [](const copc4r::ParsedVLR& vlr) -> std::string {
        if (vlr.user_id == "LASF_Projection" && vlr.record_id == 2112)
            return "WKT OGC CS";
        if (vlr.user_id == "LASF_Projection" && vlr.record_id == 34735)
            return "GeoKeyDirectoryTag";
        if (vlr.user_id == "LASF_Projection" && vlr.record_id == 34736)
            return "GeoDoubleParamsTag";
        if (vlr.user_id == "LASF_Projection" && vlr.record_id == 34737)
            return "GeoAsciiParamsTag";
        if (vlr.user_id == "copc" && vlr.record_id == 1)
            return "copc";
        if (vlr.user_id == "laszip encoded" && vlr.record_id == 22204)
            return "LASzip";
        return vlr.user_id + ":" + std::to_string(vlr.record_id);
    };

    Rcpp::List vlr_list;
    Rcpp::CharacterVector vlr_names;
    for (size_t i = 0; i < h.vlrs.size(); ++i) {
        vlr_list.push_back(format_vlr(h.vlrs[i]));
        vlr_names.push_back(vlr_name(h.vlrs[i]));
    }
    vlr_list.attr("names") = vlr_names;
    out["Variable Length Records"] = vlr_list;

    // EVLR info
    Rcpp::List evlr_list;
    Rcpp::CharacterVector evlr_names;
    for (size_t i = 0; i < h.evlrs.size(); ++i) {
        evlr_list.push_back(format_vlr(h.evlrs[i]));
        evlr_names.push_back(vlr_name(h.evlrs[i]));
    }
    evlr_list.attr("names") = evlr_names;
    out["Extended Variable Length Records"] = evlr_list;

    // COPC-specific extras
    try {
        copc4r::COPCInfo info = copc4r::parse_copc_info(h);
        Rcpp::List copc;
        copc["center_x"]    = info.center_x;
        copc["center_y"]    = info.center_y;
        copc["center_z"]    = info.center_z;
        copc["halfsize"]    = info.halfsize;
        copc["spacing"]     = info.spacing;
        copc["root_hier_offset"] = static_cast<double>(info.root_hier_offset);
        copc["root_hier_size"]   = static_cast<double>(info.root_hier_size);
        copc["gpstime_minimum"]  = info.gpstime_minimum;
        copc["gpstime_maximum"]  = info.gpstime_maximum;
        out["COPC Info"] = copc;
    } catch (...) {
        // Not all LAS files have COPC info; silently skip
    }

    return out;
}

// ═══════════════════════════════════════════════════════════════════════
// Helper: determine which columns to output based on "select" string
// Default (empty or nullptr) = all.  Lowercase letters = field groups:
//   x=XYZ, i=Intensity, t=gpstime, r=ReturnNumber, n=NumberOfReturns,
//   c=Classification, a=ScanAngle, s=PointSourceID, u=UserData,
//   p=flags (ScanDir, EdgeOfFlight, ClassFlags, ScannerChannel),
//   R=RGB (Red,Green,Blue), N=NIR, *=everything
// ═══════════════════════════════════════════════════════════════════════
struct ColumnSel {
    bool XYZ = false;
    bool intensity = false;
    bool gpstime = false;
    bool return_number = false;
    bool number_of_returns = false;
    bool classification = false;
    bool scan_angle = false;
    bool point_source_id = false;
    bool user_data = false;
    bool flags = false;    // ScanDir, Edge, ClassFlags, ScannerChannel
    bool rgb = false;
    bool nir = false;
};

static ColumnSel parse_select(const std::string& sel, uint8_t pdrf) {
    ColumnSel c;
    bool all = sel.empty() || sel == "*";
    if (all) {
        c.XYZ = c.intensity = c.gpstime = c.return_number = true;
        c.number_of_returns = c.classification = c.scan_angle = true;
        c.point_source_id = c.user_data = c.flags = true;
        c.rgb = (pdrf >= 7);
        c.nir = (pdrf >= 8);
        return c;
    }
    for (char ch : sel) {
        switch (ch) {
        case 'x': case 'y': case 'z': c.XYZ = true; break;
        case 'i': c.intensity = true; break;
        case 't': c.gpstime = true; break;
        case 'r': c.return_number = true; break;
        case 'n': c.number_of_returns = true; break;
        case 'c': c.classification = true; break;
        case 'a': c.scan_angle = true; break;
        case 's': c.point_source_id = true; break;
        case 'u': c.user_data = true; break;
        case 'p': c.flags = true; break;
        case 'R': c.rgb = (pdrf >= 7); break;
        case 'N': c.nir = (pdrf >= 8); break;
        case '*': return parse_select("", pdrf);
        default: break; // ignore unknown
        }
    }
    // Always include XYZ at minimum
    c.XYZ = true;
    return c;
}

// [[Rcpp::export]]
bool cpp_has_http_support() {
#ifdef COPC4R_WITH_CURL
    return true;
#else
    return false;
#endif
}

// [[Rcpp::export]]
Rcpp::List cpp_read_copc_header(std::string path_or_url) {
    auto reader = copc4r::make_reader(path_or_url);
    copc4r::LASHeader h = copc4r::read_las_header(*reader);
    return header_to_rlas_list(h);
}

// ═══════════════════════════════════════════════════════════════════════
// Count nodes overlapping a bbox (for density estimation)
// ═══════════════════════════════════════════════════════════════════════
// [[Rcpp::export]]
Rcpp::List cpp_count_nodes(std::string path_or_url,
                           Rcpp::NumericVector bbox) {
    auto reader = copc4r::make_reader(path_or_url);
    copc4r::LASHeader header = copc4r::read_las_header(*reader);
    copc4r::COPCInfo info = copc4r::parse_copc_info(header);

    bool has_bbox = (bbox.size() >= 4);
    double bxmin = 0, bymin = 0, bxmax = 0, bymax = 0;
    if (has_bbox) {
        bxmin = bbox[0]; bymin = bbox[1];
        bxmax = bbox[2]; bymax = bbox[3];
    }

    auto nc = copc4r::count_nodes(*reader, info,
                                   bxmin, bymin, bxmax, bymax, has_bbox);

    return Rcpp::List::create(
        Rcpp::Named("total_points") = static_cast<double>(nc.total_points),
        Rcpp::Named("num_nodes")    = static_cast<double>(nc.num_nodes)
    );
}

// ═══════════════════════════════════════════════════════════════════════
// Select hierarchy nodes (expose to R for cache orchestration)
// ═══════════════════════════════════════════════════════════════════════
// [[Rcpp::export]]
Rcpp::List cpp_select_nodes(std::string path_or_url,
                            Rcpp::NumericVector bbox,
                            Rcpp::NumericVector zrange,
                            int max_depth = -1) {
    auto reader = copc4r::make_reader(path_or_url);
    copc4r::LASHeader header = copc4r::read_las_header(*reader);
    copc4r::COPCInfo info = copc4r::parse_copc_info(header);

    bool has_bbox = (bbox.size() >= 4);
    double bxmin = 0, bymin = 0, bxmax = 0, bymax = 0;
    if (has_bbox) {
        bxmin = bbox[0]; bymin = bbox[1];
        bxmax = bbox[2]; bymax = bbox[3];
    }
    bool has_zrange = (zrange.size() >= 2);
    double zmin = 0, zmax = 0;
    if (has_zrange) {
        zmin = zrange[0]; zmax = zrange[1];
    }

    auto nodes = copc4r::select_nodes(
        *reader, info,
        bxmin, bymin, bxmax, bymax,
        zmin, zmax, has_bbox, has_zrange, max_depth);

    size_t n = nodes.size();
    Rcpp::NumericVector offsets(n), sizes(n), point_counts(n);
    Rcpp::IntegerVector levels(n);
    for (size_t i = 0; i < n; i++) {
        offsets[i]      = static_cast<double>(nodes[i].offset);
        sizes[i]        = static_cast<double>(nodes[i].byte_size);
        point_counts[i] = static_cast<double>(nodes[i].point_count);
        levels[i]       = nodes[i].key.level;
    }
    return Rcpp::List::create(
        Rcpp::Named("offset")      = offsets,
        Rcpp::Named("byte_size")   = sizes,
        Rcpp::Named("point_count") = point_counts,
        Rcpp::Named("level")       = levels
    );
}

// ═══════════════════════════════════════════════════════════════════════
// Fetch a single raw chunk (for R-level cache orchestration)
// ═══════════════════════════════════════════════════════════════════════
// [[Rcpp::export]]
Rcpp::RawVector cpp_fetch_raw_chunk(std::string path_or_url,
                                    double offset,
                                    double byte_size) {
    auto reader = copc4r::make_reader(path_or_url);
    auto data = reader->read(static_cast<uint64_t>(offset),
                             static_cast<uint64_t>(byte_size));
    Rcpp::RawVector out(data.size());
    if (!data.empty())
        std::memcpy(out.begin(), data.data(), data.size());
    return out;
}

// ═══════════════════════════════════════════════════════════════════════
// Parallel fetch of multiple raw chunks (for R-level cache orchestration)
// Creates n_threads persistent HTTP connections and batches the work.
// ═══════════════════════════════════════════════════════════════════════
// [[Rcpp::export]]
Rcpp::List cpp_fetch_raw_chunks_parallel(std::string path_or_url,
                                         Rcpp::NumericVector offsets,
                                         Rcpp::NumericVector sizes,
                                         int n_threads = 4) {
    size_t n = offsets.size();
    if (n == 0) return Rcpp::List(0);

    int actual_threads = std::min(n_threads, static_cast<int>(n));
    if (actual_threads < 1) actual_threads = 1;

    // Create persistent worker readers (one per thread for connection reuse)
    std::vector<std::unique_ptr<copc4r::RangeReader>> workers(actual_threads);
    for (int t = 0; t < actual_threads; t++) {
        workers[t] = copc4r::make_reader(path_or_url);
    }

    Rcpp::List result(n);

    for (size_t batch_start = 0; batch_start < n;
         batch_start += static_cast<size_t>(actual_threads)) {
        size_t batch_end = std::min(batch_start + static_cast<size_t>(actual_threads), n);
        std::vector<std::future<std::vector<uint8_t>>> futures;

        for (size_t i = batch_start; i < batch_end; i++) {
            size_t tid = i - batch_start;
            auto* rd = workers[tid].get();
            uint64_t off = static_cast<uint64_t>(offsets[i]);
            uint64_t sz  = static_cast<uint64_t>(sizes[i]);
            futures.push_back(std::async(std::launch::async,
                [rd, off, sz]() -> std::vector<uint8_t> {
                    return rd->read(off, sz);
                }
            ));
        }

        for (size_t i = 0; i < futures.size(); i++) {
            auto data = futures[i].get();
            Rcpp::RawVector rv(data.size());
            if (!data.empty())
                std::memcpy(rv.begin(), data.data(), data.size());
            result[batch_start + i] = rv;
        }

        Rcpp::checkUserInterrupt();
    }

    return result;
}

// ═══════════════════════════════════════════════════════════════════════
// Attribute filter: parse a lidR-style filter string and apply
// post-decompression.  Supported predicates:
//   -keep_class <c1> <c2> ...       keep only these classifications
//   -drop_class <c1> <c2> ...       drop these classifications
//   -keep_first                     ReturnNumber == 1
//   -keep_last                      ReturnNumber == NumberOfReturns
//   -keep_single                    NumberOfReturns == 1
//   -drop_withheld                  drop Withheld flag
//   -drop_overlap                   drop Overlap flag
//   -keep_intensity_above <val>     Intensity >= val
//   -keep_intensity_below <val>     Intensity <= val
// ═══════════════════════════════════════════════════════════════════════
struct PointFilter {
    bool active = false;
    std::vector<int> keep_class;
    std::vector<int> drop_class;
    bool keep_first = false;
    bool keep_last = false;
    bool keep_single = false;
    bool drop_withheld = false;
    bool drop_overlap = false;
    int  intensity_above = -1; // -1 = not set
    int  intensity_below = -1;
    // ── Additional predicates (v0.3) ─────────────────────────────
    std::vector<int> keep_return;    // specific return numbers
    std::vector<int> drop_return;
    double z_above  = -std::numeric_limits<double>::infinity();  // keep Z >= val
    double z_below  =  std::numeric_limits<double>::infinity();  // keep Z <= val
    double scan_angle_above = -std::numeric_limits<double>::infinity();
    double scan_angle_below =  std::numeric_limits<double>::infinity();
    bool   drop_noise  = false;    // drop classification 7 (low) & 18 (high)
    bool   keep_ground = false;    // keep only classification 2
    double keep_random_fraction = 1.0;    // thin: keep fraction (0,1]
    int    keep_every_nth = 0;     // keep every Nth point (0 = disabled)
    // Internal counter for keep_every_nth
    mutable uint64_t _nth_counter = 0;
    bool has_z_above() const { return z_above != -std::numeric_limits<double>::infinity(); }
    bool has_z_below() const { return z_below !=  std::numeric_limits<double>::infinity(); }
    bool has_sa_above() const { return scan_angle_above != -std::numeric_limits<double>::infinity(); }
    bool has_sa_below() const { return scan_angle_below !=  std::numeric_limits<double>::infinity(); }
};

// Helper: read consecutive integers from stream until next flag or EOF
static void read_ints(std::istringstream& iss, std::vector<int>& out) {
    std::string token;
    while (iss.peek() != '-' && !iss.eof() && iss >> token) {
        try { out.push_back(std::stoi(token)); }
        catch (...) { break; }
        // Peek past whitespace
        while (iss.peek() == ' ' || iss.peek() == '\t') iss.get();
        if (iss.peek() == '-') break;
    }
}

static PointFilter parse_filter(const std::string& fstr) {
    PointFilter f;
    if (fstr.empty()) return f;
    f.active = true;

    std::istringstream iss(fstr);
    std::string token;
    while (iss >> token) {
        if (token == "-keep_class") {
            read_ints(iss, f.keep_class);
        } else if (token == "-drop_class") {
            read_ints(iss, f.drop_class);
        } else if (token == "-keep_first") {
            f.keep_first = true;
        } else if (token == "-keep_last") {
            f.keep_last = true;
        } else if (token == "-keep_single") {
            f.keep_single = true;
        } else if (token == "-drop_withheld") {
            f.drop_withheld = true;
        } else if (token == "-drop_overlap") {
            f.drop_overlap = true;
        } else if (token == "-keep_intensity_above") {
            if (iss >> token) {
                try { f.intensity_above = std::stoi(token); }
                catch (...) {}
            }
        } else if (token == "-keep_intensity_below") {
            if (iss >> token) {
                try { f.intensity_below = std::stoi(token); }
                catch (...) {}
            }
        } else if (token == "-keep_return") {
            read_ints(iss, f.keep_return);
        } else if (token == "-drop_return") {
            read_ints(iss, f.drop_return);
        } else if (token == "-keep_z_above") {
            if (iss >> token) {
                try { f.z_above = std::stod(token); }
                catch (...) {}
            }
        } else if (token == "-keep_z_below") {
            if (iss >> token) {
                try { f.z_below = std::stod(token); }
                catch (...) {}
            }
        } else if (token == "-keep_scan_angle_above") {
            if (iss >> token) {
                try { f.scan_angle_above = std::stod(token); }
                catch (...) {}
            }
        } else if (token == "-keep_scan_angle_below") {
            if (iss >> token) {
                try { f.scan_angle_below = std::stod(token); }
                catch (...) {}
            }
        } else if (token == "-drop_noise") {
            f.drop_noise = true;
        } else if (token == "-keep_ground") {
            f.keep_ground = true;
        } else if (token == "-keep_random_fraction") {
            if (iss >> token) {
                try { f.keep_random_fraction = std::stod(token); }
                catch (...) {}
            }
        } else if (token == "-keep_every_nth") {
            if (iss >> token) {
                try { f.keep_every_nth = std::stoi(token); }
                catch (...) {}
            }
        }
    }
    return f;
}

static bool passes_filter(const copc4r::DecodedPoint& pt,
                           const PointFilter& f) {
    if (!f.active) return true;

    // ── Classification filters ────────────────────────────────────
    if (!f.keep_class.empty()) {
        bool found = false;
        for (int c : f.keep_class) {
            if (pt.classification == c) { found = true; break; }
        }
        if (!found) return false;
    }
    if (!f.drop_class.empty()) {
        for (int c : f.drop_class) {
            if (pt.classification == c) return false;
        }
    }
    if (f.keep_ground && pt.classification != 2) return false;
    if (f.drop_noise && (pt.classification == 7 || pt.classification == 18))
        return false;

    // ── Return number filters ─────────────────────────────────────
    if (f.keep_first && pt.return_number != 1) return false;
    if (f.keep_last && pt.return_number != pt.number_of_returns) return false;
    if (f.keep_single && pt.number_of_returns != 1) return false;
    if (!f.keep_return.empty()) {
        bool found = false;
        for (int r : f.keep_return) {
            if (pt.return_number == r) { found = true; break; }
        }
        if (!found) return false;
    }
    if (!f.drop_return.empty()) {
        for (int r : f.drop_return) {
            if (pt.return_number == r) return false;
        }
    }

    // ── Flag filters ──────────────────────────────────────────────
    if (f.drop_withheld && ((pt.classification_flags >> 2) & 1)) return false;
    if (f.drop_overlap && ((pt.classification_flags >> 3) & 1)) return false;

    // ── Intensity filters ─────────────────────────────────────────
    if (f.intensity_above >= 0 && pt.intensity < static_cast<uint16_t>(f.intensity_above))
        return false;
    if (f.intensity_below >= 0 && pt.intensity > static_cast<uint16_t>(f.intensity_below))
        return false;

    // ── Z-range filters ───────────────────────────────────────────
    if (f.has_z_above() && pt.Z < f.z_above) return false;
    if (f.has_z_below() && pt.Z > f.z_below) return false;

    // ── Scan-angle filters (degrees: raw * 0.006) ─────────────────
    if (f.has_sa_above() || f.has_sa_below()) {
        double sa_deg = pt.scan_angle * 0.006;
        if (f.has_sa_above() && sa_deg < f.scan_angle_above) return false;
        if (f.has_sa_below() && sa_deg > f.scan_angle_below) return false;
    }

    // ── Decimation: keep every Nth point ──────────────────────────
    if (f.keep_every_nth > 0) {
        f._nth_counter++;
        if ((f._nth_counter % f.keep_every_nth) != 1) return false;
    }

    // ── Random fraction thinning ──────────────────────────────────
    if (f.keep_random_fraction < 1.0) {
        double u = R::runif(0.0, 1.0);
        if (u > f.keep_random_fraction) return false;
    }

    return true;
}

// ═══════════════════════════════════════════════════════════════════════
// Helper: decode extra bytes for a single point into doubles
// ═══════════════════════════════════════════════════════════════════════
static double decode_extra_byte_value(const uint8_t* data,
                                       const copc4r::ExtraBytesRecord& eb) {
    double raw = 0.0;
    switch (eb.data_type) {
        case 1:  raw = static_cast<double>(*data); break;
        case 2:  raw = static_cast<double>(*reinterpret_cast<const int8_t*>(data)); break;
        case 3:  { uint16_t v; std::memcpy(&v, data, 2); raw = static_cast<double>(v); } break;
        case 4:  { int16_t v;  std::memcpy(&v, data, 2); raw = static_cast<double>(v); } break;
        case 5:  { uint32_t v; std::memcpy(&v, data, 4); raw = static_cast<double>(v); } break;
        case 6:  { int32_t v;  std::memcpy(&v, data, 4); raw = static_cast<double>(v); } break;
        case 7:  { uint64_t v; std::memcpy(&v, data, 8); raw = static_cast<double>(v); } break;
        case 8:  { int64_t v;  std::memcpy(&v, data, 8); raw = static_cast<double>(v); } break;
        case 9:  { float v;    std::memcpy(&v, data, 4); raw = static_cast<double>(v); } break;
        case 10: { double v;   std::memcpy(&v, data, 8); raw = v; } break;
        default: return 0.0;
    }
    // Apply scale and offset
    if (eb.has_scale) raw *= eb.scale;
    if (eb.has_offset) raw += eb.offset;
    return raw;
}

// [[Rcpp::export]]
Rcpp::List cpp_read_copc(std::string path_or_url,
                         Rcpp::NumericVector bbox,
                         Rcpp::NumericVector zrange,
                         std::string select,
                         double max_points_dbl,
                         bool progress,
                         int max_depth = -1,
                         std::string filter = "",
                         int n_threads = 4,
                         Rcpp::Nullable<Rcpp::List> prefetched = R_NilValue) {
    // ── 1. Open reader and parse header + COPC info ───────────────────
    auto reader = copc4r::make_reader(path_or_url);
    copc4r::LASHeader header = copc4r::read_las_header(*reader);
    copc4r::COPCInfo info = copc4r::parse_copc_info(header);

    // Strip LASzip compression bit (bit 7)
    uint8_t  pdrf = header.pub.point_data_format & 0x7F;
    uint16_t prl  = header.pub.point_data_record_length;

    // ── 1b. Extract LASzip VLR data (needed for decompression) ────────
    const copc4r::ParsedVLR* lz_vlr =
        copc4r::find_vlr(header, "laszip encoded", 22204);
    if (!lz_vlr)
        Rcpp::stop("No LASzip VLR found – file may not be LAZ-compressed");
    const std::vector<uint8_t>& laszip_vlr_data = lz_vlr->data;

    // ── 1c. Parse Extra Bytes VLR ─────────────────────────────────────
    auto extra_bytes_defs = copc4r::parse_extra_bytes_vlr(header);
    bool has_extra_bytes = !extra_bytes_defs.empty();

    // ── 1d. Parse point filter ────────────────────────────────────────
    PointFilter pf = parse_filter(filter);

    // ── 2. Parse query parameters ─────────────────────────────────────
    bool has_bbox = (bbox.size() >= 4);
    double bxmin = 0, bymin = 0, bxmax = 0, bymax = 0;
    if (has_bbox) {
        bxmin = bbox[0]; bymin = bbox[1];
        bxmax = bbox[2]; bymax = bbox[3];
    }
    bool has_zrange = (zrange.size() >= 2);
    double zmin = 0, zmax = 0;
    if (has_zrange) {
        zmin = zrange[0]; zmax = zrange[1];
    }

    uint64_t max_pts = (std::isinf(max_points_dbl) || max_points_dbl <= 0)
        ? UINT64_MAX
        : static_cast<uint64_t>(max_points_dbl);

    // ── 3. Select nodes from hierarchy ────────────────────────────────
    auto nodes = copc4r::select_nodes(
        *reader, info,
        bxmin, bymin, bxmax, bymax,
        zmin, zmax, has_bbox, has_zrange,
        max_depth);

    if (progress) {
        uint64_t total_chunk_bytes = 0;
        uint64_t total_node_points = 0;
        for (auto& nc : nodes) {
            total_chunk_bytes += nc.byte_size;
            total_node_points += nc.point_count;
        }
        Rcpp::Rcout << "Selected " << nodes.size() << " COPC node(s)"
                    << " (" << (total_chunk_bytes / 1024) << " KB compressed"
                    << ", ~" << total_node_points << " points before point-level filter)\n";
    }

    // ── 4. Determine columns ──────────────────────────────────────────
    ColumnSel sel = parse_select(select, pdrf);

    // ── 5. Decompress all selected chunks and collect points ──────────
    // Pre-allocate vectors
    uint64_t total_pts = 0;
    for (auto& nc : nodes) total_pts += nc.point_count;
    if (total_pts > max_pts) total_pts = max_pts;

    // We'll build column vectors directly for efficiency
    Rcpp::NumericVector  vX, vY, vZ, vGpstime, vScanAngle;
    Rcpp::IntegerVector  vIntensity, vReturnNumber, vNumberOfReturns;
    Rcpp::IntegerVector  vClassification, vUserData, vPointSourceID;
    Rcpp::IntegerVector  vScanDirFlag, vEdgeOfFlight, vClassFlags, vScannerChannel;
    Rcpp::IntegerVector  vR, vG, vB, vNIR;

    // Extra Bytes columns (one NumericVector per dimension)
    std::vector<Rcpp::NumericVector> vExtraCols(extra_bytes_defs.size());

    auto reserve = [&](uint64_t n) {
        if (sel.XYZ)               { vX = Rcpp::NumericVector(n); vY = Rcpp::NumericVector(n); vZ = Rcpp::NumericVector(n); }
        if (sel.intensity)         vIntensity = Rcpp::IntegerVector(n);
        if (sel.gpstime)           vGpstime = Rcpp::NumericVector(n);
        if (sel.return_number)     vReturnNumber = Rcpp::IntegerVector(n);
        if (sel.number_of_returns) vNumberOfReturns = Rcpp::IntegerVector(n);
        if (sel.classification)    vClassification = Rcpp::IntegerVector(n);
        if (sel.scan_angle)        vScanAngle = Rcpp::NumericVector(n);
        if (sel.point_source_id)   vPointSourceID = Rcpp::IntegerVector(n);
        if (sel.user_data)         vUserData = Rcpp::IntegerVector(n);
        if (sel.flags)  {
            vScanDirFlag = Rcpp::IntegerVector(n);
            vEdgeOfFlight = Rcpp::IntegerVector(n);
            vClassFlags = Rcpp::IntegerVector(n);
            vScannerChannel = Rcpp::IntegerVector(n);
        }
        if (sel.rgb)    { vR = Rcpp::IntegerVector(n); vG = Rcpp::IntegerVector(n); vB = Rcpp::IntegerVector(n); }
        if (sel.nir)    vNIR = Rcpp::IntegerVector(n);
        for (size_t e = 0; e < extra_bytes_defs.size(); e++) {
            vExtraCols[e] = Rcpp::NumericVector(n);
        }
    };
    reserve(total_pts);

    // ── 5a. Fetch compressed chunks ────────────────────────────────────
    //     Supports three modes:
    //     (1) prefetched: R-level cache provided raw bytes
    //     (2) parallel:   multi-threaded HTTP fetch via std::async
    //     (3) sequential: default single-reader path
    std::vector<std::vector<uint8_t>> compressed_chunks(nodes.size());

    bool use_prefetched = prefetched.isNotNull();
    if (use_prefetched) {
        Rcpp::List pf_list(prefetched.get());
        if (static_cast<size_t>(pf_list.size()) != nodes.size()) {
            Rcpp::warning("prefetched chunk count mismatch; falling back to fetch");
            use_prefetched = false;
        } else {
            for (size_t i = 0; i < nodes.size(); i++) {
                Rcpp::RawVector rv = pf_list[i];
                compressed_chunks[i].assign(rv.begin(), rv.end());
            }
        }
    }

    if (!use_prefetched) {
        bool is_http = (path_or_url.find("http://") == 0 ||
                        path_or_url.find("https://") == 0);

        if (n_threads > 1 && is_http && nodes.size() > 1) {
            // ── Parallel HTTP fetch ───────────────────────────────────
            int actual_threads = std::min(n_threads,
                                          static_cast<int>(nodes.size()));

            // Persistent readers (one per thread for connection reuse)
            std::vector<std::unique_ptr<copc4r::RangeReader>> workers(
                actual_threads);
            for (int t = 0; t < actual_threads; t++) {
                workers[t] = copc4r::make_reader(path_or_url);
            }

            if (progress) {
                Rcpp::Rcout << "Parallel fetch: " << nodes.size()
                            << " chunks, " << actual_threads
                            << " threads\n";
            }

            for (size_t batch_start = 0; batch_start < nodes.size();
                 batch_start += static_cast<size_t>(actual_threads)) {
                size_t batch_end = std::min(
                    batch_start + static_cast<size_t>(actual_threads),
                    nodes.size());
                std::vector<std::future<std::vector<uint8_t>>> futures;

                for (size_t i = batch_start; i < batch_end; i++) {
                    size_t tid = i - batch_start;
                    auto* rd = workers[tid].get();
                    auto off = nodes[i].offset;
                    auto sz  = nodes[i].byte_size;
                    futures.push_back(std::async(std::launch::async,
                        [rd, off, sz]() -> std::vector<uint8_t> {
                            return rd->read(off, sz);
                        }
                    ));
                }

                for (size_t i = 0; i < futures.size(); i++) {
                    compressed_chunks[batch_start + i] = futures[i].get();
                }

                Rcpp::checkUserInterrupt();
            }
        } else {
            // ── Sequential fetch ──────────────────────────────────────
            for (size_t i = 0; i < nodes.size(); i++) {
                compressed_chunks[i] = reader->read(
                    nodes[i].offset, nodes[i].byte_size);
            }
        }
    }

    // ── 5b. Decompress and filter ─────────────────────────────────────
    uint64_t idx = 0;
    for (size_t ni = 0; ni < nodes.size(); ++ni) {
        auto& nc = nodes[ni];

        if (progress && ni % 50 == 0 && ni > 0) {
            Rcpp::Rcout << "  Decompressing chunk " << ni << "/" << nodes.size() << "\n";
        }

        // Decompress (from pre-fetched buffer)
        std::vector<copc4r::DecodedPoint> pts;
        try {
            pts = copc4r::decompress_chunk(
                compressed_chunks[ni], laszip_vlr_data, pdrf, prl, nc.point_count,
                header.pub.x_scale, header.pub.y_scale, header.pub.z_scale,
                header.pub.x_offset, header.pub.y_offset, header.pub.z_offset);
        } catch (std::exception& e) {
            Rcpp::warning(std::string("Skipping chunk: ") + e.what());
            continue;
        }

        for (auto& pt : pts) {
            if (idx >= total_pts) break;

            // Optional: fine-grained bbox+zrange filtering at point level
            if (has_bbox) {
                if (pt.X < bxmin || pt.X > bxmax ||
                    pt.Y < bymin || pt.Y > bymax)
                    continue;
            }
            if (has_zrange) {
                if (pt.Z < zmin || pt.Z > zmax)
                    continue;
            }

            // Attribute filter
            if (!passes_filter(pt, pf))
                continue;

            if (sel.XYZ) { vX[idx] = pt.X; vY[idx] = pt.Y; vZ[idx] = pt.Z; }
            if (sel.intensity)         vIntensity[idx] = pt.intensity;
            if (sel.gpstime)           vGpstime[idx]   = pt.gpstime;
            if (sel.return_number)     vReturnNumber[idx] = pt.return_number;
            if (sel.number_of_returns) vNumberOfReturns[idx] = pt.number_of_returns;
            if (sel.classification)    vClassification[idx] = pt.classification;
            if (sel.scan_angle)        vScanAngle[idx] = pt.scan_angle * 0.006;
            if (sel.point_source_id)   vPointSourceID[idx] = pt.point_source_id;
            if (sel.user_data)         vUserData[idx] = pt.user_data;
            if (sel.flags) {
                vScanDirFlag[idx]    = pt.scan_direction_flag;
                vEdgeOfFlight[idx]   = pt.edge_of_flight_line;
                vClassFlags[idx]     = pt.classification_flags;
                vScannerChannel[idx] = pt.scanner_channel;
            }
            if (sel.rgb) { vR[idx] = pt.red; vG[idx] = pt.green; vB[idx] = pt.blue; }
            if (sel.nir)  vNIR[idx] = pt.nir;

            // Extra Bytes columns
            if (has_extra_bytes && !pt.extra_bytes.empty()) {
                for (size_t e = 0; e < extra_bytes_defs.size(); e++) {
                    auto& eb = extra_bytes_defs[e];
                    if (eb.byte_offset + eb.byte_size <= pt.extra_bytes.size()) {
                        vExtraCols[e][idx] = decode_extra_byte_value(
                            pt.extra_bytes.data() + eb.byte_offset, eb);
                    }
                }
            }

            idx++;
        }

        if (idx >= max_pts) break;

        // Allow R interrupt checking
        Rcpp::checkUserInterrupt();
    }

    // Build named list for data.table columns with the right length.
    Rcpp::List dt_cols;
    if (sel.XYZ) {
        dt_cols["X"] = Rcpp::NumericVector(vX.begin(), vX.begin() + idx);
        dt_cols["Y"] = Rcpp::NumericVector(vY.begin(), vY.begin() + idx);
        dt_cols["Z"] = Rcpp::NumericVector(vZ.begin(), vZ.begin() + idx);
    }
    if (sel.gpstime)           dt_cols["gpstime"]          = Rcpp::NumericVector(vGpstime.begin(), vGpstime.begin() + idx);
    if (sel.intensity)         dt_cols["Intensity"]        = Rcpp::IntegerVector(vIntensity.begin(), vIntensity.begin() + idx);
    if (sel.return_number)     dt_cols["ReturnNumber"]     = Rcpp::IntegerVector(vReturnNumber.begin(), vReturnNumber.begin() + idx);
    if (sel.number_of_returns) dt_cols["NumberOfReturns"]  = Rcpp::IntegerVector(vNumberOfReturns.begin(), vNumberOfReturns.begin() + idx);
    if (sel.classification)    dt_cols["Classification"]   = Rcpp::IntegerVector(vClassification.begin(), vClassification.begin() + idx);
    if (sel.scan_angle)        dt_cols["ScanAngleRank"]    = Rcpp::NumericVector(vScanAngle.begin(), vScanAngle.begin() + idx);
    if (sel.user_data)         dt_cols["UserData"]         = Rcpp::IntegerVector(vUserData.begin(), vUserData.begin() + idx);
    if (sel.point_source_id)   dt_cols["PointSourceID"]    = Rcpp::IntegerVector(vPointSourceID.begin(), vPointSourceID.begin() + idx);
    if (sel.flags) {
        dt_cols["ScanDirectionFlag"]  = Rcpp::IntegerVector(vScanDirFlag.begin(), vScanDirFlag.begin() + idx);
        dt_cols["EdgeOfFlightline"]   = Rcpp::IntegerVector(vEdgeOfFlight.begin(), vEdgeOfFlight.begin() + idx);
        dt_cols["Synthetic_flag"]     = Rcpp::IntegerVector(idx, 0); // derive from classification_flags bit 0
        dt_cols["Keypoint_flag"]      = Rcpp::IntegerVector(idx, 0); // bit 1
        dt_cols["Withheld_flag"]      = Rcpp::IntegerVector(idx, 0); // bit 2
        dt_cols["Overlap_flag"]       = Rcpp::IntegerVector(idx, 0); // bit 3
        dt_cols["ScannerChannel"]     = Rcpp::IntegerVector(vScannerChannel.begin(), vScannerChannel.begin() + idx);

        // Populate flag columns from classification_flags
        Rcpp::IntegerVector synv = dt_cols["Synthetic_flag"];
        Rcpp::IntegerVector keyv = dt_cols["Keypoint_flag"];
        Rcpp::IntegerVector witv = dt_cols["Withheld_flag"];
        Rcpp::IntegerVector ovlv = dt_cols["Overlap_flag"];
        for (uint64_t j = 0; j < idx; j++) {
            int cf = vClassFlags[j];
            synv[j] = (cf >> 0) & 1;
            keyv[j] = (cf >> 1) & 1;
            witv[j] = (cf >> 2) & 1;
            ovlv[j] = (cf >> 3) & 1;
        }
    }
    if (sel.rgb) {
        dt_cols["R"] = Rcpp::IntegerVector(vR.begin(), vR.begin() + idx);
        dt_cols["G"] = Rcpp::IntegerVector(vG.begin(), vG.begin() + idx);
        dt_cols["B"] = Rcpp::IntegerVector(vB.begin(), vB.begin() + idx);
    }
    if (sel.nir) {
        dt_cols["NIR"] = Rcpp::IntegerVector(vNIR.begin(), vNIR.begin() + idx);
    }

    // Extra Bytes columns
    for (size_t e = 0; e < extra_bytes_defs.size(); e++) {
        dt_cols[extra_bytes_defs[e].name] =
            Rcpp::NumericVector(vExtraCols[e].begin(),
                                vExtraCols[e].begin() + idx);
    }

    // ── 6. Build return list ──────────────────────────────────────────
    dt_cols.attr("class") = Rcpp::CharacterVector::create("data.table", "data.frame");
    dt_cols.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER, -(int)idx);

    Rcpp::List result;
    result["header"] = header_to_rlas_list(header);
    result["data"]   = dt_cols;

    return result;
}
