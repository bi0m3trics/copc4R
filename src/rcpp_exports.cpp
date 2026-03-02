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

// [[Rcpp::export]]
Rcpp::List cpp_read_copc(std::string path_or_url,
                         Rcpp::NumericVector bbox,
                         Rcpp::NumericVector zrange,
                         std::string select,
                         double max_points_dbl,
                         bool progress) {
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
        zmin, zmax, has_bbox, has_zrange);

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
    };
    reserve(total_pts);

    uint64_t idx = 0;
    for (size_t ni = 0; ni < nodes.size(); ++ni) {
        auto& nc = nodes[ni];

        if (progress && ni % 50 == 0 && ni > 0) {
            Rcpp::Rcout << "  Decompressing chunk " << ni << "/" << nodes.size() << "\n";
        }

        // Read compressed chunk from file/URL
        auto compressed = reader->read(nc.offset, nc.byte_size);

        // Decompress
        std::vector<copc4r::DecodedPoint> pts;
        try {
            pts = copc4r::decompress_chunk(
                compressed, laszip_vlr_data, pdrf, prl, nc.point_count,
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

    // ── 6. Build return list ──────────────────────────────────────────
    dt_cols.attr("class") = Rcpp::CharacterVector::create("data.table", "data.frame");
    dt_cols.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER, -(int)idx);

    Rcpp::List result;
    result["header"] = header_to_rlas_list(header);
    result["data"]   = dt_cols;

    return result;
}
