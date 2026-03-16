// copc_writer.cpp  –  Write COPC LAZ files from R data
//
// Creates a valid COPC (Cloud Optimized Point Cloud) file from a
// data.table of points and an rlas-style header.  Points are organized
// into an octree whose compressed chunks are seekable by any COPC reader.
//
// The point data section is written as a SINGLE standard LAZ stream
// (8-byte chunk table header + chunks + chunk table), with each COPC
// octree node occupying exactly one LAZ chunk.  This makes the output
// compatible with all LAZ/COPC readers (CloudCompare, PDAL, LAStools,
// QGIS, etc.) — non-COPC readers see a valid LAZ file; COPC-aware
// readers additionally use the hierarchy EVLR for spatial indexing.
//
// Uses the bundled LASzip library for LAZ compression.

#include <Rcpp.h>
#include "las14_header.h"

#include <cstring>
#include <cmath>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <cstdio>
#include <ctime>

// LASzip internal headers
#include "laszip.hpp"
#include "laswritepoint.hpp"
#include "bytestreamout_file.hpp"   // ByteStreamOutFileLE (file-based)
#include "laswriteitemraw.hpp"      // LAStempWritePoint10
#include "mydefs.hpp"               // ftell_las, fseek_las

using copc4r::LAS14PublicHeader;
using copc4r::VLRHeader;
using copc4r::EVLRHeader;

// ═══════════════════════════════════════════════════════════════════════
// Octree helpers
// ═══════════════════════════════════════════════════════════════════════
namespace {

struct VoxelKey {
    int32_t d, x, y, z;
    bool operator==(const VoxelKey& o) const {
        return d == o.d && x == o.x && y == o.y && z == o.z;
    }
};

struct VoxelKeyHash {
    size_t operator()(const VoxelKey& k) const {
        size_t h = std::hash<int32_t>{}(k.d);
        h ^= std::hash<int32_t>{}(k.x) + 0x9e3779b9 + (h << 6) + (h >> 2);
        h ^= std::hash<int32_t>{}(k.y) + 0x9e3779b9 + (h << 6) + (h >> 2);
        h ^= std::hash<int32_t>{}(k.z) + 0x9e3779b9 + (h << 6) + (h >> 2);
        return h;
    }
};

struct HierEntry {
    VoxelKey key;
    int64_t  offset;
    int32_t  byte_count;
    int32_t  point_count;
};

VoxelKey point_to_key(double px, double py, double pz,
                      double cx, double cy, double cz,
                      double halfsize, int depth) {
    if (depth <= 0) return {0, 0, 0, 0};
    double step = (2.0 * halfsize) / (1 << depth);
    int mx = (1 << depth) - 1;
    int ix = std::max(0, std::min((int)std::floor((px - (cx - halfsize)) / step), mx));
    int iy = std::max(0, std::min((int)std::floor((py - (cy - halfsize)) / step), mx));
    int iz = std::max(0, std::min((int)std::floor((pz - (cz - halfsize)) / step), mx));
    return {depth, ix, iy, iz};
}

bool key_bfs_cmp(const VoxelKey& a, const VoxelKey& b) {
    if (a.d != b.d) return a.d < b.d;
    if (a.x != b.x) return a.x < b.x;
    if (a.y != b.y) return a.y < b.y;
    return a.z < b.z;
}

// Safe column getter (returns n-length zero vector if missing)
Rcpp::IntegerVector get_int_col(const Rcpp::DataFrame& df,
                                const char* name, int n) {
    if (df.containsElementNamed(name))
        return Rcpp::as<Rcpp::IntegerVector>(df[name]);
    return Rcpp::IntegerVector(n, 0);
}
Rcpp::NumericVector get_dbl_col(const Rcpp::DataFrame& df,
                                const char* name, int n) {
    if (df.containsElementNamed(name))
        return Rcpp::as<Rcpp::NumericVector>(df[name]);
    return Rcpp::NumericVector(n, 0.0);
}

// Safe header getter
double hdr_dbl(const Rcpp::List& h, const char* name, double def) {
    if (h.containsElementNamed(name)) return Rcpp::as<double>(h[name]);
    return def;
}
int hdr_int(const Rcpp::List& h, const char* name, int def) {
    if (h.containsElementNamed(name)) return Rcpp::as<int>(h[name]);
    return def;
}

} // anonymous namespace

// ═══════════════════════════════════════════════════════════════════════
// cpp_write_copc  –  Rcpp-exported function
// ═══════════════════════════════════════════════════════════════════════

// [[Rcpp::export]]
void cpp_write_copc(std::string filename,
                    Rcpp::DataFrame data,
                    Rcpp::List header,
                    int max_depth_input,
                    bool progress) {

    // ── 1. Extract point data ───────────────────────────────────────
    Rcpp::NumericVector X = data["X"];
    Rcpp::NumericVector Y = data["Y"];
    Rcpp::NumericVector Z = data["Z"];
    int n = X.size();
    if (n == 0) Rcpp::stop("No points to write");

    Rcpp::IntegerVector vIntensity      = get_int_col(data, "Intensity", n);
    Rcpp::IntegerVector vReturnNumber   = get_int_col(data, "ReturnNumber", n);
    Rcpp::IntegerVector vNumReturns     = get_int_col(data, "NumberOfReturns", n);
    Rcpp::IntegerVector vClassification = get_int_col(data, "Classification", n);
    Rcpp::IntegerVector vUserData       = get_int_col(data, "UserData", n);
    Rcpp::IntegerVector vPointSourceID  = get_int_col(data, "PointSourceID", n);
    Rcpp::IntegerVector vScanDirFlag    = get_int_col(data, "ScanDirectionFlag", n);
    Rcpp::IntegerVector vEdgeOfFlight   = get_int_col(data, "EdgeOfFlightline", n);
    Rcpp::IntegerVector vScannerChannel = get_int_col(data, "ScannerChannel", n);
    Rcpp::NumericVector vScanAngle      = get_dbl_col(data, "ScanAngleRank", n);
    Rcpp::NumericVector vGpstime        = get_dbl_col(data, "gpstime", n);

    // Classification flags from individual flag columns
    Rcpp::IntegerVector vSynthetic = get_int_col(data, "Synthetic_flag", n);
    Rcpp::IntegerVector vKeypoint  = get_int_col(data, "Keypoint_flag", n);
    Rcpp::IntegerVector vWithheld  = get_int_col(data, "Withheld_flag", n);
    Rcpp::IntegerVector vOverlap   = get_int_col(data, "Overlap_flag", n);

    // RGB / NIR
    bool has_rgb = data.containsElementNamed("R") &&
                   data.containsElementNamed("G") &&
                   data.containsElementNamed("B");
    bool has_nir = data.containsElementNamed("NIR");

    Rcpp::IntegerVector vR   = has_rgb ? Rcpp::as<Rcpp::IntegerVector>(data["R"])   : Rcpp::IntegerVector(n, 0);
    Rcpp::IntegerVector vG   = has_rgb ? Rcpp::as<Rcpp::IntegerVector>(data["G"])   : Rcpp::IntegerVector(n, 0);
    Rcpp::IntegerVector vB   = has_rgb ? Rcpp::as<Rcpp::IntegerVector>(data["B"])   : Rcpp::IntegerVector(n, 0);
    Rcpp::IntegerVector vNIR = has_nir ? Rcpp::as<Rcpp::IntegerVector>(data["NIR"]) : Rcpp::IntegerVector(n, 0);

    // ── 2. Header metadata ──────────────────────────────────────────
    uint8_t point_format = (uint8_t)hdr_int(header, "Point Data Format ID",
                                             has_nir ? 8 : (has_rgb ? 7 : 6));
    if (point_format < 6 || point_format > 8)
        Rcpp::stop("COPC requires point format 6, 7, or 8 (got %d)", (int)point_format);

    uint16_t point_record_length;
    switch (point_format) {
        case 6:  point_record_length = 30; break;
        case 7:  point_record_length = 36; break;
        default: point_record_length = 38; break;
    }

    double x_scale  = hdr_dbl(header, "X scale factor", 0.001);
    double y_scale  = hdr_dbl(header, "Y scale factor", 0.001);
    double z_scale  = hdr_dbl(header, "Z scale factor", 0.001);
    double x_offset = hdr_dbl(header, "X offset", 0.0);
    double y_offset = hdr_dbl(header, "Y offset", 0.0);
    double z_offset = hdr_dbl(header, "Z offset", 0.0);

    // ── 3. Compute bounds ───────────────────────────────────────────
    double xmin = X[0], xmax = X[0];
    double ymin = Y[0], ymax = Y[0];
    double zmin = Z[0], zmax = Z[0];
    double tmin = vGpstime[0], tmax = vGpstime[0];

    for (int i = 1; i < n; i++) {
        if (X[i] < xmin) xmin = X[i];
        if (X[i] > xmax) xmax = X[i];
        if (Y[i] < ymin) ymin = Y[i];
        if (Y[i] > ymax) ymax = Y[i];
        if (Z[i] < zmin) zmin = Z[i];
        if (Z[i] > zmax) zmax = Z[i];
        if (vGpstime[i] < tmin) tmin = vGpstime[i];
        if (vGpstime[i] > tmax) tmax = vGpstime[i];
    }

    uint64_t pts_by_ret[15] = {0};
    for (int i = 0; i < n; i++) {
        int rn = vReturnNumber[i];
        if (rn >= 1 && rn <= 15) pts_by_ret[rn - 1]++;
    }

    // ── 4. Build octree ─────────────────────────────────────────────
    double cx = (xmin + xmax) / 2.0;
    double cy = (ymin + ymax) / 2.0;
    double cz = (zmin + zmax) / 2.0;
    double halfsize = std::max({xmax - xmin, ymax - ymin, zmax - zmin}) / 2.0;
    if (halfsize < 1e-6) halfsize = 1.0;
    halfsize *= 1.001;  // small margin to prevent boundary issues

    // Auto-compute depth (target ~10 000 points per leaf)
    int max_depth;
    if (max_depth_input >= 0) {
        max_depth = max_depth_input;
    } else {
        if (n <= 10000) {
            max_depth = 0;
        } else {
            max_depth = (int)std::ceil(std::log((double)n / 10000.0) / std::log(8.0));
        }
        max_depth = std::min(max_depth, 8);
    }

    double spacing = 2.0 * halfsize;
    if (header.containsElementNamed("COPC Info")) {
        Rcpp::List ci = Rcpp::as<Rcpp::List>(header["COPC Info"]);
        if (ci.containsElementNamed("spacing")) {
            double s = Rcpp::as<double>(ci["spacing"]);
            if (s > 0) spacing = s;
        }
    }

    // Assign points to leaf nodes
    std::unordered_map<VoxelKey, std::vector<int>, VoxelKeyHash> node_pts;
    for (int i = 0; i < n; i++) {
        VoxelKey k = point_to_key(X[i], Y[i], Z[i], cx, cy, cz,
                                  halfsize, max_depth);
        node_pts[k].push_back(i);
    }

    // Collect keys in breadth-first order
    std::vector<VoxelKey> sorted_keys;
    sorted_keys.reserve(node_pts.size());
    for (auto& kv : node_pts) sorted_keys.push_back(kv.first);
    std::sort(sorted_keys.begin(), sorted_keys.end(), key_bfs_cmp);

    if (progress)
        Rprintf("write_copc: %d points, %d nodes (depth %d)\n",
                n, (int)sorted_keys.size(), max_depth);

    // ── 5. Configure LASzip compressor ──────────────────────────────
    LASzip laszip;
    if (!laszip.setup(point_format, point_record_length,
                      LASZIP_COMPRESSOR_LAYERED_CHUNKED)) {
        Rcpp::stop("LASzip setup failed: %s",
                   laszip.get_error() ? laszip.get_error() : "unknown");
    }
    // Request v2 which clamps POINT14/RGB14/etc. to max(3,2)=3 — the
    // standard version used by PDAL's COPC writer and supported by all
    // LASzip builds (v4 requires LASzip ≥ 3.5 which rlas/lidR may lack).
    if (!laszip.request_version(2)) {
        Rcpp::stop("LASzip request_version(2) failed: %s",
                   laszip.get_error() ? laszip.get_error() : "unknown");
    }
    // chunk_size = 0 leaves the internal default at U32_MAX, which allows
    // us to call writer.chunk() manually to demarcate COPC nodes.
    laszip.chunk_size = 0;

    unsigned char* lz_vlr = nullptr;
    int lz_vlr_sz = 0;
    if (!laszip.pack(lz_vlr, lz_vlr_sz))
        Rcpp::stop("LASzip pack failed: %s",
                   laszip.get_error() ? laszip.get_error() : "unknown");

    // ── 6. Extract WKT CRS from header VLRs ────────────────────────
    std::string wkt;
    if (header.containsElementNamed("Variable Length Records")) {
        Rcpp::List vlrs = Rcpp::as<Rcpp::List>(header["Variable Length Records"]);
        Rcpp::CharacterVector vlr_names;
        if (vlrs.hasAttribute("names"))
            vlr_names = vlrs.attr("names");
        for (int i = 0; i < vlrs.size(); i++) {
            Rcpp::List v = Rcpp::as<Rcpp::List>(vlrs[i]);
            if (v.containsElementNamed("WKT OGC COORDINATE SYSTEM")) {
                wkt = Rcpp::as<std::string>(v["WKT OGC COORDINATE SYSTEM"]);
                break;
            }
        }
    }

    // ── 7. Compute VLR layout ───────────────────────────────────────
    // Order: COPC info → WKT CRS → LASzip (matching PDAL's output order)
    uint32_t num_vlrs = 2;  // COPC info + LASzip
    uint32_t vlr_bytes = (54 + 160) + (54 + (uint32_t)lz_vlr_sz);
    if (!wkt.empty()) {
        num_vlrs++;
        vlr_bytes += 54 + (uint32_t)wkt.size() + 1;
    }
    uint32_t offset_to_pts = 375 + vlr_bytes;

    // ── 8. Open file & write LAS 1.4 header ────────────────────────
    FILE* fp = std::fopen(filename.c_str(), "wb");
    if (!fp) Rcpp::stop("Cannot open '%s' for writing", filename.c_str());

    LAS14PublicHeader hdr;
    std::memset(&hdr, 0, sizeof(hdr));
    std::memcpy(hdr.file_signature, "LASF", 4);
    hdr.file_source_id   = (uint16_t)hdr_int(header, "File Source ID", 0);
    hdr.global_encoding  = 0x11;  // GPS time type (bit 0) + WKT (bit 4)
    hdr.version_major    = 1;
    hdr.version_minor    = 4;
    std::strncpy(hdr.system_identifier, "copc4R", 31);
    std::strncpy(hdr.generating_software, "copc4R (R)", 31);
    {
        std::time_t now = std::time(nullptr);
        std::tm* t = std::localtime(&now);
        hdr.creation_day  = (uint16_t)(t->tm_yday + 1);
        hdr.creation_year = (uint16_t)(t->tm_year + 1900);
    }
    hdr.header_size            = 375;
    hdr.offset_to_point_data   = offset_to_pts;
    hdr.num_vlrs               = num_vlrs;
    hdr.point_data_format      = point_format | 0x80;  // 0x80 = LAZ compressed
    hdr.point_data_record_length = point_record_length;
    hdr.legacy_point_count     = 0;
    hdr.x_scale = x_scale;  hdr.y_scale = y_scale;  hdr.z_scale = z_scale;
    hdr.x_offset = x_offset; hdr.y_offset = y_offset; hdr.z_offset = z_offset;
    hdr.max_x = xmax; hdr.min_x = xmin;
    hdr.max_y = ymax; hdr.min_y = ymin;
    hdr.max_z = zmax; hdr.min_z = zmin;
    hdr.waveform_offset = 0;
    hdr.evlr_start  = 0;
    hdr.num_evlrs   = 0;  // LASlib's EPToctree COPC handling has a
                           // bug that silently drops all points when
                           // it finds a COPC hierarchy EVLR.  We set
                           // num_evlrs = 0 so LASlib skips EVLRs.
                           // The hierarchy data is still written at
                           // the end of the file and discoverable via
                           // root_hier_offset in the COPC info VLR.
    hdr.point_count_14 = (uint64_t)n;
    std::memcpy(hdr.points_by_return_14, pts_by_ret, sizeof(pts_by_ret));

    std::fwrite(&hdr, sizeof(hdr), 1, fp);

    // ── 9. Write VLRs ───────────────────────────────────────────────
    // Track COPC info data offset for later patching
    int64_t copc_info_data_pos = ftell_las(fp) + 54;

    // a) COPC info VLR (160-byte payload, root_hier_* patched later)
    {
        VLRHeader vh;
        std::memset(&vh, 0, sizeof(vh));
        std::strncpy(vh.user_id, "copc", 15);
        vh.record_id     = 1;
        vh.record_length = 160;
        std::strncpy(vh.description, "COPC info", 31);
        std::fwrite(&vh, sizeof(vh), 1, fp);

        uint8_t buf[160];
        std::memset(buf, 0, 160);
        std::memcpy(buf +  0, &cx, 8);
        std::memcpy(buf +  8, &cy, 8);
        std::memcpy(buf + 16, &cz, 8);
        std::memcpy(buf + 24, &halfsize, 8);
        std::memcpy(buf + 32, &spacing, 8);
        // root_hier_offset (40) and root_hier_size (48): patched later
        std::memcpy(buf + 56, &tmin, 8);
        std::memcpy(buf + 64, &tmax, 8);
        std::fwrite(buf, 160, 1, fp);
    }

    // b) CRS WKT VLR (optional — written before LASzip to match PDAL order)
    if (!wkt.empty()) {
        uint16_t wkt_len = (uint16_t)(wkt.size() + 1);
        VLRHeader vh;
        std::memset(&vh, 0, sizeof(vh));
        std::memcpy(vh.user_id, "LASF_Projection", 15);
        vh.record_id     = 2112;
        vh.record_length = wkt_len;
        std::strncpy(vh.description, "OGC WKT", 31);
        std::fwrite(&vh, sizeof(vh), 1, fp);
        std::fwrite(wkt.c_str(), wkt_len, 1, fp);
    }

    // c) LASzip VLR (must be last — some readers scan VLRs sequentially)
    {
        VLRHeader vh;
        std::memset(&vh, 0, sizeof(vh));
        std::memcpy(vh.user_id, "laszip encoded", 15);
        vh.record_id     = 22204;
        vh.record_length = (uint16_t)lz_vlr_sz;
        std::strncpy(vh.description, "http://laszip.org", 31);
        std::fwrite(&vh, sizeof(vh), 1, fp);
        std::fwrite(lz_vlr, lz_vlr_sz, 1, fp);
    }

    // Verify offset
    if (ftell_las(fp) != (I64)offset_to_pts) {
        std::fclose(fp);
        Rcpp::stop("Internal error: VLR offset mismatch (%lld vs %u)",
                   (long long)ftell_las(fp), (unsigned)offset_to_pts);
    }

    // ── 10. Write compressed point data as a single LAZ stream ──────
    //
    // We write ALL nodes through one LASwritePoint instance using a
    // file-based ByteStreamOut.  Each COPC node becomes one LAZ chunk,
    // delimited by manual writer.chunk() calls.  This produces:
    //
    //   [8-byte chunk-table offset] [chunk₀] [chunk₁] … [chunkₙ] [chunk table]
    //
    // which is a standard LAZ stream that any reader (CloudCompare,
    // PDAL, LAStools, QGIS) can decompress.  The COPC hierarchy EVLR
    // additionally records per-node offsets for spatial indexing.

    ByteStreamOutFileLE* file_stream = new ByteStreamOutFileLE(fp);

    LASwritePoint writer;
    if (!writer.setup(laszip.num_items, laszip.items, &laszip)) {
        delete file_stream;
        std::fclose(fp);
        Rcpp::stop("LASwritePoint::setup() failed");
    }
    // Do NOT call disable_chunk_table() — we want the standard 8-byte
    // header + chunk table so that non-COPC readers can parse the file.
    if (!writer.init(file_stream)) {
        delete file_stream;
        std::fclose(fp);
        Rcpp::stop("LASwritePoint::init() failed");
    }

    // Allocate combo-format point buffer (reused across all nodes)
    uint32_t combo_sz = 0;
    for (U16 i = 0; i < laszip.num_items; i++)
        combo_sz += 2 * laszip.items[i].size;

    std::vector<uint8_t> pbuf(combo_sz, 0);
    std::vector<U8*> iptrs(laszip.num_items);
    {
        uint32_t off = 0;
        for (U16 i = 0; i < laszip.num_items; i++) {
            iptrs[i] = pbuf.data() + off;
            off += 2 * laszip.items[i].size;
        }
    }

    std::vector<HierEntry> hierarchy;
    hierarchy.reserve(sorted_keys.size());

    for (size_t ni = 0; ni < sorted_keys.size(); ni++) {
        const VoxelKey& key = sorted_keys[ni];
        const std::vector<int>& pts = node_pts[key];
        int nc = (int)pts.size();

        // The chunk data starts at the current stream position
        I64 node_start = file_stream->tell();

        for (int pi = 0; pi < nc; pi++) {
            int idx = pts[pi];

            // Zero the buffer
            std::memset(pbuf.data(), 0, combo_sz);

            // ── POINT14 combo (LAStempWritePoint10 layout) ──────────
            LAStempWritePoint10* pt =
                reinterpret_cast<LAStempWritePoint10*>(iptrs[0]);

            pt->X = (I32)std::round((X[idx] - x_offset) / x_scale);
            pt->Y = (I32)std::round((Y[idx] - y_offset) / y_scale);
            pt->Z = (I32)std::round((Z[idx] - z_offset) / z_scale);
            pt->intensity = (U16)vIntensity[idx];

            // Legacy bitfield (byte 14)
            pt->return_number      = (U8)(vReturnNumber[idx] & 0x07);
            pt->number_of_returns  = (U8)(vNumReturns[idx]   & 0x07);
            pt->scan_direction_flag= (U8)(vScanDirFlag[idx]  & 0x01);
            pt->edge_of_flight_line= (U8)(vEdgeOfFlight[idx] & 0x01);

            pt->classification  = (U8)vClassification[idx] |
                                  ((U8)(vSynthetic[idx] & 1) << 5) |
                                  ((U8)(vKeypoint[idx]  & 1) << 6) |
                                  ((U8)(vWithheld[idx]  & 1) << 7);
            // ScanAngleRank is in degrees; legacy rank is I8 degrees
            pt->scan_angle_rank = (I8)std::round(vScanAngle[idx]);
            pt->user_data       = (U8)vUserData[idx];
            pt->point_source_ID = (U16)vPointSourceID[idx];

            // Extended LAS 1.4 fields
            // Scan angle: ScanAngleRank (degrees) → raw I16 = deg / 0.006
            pt->extended_scan_angle  =
                (I16)std::round(vScanAngle[idx] / 0.006);
            pt->extended_point_type  = 1;  // must be 1 for POINT14 wire format
            pt->extended_scanner_channel = (U8)(vScannerChannel[idx] & 0x03);
            pt->extended_classification_flags =
                (U8)(((vSynthetic[idx] & 1))      |
                     ((vKeypoint[idx]  & 1) << 1)  |
                     ((vWithheld[idx]  & 1) << 2)  |
                     ((vOverlap[idx]   & 1) << 3));
            pt->extended_classification = (U8)vClassification[idx];
            pt->extended_return_number  = (U8)(vReturnNumber[idx] & 0x0F);
            pt->extended_number_of_returns = (U8)(vNumReturns[idx] & 0x0F);
            pt->gps_time = vGpstime[idx];

            // ── RGB / RGBNIR item ───────────────────────────────────
            if (point_format >= 7 && laszip.num_items >= 2) {
                U8* rgb = iptrs[1];
                U16 r = (U16)vR[idx], g = (U16)vG[idx], b = (U16)vB[idx];
                std::memcpy(rgb + 0, &r, 2);
                std::memcpy(rgb + 2, &g, 2);
                std::memcpy(rgb + 4, &b, 2);
                if (point_format >= 8) {
                    U16 nir = (U16)vNIR[idx];
                    std::memcpy(rgb + 6, &nir, 2);
                }
            }

            writer.write((const U8* const*)iptrs.data());
        }

        // Force a chunk boundary so the next node starts a new chunk.
        // chunk() flushes the current chunk's layered data and records
        // its byte count in the internal chunk table.
        writer.chunk();

        I64 node_end = file_stream->tell();
        int32_t node_bytes = (int32_t)(node_end - node_start);

        hierarchy.push_back({key, node_start, node_bytes, nc});

        if (progress && (ni % 50 == 0 || ni + 1 == sorted_keys.size())) {
            Rprintf("\rwrite_copc: node %d/%d",
                    (int)(ni + 1), (int)sorted_keys.size());
        }
    }

    // Finalize: flush last encoder state + write the master chunk table.
    // Because we called chunk() after the last node, the writer is in a
    // clean state with 0 pending points — done() just writes the table.
    writer.done();

    if (progress) Rprintf("\n");

    delete file_stream;
    file_stream = nullptr;

    // ── 11. Write hierarchy EVLR ────────────────────────────────────

    EVLRHeader ev;
    std::memset(&ev, 0, sizeof(ev));
    std::strncpy(ev.user_id, "copc", 15);
    ev.record_id     = 1000;
    ev.record_length = (uint64_t)(hierarchy.size() * 32);
    std::strncpy(ev.description, "EPT Hierarchy", 31);
    std::fwrite(&ev, sizeof(ev), 1, fp);

    I64 hier_data_pos = ftell_las(fp);

    for (auto& he : hierarchy) {
        std::fwrite(&he.key.d, 4, 1, fp);
        std::fwrite(&he.key.x, 4, 1, fp);
        std::fwrite(&he.key.y, 4, 1, fp);
        std::fwrite(&he.key.z, 4, 1, fp);
        std::fwrite(&he.offset, 8, 1, fp);
        std::fwrite(&he.byte_count, 4, 1, fp);
        std::fwrite(&he.point_count, 4, 1, fp);
    }

    int64_t hier_data_size = (int64_t)(hierarchy.size() * 32);

    // ── 12. Patch COPC info VLR: root_hier_offset + root_hier_size ──
    fseek_las(fp, copc_info_data_pos + 40, SEEK_SET);
    std::fwrite(&hier_data_pos, 8, 1, fp);
    std::fwrite(&hier_data_size, 8, 1, fp);

    std::fclose(fp);

    // Note: lz_vlr is owned by the laszip object (laszip.bytes) and
    // will be freed by its destructor — do NOT delete[] lz_vlr here.

    if (progress)
        Rprintf("write_copc: wrote %d points (%d nodes) to %s\n",
                n, (int)hierarchy.size(), filename.c_str());
}
