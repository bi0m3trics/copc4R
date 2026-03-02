// laszip_decompress.cpp  –  Decompress LAZ chunks from memory using LASzip
//
// In COPC, each hierarchy-selected node stores a self-contained LAZ chunk.
// The compression scheme is defined by the LASzip VLR in the LAS file
// header (user_id="laszip encoded", record_id=22204).
//
// Strategy:
//   1. Unpack the LASzip VLR to recover the num_items / LASitem array.
//   2. Create a ByteStreamInArray over the compressed chunk bytes.
//   3. Use LASreadPoint::setup(num_items, items, &laszip) then ::init()
//      to start decompression.
//   4. Read point by point; each read fills item buffers.
//   5. Reassemble item buffers into a flat raw point record and decode.
//
// LASzip internal headers used:
//   laszip.hpp          – LASzip class + LASitem
//   lasreadpoint.hpp    – LASreadPoint decompressor
//   bytestreamin_array.hpp – in-memory byte stream

#include "laszip_decompress.h"
#include <cstring>
#include <stdexcept>
#include <string>
#include <memory>

// LASzip internal headers (available because we bundle the source)
#include "laszip.hpp"
#include "lasreadpoint.hpp"
#include "bytestreamin_array.hpp"

namespace copc4r {

// ═══════════════════════════════════════════════════════════════════════
// decode_raw_point()  –  parse fields from a flat point record buffer
// ═══════════════════════════════════════════════════════════════════════
// decode_combo_point()  –  decode from LASzip's internal "combo" format
//
// LASzip internally converts POINT14 (PDRF 6/7/8) fields into a
// "LAStempReadPoint10" combo struct that merges legacy POINT10 fields
// with extended LAS 1.4 fields.  The output item buffer is TWICE the
// item size, and the layout is:
//
//   Offset  Size  Field
//   ──────  ────  ────────────────────────────────────────────────────
//     0      4    X  (I32)
//     4      4    Y  (I32)
//     8      4    Z  (I32)
//    12      2    intensity (U16)
//    14      1    legacy return/num_returns/scan_dir/edge (bitfield)
//    15      1    legacy classification (combo: flags<<5 | class)
//    16      1    scan_angle_rank (I8, already converted from I16)
//    17      1    user_data (U8)
//    18      2    point_source_ID (U16)
//    20      2    extended_scan_angle (I16, original PDRF 6 value)
//    22      1    ext_point_type(2)|scanner_channel(2)|class_flags(4)
//    23      1    extended_classification (U8, full 0-255)
//    24      1    ext_return_number(4)|ext_number_of_returns(4)
//    25      3    dummy padding
//    28      4    deleted_flag (unused)
//    32      8    gps_time (F64)
//   ──────────────────────────────────────────────────────────────────
//   Total: 40 bytes  (POINT14 item size = 30, doubled region = 60)
//
// RGB14 item follows at offset 2*POINT14.size = 60 in the combo buffer.
// RGB layout is raw 6 bytes: R(U16), G(U16), B(U16).
// RGBNIR14: 8 bytes: R, G, B, NIR (each U16).
// ═══════════════════════════════════════════════════════════════════════
DecodedPoint decode_combo_point(
    const uint8_t* combo_buf,
    const LASzip& laszip,
    uint8_t  point_format,
    double x_scale, double y_scale, double z_scale,
    double x_off,   double y_off,   double z_off)
{
    DecodedPoint p{};

    // ── POINT14 item (always item 0) → LAStempReadPoint10 combo ──────
    const uint8_t* pt = combo_buf;  // item 0 starts at offset 0

    int32_t xi, yi, zi;
    std::memcpy(&xi, pt + 0,  4);
    std::memcpy(&yi, pt + 4,  4);
    std::memcpy(&zi, pt + 8,  4);
    p.X = xi * x_scale + x_off;
    p.Y = yi * y_scale + y_off;
    p.Z = zi * z_scale + z_off;

    std::memcpy(&p.intensity, pt + 12, 2);

    // Extended fields at offsets 20+ in the combo struct
    std::memcpy(&p.scan_angle, pt + 20, 2);        // extended_scan_angle (I16)

    uint8_t byte22 = pt[22]; // ext_point_type(2) | scanner_channel(2) | class_flags(4)
    p.classification_flags = (byte22 >> 4) & 0x0F;
    p.scanner_channel      = (byte22 >> 2) & 0x03;
    // ext_point_type in bits 0-1 (not needed)

    p.classification = pt[23]; // extended_classification (full 0-255)

    uint8_t byte24 = pt[24]; // ext_return_number(4) | ext_number_of_returns(4)
    p.return_number     = byte24 & 0x0F;
    p.number_of_returns = (byte24 >> 4) & 0x0F;

    p.user_data = pt[17];
    std::memcpy(&p.point_source_id, pt + 18, 2);

    // scan_direction_flag and edge_of_flight_line from legacy byte 14
    uint8_t byte14 = pt[14];
    p.scan_direction_flag = (byte14 >> 6) & 0x01;
    p.edge_of_flight_line = (byte14 >> 7) & 0x01;

    // GPS time at offset 32 in the combo struct
    std::memcpy(&p.gpstime, pt + 32, 8);

    // ── RGB14 / RGBNIR14 items ───────────────────────────────────────
    // Items after POINT14 start at 2 * items[0].size from the base
    U16 rgb_offset = 2 * laszip.items[0].size;

    if (point_format >= 7 && laszip.num_items >= 2) {
        const uint8_t* rgb = combo_buf + rgb_offset;
        std::memcpy(&p.red,   rgb + 0, 2);
        std::memcpy(&p.green, rgb + 2, 2);
        std::memcpy(&p.blue,  rgb + 4, 2);

        // NIR for RGBNIR14
        if (point_format >= 8) {
            // Check if this item is RGBNIR14 (size 8) vs RGB14 (size 6)
            if (laszip.items[1].size >= 8) {
                std::memcpy(&p.nir, rgb + 6, 2);
            }
        }
    }

    return p;
}

// ═══════════════════════════════════════════════════════════════════════
// decompress_chunk()
// ═══════════════════════════════════════════════════════════════════════
std::vector<DecodedPoint> decompress_chunk(
    const std::vector<uint8_t>& compressed,
    const std::vector<uint8_t>& laszip_vlr_data,
    uint8_t  point_format,
    uint16_t point_record_length,
    uint64_t point_count,
    double x_scale, double y_scale, double z_scale,
    double x_off,   double y_off,   double z_off)
{
    if (point_format < 6 || point_format > 8)
        throw std::runtime_error("Unsupported point format " +
                                 std::to_string(point_format) +
                                 " (COPC requires 6/7/8)");

    if (laszip_vlr_data.empty())
        throw std::runtime_error("Missing LASzip VLR data – cannot decompress");

    // ── 1. Unpack the LASzip VLR to get compression item layout ───────
    LASzip laszip;
    if (!laszip.unpack(laszip_vlr_data.data(),
                       static_cast<int>(laszip_vlr_data.size()))) {
        throw std::runtime_error(
            std::string("LASzip::unpack() failed: ") +
            (laszip.get_error() ? laszip.get_error() : "unknown error"));
    }

    // ── 2. Create in-memory byte stream from the compressed chunk ─────
    ByteStreamInArrayLE* stream = new ByteStreamInArrayLE(
        compressed.data(), static_cast<I64>(compressed.size()));

    // ── 3. Create and initialize the point reader ─────────────────────
    LASreadPoint reader;
    if (!reader.setup(laszip.num_items, laszip.items, &laszip)) {
        delete stream;
        throw std::runtime_error(
            std::string("LASreadPoint::setup() failed: ") +
            (reader.error() ? reader.error() : "unknown"));
    }
    if (!reader.init(stream)) {
        delete stream;
        throw std::runtime_error(
            std::string("LASreadPoint::init() failed: ") +
            (reader.error() ? reader.error() : "unknown"));
    }

    // COPC: skip the chunk-table mechanism.
    reader.prepare_single_chunk();

    // ── 4. Prepare item buffers ───────────────────────────────────────
    // For layered LAS 1.4 compression, LASzip uses a "combo" format
    // where each item occupies TWICE its nominal size (the extra space
    // holds both legacy POINT10 and extended POINT14 fields).
    // We must allocate doubled buffers and use doubled offsets.
    uint32_t combo_size = 0;
    for (U16 i = 0; i < laszip.num_items; i++)
        combo_size += 2 * laszip.items[i].size;

    std::vector<uint8_t> raw_point(combo_size, 0);

    std::vector<U8*> item_ptrs(laszip.num_items);
    U32 offset = 0;
    for (U16 i = 0; i < laszip.num_items; i++) {
        item_ptrs[i] = raw_point.data() + offset;
        offset += 2 * laszip.items[i].size;
    }

    // ── 5. Read and decode each point ─────────────────────────────────
    std::vector<DecodedPoint> result;
    result.reserve(point_count);

    for (uint64_t i = 0; i < point_count; i++) {
        if (!reader.read(item_ptrs.data())) {
            break;
        }

        result.push_back(decode_combo_point(
            raw_point.data(), laszip, point_format,
            x_scale, y_scale, z_scale,
            x_off, y_off, z_off));
    }

    reader.done();
    delete stream;

    return result;
}

} // namespace copc4r
