// laszip_decompress.h  –  Decompress LAZ-compressed chunks from memory
#ifndef COPC4R_LASZIP_DECOMPRESS_H
#define COPC4R_LASZIP_DECOMPRESS_H

#include <cstdint>
#include <vector>
#include <string>
#include "las14_header.h"

namespace copc4r {

// ═══════════════════════════════════════════════════════════════════════
// Decoded point record  –  covers all fields for PDRF 6/7/8.
// Unused fields stay zero.
// ═══════════════════════════════════════════════════════════════════════
struct DecodedPoint {
    double   X, Y, Z;
    uint16_t intensity;
    uint8_t  return_number;
    uint8_t  number_of_returns;
    uint8_t  classification_flags;
    uint8_t  scanner_channel;
    uint8_t  scan_direction_flag;
    uint8_t  edge_of_flight_line;
    uint8_t  classification;
    uint8_t  user_data;
    int16_t  scan_angle;          // raw; multiply by 0.006 for degrees
    uint16_t point_source_id;
    double   gpstime;
    uint16_t red, green, blue;
    uint16_t nir;
};

// ═══════════════════════════════════════════════════════════════════════
// decompress_chunk()
//
// Takes a raw LAZ-compressed chunk and decompresses it into points.
//
// The LASzip VLR data (from the file header) is needed to initialize
// the decompressor with correct item configuration.
// ═══════════════════════════════════════════════════════════════════════
std::vector<DecodedPoint> decompress_chunk(
    const std::vector<uint8_t>& compressed,
    const std::vector<uint8_t>& laszip_vlr_data,
    uint8_t  point_format,
    uint16_t point_record_length,
    uint64_t point_count,
    double x_scale, double y_scale, double z_scale,
    double x_off,   double y_off,   double z_off);

// ═══════════════════════════════════════════════════════════════════════
// decode_raw_point()  –  (unused legacy; kept for potential non-LAS1.4 use)
// ═══════════════════════════════════════════════════════════════════════

} // namespace copc4r
#endif
