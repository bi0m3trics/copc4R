// las14_header.h  –  Minimal LAS 1.4 header + VLR parsing
#ifndef COPC4R_LAS14_HEADER_H
#define COPC4R_LAS14_HEADER_H

#include <cstdint>
#include <string>
#include <vector>
#include <cstring>

namespace copc4r {

// ── Fixed-size LAS 1.4 public header block (375 bytes) ────────────────
#pragma pack(push, 1)
struct LAS14PublicHeader {
    char     file_signature[4];       //  0  "LASF"
    uint16_t file_source_id;          //  4
    uint16_t global_encoding;         //  6
    uint32_t guid_data1;              //  8
    uint16_t guid_data2;              // 12
    uint16_t guid_data3;              // 14
    uint8_t  guid_data4[8];           // 16
    uint8_t  version_major;           // 24
    uint8_t  version_minor;           // 25
    char     system_identifier[32];   // 26
    char     generating_software[32]; // 58
    uint16_t creation_day;            // 90
    uint16_t creation_year;           // 92
    uint16_t header_size;             // 94
    uint32_t offset_to_point_data;    // 96
    uint32_t num_vlrs;                //100
    uint8_t  point_data_format;       //104
    uint16_t point_data_record_length;//105
    uint32_t legacy_point_count;      //107
    uint32_t legacy_points_by_return[5]; //111
    double   x_scale;                 //131
    double   y_scale;                 //139
    double   z_scale;                 //147
    double   x_offset;               //155
    double   y_offset;               //163
    double   z_offset;               //171
    double   max_x;                  //179
    double   min_x;                  //187
    double   max_y;                  //195
    double   min_y;                  //203
    double   max_z;                  //211
    double   min_z;                  //219
    // LAS 1.3+ : waveform data packet record start
    uint64_t waveform_offset;        //227
    // LAS 1.4 :
    uint64_t evlr_start;             //235
    uint32_t num_evlrs;              //243
    uint64_t point_count_14;         //247
    uint64_t points_by_return_14[15];//255  (15 * 8 = 120 bytes => ends at 375)
};
#pragma pack(pop)

static_assert(sizeof(LAS14PublicHeader) == 375,
              "LAS 1.4 header struct must be exactly 375 bytes");

// ── VLR / EVLR record header ──────────────────────────────────────────
#pragma pack(push, 1)
struct VLRHeader {
    uint16_t reserved;
    char     user_id[16];
    uint16_t record_id;
    uint16_t record_length;    // VLR: payload length (max 65535)
    char     description[32];
};
struct EVLRHeader {
    uint16_t reserved;
    char     user_id[16];
    uint16_t record_id;
    uint64_t record_length;    // EVLR: 8-byte length
    char     description[32];
};
#pragma pack(pop)

static_assert(sizeof(VLRHeader)  == 54, "VLR header must be 54 bytes");
static_assert(sizeof(EVLRHeader) == 60, "EVLR header must be 60 bytes");

// ── Parsed VLR (we store header + raw payload) ────────────────────────
struct ParsedVLR {
    uint16_t    reserved;
    std::string user_id;
    uint16_t    record_id;
    std::string description;
    std::vector<uint8_t> data;   // raw payload bytes
};

// ── High-level parsed header ──────────────────────────────────────────
struct LASHeader {
    LAS14PublicHeader pub;
    std::vector<ParsedVLR> vlrs;
    std::vector<ParsedVLR> evlrs;

    uint64_t point_count() const {
        return (pub.version_minor >= 4 && pub.point_count_14 > 0)
            ? pub.point_count_14
            : static_cast<uint64_t>(pub.legacy_point_count);
    }
};

// ── Parse helpers ─────────────────────────────────────────────────────
/// Read the LAS header + all VLRs from a RangeReader.
/// Also reads EVLRs if the header indicates them.
class RangeReader; // forward
LASHeader read_las_header(RangeReader& reader);

/// Find a VLR by user_id + record_id; returns nullptr if not found.
const ParsedVLR* find_vlr(const LASHeader& h,
                           const std::string& user_id,
                           uint16_t record_id);

} // namespace copc4r
#endif
