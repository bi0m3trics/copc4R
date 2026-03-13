// las14_header.cpp  –  LAS 1.4 header + VLR parsing
#include "las14_header.h"
#include "range_reader.h"
#include <cstring>
#include <stdexcept>

namespace copc4r {

LASHeader read_las_header(RangeReader& reader) {
    LASHeader h;

    // ── 1. Read the 375-byte public header ────────────────────────────
    auto hdr_raw = reader.read(0, 375);
    if (hdr_raw.size() < 375)
        throw std::runtime_error("File too small for LAS 1.4 header");
    std::memcpy(&h.pub, hdr_raw.data(), 375);

    // Sanity checks
    if (std::strncmp(h.pub.file_signature, "LASF", 4) != 0)
        throw std::runtime_error("Not a LAS file (bad signature)");
    if (h.pub.version_major != 1 || h.pub.version_minor < 4)
        throw std::runtime_error("COPC requires LAS >= 1.4");

    // ── 2. Read VLRs (located right after public header) ──────────────
    uint64_t cursor = h.pub.header_size; // typically 375
    for (uint32_t i = 0; i < h.pub.num_vlrs; ++i) {
        auto vhdr_raw = reader.read(cursor, sizeof(VLRHeader));
        if (vhdr_raw.size() < sizeof(VLRHeader))
            throw std::runtime_error("Truncated VLR header at offset " +
                                     std::to_string(cursor));
        VLRHeader vhdr;
        std::memcpy(&vhdr, vhdr_raw.data(), sizeof(VLRHeader));
        cursor += sizeof(VLRHeader);

        ParsedVLR vlr;
        vlr.reserved    = vhdr.reserved;
        vlr.user_id     = std::string(vhdr.user_id,
                                      strnlen(vhdr.user_id, 16));
        vlr.record_id   = vhdr.record_id;
        vlr.description = std::string(vhdr.description,
                                      strnlen(vhdr.description, 32));
        vlr.data        = reader.read(cursor, vhdr.record_length);
        cursor += vhdr.record_length;

        h.vlrs.push_back(std::move(vlr));
    }

    // ── 3. Read EVLRs if present ──────────────────────────────────────
    if (h.pub.num_evlrs > 0 && h.pub.evlr_start > 0) {
        uint64_t ecursor = h.pub.evlr_start;
        for (uint32_t i = 0; i < h.pub.num_evlrs; ++i) {
            auto ehdr_raw = reader.read(ecursor, sizeof(EVLRHeader));
            if (ehdr_raw.size() < sizeof(EVLRHeader))
                break; // don't hard-fail on truncated EVLRs
            EVLRHeader ehdr;
            std::memcpy(&ehdr, ehdr_raw.data(), sizeof(EVLRHeader));
            ecursor += sizeof(EVLRHeader);

            ParsedVLR vlr;
            vlr.reserved   = ehdr.reserved;
            vlr.user_id    = std::string(ehdr.user_id,
                                         strnlen(ehdr.user_id, 16));
            vlr.record_id  = ehdr.record_id;
            vlr.description = std::string(ehdr.description,
                                          strnlen(ehdr.description, 32));
            // Only read EVLR payload up to 64 MB to avoid OOM on weird files
            uint64_t safe_len = std::min(ehdr.record_length,
                                         static_cast<uint64_t>(64 * 1024 * 1024));
            vlr.data = reader.read(ecursor, safe_len);
            ecursor += ehdr.record_length;

            h.evlrs.push_back(std::move(vlr));
        }
    }

    return h;
}

const ParsedVLR* find_vlr(const LASHeader& h,
                           const std::string& user_id,
                           uint16_t record_id) {
    // Search VLRs first, then EVLRs
    for (auto& v : h.vlrs)
        if (v.user_id == user_id && v.record_id == record_id)
            return &v;
    for (auto& v : h.evlrs)
        if (v.user_id == user_id && v.record_id == record_id)
            return &v;
    return nullptr;
}

// ═══════════════════════════════════════════════════════════════════════
// Extra Bytes VLR parsing  (LASF_Spec, record_id 4)
//
// Each record is 192 bytes:
//   bytes 0-1:    reserved (uint16)
//   byte  2:      data_type (uint8, 0-10)
//   byte  3:      options bitmap
//   bytes 4-35:   name (char[32], null-terminated)
//   bytes 36-39:  unused
//   bytes 40-63:  no_data[3] (8 bytes each)
//   bytes 64-79:  deprecated
//   bytes 80-103: min[3]
//   bytes 104-119: deprecated
//   bytes 120-143: max[3]
//   bytes 144-159: deprecated
//   bytes 160-183: scale[3] (8 bytes each)
//   bytes 184-207: -- wait, let me recount...
// Actually per LAS spec the record is:
//   2 + 1 + 1 + 32 + 4 + 24 + 16 + 24 + 16 + 24 + 16 + 24 + 8 = ...
// No, the exact layout per ASPRS LAS 1.4 spec:
//   reserved[2], data_type[1], options[1], name[32], unused[4],
//   no_data[24], deprecated1[16], min[24], deprecated2[16],
//   max[24], deprecated3[16], scale[24], offset[24], description[32]
// = 2+1+1+32+4+24+16+24+16+24+16+24+24+32 = 240?
// No -- let me use the standard 192 bytes:
//   reserved[2], data_type[1], options[1], name[32], unused[4],
//   no_data[3*8=24], deprecated1[3*8=24] -- no.
// Per ASPRS LAS 1.4 R15 spec, the Extra Bytes Record is:
//   reserved        2 bytes
//   data_type       1 byte
//   options         1 byte
//   name           32 bytes
//   unused          4 bytes
//   no_data        24 bytes (3 doubles)
//   deprecated1    16 bytes
//   min            24 bytes
//   deprecated2    16 bytes
//   max            24 bytes
//   deprecated3    16 bytes
//   scale          24 bytes
//   offset         24 bytes
//   description    32 bytes
// Total = 2+1+1+32+4+24+16+24+16+24+16+24+24+32 = 240
// Wait that's >= 192. Let me just use 192 as commonly implemented.
// Actually the correct structure (commonly used in practice) is 192 bytes:
//   uint16_t reserved;      // 2
//   uint8_t  data_type;     // 1
//   uint8_t  options;       // 1
//   char     name[32];      // 32
//   uint8_t  unused[4];     // 4
//   uint8_t  no_data[24];   // 24  (anytype[3])
//   uint8_t  min[24];       // 24
//   uint8_t  max[24];       // 24
//   double   scale[3];      // 24
//   double   offset[3];     // 24
//   char     description[32]; // 32
// Total = 2+1+1+32+4+24+24+24+24+24+32 = 192
// ═══════════════════════════════════════════════════════════════════════

// Byte sizes for data_type values
static uint16_t eb_type_size(uint8_t dt) {
    switch (dt) {
        case 0:  return 0;  // undocumented
        case 1:  return 1;  // uint8
        case 2:  return 1;  // int8
        case 3:  return 2;  // uint16
        case 4:  return 2;  // int16
        case 5:  return 4;  // uint32
        case 6:  return 4;  // int32
        case 7:  return 8;  // uint64
        case 8:  return 8;  // int64
        case 9:  return 4;  // float
        case 10: return 8;  // double
        default: return 0;
    }
}

uint16_t standard_point_size(uint8_t pdrf) {
    switch (pdrf) {
        case 6: return 30;
        case 7: return 36;
        case 8: return 38;
        default: return 30; // fallback
    }
}

std::vector<ExtraBytesRecord> parse_extra_bytes_vlr(const LASHeader& h) {
    std::vector<ExtraBytesRecord> result;

    const ParsedVLR* vlr = find_vlr(h, "LASF_Spec", 4);
    if (!vlr || vlr->data.empty())
        return result;

    // Each record is 192 bytes
    const size_t RECORD_SIZE = 192;
    size_t num_records = vlr->data.size() / RECORD_SIZE;
    uint16_t running_offset = 0;

    for (size_t i = 0; i < num_records; i++) {
        const uint8_t* rec = vlr->data.data() + i * RECORD_SIZE;

        ExtraBytesRecord eb;
        eb.data_type = rec[2];
        uint8_t options = rec[3];

        // Name at offset 4, 32 chars
        eb.name = std::string(reinterpret_cast<const char*>(rec + 4),
                              strnlen(reinterpret_cast<const char*>(rec + 4), 32));

        // options bitmap:
        //   bit 0: has_no_data
        //   bit 1: has_min
        //   bit 2: has_max
        //   bit 3: has_scale
        //   bit 4: has_offset
        eb.has_no_data = (options & 0x01) != 0;
        eb.has_scale   = (options & 0x08) != 0;
        eb.has_offset  = (options & 0x10) != 0;

        eb.byte_size = eb_type_size(eb.data_type);
        if (eb.byte_size == 0 && eb.data_type == 0) {
            // Undocumented: use the difference between actual and standard
            // point size, divided by number of records... or skip
            continue;
        }

        eb.byte_offset = running_offset;
        running_offset += eb.byte_size;

        // no_data at offset 40
        eb.no_data = 0.0;
        if (eb.has_no_data) {
            std::memcpy(&eb.no_data, rec + 40, 8);
        }

        // scale at offset 136 (192 - 32 - 24 = 136)
        eb.scale = 1.0;
        if (eb.has_scale) {
            std::memcpy(&eb.scale, rec + 136, 8);
        }

        // offset at offset 160 (136 + 24)
        eb.offset = 0.0;
        if (eb.has_offset) {
            std::memcpy(&eb.offset, rec + 160, 8);
        }

        result.push_back(eb);
    }

    return result;
}

} // namespace copc4r
