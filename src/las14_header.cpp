// las14_header.cpp  –  LAS 1.4 header + VLR parsing
#include "las14_header.h"
#include "range_reader.h"
#include <cstring>
#include <stdexcept>
#include <cstdio>

namespace {
// Configuration constants
constexpr uint64_t MAX_EVLR_PAYLOAD = 64u * 1024u * 1024u;  // 64 MB EVLR size limit
} // anonymous namespace

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
    
    // Validate header size
    if (h.pub.header_size < 375) {
        throw std::runtime_error("Invalid header_size " + 
                                 std::to_string(h.pub.header_size) + 
                                 " (must be >= 375 for LAS 1.4)");
    }
    if (h.pub.header_size > 10000) {
        throw std::runtime_error("Suspicious header_size " + 
                                 std::to_string(h.pub.header_size) + 
                                 " (possibly corrupted file)");
    }

    // ── 2. Read VLRs (located right after public header) ──────────────
    uint64_t cursor = h.pub.header_size;
    for (uint32_t i = 0; i < h.pub.num_vlrs; ++i) {
        auto vhdr_raw = reader.read(cursor, sizeof(VLRHeader));
        if (vhdr_raw.size() < sizeof(VLRHeader))
            throw std::runtime_error("Truncated VLR header at offset " +
                                     std::to_string(cursor));
        VLRHeader vhdr;
        std::memcpy(&vhdr, vhdr_raw.data(), sizeof(VLRHeader));
        cursor += sizeof(VLRHeader);

        if (vhdr.reserved != 0x0000) {
            std::fprintf(stderr, "Warning: VLR reserved field is non-zero (0x%04X) at offset %lld\n",
                         (unsigned)vhdr.reserved, (long long)(cursor - sizeof(VLRHeader)));
        }

        ParsedVLR vlr;
        vlr.reserved    = vhdr.reserved;
        vlr.user_id     = std::string(vhdr.user_id,
                                      strnlen(vhdr.user_id, sizeof(vhdr.user_id)));
        vlr.record_id   = vhdr.record_id;
        vlr.description = std::string(vhdr.description,
                                      strnlen(vhdr.description, sizeof(vhdr.description)));
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

            if (ehdr.reserved != 0x0000) {
                std::fprintf(stderr, "Warning: EVLR reserved field is non-zero (0x%04X) at offset %lld\n",
                             (unsigned)ehdr.reserved, (long long)(ecursor - sizeof(EVLRHeader)));
            }

            ParsedVLR vlr;
            vlr.reserved   = ehdr.reserved;
            vlr.user_id    = std::string(ehdr.user_id,
                                         strnlen(ehdr.user_id, sizeof(ehdr.user_id)));
            vlr.record_id  = ehdr.record_id;
            vlr.description = std::string(ehdr.description,
                                          strnlen(ehdr.description, sizeof(ehdr.description)));
            // Only read EVLR payload up to 64 MB to avoid OOM on weird files
            uint64_t safe_len = std::min(ehdr.record_length, MAX_EVLR_PAYLOAD);
            vlr.data = reader.read(ecursor, safe_len);
            // Use safe_len for cursor advancement to prevent overflow from corrupted headers
            ecursor += safe_len;

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
//   bytes 64-87:  min[3] (8 bytes each)
//   bytes 88-111: max[3] (8 bytes each)
//   bytes 112-135: scale[3] (8 bytes each)
//   bytes 136-159: offset[3] (8 bytes each)
//   bytes 160-191: description[32]
//
// Per ASPRS LAS 1.4 spec (and confirmed in LAStools lasattributer.hpp),
// the Extra Bytes VLR structure is exactly 192 bytes with NO deprecated fields:
//   uint16_t reserved;        // 2 bytes
//   uint8_t  data_type;       // 1 byte
//   uint8_t  options;         // 1 byte
//   char     name[32];        // 32 bytes
//   uint8_t  unused[4];       // 4 bytes
//   uint8_t  no_data[24];     // 24 bytes (anytype[3])
//   uint8_t  min[24];         // 24 bytes (anytype[3])
//   uint8_t  max[24];         // 24 bytes (anytype[3])
//   double   scale[3];        // 24 bytes
//   double   offset[3];       // 24 bytes
//   char     description[32]; // 32 bytes
// Total = 2+1+1+32+4+24+24+24+24+24+32 = 192 bytes
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
        eb.has_min     = (options & 0x02) != 0;
        eb.has_max     = (options & 0x04) != 0;
        eb.has_scale   = (options & 0x08) != 0;
        eb.has_offset  = (options & 0x10) != 0;

        eb.byte_size = eb_type_size(eb.data_type);
        if (eb.byte_size == 0 && eb.data_type == 0) {
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
