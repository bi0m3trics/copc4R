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

} // namespace copc4r
