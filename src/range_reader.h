// range_reader.h  –  I/O abstraction for local files and HTTP range-reads
#ifndef COPC4R_RANGE_READER_H
#define COPC4R_RANGE_READER_H

#include <cstdint>
#include <string>
#include <vector>
#include <memory>
#include <fstream>
#include <stdexcept>

namespace copc4r {

// ── abstract base ──────────────────────────────────────────────────────
class RangeReader {
public:
    virtual ~RangeReader() = default;
    /// Read `len` bytes starting at absolute file offset `off`.
    virtual std::vector<uint8_t> read(uint64_t off, uint64_t len) = 0;
};

// ── local file ─────────────────────────────────────────────────────────
class LocalFileRangeReader : public RangeReader {
public:
    explicit LocalFileRangeReader(const std::string& path);
    std::vector<uint8_t> read(uint64_t off, uint64_t len) override;
private:
    std::string path_;
    std::ifstream ifs_;
};

// ── HTTP range reader (compiled only when COPC4R_WITH_CURL is defined) ─
#ifdef COPC4R_WITH_CURL
class HttpRangeReader : public RangeReader {
public:
    explicit HttpRangeReader(const std::string& url,
                             int max_retries = 3);
    ~HttpRangeReader() override;
    std::vector<uint8_t> read(uint64_t off, uint64_t len) override;
private:
    std::string url_;
    int max_retries_;
    void* curl_; // CURL* (opaque to avoid header leak)
};
#endif

// ── factory ────────────────────────────────────────────────────────────
/// Automatically picks Local vs HTTP based on scheme prefix.
std::unique_ptr<RangeReader> make_reader(const std::string& path_or_url);

} // namespace copc4r
#endif
