// range_reader.cpp  –  local + optional HTTP range-read implementation
#include "range_reader.h"
#include <cstring>
#include <algorithm>
#include <sstream>

#ifdef COPC4R_WITH_CURL
#include <curl/curl.h>
#ifdef _WIN32
#include <windows.h>
#else
#include <unistd.h>
#endif
#endif

namespace copc4r {

// ═══════════════════════════════════════════════════════════════════════
// LocalFileRangeReader
// ═══════════════════════════════════════════════════════════════════════
LocalFileRangeReader::LocalFileRangeReader(const std::string& path)
    : path_(path)
{
    ifs_.open(path, std::ios::binary);
    if (!ifs_.is_open())
        throw std::runtime_error("Cannot open file: " + path);
}

std::vector<uint8_t> LocalFileRangeReader::read(uint64_t off, uint64_t len) {
    std::vector<uint8_t> buf(len);
    ifs_.seekg(static_cast<std::streamoff>(off), std::ios::beg);
    if (!ifs_)
        throw std::runtime_error("Seek failed at offset " + std::to_string(off));
    ifs_.read(reinterpret_cast<char*>(buf.data()), static_cast<std::streamsize>(len));
    auto got = static_cast<uint64_t>(ifs_.gcount());
    if (got < len) {
        buf.resize(got);
    }
    return buf;
}

// ═══════════════════════════════════════════════════════════════════════
// HttpRangeReader  (only compiled with COPC4R_WITH_CURL)
// ═══════════════════════════════════════════════════════════════════════
#ifdef COPC4R_WITH_CURL

// libcurl write-callback: appends data to std::vector<uint8_t>*
static size_t curl_write_cb(char* ptr, size_t size, size_t nmemb,
                            void* userdata) {
    auto* vec = static_cast<std::vector<uint8_t>*>(userdata);
    size_t total = size * nmemb;
    vec->insert(vec->end(),
                reinterpret_cast<uint8_t*>(ptr),
                reinterpret_cast<uint8_t*>(ptr) + total);
    return total;
}

HttpRangeReader::HttpRangeReader(const std::string& url, int max_retries)
    : url_(url), max_retries_(max_retries), curl_(nullptr)
{
    curl_ = curl_easy_init();
    if (!curl_)
        throw std::runtime_error("curl_easy_init() failed");

    // Keep-alive and common options
    auto* c = static_cast<CURL*>(curl_);
    curl_easy_setopt(c, CURLOPT_TCP_KEEPALIVE, 1L);
    curl_easy_setopt(c, CURLOPT_FOLLOWLOCATION, 1L);
    curl_easy_setopt(c, CURLOPT_USERAGENT, "copc4R/0.1");
    curl_easy_setopt(c, CURLOPT_WRITEFUNCTION, curl_write_cb);
}

HttpRangeReader::~HttpRangeReader() {
    if (curl_)
        curl_easy_cleanup(static_cast<CURL*>(curl_));
}

std::vector<uint8_t> HttpRangeReader::read(uint64_t off, uint64_t len) {
    auto* c = static_cast<CURL*>(curl_);
    std::vector<uint8_t> buf;
    buf.reserve(len);

    // Build Range header: "bytes=off-(off+len-1)"
    // Check for overflow in range calculation
    if (len > 0 && off > UINT64_MAX - (len - 1)) {
        throw std::runtime_error("HTTP range calculation overflow");
    }
    std::ostringstream range;
    range << off << "-" << (off + len - 1);
    std::string range_str = range.str();

    curl_easy_setopt(c, CURLOPT_URL, url_.c_str());
    curl_easy_setopt(c, CURLOPT_RANGE, range_str.c_str());
    curl_easy_setopt(c, CURLOPT_WRITEDATA, &buf);

    for (int attempt = 0; attempt < max_retries_; ++attempt) {
        buf.clear();
        CURLcode res = curl_easy_perform(c);
        if (res == CURLE_OK) {
            long http_code = 0;
            curl_easy_getinfo(c, CURLINFO_RESPONSE_CODE, &http_code);
            if (http_code == 206 || http_code == 200)
                return buf;
        }
        // Retry on transient errors
        if (attempt + 1 < max_retries_) {
            // Brief pause – not ideal in R, but acceptable for retries
#ifdef _WIN32
            Sleep(500);
#else
            usleep(500000);
#endif
        }
    }
    throw std::runtime_error("HTTP range read failed after " +
                             std::to_string(max_retries_) + " attempts for " + url_);
}

#endif // COPC4R_WITH_CURL

// ═══════════════════════════════════════════════════════════════════════
// Factory
// ═══════════════════════════════════════════════════════════════════════
std::unique_ptr<RangeReader> make_reader(const std::string& path_or_url) {
    // Detect HTTP(S) URLs
    if (path_or_url.size() > 8 &&
        (path_or_url.substr(0, 7) == "http://" ||
         path_or_url.substr(0, 8) == "https://")) {
#ifdef COPC4R_WITH_CURL
        return std::make_unique<HttpRangeReader>(path_or_url);
#else
        throw std::runtime_error(
            "HTTP support not compiled. Rebuild copc4R with "
            "COPC4R_WITH_CURL defined and libcurl available.");
#endif
    }
    return std::make_unique<LocalFileRangeReader>(path_or_url);
}

} // namespace copc4r
