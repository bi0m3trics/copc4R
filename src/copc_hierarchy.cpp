// copc_hierarchy.cpp  –  COPC Info VLR + hierarchy traversal + node selection
#include "copc_hierarchy.h"
#include "range_reader.h"
#include <cstring>
#include <stdexcept>
#include <queue>
#include <algorithm>
#include <cmath>
#include <limits>

namespace copc4r {

// ═══════════════════════════════════════════════════════════════════════
// parse_copc_info
// ═══════════════════════════════════════════════════════════════════════
COPCInfo parse_copc_info(const LASHeader& header) {
    // The COPC Info VLR uses user_id = "copc" and record_id = 1.
    const ParsedVLR* vlr = find_vlr(header, "copc", 1);
    if (!vlr)
        throw std::runtime_error(
            "This file is not a valid COPC file. "
            "The mandatory COPC Info VLR (user_id='copc', record_id=1) "
            "was not found. copc4R only reads Cloud Optimized Point Cloud "
            "(.copc.laz) files. Use a tool such as 'pdal' or 'untwine' to "
            "convert standard LAS/LAZ files to COPC format.");
    if (vlr->data.size() < sizeof(COPCInfo))
        throw std::runtime_error("COPC Info VLR payload too small");

    COPCInfo info;
    std::memcpy(&info, vlr->data.data(), sizeof(COPCInfo));
    return info;
}

// ═══════════════════════════════════════════════════════════════════════
// node_bounds  –  derive AABB from octree key + COPC center/halfsize
//
// Octree traversal: root at level=0, each level subdivides by 2 per axis.
// Node (level,x,y,z) has side = 2*halfsize / 2^level
// ═══════════════════════════════════════════════════════════════════════
AABB node_bounds(const COPCInfo& info, const VoxelKey& key) {
    // Validate inputs to prevent NaN/Inf propagation
    if (!std::isfinite(info.halfsize) || !std::isfinite(info.center_x) ||
        !std::isfinite(info.center_y) || !std::isfinite(info.center_z)) {
        throw std::runtime_error("COPC Info contains non-finite values");
    }
    if (info.halfsize <= 0.0) {
        throw std::runtime_error("COPC Info halfsize must be positive");
    }
    if (key.level < 0 || key.level > 30) {
        throw std::runtime_error("COPC octree level out of valid range [0,30]: " +
                                 std::to_string(key.level));
    }
    
    // Use ldexp for safe power-of-2 calculation (avoids bit-shift overflow)
    double denominator = std::ldexp(1.0, key.level);  // 2^level
    double side = 2.0 * info.halfsize / denominator;
    AABB b;
    b.xmin = info.center_x - info.halfsize + key.x * side;
    b.ymin = info.center_y - info.halfsize + key.y * side;
    b.zmin = info.center_z - info.halfsize + key.z * side;
    b.xmax = b.xmin + side;
    b.ymax = b.ymin + side;
    b.zmax = b.zmin + side;
    return b;
}

// ═══════════════════════════════════════════════════════════════════════
// parse_hierarchy_page  –  read one page of hierarchy entries
// ═══════════════════════════════════════════════════════════════════════
static std::vector<HierarchyEntry> parse_hierarchy_page(
    RangeReader& reader, uint64_t offset, uint64_t size)
{
    auto data = reader.read(offset, size);
    size_t n = data.size() / sizeof(HierarchyEntry);
    std::vector<HierarchyEntry> entries(n);
    std::memcpy(entries.data(), data.data(), n * sizeof(HierarchyEntry));
    return entries;
}

// ═══════════════════════════════════════════════════════════════════════
// select_nodes  –  BFS through hierarchy pages, prune by bbox / zrange
// ═══════════════════════════════════════════════════════════════════════
std::vector<NodeChunk> select_nodes(
    RangeReader& reader,
    const COPCInfo& info,
    double bbox_xmin, double bbox_ymin,
    double bbox_xmax, double bbox_ymax,
    double zmin, double zmax,
    bool   has_bbox,
    bool   has_zrange,
    int    max_depth)
{
    std::vector<NodeChunk> result;

    // BFS queue of hierarchy pages to read: (offset, size)
    struct PageRef { uint64_t offset; uint64_t size; };
    std::queue<PageRef> page_queue;
    page_queue.push({info.root_hier_offset, info.root_hier_size});

    while (!page_queue.empty()) {
        auto [pg_off, pg_sz] = page_queue.front();
        page_queue.pop();

        auto entries = parse_hierarchy_page(reader, pg_off, pg_sz);

        for (auto& e : entries) {
            // Skip empty nodes
            if (e.point_count == 0)
                continue;

            // LOD depth limit: skip nodes deeper than max_depth
            if (max_depth >= 0 && e.key.level > max_depth)
                continue;

            // Spatial filtering: derive this node's bounding box
            AABB bounds = node_bounds(info, e.key);
            if (has_bbox &&
                !bounds.intersects_xy(bbox_xmin, bbox_ymin,
                                      bbox_xmax, bbox_ymax))
                continue;
            if (has_zrange && !bounds.intersects_z(zmin, zmax))
                continue;

            if (e.point_count == -1) {
                // This entry points to a child hierarchy page
                page_queue.push({e.offset,
                                 static_cast<uint64_t>(e.byte_size)});
            } else {
                // Data chunk
                NodeChunk nc;
                nc.key        = e.key;
                nc.offset     = e.offset;
                nc.byte_size  = static_cast<uint64_t>(e.byte_size);
                nc.point_count= static_cast<uint64_t>(e.point_count);
                result.push_back(nc);
            }
        }
    }

    // Sort by level (coarsest first) for consistent ordering
    std::sort(result.begin(), result.end(),
              [](const NodeChunk& a, const NodeChunk& b) {
                  return a.key.level < b.key.level;
              });

    return result;
}

// ═══════════════════════════════════════════════════════════════════════
// count_nodes  –  count points without decompressing
// ═══════════════════════════════════════════════════════════════════════
NodeCount count_nodes(
    RangeReader& reader,
    const COPCInfo& info,
    double bbox_xmin, double bbox_ymin,
    double bbox_xmax, double bbox_ymax,
    bool   has_bbox)
{
    NodeCount nc{0, 0};
    struct PageRef { uint64_t offset; uint64_t size; };
    std::queue<PageRef> page_queue;
    page_queue.push({info.root_hier_offset, info.root_hier_size});

    while (!page_queue.empty()) {
        auto [pg_off, pg_sz] = page_queue.front();
        page_queue.pop();

        auto entries = parse_hierarchy_page(reader, pg_off, pg_sz);
        for (auto& e : entries) {
            if (e.point_count == 0) continue;

            AABB bounds = node_bounds(info, e.key);
            if (has_bbox &&
                !bounds.intersects_xy(bbox_xmin, bbox_ymin,
                                      bbox_xmax, bbox_ymax))
                continue;

            if (e.point_count == -1) {
                page_queue.push({e.offset,
                                 static_cast<uint64_t>(e.byte_size)});
            } else {
                nc.total_points += static_cast<uint64_t>(e.point_count);
                nc.num_nodes++;
            }
        }
    }
    return nc;
}

} // namespace copc4r
