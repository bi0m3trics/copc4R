// copc_hierarchy.cpp  –  COPC Info VLR + hierarchy traversal + node selection
#include "copc_hierarchy.h"
#include "range_reader.h"
#include <cstring>
#include <stdexcept>
#include <queue>
#include <algorithm>

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
// The octree root node (level=0, x=0, y=0, z=0) has a cube of
//   [center_x - halfsize, center_x + halfsize] × same for y,z.
// Each level subdivides by 2 in each axis.  Node (l,x,y,z) at level l
// has a cube side = 2*halfsize / 2^l.  Its min corner is:
//   comp_min = center_comp - halfsize + x * side   (for x component)
// ═══════════════════════════════════════════════════════════════════════
AABB node_bounds(const COPCInfo& info, const VoxelKey& key) {
    double side = 2.0 * info.halfsize / (1 << key.level); // could overflow for huge levels
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
    bool   has_zrange)
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

} // namespace copc4r
