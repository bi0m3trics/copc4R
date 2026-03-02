// copc_hierarchy.h  –  COPC Info VLR + hierarchy pages + node selection
#ifndef COPC4R_COPC_HIERARCHY_H
#define COPC4R_COPC_HIERARCHY_H

#include <cstdint>
#include <vector>
#include <unordered_map>
#include <string>
#include <functional>
#include "las14_header.h"

namespace copc4r {

// ═══════════════════════════════════════════════════════════════════════
// COPC Info VLR  (record_id = 1, user_id = "copc")
// Spec: https://copc.io/  –  160 bytes total
// ═══════════════════════════════════════════════════════════════════════
#pragma pack(push, 1)
struct COPCInfo {
    // Center of the octree (X, Y, Z) – doubles
    double center_x;       //  0
    double center_y;       //  8
    double center_z;       // 16
    // Halfsize (== spacing * 2^(root_level) / 2 ... but actually
    // the spec stores "halfsize" as the half of the root cube side).
    double halfsize;       // 24
    // Spacing at the deepest level? No – spacing at level 0 (the root
    // node spacing – each lower level halves it).
    double spacing;        // 32
    // Root hierarchy page –  offset, size in file
    uint64_t root_hier_offset; // 40
    uint64_t root_hier_size;   // 48
    // GPS time range
    double gpstime_minimum;    // 56
    double gpstime_maximum;    // 64
    // Reserved / padding to 160
    uint8_t reserved[88];     // 72  (160 - 72 = 88)
};
#pragma pack(pop)
static_assert(sizeof(COPCInfo) == 160, "COPCInfo must be 160 bytes");

// ═══════════════════════════════════════════════════════════════════════
// VoxelKey  –  (level, x, y, z) identifying a node in the octree
// ═══════════════════════════════════════════════════════════════════════
struct VoxelKey {
    int32_t level;
    int32_t x;
    int32_t y;
    int32_t z;

    bool operator==(const VoxelKey& o) const {
        return level == o.level && x == o.x && y == o.y && z == o.z;
    }
};

struct VoxelKeyHash {
    size_t operator()(const VoxelKey& k) const {
        // Simple hash combining
        size_t h = std::hash<int32_t>()(k.level);
        h ^= std::hash<int32_t>()(k.x) + 0x9e3779b9 + (h << 6) + (h >> 2);
        h ^= std::hash<int32_t>()(k.y) + 0x9e3779b9 + (h << 6) + (h >> 2);
        h ^= std::hash<int32_t>()(k.z) + 0x9e3779b9 + (h << 6) + (h >> 2);
        return h;
    }
};

// ═══════════════════════════════════════════════════════════════════════
// Hierarchy entry  –  32 bytes each in a hierarchy page
//   key (16 bytes: 4×int32) + offset(8) + byteSize(4) + pointCount(4)
// But actually in COPC spec the entry is:
//   VoxelKey (16 bytes) | int64 offset | int32 byteSize | int32 pointCount
// Total = 32 bytes per entry.
//
// Interpretation:
//   - If pointCount == -1: this entry points to another hierarchy page
//     at (offset, byteSize) in the file.
//   - If pointCount == 0:  empty node (no data).
//   - If pointCount >  0:  data chunk at (offset, byteSize) with
//     pointCount compressed points.
// ═══════════════════════════════════════════════════════════════════════
#pragma pack(push, 1)
struct HierarchyEntry {
    VoxelKey key;              // 16
    uint64_t offset;           //  8
    int32_t  byte_size;        //  4
    int32_t  point_count;      //  4  (-1 = sub-page, 0 = empty, >0 = chunk)
};
#pragma pack(pop)
static_assert(sizeof(HierarchyEntry) == 32, "HierarchyEntry must be 32 bytes");

// ═══════════════════════════════════════════════════════════════════════
// Node data reference  –  what we collect for chunks to read
// ═══════════════════════════════════════════════════════════════════════
struct NodeChunk {
    VoxelKey key;
    uint64_t offset;
    uint64_t byte_size;
    uint64_t point_count;
};

// ═══════════════════════════════════════════════════════════════════════
// AABB for intersection tests
// ═══════════════════════════════════════════════════════════════════════
struct AABB {
    double xmin, ymin, zmin;
    double xmax, ymax, zmax;

    bool intersects_xy(double qxmin, double qymin,
                       double qxmax, double qymax) const {
        return xmax >= qxmin && xmin <= qxmax &&
               ymax >= qymin && ymin <= qymax;
    }
    bool intersects_z(double qzmin, double qzmax) const {
        return zmax >= qzmin && zmin <= qzmax;
    }
};

// Forward
class RangeReader;

// ═══════════════════════════════════════════════════════════════════════
// Public API
// ═══════════════════════════════════════════════════════════════════════

/// Parse the COPC Info VLR from an already-read LAS header.
COPCInfo parse_copc_info(const LASHeader& header);

/// Derive the axis-aligned bounding box of an octree node from COPC
/// center/halfsize and the node's VoxelKey.
AABB node_bounds(const COPCInfo& info, const VoxelKey& key);

/// Recursively load hierarchy pages and collect NodeChunks whose AABB
/// intersects the query bbox (and optional zrange).
/// Returns ALL matching nodes sorted by level (coarsest first).
/// If bbox is empty (all zeros), all nodes are selected.
std::vector<NodeChunk> select_nodes(
    RangeReader& reader,
    const COPCInfo& info,
    double bbox_xmin, double bbox_ymin,
    double bbox_xmax, double bbox_ymax,
    double zmin, double zmax,
    bool   has_bbox,
    bool   has_zrange);

} // namespace copc4r
#endif
