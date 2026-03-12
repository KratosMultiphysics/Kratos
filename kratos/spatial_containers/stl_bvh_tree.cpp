//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Altair Engineering
//

// System includes
#include <cstring>    // memset
#include <cfloat>     // FLT_MAX

// Project includes
#include "includes/geometrical_object.h"
#include "spatial_containers/stl_bvh_tree.h"

namespace Kratos
{


///@name StlBvhTree Private Methods
///@{

StlBvhTree::BuildResult StlBvhTree::BuildNodeInternal(int pack_begin, int pack_end)
{
    // Base case: single pack becomes a leaf child reference (no node created)
    if (pack_end - pack_begin == 1) {
        BuildResult result;
        result.ref = PackRef(pack_begin);
        // Compute pack AABB inline
        for (int k = 0; k < 3; ++k) { result.lo[k] =  FLT_MAX; result.hi[k] = -FLT_MAX; }
        const TrianglePack4& p = mPacks[pack_begin];
        for (int j = 0; j < p.count; ++j) {
            const float vx[3] = {p.v0x[j], p.v0x[j]+p.e1x[j], p.v0x[j]+p.e2x[j]};
            const float vy[3] = {p.v0y[j], p.v0y[j]+p.e1y[j], p.v0y[j]+p.e2y[j]};
            const float vz[3] = {p.v0z[j], p.v0z[j]+p.e1z[j], p.v0z[j]+p.e2z[j]};
            for (int i = 0; i < 3; ++i) {
                result.lo[0] = std::min(result.lo[0], vx[i]);
                result.hi[0] = std::max(result.hi[0], vx[i]);
                result.lo[1] = std::min(result.lo[1], vy[i]);
                result.hi[1] = std::max(result.hi[1], vy[i]);
                result.lo[2] = std::min(result.lo[2], vz[i]);
                result.hi[2] = std::max(result.hi[2], vz[i]);
            }
        }
        return result;
    }

    // Allocate this node before recursing so that root is always at index 0
    const int node_idx = static_cast<int>(mNodes.size());
    mNodes.emplace_back(BvhNode{});
    // Reserve has been called in Build(), so emplace_back won't reallocate

    const int mid = (pack_begin + pack_end) / 2;
    const BuildResult left  = BuildNodeInternal(pack_begin, mid);
    const BuildResult right = BuildNodeInternal(mid, pack_end);

    // Parent AABB = union of children AABBs
    float lo[3], hi[3];
    for (int k = 0; k < 3; ++k) {
        lo[k] = std::min(left.lo[k], right.lo[k]);
        hi[k] = std::max(left.hi[k], right.hi[k]);
    }

    // Fill node (safe because we reserved capacity before any recursive calls)
    BvhNode& node = mNodes[node_idx];
    for (int k = 0; k < 3; ++k) {
        node.base[k] = lo[k];
        node.extent[k] = std::max((hi[k] - lo[k]) / 255.0f, 1e-10f);
    }

    // Quantize each child's AABB relative to the parent
    const float* child_lo[2] = {left.lo,  right.lo};
    const float* child_hi[2] = {left.hi, right.hi};
    for (int c = 0; c < 2; ++c) {
        for (int k = 0; k < 3; ++k) {
            node.lo[c][k] = static_cast<uint8_t>(
                std::clamp((child_lo[c][k] - lo[k]) / node.extent[k], 0.0f, 255.0f));
            node.hi[c][k] = static_cast<uint8_t>(
                std::clamp((child_hi[c][k] - lo[k]) / node.extent[k], 0.0f, 255.0f));
        }
    }
    node.child[0] = left.ref;
    node.child[1] = right.ref;

    BuildResult result;
    result.ref = node_idx;
    for (int k = 0; k < 3; ++k) { result.lo[k] = lo[k]; result.hi[k] = hi[k]; }
    return result;
}

void StlBvhTree::Build(std::vector<TriangleData>& rRawTris)
{
    // --- 1. Compute scene centroid AABB for Morton normalization ---
    float scene_lo[3] = { FLT_MAX,  FLT_MAX,  FLT_MAX};
    float scene_hi[3] = {-FLT_MAX, -FLT_MAX, -FLT_MAX};
    for (const auto& t : rRawTris) {
        for (int k = 0; k < 3; ++k) {
            scene_lo[k] = std::min(scene_lo[k], t.centroid[k]);
            scene_hi[k] = std::max(scene_hi[k], t.centroid[k]);
        }
    }

    const float extent[3] = {
        scene_hi[0] - scene_lo[0] + 1e-10f,
        scene_hi[1] - scene_lo[1] + 1e-10f,
        scene_hi[2] - scene_lo[2] + 1e-10f
    };

    // --- 2. Sort triangles by Morton code for spatial coherence ---
    auto morton_key = [&](const TriangleData& t) -> uint32_t {
        const uint32_t x = static_cast<uint32_t>(
            std::clamp((t.centroid[0] - scene_lo[0]) / extent[0], 0.0f, 1.0f) * 1023.0f);
        const uint32_t y = static_cast<uint32_t>(
            std::clamp((t.centroid[1] - scene_lo[1]) / extent[1], 0.0f, 1.0f) * 1023.0f);
        const uint32_t z = static_cast<uint32_t>(
            std::clamp((t.centroid[2] - scene_lo[2]) / extent[2], 0.0f, 1.0f) * 1023.0f);
        return EncodeMorton3(x, y, z);
    };

    std::sort(rRawTris.begin(), rRawTris.end(),
        [&](const TriangleData& a, const TriangleData& b) {
            return morton_key(a) < morton_key(b);
        });

    // --- 3. Pack triangles into groups of 4 (SOA layout) ---
    const int n_tris  = static_cast<int>(rRawTris.size());
    const int n_packs = (n_tris + 3) / 4;
    mPacks.resize(n_packs);

    for (int i = 0; i < n_packs; ++i) {
        TrianglePack4& pack     = mPacks[i];
        const int count    = std::min(4, n_tris - i * 4);
        pack.count         = count;

        for (int j = 0; j < 4; ++j) {
            // Duplicate last triangle into unused slots so SIMD lanes are safe
            const int ti = i * 4 + (j < count ? j : count - 1);
            const TriangleData& t = rRawTris[ti];
            pack.v0x[j] = t.v0[0]; pack.v0y[j] = t.v0[1]; pack.v0z[j] = t.v0[2];
            pack.e1x[j] = t.e1[0]; pack.e1y[j] = t.e1[1]; pack.e1z[j] = t.e1[2];
            pack.e2x[j] = t.e2[0]; pack.e2y[j] = t.e2[1]; pack.e2z[j] = t.e2[2];
            pack.objects[j] = (j < count) ? t.p_object : nullptr;
        }
    }

    // --- 4. Build BVH2 ---
    // A binary tree over n_packs leaves has at most 2*n_packs - 1 inner nodes.
    // +1 for the potential single-pack wrapper root.
    mNodes.reserve(2 * n_packs + 1);

    if (n_packs == 1) {
        // Special case: wrap the single pack in a minimal root node
        mNodes.emplace_back(BvhNode{});
        BvhNode& root = mNodes[0];
        float lo[3] = { FLT_MAX,  FLT_MAX,  FLT_MAX};
        float hi[3] = {-FLT_MAX, -FLT_MAX, -FLT_MAX};
        const TrianglePack4& p0 = mPacks[0];
        for (int j = 0; j < p0.count; ++j) {
            const float vx[3] = {p0.v0x[j], p0.v0x[j]+p0.e1x[j], p0.v0x[j]+p0.e2x[j]};
            const float vy[3] = {p0.v0y[j], p0.v0y[j]+p0.e1y[j], p0.v0y[j]+p0.e2y[j]};
            const float vz[3] = {p0.v0z[j], p0.v0z[j]+p0.e1z[j], p0.v0z[j]+p0.e2z[j]};
            for (int i = 0; i < 3; ++i) {
                lo[0]=std::min(lo[0],vx[i]); hi[0]=std::max(hi[0],vx[i]);
                lo[1]=std::min(lo[1],vy[i]); hi[1]=std::max(hi[1],vy[i]);
                lo[2]=std::min(lo[2],vz[i]); hi[2]=std::max(hi[2],vz[i]);
            }
        }
        for (int k = 0; k < 3; ++k) {
            root.base[k] = lo[k];
            root.extent[k] = std::max((hi[k] - lo[k]) / 255.0f, 1e-10f);
            root.lo[0][k] = 0;
            root.hi[0][k] = 255;
        }
        root.child[0] = PackRef(0);
        root.child[1] = kEmptySlot;
        return;
    }

    BuildNodeInternal(0, n_packs);
    // Root is always at mNodes[0] (first allocation in outer-most recursive call)
}

///@}

} // namespace Kratos
