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

#pragma once

// System includes
#include <vector>
#include <cstdint>
#include <algorithm>
#include <limits>
#include <cmath>

// Project includes
#include "includes/define.h"
#include "geometries/geometry_data.h"

namespace Kratos
{

class GeometricalObject; // forward declaration to reduce compile time

///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

/**
 * @class StlBvhTree
 * @ingroup KratosCore
 * @brief A BVH (Bounding Volume Hierarchy) tree for fast radius searches over triangle meshes.
 * @details Constructs from an iterator range. Only Triangle3D3 geometries are stored.
 *          Triangles are grouped in packs of 4 (SOA layout) for vectorizable leaf tests.
 *          Nodes are quantized to fit in one 64-byte cache line.
 *          Morton codes are used to sort primitives before building for spatial coherence.
 *          The tree is immutable after construction.
 */
class KRATOS_API(KRATOS_CORE) StlBvhTree
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(StlBvhTree);

    using PointType = array_1d<double, 3>;

    /// Result returned by SearchNearest
    struct NearestResult {
        GeometricalObject* p_object = nullptr; ///< Nearest object; nullptr if tree is empty
        double Distance = std::numeric_limits<double>::max(); ///< Distance to the nearest triangle surface
    };

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Constructor from an iterator range of objects with GetGeometry().
     * @details Iterates the range, stores Triangle3D3 geometries, and builds the BVH.
     * @param rBegin Begin iterator
     * @param rEnd   End iterator
     * @tparam TIteratorType Iterator whose dereferenced type exposes GetGeometry()
     */
    template<typename TIteratorType>
    StlBvhTree(TIteratorType rBegin, TIteratorType rEnd)
    {
        std::vector<TriangleData> raw_triangles;
        for (auto it = rBegin; it != rEnd; ++it) {
            const auto& r_geom = it->GetGeometry();
            if (r_geom.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Triangle3D3) {
                raw_triangles.push_back(ExtractTriangle(r_geom, &(*it)));
            }
        }
        if (!raw_triangles.empty()) {
            Build(raw_triangles);
        }
    }

    ~StlBvhTree() = default;

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Applies a function to every triangle within radius of a center point.
     * @details Traversal stops early if rFunction returns false.
     * @param rCenter   The query point
     * @param Radius    The search radius
     * @param rFunction Callable(GeometricalObject&) -> bool; return false to stop
     * @tparam TFunction Callable type
     */
    template<typename TFunction>
    void ApplyInRadius(
        const PointType& rCenter,
        double Radius,
        TFunction&& rFunction) const
    {
        if (mNodes.empty()) return;

        const float cx = static_cast<float>(rCenter[0]);
        const float cy = static_cast<float>(rCenter[1]);
        const float cz = static_cast<float>(rCenter[2]);
        const float r2 = static_cast<float>(Radius * Radius);

        // Stack-based traversal avoids heap allocation per query
        int stack[256];
        int top = 0;
        stack[top++] = 0; // root node index

        while (top > 0) {
            const BvhNode& node = mNodes[stack[--top]];

            for (int child_index = 0; child_index < 2; ++child_index) {
                const int32_t slot = node.child[child_index];
                if (slot == kEmptySlot) continue;

                // Dequantize child AABB
                const float lx = node.base[0] + node.lo[child_index][0] * node.extent[0];
                const float ly = node.base[1] + node.lo[child_index][1] * node.extent[1];
                const float lz = node.base[2] + node.lo[child_index][2] * node.extent[2];
                const float ux = node.base[0] + node.hi[child_index][0] * node.extent[0];
                const float uy = node.base[1] + node.hi[child_index][1] * node.extent[1];
                const float uz = node.base[2] + node.hi[child_index][2] * node.extent[2];

                if (AABBSphereDistance2(cx, cy, cz, lx, ly, lz, ux, uy, uz) > r2) continue;

                if (slot < 0) {
                    // Leaf: test exact point-triangle distance for each triangle in pack
                    const TrianglePack4& pack = mPacks[PackIndex(slot)];
                    for (int j = 0; j < pack.count; ++j) {
                        const float v0[3] = {pack.v0x[j], pack.v0y[j], pack.v0z[j]};
                        const float e1[3] = {pack.e1x[j], pack.e1y[j], pack.e1z[j]};
                        const float e2[3] = {pack.e2x[j], pack.e2y[j], pack.e2z[j]};
                        if (PointTriangleDistance2(cx, cy, cz, v0, e1, e2) <= r2) {
                            if (!rFunction(*pack.objects[j])) return;
                        }
                    }
                } else {
                    // Inner node: push onto stack
                    stack[top++] = slot;
                }
            }
        }
    }

    /**
     * @brief Finds the nearest triangle to a query point.
     * @details Two-step algorithm sharing a single traversal stack:
     *   1. Greedy descent: at each inner node follow the nearest child, pushing the
     *      farther sibling onto the stack. On reaching a leaf, evaluate its pack and
     *      record best_dist2.  The stack now holds all unvisited siblings.
     *   2. Verification: continue consuming the same stack.  Any node whose AABB
     *      min-distance2 >= best_dist2 is pruned; leaves are tested directly using
     *      pack data (no geometry indirection).  best_dist2 shrinks as closer
     *      triangles are found, pruning progressively more of the remaining nodes.
     * @param rCenter The query point.
     * @return NearestResult with the nearest GeometricalObject pointer and its distance.
     *         If the tree is empty, p_object is nullptr and Distance is max double.
     */
    NearestResult SearchNearest(const PointType& rCenter) const;

    ///@}
    ///@name Input and output
    ///@{

    std::string Info() const { return "StlBvhTree"; }

    ///@}

private:
    ///@name Private Types
    ///@{

    /// Intermediate data per triangle, used during construction only
    struct TriangleData {
        float v0[3];
        float e1[3]; /// edge1 = v1 - v0
        float e2[3]; /// edge2 = v2 - v0
        float centroid[3];
        GeometricalObject* p_object;
    };

    /// Four triangles packed in SOA layout for vectorized leaf queries
    struct alignas(64) TrianglePack4 {
        float v0x[4], v0y[4], v0z[4];
        float e1x[4], e1y[4], e1z[4];
        float e2x[4], e2y[4], e2z[4];
        GeometricalObject* objects[4];
        int count; /// number of valid triangles (1-4)
    };

    /// Quantized BVH2 node: stores AABB bounds for 2 children, fits in one 64-byte cache line
    struct alignas(64) BvhNode {
        float   base[3];    /// 12 bytes – origin for child AABB dequantization
        float   extent[3];    /// 12 bytes – (aabb_extent / 255) per axis
        uint8_t lo[2][3];   ///  6 bytes – quantized lower bounds for children 0 and 1
        uint8_t hi[2][3];   ///  6 bytes – quantized upper bounds for children 0 and 1
        int32_t child[2];   ///  8 bytes – >=0: inner node index; <0: leaf pack; kEmptySlot: unused
        uint8_t pad[20];    /// 20 bytes – padding to reach 64 bytes
        //                                 12+12+6+6+8+20 = 64 bytes
    };

    static constexpr int32_t kEmptySlot = 0x7fffffff;

    /// Encode a pack index as a negative child slot value
    static int32_t PackRef(int pack_idx) { return -(pack_idx + 1); }

    /// Decode a pack index from a negative child slot value
    static int PackIndex(int32_t ref) { return -(ref + 1); }

    ///@}
    ///@name Private Members
    ///@{

    std::vector<BvhNode>  mNodes;
    std::vector<TrianglePack4> mPacks;

    ///@}
    ///@name Private Helpers (inline)
    ///@{

    static TriangleData ExtractTriangle(
        const Geometry<Node>& rGeom,
        GeometricalObject* pObj)
    {
        const auto& p0 = rGeom[0].Coordinates();
        const auto& p1 = rGeom[1].Coordinates();
        const auto& p2 = rGeom[2].Coordinates();
        TriangleData t;
        t.v0[0] = static_cast<float>(p0[0]);
        t.v0[1] = static_cast<float>(p0[1]);
        t.v0[2] = static_cast<float>(p0[2]);
        t.e1[0] = static_cast<float>(p1[0] - p0[0]);
        t.e1[1] = static_cast<float>(p1[1] - p0[1]);
        t.e1[2] = static_cast<float>(p1[2] - p0[2]);
        t.e2[0] = static_cast<float>(p2[0] - p0[0]);
        t.e2[1] = static_cast<float>(p2[1] - p0[1]);
        t.e2[2] = static_cast<float>(p2[2] - p0[2]);
        t.centroid[0] = t.v0[0] + (t.e1[0] + t.e2[0]) / 3.0f;
        t.centroid[1] = t.v0[1] + (t.e1[1] + t.e2[1]) / 3.0f;
        t.centroid[2] = t.v0[2] + (t.e1[2] + t.e2[2]) / 3.0f;
        t.p_object = pObj;
        return t;
    }

    static uint32_t EncodeMorton3(uint32_t x, uint32_t y, uint32_t z)
    {
        auto spread = [](uint32_t v) -> uint32_t {
            v &= 0x000003ffu;
            v = (v | (v << 16u)) & 0x030000ffu;
            v = (v | (v <<  8u)) & 0x0300f00fu;
            v = (v | (v <<  4u)) & 0x030c30c3u;
            v = (v | (v <<  2u)) & 0x09249249u;
            return v;
        };
        return spread(x) | (spread(y) << 1u) | (spread(z) << 2u);
    }

    /// Squared distance from point (px,py,pz) to nearest point on triangle V0, V0+e1, V0+e2.
    /// Implements the Ericson closest-point-on-triangle method.
    static float PointTriangleDistance2(
        float px, float py, float pz,
        const float v0[3],
        const float e1[3],
        const float e2[3])
    {
        // d = v0 - P
        const float dx = v0[0] - px, dy = v0[1] - py, dz = v0[2] - pz;
        const float a  = e1[0]*e1[0] + e1[1]*e1[1] + e1[2]*e1[2]; // dot(e1,e1)
        const float b  = e1[0]*e2[0] + e1[1]*e2[1] + e1[2]*e2[2]; // dot(e1,e2)
        const float c  = e2[0]*e2[0] + e2[1]*e2[1] + e2[2]*e2[2]; // dot(e2,e2)
        const float dd = e1[0]*dx    + e1[1]*dy    + e1[2]*dz;     // dot(e1,d)
        const float e  = e2[0]*dx    + e2[1]*dy    + e2[2]*dz;     // dot(e2,d)
        const float f  = dx*dx       + dy*dy        + dz*dz;       // dot(d,d)

        const float det = a*c - b*b;
        float s = b*e - c*dd;
        float t = b*dd - a*e;

        auto dist2_at = [&](float s_, float t_) {
            const float qx = v0[0] + s_*e1[0] + t_*e2[0] - px;
            const float qy = v0[1] + s_*e1[1] + t_*e2[1] - py;
            const float qz = v0[2] + s_*e1[2] + t_*e2[2] - pz;
            return qx*qx + qy*qy + qz*qz;
        };

        if (s + t <= det) {
            if (s < 0.0f) {
                if (t < 0.0f) {
                    // Region 4: closest to v0
                    return f;
                } else {
                    // Region 3: clamp to edge e2
                    s = 0.0f;
                    t = (e >= 0.0f) ? 0.0f : ((-e >= c) ? 1.0f : -e/c);
                }
            } else if (t < 0.0f) {
                // Region 5: clamp to edge e1
                t = 0.0f;
                s = (dd >= 0.0f) ? 0.0f : ((-dd >= a) ? 1.0f : -dd/a);
            } else {
                // Region 0: interior, use unconstrained minimizer
                const float inv_det = (det > 1e-10f) ? (1.0f / det) : 0.0f;
                s *= inv_det;
                t *= inv_det;
            }
        } else {
            if (s < 0.0f) {
                // Region 2
                const float tmp0 = b + dd, tmp1 = c + e;
                if (tmp1 > tmp0) {
                    const float numer = tmp1 - tmp0;
                    const float denom = a - 2.0f*b + c;
                    s = (numer >= denom) ? 1.0f : numer/denom;
                    t = 1.0f - s;
                } else {
                    s = 0.0f;
                    t = (tmp1 <= 0.0f) ? 1.0f : ((e >= 0.0f) ? 0.0f : -e/c);
                }
            } else if (t < 0.0f) {
                // Region 6
                const float tmp0 = b + e, tmp1 = a + dd;
                if (tmp1 > tmp0) {
                    const float numer = tmp1 - tmp0;
                    const float denom = a - 2.0f*b + c;
                    t = (numer >= denom) ? 1.0f : numer/denom;
                    s = 1.0f - t;
                } else {
                    t = 0.0f;
                    s = (tmp1 <= 0.0f) ? 1.0f : ((dd >= 0.0f) ? 0.0f : -dd/a);
                }
            } else {
                // Region 1: clamp to hypotenuse s+t=1
                const float numer = c + e - b - dd;
                if (numer <= 0.0f) {
                    s = 0.0f;
                } else {
                    const float denom = a - 2.0f*b + c;
                    s = (numer >= denom) ? 1.0f : numer/denom;
                }
                t = 1.0f - s;
            }
        }
        return dist2_at(s, t);
    }

    /// Squared distance from point to an AABB (0 if point is inside)
    static float AABBSphereDistance2(
        float cx, float cy, float cz,
        float lx, float ly, float lz,
        float ux, float uy, float uz)
    {
        const float dx = std::max(0.0f, std::max(lx - cx, cx - ux));
        const float dy = std::max(0.0f, std::max(ly - cy, cy - uy));
        const float dz = std::max(0.0f, std::max(lz - cz, cz - uz));
        return dx*dx + dy*dy + dz*dz;
    }

    ///@}
    ///@name Private Methods (defined in .cpp)
    ///@{

    void Build(std::vector<TriangleData>& rRawTris);

    struct BuildResult { int ref; float lo[3], hi[3]; };
    BuildResult BuildNodeInternal(int pack_begin, int pack_end);

    ///@}
};

///@}
///@}

} // namespace Kratos
