//  ██████   ██████ ██████████  █████████  █████   █████ █████    ███████
// ░░██████ ██████ ░░███░░░░░█ ███░░░░░███░░███   ░░███ ░░███   ███░░░░░███      ███         ███
//  ░███░█████░███  ░███  █ ░ ░███    ░░░  ░███    ░███  ░███  ███     ░░███    ░███        ░███
//  ░███░░███ ░███  ░██████   ░░█████████  ░███████████  ░███ ░███      ░███ ███████████ ███████████
//  ░███ ░░░  ░███  ░███░░█    ░░░░░░░░███ ░███░░░░░███  ░███ ░███      ░███░░░░░███░░░ ░░░░░███░░░
//  ░███      ░███  ░███ ░   █ ███    ░███ ░███    ░███  ░███ ░░███     ███     ░███        ░███
//  █████     █████ ██████████░░█████████  █████   █████ █████ ░░░███████░      ░░░         ░░░
// ░░░░░     ░░░░░ ░░░░░░░░░░  ░░░░░░░░░  ░░░░░   ░░░░░ ░░░░░    ░░░░░░░
//
//
//  License:         MIT License
//                   meshio++ default license: LICENSE
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//
// Boundary-skin extraction, following the face-hashing algorithm of Kratos
// Multiphysics' SkinDetectionProcess (kratos/processes/skin_detection_process):
// a face whose sorted corner-node key occurs exactly once across all volume
// cells is a boundary face; faces seen twice are interior and cancel out.

// System includes
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstring>
#include <functional>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

// Project includes
#include "meshioplusplus/skin.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/cell_faces.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/log.hpp"

namespace meshioplusplus {

namespace {

// Sorted corner ids of one face, padded with -1 up to 4 entries.
using SkinFaceKey = std::array<std::int64_t, 4>;

struct SkinFaceKeyHash {
    std::size_t operator()(const SkinFaceKey& rKey) const {
        std::size_t h = 0;
        for (std::int64_t v : rKey)
            h ^= std::hash<std::int64_t>{}(v) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        return h;
    }
};

bool skin_block_supported(const CellType type, const bool is_ragged) {
    return cell_type_dimension(type) == 3 && !is_ragged && detail::skin_supported(type);
}

SkinFaceKey skin_face_key(const NDArray& rConn, std::size_t row_offset,
                          const detail::CellFaceDef& rFace) {
    SkinFaceKey key = {-1, -1, -1, -1};
    for (std::uint8_t k = 0; k < rFace.mNumCorners; ++k)
        key[k] = detail::read_int(rConn, row_offset + rFace.mNodes[k]);
    std::sort(key.begin(), key.end());
    return key;
}

// Invoke f(conn, row_offset, face_def) for every face of every cell of every
// supported volume block, in block-major / cell-major / face-slot order.
template <class F>
void skin_for_each_face(const Mesh& rMesh, F&& f) {
    for (const auto cb : rMesh.CellRange()) {
        const CellType ct = cell_type_from_name(cb.Type());
        if (!skin_block_supported(ct, cb.IsRagged()))
            continue;
        const std::vector<detail::CellFaceDef>& faces = detail::cell_faces(ct);
        const NDArray& conn = cb.Conn();
        const std::size_t npc = cb.NodesPerCell();
        const std::size_t nc = cb.NumCells();
        for (std::size_t r = 0; r < nc; ++r)
            for (const detail::CellFaceDef& fd : faces)
                f(conn, r * npc, fd);
    }
}

// Canonical output block order: triangle, triangle6, quad, quad8, quad9.
constexpr std::size_t SKIN_NUM_OUT_TYPES = 5;

std::size_t skin_out_type_index(CellType type) {
    switch (type) {
        case CellType::Triangle:
            return 0;
        case CellType::Triangle6:
            return 1;
        case CellType::Quad:
            return 2;
        case CellType::Quad8:
            return 3;
        default:  // Quad9
            return 4;
    }
}

const char* skin_out_type_name(std::size_t index) {
    static const char* names[SKIN_NUM_OUT_TYPES] = {"triangle", "triangle6", "quad", "quad8",
                                                    "quad9"};
    return names[index];
}

std::size_t skin_out_type_nodes(std::size_t index) {
    static const std::size_t counts[SKIN_NUM_OUT_TYPES] = {3, 6, 4, 8, 9};
    return counts[index];
}

}  // namespace

bool has_skinnable_cells(const Mesh& rMesh) {
    for (const auto cb : rMesh.CellRange()) {
        if (skin_block_supported(cell_type_from_name(cb.Type()), cb.IsRagged()))
            return true;
    }
    return false;
}

Mesh extract_skin(const Mesh& rMesh, bool linearize) {
    // Pre-scan: warn once per unsupported 3D block, require at least one
    // supported one.
    bool any_supported = false;
    for (const auto cb : rMesh.CellRange()) {
        const CellType ct = cell_type_from_name(cb.Type());
        const bool is_volume = cell_type_dimension(ct) == 3 || cb.IsPolyhedron();
        if (!is_volume)
            continue;
        if (skin_block_supported(ct, cb.IsRagged())) {
            any_supported = true;
        } else {
            log::warn("extract_skin: volume cell block '{}' is not supported; skipping it.",
                      cb.Type());
        }
    }
    if (!any_supported)
        throw std::invalid_argument(
            "extract_skin: mesh contains no supported 3D volume cell block");

    // Pass 1: count occurrences of every face key.
    std::unordered_map<SkinFaceKey, std::uint32_t, SkinFaceKeyHash> face_count;
    skin_for_each_face(
        rMesh, [&](const NDArray& rConn, std::size_t row_offset, const detail::CellFaceDef& rFace) {
            ++face_count[skin_face_key(rConn, row_offset, rFace)];
        });

    // Pass 2: collect boundary faces (count == 1) in enumeration order,
    // partitioned by output face type.
    std::array<std::vector<std::int64_t>, SKIN_NUM_OUT_TYPES> out_conn;
    skin_for_each_face(
        rMesh, [&](const NDArray& rConn, std::size_t row_offset, const detail::CellFaceDef& rFace) {
            if (face_count[skin_face_key(rConn, row_offset, rFace)] != 1)
                return;
            const std::uint8_t n_out = linearize ? rFace.mNumCorners : rFace.mNumNodes;
            const CellType out_type =
                linearize ? (rFace.mNumCorners == 3 ? CellType::Triangle : CellType::Quad)
                          : rFace.mFaceType;
            std::vector<std::int64_t>& dst = out_conn[skin_out_type_index(out_type)];
            for (std::uint8_t k = 0; k < n_out; ++k)
                dst.push_back(detail::read_int(rConn, row_offset + rFace.mNodes[k]));
        });

    // Compaction: keep only referenced points, in ascending original-id order
    // (the numpy fallback relies on np.unique's sortedness for the same order).
    const std::size_t num_points = rMesh.NumPoints();
    std::vector<char> used(num_points, 0);
    for (const std::vector<std::int64_t>& conn : out_conn)
        for (std::int64_t id : conn)
            used[static_cast<std::size_t>(id)] = 1;
    std::vector<std::int64_t> remap(num_points, -1);
    std::size_t num_used = 0;
    for (std::size_t i = 0; i < num_points; ++i)
        if (used[i])
            remap[i] = static_cast<std::int64_t>(num_used++);

    Mesh skin;

    // Points: byte-row subset preserving the source dtype.
    {
        const NDArray& points = rMesh.Points();
        const std::size_t dim = detail::cols(points);
        const std::size_t row_bytes = dim * dtype_size(points.Dtype());
        NDArray new_points = NDArray::Uninit(points.Dtype(), {num_used, dim});
        const std::byte* src = points.Data();
        std::byte* dst = new_points.Data();
        std::size_t w = 0;
        for (std::size_t i = 0; i < num_points; ++i)
            if (used[i])
                std::memcpy(dst + (w++) * row_bytes, src + i * row_bytes, row_bytes);
        skin.AssignPoints(std::move(new_points));
    }

    // Cell blocks, canonical order, remapped connectivity.
    for (std::size_t t = 0; t < SKIN_NUM_OUT_TYPES; ++t) {
        const std::vector<std::int64_t>& conn = out_conn[t];
        if (conn.empty())
            continue;
        const std::size_t npc = skin_out_type_nodes(t);
        NDArray block = NDArray::Uninit(DType::Int64, {conn.size() / npc, npc});
        std::int64_t* dst = block.As<std::int64_t>();
        for (std::size_t i = 0; i < conn.size(); ++i)
            dst[i] = remap[static_cast<std::size_t>(conn[i])];
        skin.AddCellBlock(skin_out_type_name(t), std::move(block));
    }

    // point_data: subset rows through the same remap; other data maps are
    // dropped (no 1:1 face -> source cell mapping).
    for (const std::string& name : rMesh.PointDataNames()) {
        const NDArray& a = rMesh.PointData(name);
        if (detail::rows(a) != num_points)
            continue;  // shape does not match the point table; drop
        std::vector<std::size_t> new_shape = a.Shape();
        new_shape[0] = num_used;
        const std::size_t row_bytes = num_points == 0 ? 0 : a.Nbytes() / num_points;
        NDArray b = NDArray::Uninit(a.Dtype(), std::move(new_shape));
        const std::byte* src = a.Data();
        std::byte* dst = b.Data();
        std::size_t w = 0;
        for (std::size_t i = 0; i < num_points; ++i)
            if (used[i])
                std::memcpy(dst + (w++) * row_bytes, src + i * row_bytes, row_bytes);
        skin.AddPointData(name, std::move(b));
    }

    return skin;
}

}  // namespace meshioplusplus
