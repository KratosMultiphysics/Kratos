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
#pragma once

/**
 * @file meshio_mesh.hpp
 * @brief The MESHIO mesh backend: `meshioplusplus::Mesh` and
 * `meshioplusplus::CellBlock`, the meshio-mirroring in-memory representation.
 *
 * This is the default mesh backend (see `mesh.hpp` for the compile-time
 * dispatch and `mesh_api.hpp` for the uniform format-facing API it
 * implements), and the only one compatible with the pybind11 extension —
 * `bindings/np_conversions.hpp` is written against these exact members.
 *
 * It mirrors the fields of the pure-Python `meshio.Mesh`: it is the type
 * every C++ format reader produces and every C++ format writer consumes.
 * The pybind11 binding layer (`bindings/np_conversions.hpp`) converts between
 * this type and the pure-Python `meshio.Mesh` at the I/O boundary, following
 * a "zero-copy at the boundary" strategy: `py_to_mesh` builds non-owning
 * `NDArray` *views* over the caller's numpy buffers (write path, no input
 * copy) and `mesh_to_py` moves each `NDArray`'s owned buffer into a capsule
 * backing a writeable numpy array (read path, no output copy).
 *
 * The conversion layer carries `points`, `cells`, `point_data`, `cell_data`,
 * and `field_data` — but deliberately **not** `mesh.info`, `cell_sets`, or
 * `point_sets`, which are custom attributes that live only on the Python
 * `Mesh`. Formats that need those either defer entirely to the Python
 * fallback or carry the extra data out-of-band via a side-channel struct
 * that the binding `setattr`s onto the Python `Mesh` object after
 * conversion (e.g. `MedInfo`/`AnsysInfo` for `point_sets`/`cell_sets`,
 * `OpenFoamInfo` for `cell_tags`).
 */

// System includes
#include <cstdint>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/detail/map_order.hpp"
#include "meshioplusplus/mesh_api.hpp"
#include "meshioplusplus/ndarray.hpp"

namespace meshioplusplus {

/**
 * @brief One homogeneous block of cells of a single meshio cell type.
 *
 * Mirrors a Python `meshio.CellBlock`. Most blocks are *rectangular*: `mData`
 * is a `(num_cells, nodes_per_cell)` integer `NDArray` of node indices.
 * Some formats, however, produce cells that cannot be described by a fixed
 * nodes-per-cell count, so `CellBlock` also carries two optional *ragged*
 * (jagged) representations:
 *
 *  - `mPolygonRows` — 1-level ragged: a `"polygon"` block whose cells have
 *    varying node counts (e.g. MED POG Voronoi meshes). Row `i` is the list
 *    of node ids for cell `i`.
 *  - `mPolyhedronRows` — 2-level ragged: a `"polyhedron"` block. Cell `i` is
 *    a list of faces, each face itself a list of node ids.
 *
 * Exactly one of `mData`, `mPolygonRows`, `mPolyhedronRows` is populated per
 * block (see `IsRagged()`); the unused members are left empty, so ordinary
 * rectangular blocks (the overwhelming majority) are unaffected. Zero-copy
 * numpy conversion at the binding boundary only applies to the rectangular
 * `mData` case — ragged blocks are always *copied* across the boundary, and
 * `py_to_mesh`'s `allow_ragged` flag is off by default so a rectangular-only
 * writer given a ragged mesh safely throws and triggers the Python fallback;
 * only ragged-aware bindings (e.g. MED write) opt in.
 */
struct CellBlock {
    std::string mType;  // meshio cell type, e.g. "triangle"
    NDArray mData;      // (num_cells, nodes_per_cell), integer dtype
    std::vector<std::string> mTags;

    // Ragged (jagged) representations, used only for cell types whose rows do
    // not fit a rectangular buffer. Exactly one of `mData` / `mPolygonRows` /
    // `mPolyhedronRows` is populated per block; the two ragged members are
    // empty for every rectangular block (all rectangular formats unaffected).
    //
    //  * mPolygonRows    — 1-level ragged: a "polygon" block whose cells have
    //                      varying node counts (e.g. MED POG Voronoi meshes).
    //                      Row i = mPolygonRows[i] = node ids of cell i.
    //  * mPolyhedronRows — 2-level ragged: a "polyhedron" block. Cell i is a
    //                      list of faces; each face is a list of node ids.
    std::vector<std::vector<std::int64_t>> mPolygonRows;
    std::vector<std::vector<std::vector<std::int64_t>>> mPolyhedronRows;

    CellBlock() = default;
    CellBlock(std::string t, NDArray d) : mType(std::move(t)), mData(std::move(d)) {}

    /**
     * @brief Whether this block uses one of the ragged representations.
     * @return `true` iff `mPolygonRows` or `mPolyhedronRows` is non-empty.
     */
    bool IsRagged() const { return !mPolygonRows.empty() || !mPolyhedronRows.empty(); }

    /**
     * @brief Number of cells in this block, whichever representation is active.
     * @return `mPolygonRows.size()`, else `mPolyhedronRows.size()`, else the
     *         first dimension of `mData` (0 if `mData` has no shape).
     */
    std::size_t NumCells() const {
        if (!mPolygonRows.empty())
            return mPolygonRows.size();
        if (!mPolyhedronRows.empty())
            return mPolyhedronRows.size();
        return mData.Shape().empty() ? 0 : mData.Shape()[0];
    }
};

/**
 * @brief The C++ in-memory mesh: points, cell blocks, and field data.
 *
 * Produced by every C++ format reader and consumed by every C++ format
 * writer; see the file-level comment for how this maps to/from the
 * pure-Python `meshio.Mesh` at the pybind11 boundary. Note what is
 * deliberately absent from this struct: `mesh.info`, `point_sets`, and
 * `cell_sets` are Python-only attributes not represented here (they travel,
 * when needed, through a per-format side-channel struct instead).
 */
struct Mesh {
    NDArray mPoints;  // (num_points, dim)
    std::vector<CellBlock> mCells;

    // Field data. mCellData holds one NDArray per cell block, in mCells order.
    // These are unordered_map for O(1) name lookup; where key *order* is
    // observable (Python dict order, on-disk field order) call
    // detail::sorted_keys (map_order.hpp) at the consumption site.
    std::unordered_map<std::string, NDArray> mPointData;
    std::unordered_map<std::string, std::vector<NDArray>> mCellData;
    std::unordered_map<std::string, NDArray> mFieldData;

    /**
     * @brief Number of points in the mesh.
     * @return The first dimension of `mPoints.Shape()`, or 0 if unset.
     */
    std::size_t NumPoints() const { return mPoints.Shape().empty() ? 0 : mPoints.Shape()[0]; }

    // -----------------------------------------------------------------
    // Uniform format-facing API (the compile-time contract shared by all
    // mesh backends — see mesh_api.hpp). Format code must go through these
    // methods, never the public members above; on this backend every method
    // is a trivial inline forward, so there is zero cost over direct access.
    // -----------------------------------------------------------------

    /**
     * @brief Cheap, copyable view over one cell block (see `mesh_api.hpp`).
     *
     * On this backend it simply wraps a `const CellBlock*`; it stays valid
     * only while the underlying `Mesh` is alive and its `mCells` vector is
     * not resized.
     */
    class CellView {
    public:
        explicit CellView(const CellBlock& rBlock) : mpBlock(&rBlock) {}
        /** @brief The meshio cell-type name (e.g. `"triangle"`). */
        const std::string& Type() const { return mpBlock->mType; }
        /** @brief Number of cells in the block (any representation). */
        std::size_t NumCells() const { return mpBlock->NumCells(); }
        /** @brief Nodes per cell for rectangular blocks; 0 for ragged ones. */
        std::size_t NodesPerCell() const {
            return mpBlock->mData.Ndim() >= 2 ? mpBlock->mData.Shape()[1] : 0;
        }
        /** @brief Whether the block uses a ragged representation. */
        bool IsRagged() const { return mpBlock->IsRagged(); }
        /** @brief Whether the block is 2-level ragged (list of faces per cell). */
        bool IsPolyhedron() const { return !mpBlock->mPolyhedronRows.empty(); }
        /** @brief Rectangular `(num_cells, nodes_per_cell)` connectivity (empty if ragged). */
        const NDArray& Conn() const { return mpBlock->mData; }
        /** @brief Node count of polygon cell @p cell (1-level ragged blocks). */
        std::size_t RowSize(std::size_t cell) const { return mpBlock->mPolygonRows[cell].size(); }
        /** @brief Node ids of polygon cell @p cell (1-level ragged blocks). */
        const std::int64_t* Row(std::size_t cell) const {
            return mpBlock->mPolygonRows[cell].data();
        }
        /** @brief Face count of polyhedron cell @p cell (2-level ragged blocks). */
        std::size_t NumFaces(std::size_t cell) const {
            return mpBlock->mPolyhedronRows[cell].size();
        }
        /** @brief `{node ids, count}` of face @p face of polyhedron cell @p cell. */
        std::pair<const std::int64_t*, std::size_t> Face(std::size_t cell, std::size_t face) const {
            const auto& r_face = mpBlock->mPolyhedronRows[cell][face];
            return {r_face.data(), r_face.size()};
        }

    private:
        const CellBlock* mpBlock;
    };

    // --- reader-side ingestion ---

    /** @brief Takes ownership of the point array (float dtype, shape `(n, dim)`). */
    void AssignPoints(NDArray points) { mPoints = std::move(points); }
    /** @brief Appends a rectangular cell block (integer dtype, shape `(n, npc)`). */
    void AddCellBlock(std::string type, NDArray conn) {
        mCells.emplace_back(std::move(type), std::move(conn));
    }
    /** @brief Appends a 1-level ragged (polygon) cell block. */
    void AddPolygonBlock(std::string type, std::vector<std::vector<std::int64_t>> rows) {
        CellBlock cb;
        cb.mType = std::move(type);
        cb.mPolygonRows = std::move(rows);
        mCells.push_back(std::move(cb));
    }
    /** @brief Appends a 2-level ragged (polyhedron) cell block. */
    void AddPolyhedronBlock(std::string type,
                            std::vector<std::vector<std::vector<std::int64_t>>> cells) {
        CellBlock cb;
        cb.mType = std::move(type);
        cb.mPolyhedronRows = std::move(cells);
        mCells.push_back(std::move(cb));
    }
    /** @brief Inserts or replaces a named per-point data array. */
    void AddPointData(std::string name, NDArray data) {
        mPointData[std::move(name)] = std::move(data);
    }
    /** @brief Inserts or replaces a named per-cell data array list (one per block). */
    void AddCellData(std::string name, std::vector<NDArray> blocks) {
        mCellData[std::move(name)] = std::move(blocks);
    }
    /** @brief Appends one block's array to a named cell-data list (creating it if new). */
    void AppendCellData(const std::string& rName, NDArray block) {
        mCellData[rName].push_back(std::move(block));
    }
    /** @brief Inserts or replaces a named field-data array. */
    void AddFieldData(std::string name, NDArray data) {
        mFieldData[std::move(name)] = std::move(data);
    }

    // --- writer-side accessors ---

    /** @brief Spatial dimension of the points (second shape entry), or 0 if unset. */
    std::size_t PointDim() const { return mPoints.Ndim() >= 2 ? mPoints.Shape()[1] : 0; }
    /** @brief The `(num_points, dim)` point array. */
    const NDArray& Points() const { return mPoints; }
    /** @brief Number of cell blocks. */
    std::size_t NumCellBlocks() const { return mCells.size(); }
    /** @brief View over cell block @p i (in insertion order). */
    CellView Cells(std::size_t i) const { return CellView(mCells[i]); }
    /** @brief Range over all cell blocks: `for (const auto cb : mesh.CellRange())`. */
    detail::CellBlockRange<Mesh> CellRange() const { return detail::CellBlockRange<Mesh>(*this); }

    /** @brief Point-data names in sorted order (drives on-disk field order). */
    std::vector<std::string> PointDataNames() const { return detail::sorted_keys(mPointData); }
    /** @brief Number of named point-data arrays. */
    std::size_t NumPointData() const { return mPointData.size(); }
    /** @brief Whether a point-data array named @p rName exists. */
    bool HasPointData(const std::string& rName) const { return mPointData.count(rName) > 0; }
    /** @brief The point-data array named @p rName (throws if absent). */
    const NDArray& PointData(const std::string& rName) const { return mPointData.at(rName); }

    /** @brief Cell-data names in sorted order (drives on-disk field order). */
    std::vector<std::string> CellDataNames() const { return detail::sorted_keys(mCellData); }
    /** @brief Number of named cell-data array lists. */
    std::size_t NumCellData() const { return mCellData.size(); }
    /** @brief Whether a cell-data list named @p rName exists. */
    bool HasCellData(const std::string& rName) const { return mCellData.count(rName) > 0; }
    /** @brief Block @p block of the cell-data list named @p rName (throws if absent). */
    const NDArray& CellData(const std::string& rName, std::size_t block) const {
        return mCellData.at(rName)[block];
    }
    /** @brief Number of blocks in the cell-data list named @p rName (throws if absent). */
    std::size_t CellDataNumBlocks(const std::string& rName) const {
        return mCellData.at(rName).size();
    }

    /** @brief Field-data names in sorted order (drives on-disk field order). */
    std::vector<std::string> FieldDataNames() const { return detail::sorted_keys(mFieldData); }
    /** @brief Number of named field-data arrays. */
    std::size_t NumFieldData() const { return mFieldData.size(); }
    /** @brief Whether a field-data array named @p rName exists. */
    bool HasFieldData(const std::string& rName) const { return mFieldData.count(rName) > 0; }
    /** @brief The field-data array named @p rName (throws if absent). */
    const NDArray& FieldData(const std::string& rName) const { return mFieldData.at(rName); }
};

}  // namespace meshioplusplus
