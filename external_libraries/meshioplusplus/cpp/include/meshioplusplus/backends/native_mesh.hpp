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
 * @file native_mesh.hpp
 * @brief The NATIVE mesh backend: `meshioplusplus::NativeMesh`, a canonical
 * statically-typed in-memory mesh built for downstream C++ consumers.
 *
 * Selected by `MESHIOPLUSPLUS_MESH_BACKEND=NATIVE` (see `mesh.hpp`); the
 * WebAssembly build uses it. Where the MESHIO backend stores every array
 * with whatever dtype the file supplied (so the Python boundary can be
 * zero-copy), NATIVE canonicalizes at ingest — points are always contiguous
 * `double`, connectivity always `std::int64_t`, data arrays Float64/Int64
 * (never int -> float, so integer tag conventions survive) — and identifies
 * cell types by the `CellType` enum instead of a string. An owning array
 * that is already canonical is *moved* in, not copied, and format readers
 * produce canonical dtypes almost everywhere, so ingest is near-free.
 *
 * Ragged (polygon/polyhedron) blocks are stored CSR-style — one flat node
 * buffer plus offset arrays — rather than nested vectors: one allocation
 * per level, cache-friendly iteration, and the natural shape a FEM/graphics
 * consumer wants. On top of the uniform format-facing API (`mesh_api.hpp`)
 * it adds a fast-consumer surface: `PointsData()`, `ConnSpan()`,
 * `BlockType()`, and a lazily-built whole-mesh CSR (`GlobalConnectivity()`).
 */

// System includes
#include <cstdint>
#include <optional>
#include <span>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/named_arrays.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/mesh_api.hpp"
#include "meshioplusplus/ndarray.hpp"

namespace meshioplusplus {

namespace detail {

/**
 * @brief Canonicalize an array *within kind*: float dtypes -> Float64,
 * integer dtypes -> Int64.
 *
 * Already-canonical owning arrays are moved through untouched (the fast
 * path); views are copied into owned canonical storage so the mesh always
 * owns its memory.
 * @param a The array to canonicalize (consumed).
 * @return An owning canonical array.
 */
inline NDArray canonicalize_array(NDArray a) {
    const DType target = is_float_dtype(a.Dtype()) ? DType::Float64 : DType::Int64;
    if (a.Dtype() == target) {
        a.MakeOwned();  // no-op when already owning
        return a;
    }
    NDArray out = NDArray::Uninit(target, a.Shape());
    const std::size_t n = a.Size();
    dispatch_dtype(a.Dtype(), [&]<class T>() {
        const T* src = a.As<T>();
        if (target == DType::Float64) {
            double* dst = out.As<double>();
            for (std::size_t i = 0; i < n; ++i)
                dst[i] = static_cast<double>(src[i]);
        } else {
            std::int64_t* dst = out.As<std::int64_t>();
            for (std::size_t i = 0; i < n; ++i)
                dst[i] = static_cast<std::int64_t>(src[i]);
        }
    });
    return out;
}

/** @brief Canonicalize to Float64 regardless of kind (for point arrays). */
inline NDArray canonicalize_float64(NDArray a) {
    if (a.Dtype() == DType::Float64) {
        a.MakeOwned();
        return a;
    }
    NDArray out = NDArray::Uninit(DType::Float64, a.Shape());
    const std::size_t n = a.Size();
    dispatch_dtype(a.Dtype(), [&]<class T>() {
        const T* src = a.As<T>();
        double* dst = out.As<double>();
        for (std::size_t i = 0; i < n; ++i)
            dst[i] = static_cast<double>(src[i]);
    });
    return out;
}

/** @brief Canonicalize to Int64 regardless of kind (for connectivity). */
inline NDArray canonicalize_int64(NDArray a) {
    if (a.Dtype() == DType::Int64) {
        a.MakeOwned();
        return a;
    }
    NDArray out = NDArray::Uninit(DType::Int64, a.Shape());
    const std::size_t n = a.Size();
    dispatch_dtype(a.Dtype(), [&]<class T>() {
        const T* src = a.As<T>();
        std::int64_t* dst = out.As<std::int64_t>();
        for (std::size_t i = 0; i < n; ++i)
            dst[i] = static_cast<std::int64_t>(src[i]);
    });
    return out;
}

}  // namespace detail

/**
 * @brief One cell block of a `NativeMesh`.
 *
 * Rectangular blocks live in `mConn` (Int64, `(n, nodes_per_cell)`).
 * Ragged blocks are CSR-shaped in `mFlat`/`mRowOffsets` (polygon: row `i`'s
 * nodes are `mFlat[mRowOffsets[i] .. mRowOffsets[i+1])`) with the extra
 * `mFaceOffsets` level for polyhedra (cell `c`'s faces are rows
 * `mFaceOffsets[c] .. mFaceOffsets[c+1])` of `mRowOffsets`).
 */
struct NativeCellBlock {
    CellType mType = CellType::Custom;
    std::string mTypeName;                   // canonical meshio name; preserves spellings the
                                             // enum can't represent (e.g. "polyhedron12")
    NDArray mConn;                           // Int64 (n, nodes_per_cell); empty for ragged blocks
    std::vector<std::int64_t> mFlat;         // ragged: all node ids, row-major
    std::vector<std::int64_t> mRowOffsets;   // ragged: nrows+1 offsets into mFlat
    std::vector<std::int64_t> mFaceOffsets;  // polyhedron only: ncells+1 offsets
                                             // into mRowOffsets' rows

    bool IsRagged() const { return !mRowOffsets.empty(); }
    bool IsPolyhedron() const { return !mFaceOffsets.empty(); }
    std::size_t NumCells() const {
        if (IsPolyhedron())
            return mFaceOffsets.size() - 1;
        if (IsRagged())
            return mRowOffsets.size() - 1;
        return mConn.Shape().empty() ? 0 : mConn.Shape()[0];
    }
};

/**
 * @brief The NATIVE mesh backend (aliased to `meshioplusplus::Mesh` when
 * `MESHIOPLUSPLUS_MESH_BACKEND_NATIVE` is defined).
 *
 * Implements the uniform format-facing API (`mesh_api.hpp`) over canonical
 * statically-typed storage, plus a fast-consumer surface for downstream C++
 * users. Data-array names are stored insertion-ordered with O(1) lookup;
 * the API's observable order is sorted, like every backend.
 */
class NativeMesh {
public:
    // --- uniform API: reader-side ingestion -------------------------------

    /** @brief Takes ownership of the point array, canonicalized to Float64. */
    void AssignPoints(NDArray points) { mPoints = detail::canonicalize_float64(std::move(points)); }
    /** @brief Appends a rectangular cell block, connectivity canonicalized to Int64. */
    void AddCellBlock(std::string type, NDArray conn) {
        NativeCellBlock b;
        b.mType = cell_type_from_name(type);
        b.mTypeName = std::move(type);
        b.mConn = detail::canonicalize_int64(std::move(conn));
        mBlocks.push_back(std::move(b));
        mGlobalCsr.reset();
    }
    /** @brief Appends a 1-level ragged (polygon) block, stored CSR-style. */
    void AddPolygonBlock(std::string type, std::vector<std::vector<std::int64_t>> rows) {
        NativeCellBlock b;
        b.mType = cell_type_from_name(type);
        b.mTypeName = std::move(type);
        std::size_t total = 0;
        for (const auto& r_row : rows)
            total += r_row.size();
        b.mFlat.reserve(total);
        b.mRowOffsets.reserve(rows.size() + 1);
        b.mRowOffsets.push_back(0);
        for (const auto& r_row : rows) {
            b.mFlat.insert(b.mFlat.end(), r_row.begin(), r_row.end());
            b.mRowOffsets.push_back(static_cast<std::int64_t>(b.mFlat.size()));
        }
        mBlocks.push_back(std::move(b));
        mGlobalCsr.reset();
    }
    /** @brief Appends a 2-level ragged (polyhedron) block, stored CSR-style. */
    void AddPolyhedronBlock(std::string type,
                            std::vector<std::vector<std::vector<std::int64_t>>> cells) {
        NativeCellBlock b;
        b.mType = cell_type_from_name(type);
        b.mTypeName = std::move(type);
        b.mFaceOffsets.reserve(cells.size() + 1);
        b.mFaceOffsets.push_back(0);
        std::size_t nrows = 0;
        for (const auto& r_cell : cells)
            nrows += r_cell.size();
        b.mRowOffsets.reserve(nrows + 1);
        b.mRowOffsets.push_back(0);
        for (const auto& r_cell : cells) {
            for (const auto& r_face : r_cell) {
                b.mFlat.insert(b.mFlat.end(), r_face.begin(), r_face.end());
                b.mRowOffsets.push_back(static_cast<std::int64_t>(b.mFlat.size()));
            }
            b.mFaceOffsets.push_back(static_cast<std::int64_t>(b.mRowOffsets.size() - 1));
        }
        mBlocks.push_back(std::move(b));
        mGlobalCsr.reset();
    }
    /** @brief Inserts or replaces a named per-point data array (canonicalized). */
    void AddPointData(std::string name, NDArray data) {
        mPointData.Set(std::move(name), detail::canonicalize_array(std::move(data)));
    }
    /** @brief Inserts or replaces a named per-cell data array list (canonicalized). */
    void AddCellData(std::string name, std::vector<NDArray> blocks) {
        for (auto& r_b : blocks)
            r_b = detail::canonicalize_array(std::move(r_b));
        mCellData.Set(std::move(name), std::move(blocks));
    }
    /** @brief Appends one block's array to a named cell-data list (canonicalized). */
    void AppendCellData(const std::string& rName, NDArray block) {
        mCellData.GetOrCreate(rName).push_back(detail::canonicalize_array(std::move(block)));
    }
    /** @brief Inserts or replaces a named field-data array (canonicalized). */
    void AddFieldData(std::string name, NDArray data) {
        mFieldData.Set(std::move(name), detail::canonicalize_array(std::move(data)));
    }

    // --- uniform API: writer-side accessors -------------------------------

    /** @brief Cheap, copyable view over one cell block (see `mesh_api.hpp`). */
    class CellView {
    public:
        explicit CellView(const NativeCellBlock& rBlock) : mpBlock(&rBlock) {}
        const std::string& Type() const { return mpBlock->mTypeName; }
        std::size_t NumCells() const { return mpBlock->NumCells(); }
        std::size_t NodesPerCell() const {
            return mpBlock->mConn.Ndim() >= 2 ? mpBlock->mConn.Shape()[1] : 0;
        }
        bool IsRagged() const { return mpBlock->IsRagged(); }
        bool IsPolyhedron() const { return mpBlock->IsPolyhedron(); }
        const NDArray& Conn() const { return mpBlock->mConn; }
        std::size_t RowSize(std::size_t cell) const {
            return static_cast<std::size_t>(mpBlock->mRowOffsets[cell + 1] -
                                            mpBlock->mRowOffsets[cell]);
        }
        const std::int64_t* Row(std::size_t cell) const {
            return mpBlock->mFlat.data() + mpBlock->mRowOffsets[cell];
        }
        std::size_t NumFaces(std::size_t cell) const {
            return static_cast<std::size_t>(mpBlock->mFaceOffsets[cell + 1] -
                                            mpBlock->mFaceOffsets[cell]);
        }
        std::pair<const std::int64_t*, std::size_t> Face(std::size_t cell, std::size_t face) const {
            const std::size_t row = static_cast<std::size_t>(mpBlock->mFaceOffsets[cell]) + face;
            return {mpBlock->mFlat.data() + mpBlock->mRowOffsets[row],
                    static_cast<std::size_t>(mpBlock->mRowOffsets[row + 1] -
                                             mpBlock->mRowOffsets[row])};
        }

    private:
        const NativeCellBlock* mpBlock;
    };

    std::size_t NumPoints() const { return mPoints.Shape().empty() ? 0 : mPoints.Shape()[0]; }
    std::size_t PointDim() const { return mPoints.Ndim() >= 2 ? mPoints.Shape()[1] : 0; }
    const NDArray& Points() const { return mPoints; }
    std::size_t NumCellBlocks() const { return mBlocks.size(); }
    CellView Cells(std::size_t i) const { return CellView(mBlocks[i]); }
    detail::CellBlockRange<NativeMesh> CellRange() const {
        return detail::CellBlockRange<NativeMesh>(*this);
    }

    std::vector<std::string> PointDataNames() const { return mPointData.SortedNames(); }
    std::size_t NumPointData() const { return mPointData.Size(); }
    bool HasPointData(const std::string& rName) const { return mPointData.Has(rName); }
    const NDArray& PointData(const std::string& rName) const { return mPointData.Get(rName); }

    std::vector<std::string> CellDataNames() const { return mCellData.SortedNames(); }
    std::size_t NumCellData() const { return mCellData.Size(); }
    bool HasCellData(const std::string& rName) const { return mCellData.Has(rName); }
    const NDArray& CellData(const std::string& rName, std::size_t block) const {
        return mCellData.Get(rName)[block];
    }
    std::size_t CellDataNumBlocks(const std::string& rName) const {
        return mCellData.Get(rName).size();
    }

    std::vector<std::string> FieldDataNames() const { return mFieldData.SortedNames(); }
    std::size_t NumFieldData() const { return mFieldData.Size(); }
    bool HasFieldData(const std::string& rName) const { return mFieldData.Has(rName); }
    const NDArray& FieldData(const std::string& rName) const { return mFieldData.Get(rName); }

    // --- fast-consumer surface (NATIVE-only extras) -----------------------

    /** @brief Contiguous `(NumPoints() * PointDim())` Float64 coordinate buffer. */
    const double* PointsData() const { return mPoints.As<double>(); }
    /** @brief Block @p block's rectangular Int64 connectivity as a span. */
    std::span<const std::int64_t> ConnSpan(std::size_t block) const {
        const NDArray& conn = mBlocks[block].mConn;
        return {conn.As<std::int64_t>(), conn.Size()};
    }
    /** @brief Block @p block's cell type as the compact enum. */
    CellType BlockType(std::size_t block) const { return mBlocks[block].mType; }

    /**
     * @brief Whole-mesh CSR connectivity over all *rectangular* blocks, in
     * block order: cell `i`'s nodes are `mConn[mOffsets[i] .. mOffsets[i+1])`
     * and its type `mTypes[i]`. Ragged blocks are skipped.
     */
    struct GlobalCsr {
        std::vector<std::int64_t> mOffsets;  // ncells+1
        std::vector<std::int64_t> mConn;     // flat node ids
        std::vector<CellType> mTypes;        // one per cell
    };

    /**
     * @brief The whole-mesh CSR, built lazily on first call and cached
     * (invalidated whenever a block is added).
     * @return Reference to the cached CSR (valid until the next mutation).
     */
    const GlobalCsr& GlobalConnectivity() const {
        if (!mGlobalCsr) {
            GlobalCsr csr;
            std::size_t ncells = 0, nconn = 0;
            for (const auto& r_b : mBlocks) {
                if (r_b.IsRagged())
                    continue;
                ncells += r_b.NumCells();
                nconn += r_b.mConn.Size();
            }
            csr.mOffsets.reserve(ncells + 1);
            csr.mConn.reserve(nconn);
            csr.mTypes.reserve(ncells);
            csr.mOffsets.push_back(0);
            for (const auto& r_b : mBlocks) {
                if (r_b.IsRagged())
                    continue;
                const std::size_t n = r_b.NumCells();
                const std::size_t k = r_b.mConn.Ndim() >= 2 ? r_b.mConn.Shape()[1] : 0;
                const std::int64_t* src = r_b.mConn.As<std::int64_t>();
                for (std::size_t c = 0; c < n; ++c) {
                    csr.mConn.insert(csr.mConn.end(), src + c * k, src + (c + 1) * k);
                    csr.mOffsets.push_back(static_cast<std::int64_t>(csr.mConn.size()));
                    csr.mTypes.push_back(r_b.mType);
                }
            }
            mGlobalCsr = std::move(csr);
        }
        return *mGlobalCsr;
    }

    /** @brief The per-block storage (for direct fast-path consumers). */
    const std::vector<NativeCellBlock>& Blocks() const { return mBlocks; }

private:
    NDArray mPoints;  // always Float64 (n, dim)
    std::vector<NativeCellBlock> mBlocks;
    detail::NamedArrays mPointData;     // canonical Float64/Int64 arrays
    detail::NamedArrayLists mCellData;  // one array per block, block order
    detail::NamedArrays mFieldData;
    mutable std::optional<GlobalCsr> mGlobalCsr;
};

}  // namespace meshioplusplus
