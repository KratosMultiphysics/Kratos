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
 * @file kratos_mesh.hpp
 * @brief The KRATOS mesh backend: `meshioplusplus::KratosMesh`, the uniform
 * mesh API implemented over a Kratos-style `ModelPart`.
 *
 * Selected by `MESHIOPLUSPLUS_MESH_BACKEND=KRATOS` (see `mesh.hpp`). Format
 * readers ingest into a canonical staging structure (a `NativeMesh` — same
 * canonical Float64/Int64 storage), and the `ModelPart` is **materialized
 * lazily** on the first `GetModelPart()` call:
 *
 *  - Nodes get Ids `index + 1` (z = 0-padded for 2-D points).
 *  - Cell blocks whose topological dimension equals the mesh's maximum
 *    become **Elements**; lower-dimension blocks become **Conditions** (the
 *    Kratos convention, matching `src/meshioplusplus/mdpa/_mdpa.py`), each
 *    kind Id-numbered 1..N in block order, with default Kratos names from
 *    `kratos_names.hpp`.
 *  - `point_data` becomes nodal data; `cell_data` becomes elemental /
 *    conditional data (concatenated per kind, entity order); `field_data`
 *    stays on the staging mesh (Kratos has no equivalent).
 *  - **Integer tag arrays** under well-known names (`gmsh:physical`,
 *    `su2:tag`, `medit:ref`, `cell_tags`, ...) automatically become named
 *    SubModelParts (`gmsh_physical_1`, ...) containing the tagged entities
 *    and their nodes — disable with `SetBuildSubModelPartsFromTags(false)`.
 *    The tag arrays remain as elemental/conditional data either way, so
 *    writer round-trips are unaffected.
 *  - **Ragged blocks** (polygon/polyhedron) have no ModelPart geometry;
 *    they stay in staging pass-through (round-trips keep working) and do
 *    not become entities.
 *
 * Writer accessors serve from the staging structure, so a read -> write
 * round-trip never pays for (or builds) the ModelPart at all, and output
 * bytes match the NATIVE backend exactly. After mutating the ModelPart
 * directly, call `InvalidateBlocks()`: the staging is then rebuilt from the
 * ModelPart on the next accessor use (consecutive same-type Elements are
 * grouped into blocks, then Conditions; ragged pass-through blocks and
 * SubModelPart structure are not representable back and are dropped —
 * a documented sharp edge of the mutation path).
 */

// System includes
#include <algorithm>
#include <cstdint>
#include <cstring>
#include <memory>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/backends/model_part.hpp"
#include "meshioplusplus/backends/native_mesh.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/mesh_api.hpp"
#include "meshioplusplus/ndarray.hpp"
#include "meshioplusplus/types.hpp"

namespace meshioplusplus {

/**
 * @brief The KRATOS mesh backend (aliased to `meshioplusplus::Mesh` when
 * `MESHIOPLUSPLUS_MESH_BACKEND_KRATOS` is defined). See the file-level
 * comment for the staging/materialization design.
 */
class KratosMesh {
public:
    /** @brief Cell-data names treated as entity tags for automatic SubModelParts. */
    static const std::vector<std::string>& KnownTagKeys() {
        static const std::vector<std::string> keys = {
            "cell_tags", "gmsh:physical", "su2:tag", "medit:ref",  "avsucd:material", "freefem:ref",
            "mfm:ref",   "netgen:index",  "pf3:ref", "tetgen:ref", "ugrid:ref",       "unv:pid",
        };
        return keys;
    }

    // --- uniform API: reader-side ingestion (forwarded to staging) ---------

    void AssignPoints(NDArray points) {
        ResetModelPartOnly();
        mStage.AssignPoints(std::move(points));
    }
    void AddCellBlock(std::string type, NDArray conn) {
        ResetModelPartOnly();
        mStage.AddCellBlock(std::move(type), std::move(conn));
    }
    void AddPolygonBlock(std::string type, std::vector<std::vector<std::int64_t>> rows) {
        ResetModelPartOnly();
        mStage.AddPolygonBlock(std::move(type), std::move(rows));
    }
    void AddPolyhedronBlock(std::string type,
                            std::vector<std::vector<std::vector<std::int64_t>>> cells) {
        ResetModelPartOnly();
        mStage.AddPolyhedronBlock(std::move(type), std::move(cells));
    }
    void AddPointData(std::string name, NDArray data) {
        ResetModelPartOnly();
        mStage.AddPointData(std::move(name), std::move(data));
    }
    void AddCellData(std::string name, std::vector<NDArray> blocks) {
        ResetModelPartOnly();
        mStage.AddCellData(std::move(name), std::move(blocks));
    }
    void AppendCellData(const std::string& rName, NDArray block) {
        ResetModelPartOnly();
        mStage.AppendCellData(rName, std::move(block));
    }
    void AddFieldData(std::string name, NDArray data) {
        ResetModelPartOnly();
        mStage.AddFieldData(std::move(name), std::move(data));
    }

    // --- uniform API: writer-side accessors (staging, stale-synced) --------

    using CellView = NativeMesh::CellView;

    std::size_t NumPoints() const { return Stage().NumPoints(); }
    std::size_t PointDim() const { return Stage().PointDim(); }
    const NDArray& Points() const { return Stage().Points(); }
    std::size_t NumCellBlocks() const { return Stage().NumCellBlocks(); }
    CellView Cells(std::size_t i) const { return Stage().Cells(i); }
    detail::CellBlockRange<KratosMesh> CellRange() const {
        return detail::CellBlockRange<KratosMesh>(*this);
    }

    std::vector<std::string> PointDataNames() const { return Stage().PointDataNames(); }
    std::size_t NumPointData() const { return Stage().NumPointData(); }
    bool HasPointData(const std::string& rName) const { return Stage().HasPointData(rName); }
    const NDArray& PointData(const std::string& rName) const { return Stage().PointData(rName); }

    std::vector<std::string> CellDataNames() const { return Stage().CellDataNames(); }
    std::size_t NumCellData() const { return Stage().NumCellData(); }
    bool HasCellData(const std::string& rName) const { return Stage().HasCellData(rName); }
    const NDArray& CellData(const std::string& rName, std::size_t block) const {
        return Stage().CellData(rName, block);
    }
    std::size_t CellDataNumBlocks(const std::string& rName) const {
        return Stage().CellDataNumBlocks(rName);
    }

    std::vector<std::string> FieldDataNames() const { return Stage().FieldDataNames(); }
    std::size_t NumFieldData() const { return Stage().NumFieldData(); }
    bool HasFieldData(const std::string& rName) const { return Stage().HasFieldData(rName); }
    const NDArray& FieldData(const std::string& rName) const { return Stage().FieldData(rName); }

    // --- KRATOS-specific surface -------------------------------------------

    /**
     * @brief The ModelPart view of this mesh, materialized on first call.
     * @return The root `ModelPart` (named "Main").
     */
    ModelPart& GetModelPart() {
        EnsureStage();  // a pending user mutation must be folded in first
        if (!mMaterialized)
            Materialize();
        return *mpRoot;
    }
    /** @brief Whether `GetModelPart()` has materialized the ModelPart yet. */
    bool IsMaterialized() const { return mMaterialized; }

    /**
     * @brief Declare that the ModelPart was mutated directly: the block/
     * point staging is rebuilt from the ModelPart on the next accessor use
     * (grouping consecutive same-type Elements, then Conditions; ragged
     * pass-through blocks and SubModelPart structure are dropped).
     */
    void InvalidateBlocks() {
        if (mMaterialized)
            mStale = true;
    }

    /**
     * @brief Enable/disable automatic tag -> SubModelPart creation at
     * materialization (default: enabled). Call before `GetModelPart()`.
     */
    void SetBuildSubModelPartsFromTags(bool enable) {
        mTagsToSubModelParts = enable;
        if (mMaterialized && !mStale) {
            mMaterialized = false;  // re-materialize with the new setting
            mpRoot.reset();
        }
    }
    bool BuildSubModelPartsFromTags() const { return mTagsToSubModelParts; }

private:
    /** @brief Per staged block: what it became in the ModelPart. */
    struct BlockRecord {
        enum class Kind { Element, Condition, Ragged } mKind = Kind::Ragged;
        IndexType mFirstId = 0;  // first entity Id (Element/Condition kinds)
        std::size_t mCount = 0;
    };

    void ResetModelPartOnly() {
        // Ingestion after materialization restarts the ModelPart view.
        if (mMaterialized) {
            mpRoot.reset();
            mRecords.clear();
            mMaterialized = false;
            mStale = false;
        }
    }

    const NativeMesh& Stage() const {
        EnsureStage();
        return mStage;
    }
    void EnsureStage() const {
        if (mStale) {
            RebuildStageFromModelPart();
            mStale = false;
        }
    }

    /** @brief Topological dimension of a staged block (Custom via the name table). */
    static int BlockDimension(const NativeCellBlock& rBlock) {
        const int dim = cell_type_dimension(rBlock.mType);
        if (dim >= 0)
            return dim;
        auto it = topological_dimension().find(rBlock.mTypeName);
        if (it != topological_dimension().end())
            return it->second;
        return 3;
    }

    void Materialize() {
        mpRoot = std::make_unique<ModelPart>("Main");
        mRecords.clear();
        ModelPart& r_mp = *mpRoot;

        // Nodes: Id = index + 1, z zero-padded for 2-D points.
        const NDArray& points = mStage.Points();
        const std::size_t npts = mStage.NumPoints();
        const std::size_t dim = mStage.PointDim();
        const double* p = npts ? points.As<double>() : nullptr;
        for (std::size_t i = 0; i < npts; ++i)
            r_mp.CreateNewNode(i + 1, p[i * dim], dim > 1 ? p[i * dim + 1] : 0.0,
                               dim > 2 ? p[i * dim + 2] : 0.0);

        // Elements/Conditions split: block dim == mesh max dim -> Element.
        int mesh_dim = 0;
        for (const auto& r_b : mStage.Blocks())
            mesh_dim = std::max(mesh_dim, BlockDimension(r_b));

        IndexType next_elem = 1, next_cond = 1;
        std::size_t n_elem_rows = 0, n_cond_rows = 0;
        for (const auto& r_b : mStage.Blocks()) {
            BlockRecord rec;
            rec.mCount = r_b.NumCells();
            if (r_b.IsRagged()) {
                rec.mKind = BlockRecord::Kind::Ragged;  // pass-through, no entities
            } else {
                const bool is_elem = BlockDimension(r_b) == mesh_dim;
                rec.mKind = is_elem ? BlockRecord::Kind::Element : BlockRecord::Kind::Condition;
                rec.mFirstId = is_elem ? next_elem : next_cond;
                const std::size_t n = r_b.NumCells();
                const std::size_t k = r_b.mConn.Ndim() >= 2 ? r_b.mConn.Shape()[1] : 0;
                const std::int64_t* conn = r_b.mConn.As<std::int64_t>();
                for (std::size_t c = 0; c < n; ++c) {
                    std::vector<IndexType> ids(k);
                    for (std::size_t j = 0; j < k; ++j)
                        ids[j] = static_cast<IndexType>(conn[c * k + j]) + 1;
                    if (is_elem)
                        r_mp.CreateNewElement(r_b.mType, next_elem++, std::move(ids));
                    else
                        r_mp.CreateNewCondition(r_b.mType, next_cond++, std::move(ids));
                }
                (is_elem ? n_elem_rows : n_cond_rows) += n;
            }
            mRecords.push_back(rec);
        }

        // point_data -> nodal data (shared row order: node index).
        for (const auto& r_name : mStage.PointDataNames())
            r_mp.SetNodalData(r_name, mStage.PointData(r_name));

        // cell_data -> elemental/conditional data, concatenated per kind in
        // entity order (== block order within each kind).
        for (const auto& r_name : mStage.CellDataNames()) {
            if (mStage.CellDataNumBlocks(r_name) != mStage.NumCellBlocks())
                continue;  // partial data cannot be aligned with entities
            if (n_elem_rows > 0) {
                NDArray col = ConcatKind(r_name, BlockRecord::Kind::Element, n_elem_rows);
                if (col.Size() > 0)
                    r_mp.SetElementalData(r_name, std::move(col));
            }
            if (n_cond_rows > 0) {
                NDArray col = ConcatKind(r_name, BlockRecord::Kind::Condition, n_cond_rows);
                if (col.Size() > 0)
                    r_mp.SetConditionalData(r_name, std::move(col));
            }
        }

        // Integer tags -> SubModelParts (unless disabled).
        if (mTagsToSubModelParts)
            for (const auto& r_key : KnownTagKeys())
                if (mStage.HasCellData(r_key) &&
                    mStage.CellDataNumBlocks(r_key) == mStage.NumCellBlocks())
                    BuildSubModelPartsFor(r_key);

        mMaterialized = true;
    }

    /** @brief Concatenate one cell-data name's arrays over blocks of one kind. */
    NDArray ConcatKind(const std::string& rName, BlockRecord::Kind kind,
                       std::size_t totalRows) const {
        // Determine dtype/trailing shape from the first contributing block;
        // bail out (empty result) if the blocks disagree on dtype.
        const NDArray* p_first = nullptr;
        for (std::size_t b = 0; b < mRecords.size(); ++b) {
            if (mRecords[b].mKind != kind || mRecords[b].mCount == 0)
                continue;
            const NDArray& blk = mStage.CellData(rName, b);
            if (!p_first)
                p_first = &blk;
            else if (blk.Dtype() != p_first->Dtype())
                return NDArray{};
        }
        if (!p_first || totalRows == 0)
            return NDArray{};
        std::vector<std::size_t> shape = p_first->Shape();
        if (shape.empty())
            shape = {0};
        shape[0] = totalRows;
        NDArray out = NDArray::Uninit(p_first->Dtype(), shape);
        std::size_t off = 0;
        for (std::size_t b = 0; b < mRecords.size(); ++b) {
            if (mRecords[b].mKind != kind)
                continue;
            const NDArray& blk = mStage.CellData(rName, b);
            std::memcpy(out.Data() + off, blk.Data(), blk.Nbytes());
            off += blk.Nbytes();
        }
        return out;
    }

    void BuildSubModelPartsFor(const std::string& rKey) {
        // Only integer-kind scalar-per-cell arrays qualify as tags.
        for (std::size_t b = 0; b < mRecords.size(); ++b) {
            if (mRecords[b].mKind == BlockRecord::Kind::Ragged)
                continue;
            const NDArray& a = mStage.CellData(rKey, b);
            if (detail::is_float_dtype(a.Dtype()) || a.Size() != mRecords[b].mCount)
                return;
        }
        std::string prefix = rKey;
        for (char& r_c : prefix)
            if (r_c == ':')
                r_c = '_';

        // tag value -> (element ids, condition ids)
        struct Members {
            std::vector<IndexType> mElems, mConds;
        };
        std::unordered_map<std::int64_t, Members> groups;
        std::vector<std::int64_t> order;  // first-seen order for determinism
        for (std::size_t b = 0; b < mRecords.size(); ++b) {
            const BlockRecord& rec = mRecords[b];
            if (rec.mKind == BlockRecord::Kind::Ragged)
                continue;
            const NDArray& a = mStage.CellData(rKey, b);
            for (std::size_t c = 0; c < rec.mCount; ++c) {
                const std::int64_t tag = detail::read_int(a, c);
                auto [it, inserted] = groups.try_emplace(tag);
                if (inserted)
                    order.push_back(tag);
                if (rec.mKind == BlockRecord::Kind::Element)
                    it->second.mElems.push_back(rec.mFirstId + c);
                else
                    it->second.mConds.push_back(rec.mFirstId + c);
            }
        }
        for (const std::int64_t tag : order) {
            const std::string name = prefix + "_" + std::to_string(tag);
            if (mpRoot->HasSubModelPart(name))
                continue;  // an earlier tag key already claimed the name
            ModelPart& r_smp = mpRoot->CreateSubModelPart(name);
            const Members& r_m = groups.at(tag);
            r_smp.AddElements(r_m.mElems);
            r_smp.AddConditions(r_m.mConds);
            // Kratos convention: a sub model part contains its entities' nodes.
            std::vector<IndexType> node_ids;
            detail::IdList seen;
            for (IndexType eid : r_m.mElems)
                for (IndexType nid : mpRoot->GetElement(eid).NodeIds())
                    if (seen.Add(nid))
                        node_ids.push_back(nid);
            for (IndexType cid : r_m.mConds)
                for (IndexType nid : mpRoot->GetCondition(cid).NodeIds())
                    if (seen.Add(nid))
                        node_ids.push_back(nid);
            r_smp.AddNodes(node_ids);
        }
    }

    void RebuildStageFromModelPart() const {
        const ModelPart& r_mp = *mpRoot;
        NativeMesh fresh;

        // Points from nodes in container order; node Id -> 0-based index.
        const std::size_t n = r_mp.Nodes().Size();
        NDArray pts = NDArray::Uninit(DType::Float64, {n, 3});
        double* p = pts.As<double>();
        std::size_t i = 0;
        for (const Node& r_node : r_mp.Nodes()) {
            p[i * 3 + 0] = r_node.X();
            p[i * 3 + 1] = r_node.Y();
            p[i * 3 + 2] = r_node.Z();
            ++i;
        }
        fresh.AssignPoints(std::move(pts));

        // Consecutive same-type runs -> blocks (Elements first, then
        // Conditions), matching the materialization convention.
        mRecords.clear();
        AppendEntityBlocks(r_mp, r_mp.Elements(), BlockRecord::Kind::Element, fresh);
        AppendEntityBlocks(r_mp, r_mp.Conditions(), BlockRecord::Kind::Condition, fresh);

        // Nodal data -> point_data; elemental/conditional -> per-block slices.
        for (const auto& r_name : r_mp.NodalDataNames())
            fresh.AddPointData(r_name, r_mp.GetNodalData(r_name));
        RestoreCellData(r_mp.ElementalDataNames(), BlockRecord::Kind::Element, r_mp, fresh);
        RestoreCellData(r_mp.ConditionalDataNames(), BlockRecord::Kind::Condition, r_mp, fresh);
        for (const auto& r_name : mStage.FieldDataNames())
            fresh.AddFieldData(r_name, mStage.FieldData(r_name));

        mStage = std::move(fresh);
    }

    template <class TContainer>
    void AppendEntityBlocks(const ModelPart& rMp, const TContainer& rEntities,
                            BlockRecord::Kind kind, NativeMesh& rOut) const {
        std::vector<const GeometricalEntity*> run;
        auto flush = [&]() {
            if (run.empty())
                return;
            const std::size_t nc = run.size();
            const std::size_t k = run.front()->NumberOfNodes();
            NDArray conn = NDArray::Uninit(DType::Int64, {nc, k});
            std::int64_t* c = conn.As<std::int64_t>();
            for (std::size_t r = 0; r < nc; ++r)
                for (std::size_t j = 0; j < k; ++j)
                    c[r * k + j] =
                        static_cast<std::int64_t>(rMp.Nodes().IndexOf(run[r]->NodeIds()[j]));
            const CellType type = run.front()->Type();
            BlockRecord rec;
            rec.mKind = kind;
            rec.mFirstId = run.front()->Id();
            rec.mCount = nc;
            mRecords.push_back(rec);
            rOut.AddCellBlock(cell_type_name(type), std::move(conn));
            run.clear();
        };
        for (const auto& r_e : rEntities) {
            if (!run.empty() && (run.front()->Type() != r_e.Type() ||
                                 run.front()->NumberOfNodes() != r_e.NumberOfNodes()))
                flush();
            run.push_back(&r_e);
        }
        flush();
    }

    void RestoreCellData(const std::vector<std::string>& rNames, BlockRecord::Kind kind,
                         const ModelPart& rMp, NativeMesh& rOut) const {
        for (const auto& r_name : rNames) {
            const NDArray& col = kind == BlockRecord::Kind::Element
                                     ? rMp.GetElementalData(r_name)
                                     : rMp.GetConditionalData(r_name);
            const std::size_t ncols = col.Ndim() >= 2 ? col.Shape()[1] : 1;
            const std::size_t isz = dtype_size(col.Dtype());
            std::size_t row = 0;
            for (const auto& rec : mRecords) {
                if (rec.mKind != kind)
                    continue;
                std::vector<std::size_t> shape = col.Shape();
                if (shape.empty())
                    shape = {0};
                shape[0] = rec.mCount;
                NDArray slice = NDArray::Uninit(col.Dtype(), shape);
                std::memcpy(slice.Data(), col.Data() + row * ncols * isz, rec.mCount * ncols * isz);
                row += rec.mCount;
                rOut.AppendCellData(r_name, std::move(slice));
            }
        }
    }

    mutable NativeMesh mStage;  // staging + writer-serving storage
    std::unique_ptr<ModelPart> mpRoot;
    mutable std::vector<BlockRecord> mRecords;
    bool mMaterialized = false;
    mutable bool mStale = false;
    bool mTagsToSubModelParts = true;
};

}  // namespace meshioplusplus
