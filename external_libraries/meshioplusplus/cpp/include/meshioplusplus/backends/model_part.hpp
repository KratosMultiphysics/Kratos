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
 * @file model_part.hpp
 * @brief `meshioplusplus::ModelPart`: a standalone, Kratos-Multiphysics-style
 * mesh container (Nodes / Elements / Conditions / SubModelParts).
 *
 * A clean-room implementation of the *semantics* of Kratos's `ModelPart`
 * (and of CoSimIO's simplified one) — no Kratos code is used — so meshes can
 * be exchanged with Kratos at the cost of one bulk `CreateNew*` loop, the
 * same cost Kratos's own CoSimIO bridge pays (see `kratos_bridge.hpp`).
 * Key Kratos conventions preserved:
 *
 *  - Entity **Ids are 1-based** (`Id >= 1`) and unique per entity kind
 *    within the root; duplicate creation throws.
 *  - The **root ModelPart owns all entities**; a *sub* model part is a named
 *    nested view referencing entities by Id. Creating an entity on a sub
 *    part inserts it into the root and records membership in that sub part
 *    and every ancestor (Kratos's upward propagation).
 *  - **Elements vs Conditions**: volume/bulk cells vs boundary cells, each
 *    with its own Id space and a Kratos entity name (e.g. `"Element3D4N"`,
 *    `"SurfaceCondition3D3N"` — see `kratos_names.hpp`).
 *  - **Variable data** is simplified to named per-entity columns (an
 *    `NDArray` row per entity, in container order) rather than Kratos's
 *    full `Variable<T>`/solution-step machinery.
 *
 * Containers follow the CoSimIO `IndexedVector` idea (reimplemented):
 * insertion-ordered contiguous storage plus an `Id -> index` hash map for
 * O(1) lookup. `std::deque` keeps entity references stable across growth.
 *
 * This header is backend-independent (it never includes `mesh.hpp`), so the
 * `ModelPart` type and the templated Kratos bridge are usable from *any*
 * mesh-backend build; the KRATOS backend (`kratos_mesh.hpp`) wraps it
 * behind the uniform mesh API.
 */

// System includes
#include <array>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/backends/kratos_names.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/named_arrays.hpp"
#include "meshioplusplus/ndarray.hpp"

namespace meshioplusplus {

/** @brief Entity id type; ids are 1-based (0 is never a valid id). */
using IndexType = std::size_t;

/** @brief A geometric node: 1-based Id plus 3-D coordinates. */
class Node {
public:
    Node(IndexType id, double x, double y, double z) : mId(id), mX(x), mY(y), mZ(z) {
        if (id < 1)
            throw std::invalid_argument("meshio++ ModelPart: node Id must be >= 1");
    }
    IndexType Id() const { return mId; }
    double X() const { return mX; }
    double Y() const { return mY; }
    double Z() const { return mZ; }
    std::array<double, 3> Coordinates() const { return {mX, mY, mZ}; }

private:
    IndexType mId;
    double mX, mY, mZ;
};

/**
 * @brief Shared shape of `Element` and `Condition`: 1-based Id, cell type,
 * properties id, and connectivity as node Ids (1-based, referencing the
 * root's nodes).
 */
class GeometricalEntity {
public:
    GeometricalEntity(IndexType id, CellType type, std::vector<IndexType> nodeIds,
                      IndexType propertiesId)
        : mId(id), mType(type), mPropertiesId(propertiesId), mNodeIds(std::move(nodeIds)) {
        if (id < 1)
            throw std::invalid_argument("meshio++ ModelPart: entity Id must be >= 1");
        if (mNodeIds.empty())
            throw std::invalid_argument("meshio++ ModelPart: entity needs at least one node");
        const int expected = cell_type_num_nodes(type);
        if (expected > 0 && static_cast<std::size_t>(expected) != mNodeIds.size())
            throw std::invalid_argument(
                "meshio++ ModelPart: entity of type '" + cell_type_name(type) + "' expects " +
                std::to_string(expected) + " nodes, got " + std::to_string(mNodeIds.size()));
    }
    IndexType Id() const { return mId; }
    CellType Type() const { return mType; }
    IndexType PropertiesId() const { return mPropertiesId; }
    const std::vector<IndexType>& NodeIds() const { return mNodeIds; }
    std::size_t NumberOfNodes() const { return mNodeIds.size(); }

private:
    IndexType mId;
    CellType mType;
    IndexType mPropertiesId;
    std::vector<IndexType> mNodeIds;
};

/** @brief A bulk (typically max-dimension) entity. */
class Element : public GeometricalEntity {
    using GeometricalEntity::GeometricalEntity;
};
/** @brief A boundary (typically lower-dimension) entity. */
class Condition : public GeometricalEntity {
    using GeometricalEntity::GeometricalEntity;
};

namespace detail {

/**
 * @brief Insertion-ordered entity store with O(1) lookup by 1-based Id —
 * the CoSimIO `IndexedVector` pattern, reimplemented. `std::deque` storage
 * keeps references stable while the container grows.
 * @tparam TEntity `Node`, `Element`, or `Condition`.
 */
template <class TEntity>
class EntityContainer {
public:
    /** @brief Constructs an entity in place; throws on duplicate Id. */
    template <class... TArgs>
    TEntity& Create(IndexType id, TArgs&&... rArgs) {
        if (mIdIndex.count(id))
            throw std::invalid_argument("meshio++ ModelPart: duplicate entity Id " +
                                        std::to_string(id));
        mData.emplace_back(id, std::forward<TArgs>(rArgs)...);
        mIdIndex.emplace(id, mData.size() - 1);
        return mData.back();
    }
    bool Has(IndexType id) const { return mIdIndex.count(id) > 0; }
    const TEntity& Get(IndexType id) const {
        auto it = mIdIndex.find(id);
        if (it == mIdIndex.end())
            throw std::out_of_range("meshio++ ModelPart: no entity with Id " + std::to_string(id));
        return mData[it->second];
    }
    /** @brief 0-based position of Id in insertion order (for data columns). */
    std::size_t IndexOf(IndexType id) const {
        auto it = mIdIndex.find(id);
        if (it == mIdIndex.end())
            throw std::out_of_range("meshio++ ModelPart: no entity with Id " + std::to_string(id));
        return it->second;
    }
    std::size_t Size() const { return mData.size(); }
    auto begin() const { return mData.begin(); }
    auto end() const { return mData.end(); }
    void Reserve(std::size_t n) { mIdIndex.reserve(n); }

private:
    std::deque<TEntity> mData;  // insertion order == container order
    std::unordered_map<IndexType, std::size_t> mIdIndex;
};

/** @brief Ordered id-membership list with O(1) `Has` (for sub model parts). */
class IdList {
public:
    bool Add(IndexType id) {  // returns false if already present
        if (!mSet.insert(id).second)
            return false;
        mIds.push_back(id);
        return true;
    }
    bool Has(IndexType id) const { return mSet.count(id) > 0; }
    std::size_t Size() const { return mIds.size(); }
    const std::vector<IndexType>& Ids() const { return mIds; }

private:
    std::vector<IndexType> mIds;
    std::unordered_set<IndexType> mSet;
};

}  // namespace detail

/**
 * @brief The Kratos-style mesh container: nodes, elements, conditions,
 * nested sub model parts, and simplified per-entity variable data.
 *
 * See the file-level comment for the semantics. Entity creation via a name
 * string accepts both Kratos entity names (`"Element3D4N"`,
 * `"SurfaceCondition3D3N"`, geometry names like `"Tetrahedra3D4"`) and
 * meshio cell-type names (`"tetra"`) — resolution goes through
 * `kratos_names.hpp`'s tables, forward-declared here and defined there to
 * keep this header self-contained for the bridge.
 */
class ModelPart {
public:
    explicit ModelPart(std::string name = "Main") : mName(std::move(name)) {}

    ModelPart(const ModelPart&) = delete;  // entity graph + parent pointers: move-only
    ModelPart& operator=(const ModelPart&) = delete;
    // Moves must re-point the children's parent pointers at the new address
    // (the children themselves are unique_ptr-owned, so their addresses are
    // stable and only the back-pointers need fixing).
    ModelPart(ModelPart&& rOther) noexcept
        : mName(std::move(rOther.mName)),
          mpParent(rOther.mpParent),
          mNodes(std::move(rOther.mNodes)),
          mElements(std::move(rOther.mElements)),
          mConditions(std::move(rOther.mConditions)),
          mLocalNodeIds(std::move(rOther.mLocalNodeIds)),
          mLocalElementIds(std::move(rOther.mLocalElementIds)),
          mLocalConditionIds(std::move(rOther.mLocalConditionIds)),
          mSubModelParts(std::move(rOther.mSubModelParts)),
          mSubIndex(std::move(rOther.mSubIndex)),
          mNodalData(std::move(rOther.mNodalData)),
          mElementalData(std::move(rOther.mElementalData)),
          mConditionalData(std::move(rOther.mConditionalData)) {
        for (auto& r_p : mSubModelParts)
            r_p->mpParent = this;
    }
    ModelPart& operator=(ModelPart&& rOther) noexcept {
        if (this != &rOther) {
            this->~ModelPart();
            new (this) ModelPart(std::move(rOther));
        }
        return *this;
    }

    const std::string& Name() const { return mName; }
    bool IsSubModelPart() const { return mpParent != nullptr; }
    ModelPart& GetRootModelPart() { return mpParent ? mpParent->GetRootModelPart() : *this; }
    const ModelPart& GetRootModelPart() const {
        return mpParent ? mpParent->GetRootModelPart() : *this;
    }
    /** @brief Dotted path from the root (Kratos `FullName()`). */
    std::string FullName() const { return mpParent ? mpParent->FullName() + "." + mName : mName; }

    // --- creation (Kratos semantics: entities live in the root; creating on
    // --- a sub part records membership here and in every ancestor) ---------

    Node& CreateNewNode(IndexType id, double x, double y, double z) {
        Node& r_node = GetRootModelPart().mNodes.Create(id, x, y, z);
        RecordMembership(&ModelPart::mLocalNodeIds, id);
        return r_node;
    }
    Element& CreateNewElement(const std::string& rKratosName, IndexType id,
                              std::vector<IndexType> nodeIds, IndexType propertiesId = 0) {
        Element& r_elem = GetRootModelPart().mElements.Create(
            id, ResolveEntityType(rKratosName), ValidatedNodeIds(std::move(nodeIds)), propertiesId);
        RecordMembership(&ModelPart::mLocalElementIds, id);
        return r_elem;
    }
    Condition& CreateNewCondition(const std::string& rKratosName, IndexType id,
                                  std::vector<IndexType> nodeIds, IndexType propertiesId = 0) {
        Condition& r_cond = GetRootModelPart().mConditions.Create(
            id, ResolveEntityType(rKratosName), ValidatedNodeIds(std::move(nodeIds)), propertiesId);
        RecordMembership(&ModelPart::mLocalConditionIds, id);
        return r_cond;
    }
    // CellType overloads: the bulk-ingest fast path (no per-entity name
    // resolution; connectivity validation is the caller's responsibility).
    Element& CreateNewElement(CellType type, IndexType id, std::vector<IndexType> nodeIds,
                              IndexType propertiesId = 0) {
        Element& r_elem =
            GetRootModelPart().mElements.Create(id, type, std::move(nodeIds), propertiesId);
        RecordMembership(&ModelPart::mLocalElementIds, id);
        return r_elem;
    }
    Condition& CreateNewCondition(CellType type, IndexType id, std::vector<IndexType> nodeIds,
                                  IndexType propertiesId = 0) {
        Condition& r_cond =
            GetRootModelPart().mConditions.Create(id, type, std::move(nodeIds), propertiesId);
        RecordMembership(&ModelPart::mLocalConditionIds, id);
        return r_cond;
    }

    // --- membership (add existing root entities to a sub part) -------------

    void AddNodes(const std::vector<IndexType>& rIds) {
        AddExisting(&ModelPart::mLocalNodeIds, &ModelPart::mNodes, rIds, "node");
    }
    void AddElements(const std::vector<IndexType>& rIds) {
        AddExisting(&ModelPart::mLocalElementIds, &ModelPart::mElements, rIds, "element");
    }
    void AddConditions(const std::vector<IndexType>& rIds) {
        AddExisting(&ModelPart::mLocalConditionIds, &ModelPart::mConditions, rIds, "condition");
    }

    // --- access -------------------------------------------------------------

    /** @brief The ROOT's node container (all entities live in the root). */
    const detail::EntityContainer<Node>& Nodes() const { return GetRootModelPart().mNodes; }
    const detail::EntityContainer<Element>& Elements() const {
        return GetRootModelPart().mElements;
    }
    const detail::EntityContainer<Condition>& Conditions() const {
        return GetRootModelPart().mConditions;
    }
    /** @brief This part's member ids (root: every id, in container order). */
    std::vector<IndexType> NodeIds() const { return MemberIds(&ModelPart::mLocalNodeIds, mNodes); }
    std::vector<IndexType> ElementIds() const {
        return MemberIds(&ModelPart::mLocalElementIds, mElements);
    }
    std::vector<IndexType> ConditionIds() const {
        return MemberIds(&ModelPart::mLocalConditionIds, mConditions);
    }

    bool HasNode(IndexType id) const { return mpParent ? mLocalNodeIds.Has(id) : mNodes.Has(id); }
    bool HasElement(IndexType id) const {
        return mpParent ? mLocalElementIds.Has(id) : mElements.Has(id);
    }
    bool HasCondition(IndexType id) const {
        return mpParent ? mLocalConditionIds.Has(id) : mConditions.Has(id);
    }
    const Node& GetNode(IndexType id) const { return Nodes().Get(id); }
    const Element& GetElement(IndexType id) const { return Elements().Get(id); }
    const Condition& GetCondition(IndexType id) const { return Conditions().Get(id); }

    std::size_t NumberOfNodes() const { return mpParent ? mLocalNodeIds.Size() : mNodes.Size(); }
    std::size_t NumberOfElements() const {
        return mpParent ? mLocalElementIds.Size() : mElements.Size();
    }
    std::size_t NumberOfConditions() const {
        return mpParent ? mLocalConditionIds.Size() : mConditions.Size();
    }

    // --- sub model parts ----------------------------------------------------

    ModelPart& CreateSubModelPart(const std::string& rName) {
        if (rName.empty() || rName.find('.') != std::string::npos)
            throw std::invalid_argument("meshio++ ModelPart: invalid sub model part name '" +
                                        rName + "'");
        if (HasSubModelPart(rName))
            throw std::invalid_argument("meshio++ ModelPart: sub model part '" + rName +
                                        "' already exists");
        auto p_smp = std::make_unique<ModelPart>(rName);
        p_smp->mpParent = this;
        mSubIndex.emplace(rName, mSubModelParts.size());
        mSubModelParts.push_back(std::move(p_smp));
        return *mSubModelParts.back();
    }
    bool HasSubModelPart(const std::string& rName) const { return mSubIndex.count(rName) > 0; }
    ModelPart& GetSubModelPart(const std::string& rName) {
        auto it = mSubIndex.find(rName);
        if (it == mSubIndex.end())
            throw std::out_of_range("meshio++ ModelPart: no sub model part named '" + rName + "'");
        return *mSubModelParts[it->second];
    }
    const ModelPart& GetSubModelPart(const std::string& rName) const {
        return const_cast<ModelPart*>(this)->GetSubModelPart(rName);
    }
    std::size_t NumberOfSubModelParts() const { return mSubModelParts.size(); }
    std::vector<std::string> SubModelPartNames() const {
        std::vector<std::string> names;
        names.reserve(mSubModelParts.size());
        for (const auto& r_p : mSubModelParts)
            names.push_back(r_p->mName);
        return names;
    }

    // --- simplified variable data (root-level, container order) ------------
    //
    // A "variable" is a named NDArray with one row per entity in the ROOT
    // container's insertion order — a pragmatic stand-in for Kratos's
    // Variable<T>/solution-step machinery, sufficient to round-trip
    // point_data/cell_data. Setting from a sub part forwards to the root.

    void SetNodalData(const std::string& rName, NDArray data) {
        GetRootModelPart().mNodalData.Set(rName, std::move(data));
    }
    void SetElementalData(const std::string& rName, NDArray data) {
        GetRootModelPart().mElementalData.Set(rName, std::move(data));
    }
    void SetConditionalData(const std::string& rName, NDArray data) {
        GetRootModelPart().mConditionalData.Set(rName, std::move(data));
    }
    bool HasNodalData(const std::string& rName) const {
        return GetRootModelPart().mNodalData.Has(rName);
    }
    bool HasElementalData(const std::string& rName) const {
        return GetRootModelPart().mElementalData.Has(rName);
    }
    bool HasConditionalData(const std::string& rName) const {
        return GetRootModelPart().mConditionalData.Has(rName);
    }
    const NDArray& GetNodalData(const std::string& rName) const {
        return GetRootModelPart().mNodalData.Get(rName);
    }
    const NDArray& GetElementalData(const std::string& rName) const {
        return GetRootModelPart().mElementalData.Get(rName);
    }
    const NDArray& GetConditionalData(const std::string& rName) const {
        return GetRootModelPart().mConditionalData.Get(rName);
    }
    std::vector<std::string> NodalDataNames() const {
        return GetRootModelPart().mNodalData.SortedNames();
    }
    std::vector<std::string> ElementalDataNames() const {
        return GetRootModelPart().mElementalData.SortedNames();
    }
    std::vector<std::string> ConditionalDataNames() const {
        return GetRootModelPart().mConditionalData.SortedNames();
    }

    /** @brief Scalar component of a nodal variable, addressed by node Id. */
    double GetNodalValue(const std::string& rName, IndexType nodeId,
                         std::size_t component = 0) const {
        const ModelPart& r_root = GetRootModelPart();
        const NDArray& a = r_root.mNodalData.Get(rName);
        const std::size_t ncomp = a.Ndim() >= 2 ? a.Shape()[1] : 1;
        return a.As<double>()[r_root.mNodes.IndexOf(nodeId) * ncomp + component];
    }

private:
    /** @brief Resolve a Kratos entity/geometry name or meshio name to a CellType. */
    static CellType ResolveEntityType(const std::string& rKratosName) {
        const CellType type = cell_type_from_kratos_name(rKratosName);
        // Custom is only acceptable when the spelling itself carries meaning
        // (variable-node-count meshio names like "polyhedron12"); a name
        // neither table knows and that is not a meshio spelling is an error.
        if (type == CellType::Custom && cell_type_from_name(rKratosName) == CellType::Custom &&
            rKratosName.rfind("polygon", 0) != 0 && rKratosName.rfind("polyhedron", 0) != 0)
            throw std::invalid_argument("meshio++ ModelPart: unknown entity type name '" +
                                        rKratosName + "'");
        return type;
    }

    std::vector<IndexType> ValidatedNodeIds(std::vector<IndexType> ids) {
        const ModelPart& r_root = GetRootModelPart();
        for (IndexType id : ids)
            if (!r_root.mNodes.Has(id))
                throw std::invalid_argument("meshio++ ModelPart: connectivity references " +
                                            std::string("unknown node Id ") + std::to_string(id));
        return ids;
    }
    void RecordMembership(detail::IdList ModelPart::* pList, IndexType id) {
        for (ModelPart* p = this; p->mpParent != nullptr; p = p->mpParent)
            (p->*pList).Add(id);
    }
    template <class TEntity>
    void AddExisting(detail::IdList ModelPart::* pList,
                     detail::EntityContainer<TEntity> ModelPart::* pContainer,
                     const std::vector<IndexType>& rIds, const char* pKind) {
        const ModelPart& r_root = GetRootModelPart();
        for (IndexType id : rIds)
            if (!(r_root.*pContainer).Has(id))
                throw std::invalid_argument("meshio++ ModelPart: cannot add unknown " +
                                            std::string(pKind) + " Id " + std::to_string(id));
        for (IndexType id : rIds)
            for (ModelPart* p = this; p->mpParent != nullptr; p = p->mpParent)
                (p->*pList).Add(id);
    }
    template <class TEntity>
    std::vector<IndexType> MemberIds(const detail::IdList ModelPart::* pList,
                                     const detail::EntityContainer<TEntity>& rRootContainer) const {
        if (mpParent)
            return (this->*pList).Ids();
        std::vector<IndexType> ids;
        ids.reserve(rRootContainer.Size());
        for (const auto& r_e : rRootContainer)
            ids.push_back(r_e.Id());
        return ids;
    }

    std::string mName;
    ModelPart* mpParent = nullptr;

    // Root-only entity storage (empty on sub parts).
    detail::EntityContainer<Node> mNodes;
    detail::EntityContainer<Element> mElements;
    detail::EntityContainer<Condition> mConditions;

    // Sub-part membership (unused on the root).
    detail::IdList mLocalNodeIds, mLocalElementIds, mLocalConditionIds;

    // Nested sub model parts, insertion-ordered with O(1) name lookup.
    std::vector<std::unique_ptr<ModelPart>> mSubModelParts;
    std::unordered_map<std::string, std::size_t> mSubIndex;

    // Simplified variables (root-only).
    detail::NamedArrays mNodalData, mElementalData, mConditionalData;
};

}  // namespace meshioplusplus
