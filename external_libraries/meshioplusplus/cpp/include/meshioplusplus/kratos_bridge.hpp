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
 * @file kratos_bridge.hpp
 * @brief Header-only, templated bridge between `meshioplusplus::ModelPart`
 * and any Kratos-like model part class — including the real
 * `Kratos::ModelPart` — with no Kratos build dependency.
 *
 * `to_model_part` populates a destination through nothing but the narrow
 * Kratos creation API (`CreateNewNode(id, x, y, z)`,
 * `CreateNewElement(name, id, node_ids, properties)`,
 * `CreateNewCondition(...)`, and — when the destination supports it —
 * `CreateSubModelPart(name)` / `AddNodes` / `AddElements` /
 * `AddConditions`), so the destination can be:
 *
 *  - a real `Kratos::ModelPart` — pass a properties getter that maps a
 *    properties id to a `Properties::Pointer`:
 *    @code
 *    meshioplusplus::to_model_part(source, kratos_mp, [&](auto pid) {
 *        return kratos_mp.HasProperties(pid) ? kratos_mp.pGetProperties(pid)
 *                                            : kratos_mp.CreateNewProperties(pid);
 *    });
 *    @endcode
 *  - a `CoSimIO`-style or mock model part (the overload without a getter
 *    forwards the raw properties id, matching `meshioplusplus::ModelPart`'s
 *    own signature).
 *
 * `from_model_part` walks a Kratos-like source duck-typed through
 * `bridge_traits`, whose default expects `meshioplusplus::ModelPart`'s
 * accessor shape (`Nodes()`/`Elements()`/`Conditions()` ranges of entities
 * with `Id()`/`X()`/`NodeIds()`...). For classes with a different surface
 * (real Kratos exposes connectivity via `GetGeometry()`), specialize
 * `bridge_traits<YourModelPart>` — every customization point is a static
 * function, so a specialization only overrides what differs.
 *
 * Costless in the Kratos sense: conversion is one O(n) bulk-create pass —
 * the same cost Kratos's own CoSimIO conversion utilities pay — because
 * Kratos's pointer-based entity storage cannot be aliased from outside.
 *
 * Backend-independent: usable from any `MESHIOPLUSPLUS_MESH_BACKEND` build
 * (it only needs `model_part.hpp`, never `mesh.hpp`).
 */

// System includes
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/backends/kratos_names.hpp"
#include "meshioplusplus/backends/model_part.hpp"

namespace meshioplusplus {

/**
 * @brief Customization point for `from_model_part`: how to read entities out
 * of a Kratos-like source class. The primary template matches
 * `meshioplusplus::ModelPart`'s own accessor shape; specialize for classes
 * with different spellings (e.g. real Kratos's `GetGeometry()`).
 * @tparam TModelPart The source model part class.
 */
template <class TModelPart>
struct bridge_traits {
    template <class TEntity>
    static IndexType IdOf(const TEntity& rEntity) {
        return static_cast<IndexType>(rEntity.Id());
    }
    template <class TNode>
    static double XOf(const TNode& rNode) {
        return rNode.X();
    }
    template <class TNode>
    static double YOf(const TNode& rNode) {
        return rNode.Y();
    }
    template <class TNode>
    static double ZOf(const TNode& rNode) {
        return rNode.Z();
    }
    /** @brief Connectivity as 1-based node ids. */
    template <class TEntity>
    static std::vector<IndexType> ConnectivityOf(const TEntity& rEntity) {
        const auto& r_ids = rEntity.NodeIds();
        return std::vector<IndexType>(r_ids.begin(), r_ids.end());
    }
    /** @brief The entity's cell type (real-Kratos specializations map
     * GetGeometry().GetGeometryType()). */
    template <class TEntity>
    static CellType TypeOf(const TEntity& rEntity) {
        return rEntity.Type();
    }
    /** @brief The entity's properties id (0 if the class has none). */
    template <class TEntity>
    static IndexType PropertiesIdOf(const TEntity& rEntity) {
        return rEntity.PropertiesId();
    }
};

namespace detail {

template <class TDestModelPart>
void add_sub_model_part_members(const ModelPart& rSourceSmp, TDestModelPart& rDestSmp) {
    if constexpr (requires(TDestModelPart mp, std::vector<IndexType> ids) { mp.AddNodes(ids); }) {
        rDestSmp.AddNodes(rSourceSmp.NodeIds());
        rDestSmp.AddElements(rSourceSmp.ElementIds());
        rDestSmp.AddConditions(rSourceSmp.ConditionIds());
    }
}

template <class TSourceModelPart>
void copy_sub_model_parts_from(const TSourceModelPart& rSource, ModelPart& rDest) {
    if constexpr (requires(const TSourceModelPart mp) { mp.SubModelPartNames(); }) {
        for (const auto& r_name : rSource.SubModelPartNames()) {
            const auto& r_src_smp = rSource.GetSubModelPart(r_name);
            ModelPart& r_smp = rDest.CreateSubModelPart(r_name);
            r_smp.AddNodes(r_src_smp.NodeIds());
            r_smp.AddElements(r_src_smp.ElementIds());
            r_smp.AddConditions(r_src_smp.ConditionIds());
            copy_sub_model_parts_from(r_src_smp, r_smp);  // nested sub parts
        }
    }
}

template <class TDestModelPart>
void copy_sub_model_parts(const ModelPart& rSource, TDestModelPart& rDest) {
    if constexpr (requires(TDestModelPart mp, std::string name) { mp.CreateSubModelPart(name); }) {
        for (const auto& r_name : rSource.SubModelPartNames()) {
            const ModelPart& r_src_smp = rSource.GetSubModelPart(r_name);
            auto& r_dest_smp = rDest.CreateSubModelPart(r_name);
            add_sub_model_part_members(r_src_smp, r_dest_smp);
            copy_sub_model_parts(r_src_smp, r_dest_smp);  // nested sub parts
        }
    }
}

}  // namespace detail

/**
 * @brief Populate a Kratos-like destination model part from a
 * `meshioplusplus::ModelPart` (one bulk O(n) creation pass).
 *
 * @tparam TModelPart The destination class (real Kratos, CoSimIO-like, ...).
 * @tparam TPropertiesGetter Callable `IndexType -> ` whatever the
 *         destination's `CreateNewElement` takes as its properties argument.
 * @param rSource The source model part (must be a root).
 * @param rDest The destination; expected empty (ids are created verbatim).
 * @param rGetProperties Maps a source properties id to the destination's
 *        properties handle (see the file-level real-Kratos example).
 */
template <class TModelPart, class TPropertiesGetter>
void to_model_part(const ModelPart& rSource, TModelPart& rDest,
                   TPropertiesGetter&& rGetProperties) {
    for (const Node& r_node : rSource.Nodes())
        rDest.CreateNewNode(r_node.Id(), r_node.X(), r_node.Y(), r_node.Z());
    for (const Element& r_elem : rSource.Elements())
        rDest.CreateNewElement(kratos_element_name(r_elem.Type()), r_elem.Id(), r_elem.NodeIds(),
                               rGetProperties(r_elem.PropertiesId()));
    for (const Condition& r_cond : rSource.Conditions())
        rDest.CreateNewCondition(kratos_condition_name(r_cond.Type()), r_cond.Id(),
                                 r_cond.NodeIds(), rGetProperties(r_cond.PropertiesId()));
    detail::copy_sub_model_parts(rSource, rDest);
}

/**
 * @brief `to_model_part` overload forwarding the raw properties id (matches
 * `meshioplusplus::ModelPart`'s own creation signature and integer-taking
 * mocks; real Kratos needs the getter overload).
 */
template <class TModelPart>
void to_model_part(const ModelPart& rSource, TModelPart& rDest) {
    to_model_part(rSource, rDest, [](IndexType propertiesId) { return propertiesId; });
}

/**
 * @brief Build a `meshioplusplus::ModelPart` from a Kratos-like source.
 *
 * Reads through `bridge_traits<TModelPart>` (specialize it for classes whose
 * accessors differ from `meshioplusplus::ModelPart`'s shape). Sub model
 * parts are copied when the source exposes `SubModelPartNames()` /
 * `GetSubModelPart()` / per-part id lists.
 *
 * @tparam TModelPart The source class.
 * @param rSource The source model part.
 * @param rName Name for the resulting root (default "Main").
 * @return A freshly-built `meshioplusplus::ModelPart`.
 */
template <class TModelPart>
ModelPart from_model_part(const TModelPart& rSource, std::string rName = "Main") {
    using Traits = bridge_traits<TModelPart>;
    ModelPart out(std::move(rName));
    for (const auto& r_node : rSource.Nodes())
        out.CreateNewNode(Traits::IdOf(r_node), Traits::XOf(r_node), Traits::YOf(r_node),
                          Traits::ZOf(r_node));
    for (const auto& r_elem : rSource.Elements())
        out.CreateNewElement(Traits::TypeOf(r_elem), Traits::IdOf(r_elem),
                             Traits::ConnectivityOf(r_elem), Traits::PropertiesIdOf(r_elem));
    for (const auto& r_cond : rSource.Conditions())
        out.CreateNewCondition(Traits::TypeOf(r_cond), Traits::IdOf(r_cond),
                               Traits::ConnectivityOf(r_cond), Traits::PropertiesIdOf(r_cond));
    detail::copy_sub_model_parts_from(rSource, out);
    return out;
}

}  // namespace meshioplusplus
