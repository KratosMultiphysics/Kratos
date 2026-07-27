//  ██████   ██████ ██████████  █████████  █████   █████ █████    ███████
// ░░██████ ██████ ░░███░░░░░█ ███░░░░░███░░███   ░░███ ░░███   ███░░░░░███      ███         ███
//  ░███░█████░███  ░███  █ ░ ░███    ░░░  ░███    ░███  ░███  ███     ░░███    ░███        ░███
//  ░███░░███ ░███  ░██████   ░░█████████  ░███████████  ░███ ░███      ░███ ███████████ ███████████
//  ░███ ░░░  ░███  ░███░░█    ░░░░░░░░███ ░███░░░░░███  ░███ ░███      ░███░░░░░███░░░ ░░░░░███░░░
//  ░███      ░███  ░███ ░   █ ███    ░███ ░███    ░███  ░███ ░░███     ███     ░███        ░███
//  █████     █████ ██████████░░█████████  █████   █████ █████ ░░░███████░      ░░░         ░░░
// ░░░░░     ░░░░░ ░░░░░░░░░░  ░░░░░░░░░  ░░░░░   ░░░░░ ░░░░░    ░░░░░░░                        Application
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#pragma once

// System includes
#include <string>
#include <vector>

// External includes
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/kratos_bridge.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/properties.hpp"

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "geometries/geometry_data.h"

namespace Kratos::Internals
{
///@name Kratos Classes
///@{

/**
 * @brief Maps a Kratos geometry type to the equivalent meshio++ cell type.
 * @param GeometryType The Kratos geometry type.
 * @return The matching meshio++ cell type.
 * @throws If the geometry has no meshio++ equivalent.
 */
KRATOS_API(KRATOS_MESHIOPLUSPLUS_APPLICATION)
meshioplusplus::CellType MeshioCellTypeFromKratosGeometry(const GeometryData::KratosGeometryType GeometryType);

/**
 * @brief The name an element is registered under in @ref KratosComponents.
 * @details Wraps @ref CompareElementsAndConditionsUtility::GetRegisteredName, which is a
 * linear scan over every registered component and is therefore far too costly to call once
 * per entity. The result is memoized on the pair actually compared by
 * `GeometricalObject::IsSame` — the C++ type and the geometry type — so the scan runs once
 * per distinct entity kind rather than once per entity.
 * @param rElement The element to name.
 * @return The registered name, or an empty string when the element is not registered (in
 * which case meshio++ derives a name from the cell type instead of the write failing).
 */
KRATOS_API(KRATOS_MESHIOPLUSPLUS_APPLICATION)
const std::string& RegisteredEntityName(const Element& rElement);

/**
 * @brief The name a condition is registered under in @ref KratosComponents.
 * @details See the @ref Element overload.
 * @param rCondition The condition to name.
 * @return The registered name, or an empty string when it is not registered.
 */
KRATOS_API(KRATOS_MESHIOPLUSPLUS_APPLICATION)
const std::string& RegisteredEntityName(const Condition& rCondition);

/**
 * @brief Assigns one meshio++ property value to a Kratos Properties block.
 * @details meshio++ carries material data as untyped key/value pairs, because filling a
 * real @ref Properties needs @ref KratosComponents — Kratos's own variable registry, which
 * meshio++ deliberately does not link. This is the typed assignment it hands back to the
 * consumer: the key is resolved against the registered `double`, `int`, `bool`, `Vector`
 * and `array_1d<double, 3>` variables.
 *
 * Entries that cannot be represented (a table, a non-numeric value such as a constitutive
 * law name, or a key that is not a registered variable) are skipped with a warning rather
 * than aborting the read, since the mesh itself is still perfectly usable.
 * @param pProperties The destination properties block.
 * @param rValue The meshio++ key/value pair to assign.
 */
KRATOS_API(KRATOS_MESHIOPLUSPLUS_APPLICATION)
void ApplyMeshioProperty(Properties::Pointer pProperties, const meshioplusplus::PropertyValue& rValue);

/**
 * @brief Populates a meshio++ model part from a Kratos one (one bulk O(n) pass).
 * @details meshio++'s own `from_model_part` cannot be used here: its sub-model-part copy is
 * gated on methods `Kratos::ModelPart` spells differently (`SubModelPartNames` versus
 * `GetSubModelPartNames`) and would silently drop them. The entity mapping still goes
 * through @ref meshioplusplus::bridge_traits<Kratos::ModelPart>, and the Kratos registration
 * names are carried across so a round trip does not degrade them.
 * @param rSource The Kratos model part to read.
 * @param rDestination The meshio++ model part to fill (expected empty).
 * @param WriteElements Whether elements are transferred.
 * @param WriteConditions Whether conditions are transferred.
 * @param WriteDeformedConfiguration True writes the current coordinates, false the initial ones.
 */
KRATOS_API(KRATOS_MESHIOPLUSPLUS_APPLICATION)
void FillMeshioModelPart(
    const ModelPart& rSource,
    meshioplusplus::ModelPart& rDestination,
    const bool WriteElements,
    const bool WriteConditions,
    const bool WriteDeformedConfiguration
    );

/**
 * @brief Builds a staged meshio++ mesh from a Kratos model part.
 * @details Convenience wrapper over @ref FillMeshioModelPart that also rebuilds the
 * point/cell staging the operations and format writers read from.
 * @param rSource The Kratos model part to convert.
 * @param WriteElements Whether elements are transferred.
 * @param WriteConditions Whether conditions are transferred.
 * @param WriteDeformedConfiguration True uses the current coordinates, false the initial ones.
 * @return The equivalent meshio++ mesh.
 */
KRATOS_API(KRATOS_MESHIOPLUSPLUS_APPLICATION)
meshioplusplus::Mesh ModelPartToMesh(
    const ModelPart& rSource,
    const bool WriteElements = true,
    const bool WriteConditions = true,
    const bool WriteDeformedConfiguration = false
    );

/**
 * @brief Fills a Kratos model part from a meshio++ mesh.
 * @details Creates the properties blocks referenced by the entities, transfers the material
 * data through @ref ApplyMeshioProperty, and restores sub model parts from the named regions.
 * @param rSource The meshio++ mesh to read (mutable: materializing its model part view is lazy).
 * @param rDestination The Kratos model part to fill (expected empty).
 */
KRATOS_API(KRATOS_MESHIOPLUSPLUS_APPLICATION)
void MeshToModelPart(meshioplusplus::Mesh& rSource, ModelPart& rDestination);

///@}

} // namespace Kratos::Internals

namespace meshioplusplus
{

/**
 * @brief Specialization of the meshio++ bridge customization point for the real Kratos
 * ModelPart: the single place mapping Kratos entities to meshio++ ones.
 * @details Used by `meshioplusplus::from_model_part` to read a `Kratos::ModelPart`, whose
 * accessors (`GetGeometry()`, `GetProperties()`) differ from the shape the primary template
 * expects. `to_model_part` needs no specialization; it only uses the narrow Kratos creation
 * API, which the real ModelPart already provides.
 */
template <>
struct bridge_traits<Kratos::ModelPart>
{
    template <class TEntity>
    static IndexType IdOf(const TEntity& rEntity)
    {
        return rEntity.Id();
    }

    // Current coordinates (deformed configuration); the initial ones are read
    // directly from the node when "write_deformed_configuration" is false.
    static double XOf(const Kratos::Node& rNode) { return rNode.X(); }
    static double YOf(const Kratos::Node& rNode) { return rNode.Y(); }
    static double ZOf(const Kratos::Node& rNode) { return rNode.Z(); }

    template <class TEntity>
    static std::vector<IndexType> ConnectivityOf(const TEntity& rEntity)
    {
        const auto& r_geometry = rEntity.GetGeometry();
        std::vector<IndexType> node_ids;
        node_ids.reserve(r_geometry.size());
        for (const auto& r_node : r_geometry) {
            node_ids.push_back(r_node.Id());
        }
        return node_ids;
    }

    template <class TEntity>
    static CellType TypeOf(const TEntity& rEntity)
    {
        return Kratos::Internals::MeshioCellTypeFromKratosGeometry(rEntity.GetGeometry().GetGeometryType());
    }

    template <class TEntity>
    static IndexType PropertiesIdOf(const TEntity& rEntity)
    {
        return rEntity.GetProperties().Id();
    }

    /**
     * @brief The Kratos registration name, so that a concrete
     * "SmallDisplacementElement3D4N" survives a round trip instead of degrading to the
     * generic "Element3D4N" that the cell type alone would produce.
     */
    template <class TEntity>
    static std::string NameOf(const TEntity& rEntity)
    {
        return Kratos::Internals::RegisteredEntityName(rEntity);
    }
};

} // namespace meshioplusplus
