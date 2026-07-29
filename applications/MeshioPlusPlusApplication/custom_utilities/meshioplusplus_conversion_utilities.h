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

// External includes
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/kratos_bridge.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/properties.hpp"
#include "meshioplusplus/formats/xdmf_time_series.hpp"

// Project includes
#include "includes/model_part.h"
#include "geometries/geometry_data.h"

namespace Kratos::Internals
{
///@name Type Definitions
///@{

/// A named data array (row-major, rows x components) - shared with @ref MeshioPlusPlusIO,
/// which needs the same shape for its own transient writer.
using DataArray = meshioplusplus::XdmfTimeSeriesWriter::NamedArray;

/**
 * @brief Which nodal/elemental/conditional variables, flags and ids to stage as mesh
 * point/cell data - the vocabulary @ref MeshioPlusPlusIO uses for its own settings, shared so
 * @ref MeshioPlusPlusMeshOperations can carry field data through the operations layer too.
 * @details A name that is not a registered @ref Variable (or, for
 * @ref FieldDataSelection::NodalSolutionStepVariables, not added to the model part's solution
 * step data) is skipped with a warning rather than failing the whole conversion.
 */
struct FieldDataSelection
{
    /// Historical nodal variables (validated against the model part's solution step data).
    std::vector<std::string> NodalSolutionStepVariables;
    /// Non-historical nodal variables.
    std::vector<std::string> NodalDataValueVariables;
    /// Nodal flags, staged as 1 (set) / 0 (unset) / -1 (undefined).
    std::vector<std::string> NodalFlags;
    /// Non-historical elemental variables.
    std::vector<std::string> ElementDataValueVariables;
    /// Elemental flags, staged as 1 / 0 / -1.
    std::vector<std::string> ElementFlags;
    /// Non-historical conditional variables.
    std::vector<std::string> ConditionDataValueVariables;
    /// Conditional flags, staged as 1 / 0 / -1.
    std::vector<std::string> ConditionFlags;
    /// Gauss point variables, averaged over the integration points of each entity. Matches
    /// @ref MeshioPlusPlusIO's "gauss_point_variables_in_elements" setting exactly, including
    /// its scope: applied to *both* elements and conditions, not elements only.
    std::vector<std::string> GaussPointVariables;
    /// Whether to also stage KRATOS_NODE_ID / KRATOS_ELEMENT_ID / KRATOS_CONDITION_ID and
    /// PROPERTIES_ID arrays.
    bool WriteIds = false;
};

///@}
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
 * Multi-component values are recovered from text where the format stored them that way:
 * `mdpa` only turns a single number into a numeric array, so `VOLUME_ACCELERATION 0 0 -9.81`
 * arrives verbatim as text and is tokenized back into a numeric list here.
 *
 * Entries that cannot be represented (a table, a genuinely non-numeric value such as a
 * constitutive law name, or a key that is not a registered variable) are skipped with a
 * warning rather than aborting the read, since the mesh itself is still perfectly usable.
 * @param pProperties The destination properties block.
 * @param rValue The meshio++ key/value pair to assign.
 */
KRATOS_API(KRATOS_MESHIOPLUSPLUS_APPLICATION)
void ApplyMeshioProperty(Properties::Pointer pProperties, const meshioplusplus::PropertyValue& rValue);

/**
 * @brief Exports a Kratos Properties block's material data as meshio++ property values.
 * @details The write-side counterpart of @ref ApplyMeshioProperty, and the reverse of its
 * type dispatch: each entry of the block's variable container is resolved against the
 * registered `double`, `int`, `bool`, `Vector` and `array_1d<double, 3>` variables and
 * appended as a key/value pair. Without this the entities would carry only their properties
 * *id*, and a write to a format that stores material data (`mdpa`) would emit empty blocks.
 *
 * Variables of a type meshio++ has no representation for are skipped with a warning; the
 * mesh is still written. Keep this in step with @ref ApplyMeshioProperty - the two dispatch
 * lists are deliberately adjacent so they cannot drift apart.
 * @param rProperties The source properties block.
 * @param rDestination The meshio++ property set to append to.
 */
KRATOS_API(KRATOS_MESHIOPLUSPLUS_APPLICATION)
void ExportMeshioProperties(const Properties& rProperties, meshioplusplus::PropertySet& rDestination);

/**
 * @brief Populates a meshio++ model part from a Kratos one (one bulk O(n) pass).
 * @details meshio++'s own `from_model_part` cannot be used here: its sub-model-part copy is
 * gated on methods `Kratos::ModelPart` spells differently (`SubModelPartNames` versus
 * `GetSubModelPartNames`) and would silently drop them. The entity mapping still goes
 * through @ref meshioplusplus::bridge_traits<Kratos::ModelPart>, and the Kratos registration
 * names are carried across so a round trip does not degrade them. Material data is carried
 * too, via @ref ExportMeshioProperties.
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
 * @brief Collects the selected nodal data as flat point-data arrays (node container order).
 * @param rSource The model part to read from.
 * @param rSelection Which variables/flags/ids to collect; see @ref FieldDataSelection.
 * @return One array per collected name; a name that could not be resolved is skipped (a
 * warning is logged) rather than aborting the collection.
 */
KRATOS_API(KRATOS_MESHIOPLUSPLUS_APPLICATION)
std::vector<DataArray> CollectPointData(const ModelPart& rSource, const FieldDataSelection& rSelection);

/**
 * @brief Collects the selected elemental/conditional data as flat cell-data arrays.
 * @details Element rows first, then condition rows, zero-filled for the entity kind an array
 * does not apply to - the layout @ref ModelPartToMeshWithData splits back into the per-kind
 * containers meshio++ restores per block.
 * @param rSource The model part to read from.
 * @param WriteElements Whether elements are considered.
 * @param WriteConditions Whether conditions are considered.
 * @param rSelection Which variables/flags/ids to collect; see @ref FieldDataSelection.
 * @return One combined array per collected name.
 */
KRATOS_API(KRATOS_MESHIOPLUSPLUS_APPLICATION)
std::vector<DataArray> CollectCellData(
    const ModelPart& rSource,
    const bool WriteElements,
    const bool WriteConditions,
    const FieldDataSelection& rSelection
    );

/**
 * @brief Builds a staged meshio++ mesh from a Kratos model part, carrying field data too.
 * @details The data-carrying counterpart of @ref ModelPartToMesh: builds the mesh the same
 * way, then stages the arrays @ref CollectPointData / @ref CollectCellData collect as the
 * mesh's point_data / cell_data (the same "set nodal/elemental/conditional data, then
 * InvalidateBlocks()" sequence @ref MeshioPlusPlusIO uses for its own transient writes).
 * @param rSource The Kratos model part to convert.
 * @param WriteElements Whether elements are transferred.
 * @param WriteConditions Whether conditions are transferred.
 * @param WriteDeformedConfiguration True uses the current coordinates, false the initial ones.
 * @param rSelection Which variables/flags/ids to carry; see @ref FieldDataSelection.
 * @return The equivalent meshio++ mesh, data included.
 */
KRATOS_API(KRATOS_MESHIOPLUSPLUS_APPLICATION)
meshioplusplus::Mesh ModelPartToMeshWithData(
    const ModelPart& rSource,
    const bool WriteElements,
    const bool WriteConditions,
    const bool WriteDeformedConfiguration,
    const FieldDataSelection& rSelection
    );

/**
 * @brief Fills a Kratos model part from a meshio++ mesh.
 * @details Creates the properties blocks referenced by the entities, transfers the material
 * data through @ref ApplyMeshioProperty, restores sub model parts from the named regions, and
 * carries the mesh's point_data / cell_data back as non-historical nodal / elemental /
 * conditional @ref Variable data: an array survives only when it is `Float64` or `Int64` and its
 * name matches a registered variable whose component count agrees (the same
 * `double`/`int`/`bool`/`array_1d`/`Vector` dispatch @ref ApplyMeshioProperty uses); anything
 * else is skipped with a warning, because a Kratos entity can only hold data under a name that
 * is an actual registered `Variable` - an operation's own invented array names
 * (`attach_quality`'s `"quality:scaled_jacobian"`, `refine`'s `"refine:level"`, say) never carry
 * through, and that is by design rather than a gap.
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
