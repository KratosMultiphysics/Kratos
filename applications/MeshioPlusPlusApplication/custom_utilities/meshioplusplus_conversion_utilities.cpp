//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes
#include <map>
#include <mutex>
#include <typeindex>
#include <utility>

// External includes

// Project includes
#include "custom_utilities/meshioplusplus_conversion_utilities.h"
#include "utilities/compare_elements_and_conditions_utility.h"
#include "utilities/geometry_utilities.h"
#include "includes/kratos_components.h"
#include "includes/variables.h"

namespace mio = meshioplusplus;

namespace Kratos::Internals
{
namespace
{

/// Key of the entity-name cache: exactly the pair GeometricalObject::IsSame compares.
using EntityKindKey = std::pair<std::type_index, GeometryData::KratosGeometryType>;

/**
 * @brief Memoized wrapper around CompareElementsAndConditionsUtility::GetRegisteredName.
 * @details The underlying lookup is a linear scan over every registered component, so it is
 * resolved once per distinct entity kind and reused afterwards. An unregistered entity is
 * cached as the empty string, which meshio++ reads as "derive the name from the cell type";
 * that keeps a mesh of ad-hoc entities writable instead of aborting the whole write.
 */
template <class TEntity>
const std::string& CachedRegisteredName(const TEntity& rEntity)
{
    // The IO is serial, but the cache is shared across model parts and the surrounding
    // Kratos process may be running inside an OpenMP region, so guard it.
    static std::map<EntityKindKey, std::string> cache;
    static std::mutex mutex;

    const EntityKindKey key{std::type_index(typeid(rEntity)), rEntity.GetGeometry().GetGeometryType()};

    const std::lock_guard<std::mutex> lock(mutex);
    const auto it = cache.find(key);
    if (it != cache.end()) {
        return it->second;
    }

    std::string name;
    try {
        CompareElementsAndConditionsUtility::GetRegisteredName(rEntity, name);
    } catch (const std::exception&) {
        // Not registered: fall back to the cell-type-derived name rather than failing.
        name.clear();
    }
    return cache.emplace(key, std::move(name)).first->second;
}

} // namespace

const std::string& RegisteredEntityName(const Element& rElement)
{
    return CachedRegisteredName(rElement);
}

const std::string& RegisteredEntityName(const Condition& rCondition)
{
    return CachedRegisteredName(rCondition);
}

namespace
{

/**
 * @brief Recursively copies the sub model part structure (names + memberships).
 */
void CopySubModelParts(
    const ModelPart& rSource,
    mio::ModelPart& rDestination,
    const bool WriteElements,
    const bool WriteConditions
    )
{
    for (const auto& r_name : rSource.GetSubModelPartNames()) {
        const ModelPart& r_source_smp = rSource.GetSubModelPart(r_name);
        mio::ModelPart& r_destination_smp = rDestination.CreateSubModelPart(r_name);

        std::vector<std::size_t> ids;
        ids.reserve(r_source_smp.NumberOfNodes());
        for (const auto& r_node : r_source_smp.Nodes()) {
            ids.push_back(r_node.Id());
        }
        r_destination_smp.AddNodes(ids);

        if (WriteElements) {
            ids.clear();
            ids.reserve(r_source_smp.NumberOfElements());
            for (const auto& r_element : r_source_smp.Elements()) {
                ids.push_back(r_element.Id());
            }
            r_destination_smp.AddElements(ids);
        }

        if (WriteConditions) {
            ids.clear();
            ids.reserve(r_source_smp.NumberOfConditions());
            for (const auto& r_condition : r_source_smp.Conditions()) {
                ids.push_back(r_condition.Id());
            }
            r_destination_smp.AddConditions(ids);
        }

        CopySubModelParts(r_source_smp, r_destination_smp, WriteElements, WriteConditions);
    }
}

} // namespace

/**
 * @brief Populates a meshio++ model part from a Kratos one (one O(n) pass).
 * @note meshio++'s from_model_part() cannot be used here: its sub-model-part
 * copy is gated on methods Kratos::ModelPart spells differently
 * (SubModelPartNames vs GetSubModelPartNames) and would silently drop them.
 * The entity mapping still goes through bridge_traits<Kratos::ModelPart>.
 */
void FillMeshioModelPart(
    const ModelPart& rSource,
    mio::ModelPart& rDestination,
    const bool WriteElements,
    const bool WriteConditions,
    const bool WriteDeformedConfiguration
    )
{
    using BridgeTraits = mio::bridge_traits<ModelPart>;

    if (WriteDeformedConfiguration) {
        for (const auto& r_node : rSource.Nodes()) {
            rDestination.CreateNewNode(BridgeTraits::IdOf(r_node), BridgeTraits::XOf(r_node),
                                       BridgeTraits::YOf(r_node), BridgeTraits::ZOf(r_node));
        }
    } else {
        for (const auto& r_node : rSource.Nodes()) {
            rDestination.CreateNewNode(BridgeTraits::IdOf(r_node), r_node.X0(), r_node.Y0(), r_node.Z0());
        }
    }
    // The five-argument overloads also carry the Kratos registration name, so that a
    // concrete "SmallDisplacementElement3D4N" survives instead of degrading to the generic
    // "Element3D4N" that the cell type alone would produce.
    if (WriteElements) {
        for (const auto& r_element : rSource.Elements()) {
            rDestination.CreateNewElement(BridgeTraits::TypeOf(r_element), BridgeTraits::IdOf(r_element),
                                          BridgeTraits::ConnectivityOf(r_element),
                                          BridgeTraits::PropertiesIdOf(r_element),
                                          BridgeTraits::NameOf(r_element));
        }
    }
    if (WriteConditions) {
        for (const auto& r_condition : rSource.Conditions()) {
            rDestination.CreateNewCondition(BridgeTraits::TypeOf(r_condition), BridgeTraits::IdOf(r_condition),
                                            BridgeTraits::ConnectivityOf(r_condition),
                                            BridgeTraits::PropertiesIdOf(r_condition),
                                            BridgeTraits::NameOf(r_condition));
        }
    }

    CopySubModelParts(rSource, rDestination, WriteElements, WriteConditions);
}

meshioplusplus::Mesh ModelPartToMesh(
    const ModelPart& rSource,
    const bool WriteElements,
    const bool WriteConditions,
    const bool WriteDeformedConfiguration
    )
{
    KRATOS_TRY

    mio::Mesh mesh;
    FillMeshioModelPart(rSource, mesh.GetModelPart(), WriteElements, WriteConditions,
                        WriteDeformedConfiguration);
    // The model part view was mutated directly: rebuild the point/cell staging that the
    // operations and the format writers read from.
    mesh.InvalidateBlocks();
    return mesh;

    KRATOS_CATCH("")
}

void MeshToModelPart(meshioplusplus::Mesh& rSource, ModelPart& rDestination)
{
    KRATOS_TRY

    mio::to_model_part(
        rSource.GetModelPart(), rDestination,
        [&rDestination](std::size_t PropertiesId) {
            return rDestination.HasProperties(PropertiesId)
                ? rDestination.pGetProperties(PropertiesId)
                : rDestination.CreateNewProperties(PropertiesId);
        },
        [](Properties::Pointer pProperties, const mio::PropertyValue& rValue) {
            ApplyMeshioProperty(pProperties, rValue);
        });

    KRATOS_CATCH("")
}

void ApplyMeshioProperty(Properties::Pointer pProperties, const mio::PropertyValue& rValue)
{
    KRATOS_TRY

    // A table is a curve, not a value; a text entry is something like a constitutive law
    // name, which needs the application's own registry to instantiate. Both are left to
    // the caller rather than guessed at here.
    if (rValue.mIsTable || rValue.IsText()) {
        KRATOS_WARNING("MeshioPlusPlusIO")
            << "Property \"" << rValue.mKey << "\" is "
            << (rValue.mIsTable ? "a table" : "not numeric")
            << " and was not assigned; set it from the materials file instead." << std::endl;
        return;
    }

    const double* p_values = rValue.mValues.As<double>();
    const std::size_t number_of_values = rValue.mValues.Size();

    if (KratosComponents<Variable<double>>::Has(rValue.mKey)) {
        pProperties->SetValue(KratosComponents<Variable<double>>::Get(rValue.mKey), p_values[0]);
    } else if (KratosComponents<Variable<int>>::Has(rValue.mKey)) {
        pProperties->SetValue(KratosComponents<Variable<int>>::Get(rValue.mKey),
                              static_cast<int>(p_values[0]));
    } else if (KratosComponents<Variable<bool>>::Has(rValue.mKey)) {
        pProperties->SetValue(KratosComponents<Variable<bool>>::Get(rValue.mKey), p_values[0] != 0.0);
    } else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(rValue.mKey)) {
        array_1d<double, 3> value = ZeroVector(3);
        for (std::size_t i = 0; i < std::min<std::size_t>(3, number_of_values); ++i) {
            value[i] = p_values[i];
        }
        pProperties->SetValue(KratosComponents<Variable<array_1d<double, 3>>>::Get(rValue.mKey), value);
    } else if (KratosComponents<Variable<Vector>>::Has(rValue.mKey)) {
        Vector value(number_of_values);
        for (std::size_t i = 0; i < number_of_values; ++i) {
            value[i] = p_values[i];
        }
        pProperties->SetValue(KratosComponents<Variable<Vector>>::Get(rValue.mKey), value);
    } else {
        KRATOS_WARNING("MeshioPlusPlusIO")
            << "Property \"" << rValue.mKey << "\" is not a registered variable and was skipped."
            << std::endl;
    }

    KRATOS_CATCH("")
}

mio::CellType MeshioCellTypeFromKratosGeometry(const GeometryData::KratosGeometryType GeometryType)
{
    using KGT = GeometryData::KratosGeometryType;
    switch (GeometryType) {
        case KGT::Kratos_Point2D:
        case KGT::Kratos_Point3D:              return mio::CellType::Vertex;
        case KGT::Kratos_Line2D2:
        case KGT::Kratos_Line3D2:              return mio::CellType::Line;
        case KGT::Kratos_Line2D3:
        case KGT::Kratos_Line3D3:              return mio::CellType::Line3;
        case KGT::Kratos_Line2D4:              return mio::CellType::Line4;
        case KGT::Kratos_Line2D5:              return mio::CellType::Line5;
        case KGT::Kratos_Triangle2D3:
        case KGT::Kratos_Triangle3D3:          return mio::CellType::Triangle;
        case KGT::Kratos_Triangle2D6:
        case KGT::Kratos_Triangle3D6:          return mio::CellType::Triangle6;
        case KGT::Kratos_Triangle2D10:         return mio::CellType::Triangle10;
        case KGT::Kratos_Triangle2D15:         return mio::CellType::Triangle15;
        case KGT::Kratos_Quadrilateral2D4:
        case KGT::Kratos_Quadrilateral3D4:     return mio::CellType::Quad;
        case KGT::Kratos_Quadrilateral2D8:
        case KGT::Kratos_Quadrilateral3D8:     return mio::CellType::Quad8;
        case KGT::Kratos_Quadrilateral2D9:
        case KGT::Kratos_Quadrilateral3D9:     return mio::CellType::Quad9;
        case KGT::Kratos_Tetrahedra3D4:        return mio::CellType::Tetra;
        case KGT::Kratos_Tetrahedra3D10:       return mio::CellType::Tetra10;
        case KGT::Kratos_Prism3D6:             return mio::CellType::Wedge;
        case KGT::Kratos_Prism3D15:            return mio::CellType::Wedge15;
        case KGT::Kratos_Pyramid3D5:           return mio::CellType::Pyramid;
        case KGT::Kratos_Pyramid3D13:          return mio::CellType::Pyramid13;
        case KGT::Kratos_Hexahedra3D8:         return mio::CellType::Hexahedron;
        case KGT::Kratos_Hexahedra3D20:        return mio::CellType::Hexahedron20;
        case KGT::Kratos_Hexahedra3D27:        return mio::CellType::Hexahedron27;
        default:
            KRATOS_ERROR << "Geometry type " << GeometryUtils::GetGeometryName(GeometryType)
                         << " is not supported by MeshioPlusPlusIO" << std::endl;
    }
}

} // namespace Kratos::Internals
