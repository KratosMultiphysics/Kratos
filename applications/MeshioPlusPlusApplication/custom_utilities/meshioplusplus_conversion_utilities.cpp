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

// System includes
#include <map>
#include <mutex>
#include <sstream>
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

/***********************************************************************************/
/***********************************************************************************/

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

/***********************************************************************************/
/***********************************************************************************/

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

    // Material data. The entity loops above register the properties *ids* they reference;
    // this carries the values, so a format that stores them (mdpa) writes real blocks
    // instead of empty ones. CreateNewProperties is find-or-create, so a block already
    // registered by an entity is filled in place rather than duplicated.
    for (const auto& r_properties : rSource.GetMesh().Properties()) {
        ExportMeshioProperties(r_properties,
                               rDestination.CreateNewProperties(r_properties.Id()));
    }

    CopySubModelParts(rSource, rDestination, WriteElements, WriteConditions);
}

/***********************************************************************************/
/***********************************************************************************/

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

/***********************************************************************************/
/***********************************************************************************/

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

/***********************************************************************************/
/***********************************************************************************/

void ApplyMeshioProperty(Properties::Pointer pProperties, const mio::PropertyValue& rValue)
{
    KRATOS_TRY

    // A table is a curve, not a value: it needs the application's own registry and is left
    // to the caller rather than guessed at here.
    if (rValue.mIsTable) {
        KRATOS_WARNING("MeshioPlusPlusIO")
            << "Property \"" << rValue.mKey << "\" is a table"
            << " and was not assigned; set it from the materials file instead." << std::endl;
        return;
    }

    // The values, whichever way the format stored them. mdpa only turns a *single* number
    // into a numeric array; a multi-component value such as "0 0 -9.81" arrives verbatim as
    // text, so it is tokenized here. Anything that is not a whitespace-separated list of
    // numbers - a constitutive law name, say - stays unparsed and is skipped below.
    std::vector<double> parsed_values;
    if (rValue.IsText()) {
        std::istringstream stream(rValue.mText);
        double token = 0.0;
        while (stream >> token) {
            parsed_values.push_back(token);
        }
        if (parsed_values.empty() || !stream.eof()) {
            KRATOS_WARNING("MeshioPlusPlusIO")
                << "Property \"" << rValue.mKey << "\" is not numeric"
                << " and was not assigned; set it from the materials file instead." << std::endl;
            return;
        }
    } else {
        const double* p_raw = rValue.mValues.As<double>();
        parsed_values.assign(p_raw, p_raw + rValue.mValues.Size());
    }

    const double* p_values = parsed_values.data();
    const std::size_t number_of_values = parsed_values.size();

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

/***********************************************************************************/
/***********************************************************************************/

void ExportMeshioProperties(const Properties& rProperties, mio::PropertySet& rDestination)
{
    KRATOS_TRY

    // Builds one key/value pair; the NDArray is always Float64, which is what meshio++
    // carries and what ApplyMeshioProperty reads back.
    auto append = [&rDestination](const std::string& rKey, const std::vector<double>& rValues) {
        mio::PropertyValue value;
        value.mKey = rKey;
        value.mValues = mio::NDArray::Uninit(mio::DType::Float64, {rValues.size()});
        std::copy(rValues.begin(), rValues.end(), value.mValues.As<double>());
        rDestination.mValues.push_back(std::move(value));
    };

    for (const auto& r_item : rProperties.GetData()) {
        const std::string variable_name = r_item.first->Name();

        // The reverse of ApplyMeshioProperty's dispatch, in the same order. The value is
        // read through the typed variable rather than the container's void*, so no cast is
        // involved and the component variables resolve exactly as they were registered.
        if (KratosComponents<Variable<double>>::Has(variable_name)) {
            append(variable_name, {rProperties.GetValue(KratosComponents<Variable<double>>::Get(variable_name))});
        } else if (KratosComponents<Variable<int>>::Has(variable_name)) {
            append(variable_name, {static_cast<double>(
                rProperties.GetValue(KratosComponents<Variable<int>>::Get(variable_name)))});
        } else if (KratosComponents<Variable<bool>>::Has(variable_name)) {
            append(variable_name, {rProperties.GetValue(
                KratosComponents<Variable<bool>>::Get(variable_name)) ? 1.0 : 0.0});
        } else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(variable_name)) {
            const auto& r_value = rProperties.GetValue(
                KratosComponents<Variable<array_1d<double, 3>>>::Get(variable_name));
            append(variable_name, {r_value[0], r_value[1], r_value[2]});
        } else if (KratosComponents<Variable<Vector>>::Has(variable_name)) {
            const auto& r_value = rProperties.GetValue(
                KratosComponents<Variable<Vector>>::Get(variable_name));
            append(variable_name, std::vector<double>(r_value.begin(), r_value.end()));
        } else {
            // Matrix, string, table-valued and application-specific types have no meshio++
            // representation; the geometry is still written.
            KRATOS_WARNING_ONCE("MeshioPlusPlusIO")
                << "Property \"" << variable_name << "\" of properties block "
                << rProperties.Id() << " has a type meshio++ cannot carry and was not written."
                << std::endl;
        }
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

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
