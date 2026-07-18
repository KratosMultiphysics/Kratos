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
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <set>
#include <sstream>
#include <unordered_map>

// External includes
#include "pugixml.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/registry.hpp"
#include "meshioplusplus/kratos_bridge.hpp"
#include "meshioplusplus/detail/xdmf_common.hpp"
#include "meshioplusplus/formats/ansys.hpp"
#include "meshioplusplus/formats/ensight.hpp"
#include "meshioplusplus/formats/flac3d.hpp"
#include "meshioplusplus/formats/gmsh.hpp"
#include "meshioplusplus/formats/ply.hpp"
#include "meshioplusplus/formats/stl.hpp"
#include "meshioplusplus/formats/vtk.hpp"
#include "meshioplusplus/formats/vtp.hpp"
#include "meshioplusplus/formats/vtu.hpp"
#ifdef MESHIOPLUSPLUS_HAS_HDF5
#include "meshioplusplus/detail/hdf5_util.hpp"
#endif

// Project includes
#include "input_output/meshioplusplus_io.h"
#include "includes/kratos_components.h"
#include "includes/model_part.h"
#include "includes/variables.h"
#include "processes/integration_values_extrapolation_to_nodes_process.h"
#include "utilities/geometry_utilities.h"

namespace mio = meshioplusplus;

namespace Kratos::Internals
{

/**
 * @brief Maps a Kratos geometry type to the equivalent meshio++ cell type.
 */
inline mio::CellType MeshioCellTypeFromKratosGeometry(const GeometryData::KratosGeometryType GeometryType)
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

namespace meshioplusplus
{

// Specialization of the meshio++ bridge customization point for the real
// Kratos ModelPart: the single place mapping Kratos entities to meshio++ ones.
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
};

} // namespace meshioplusplus

namespace Kratos
{
namespace
{

using DataArray = XdmfTimeSeriesWriter::DataArray;

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
    if (WriteElements) {
        for (const auto& r_element : rSource.Elements()) {
            rDestination.CreateNewElement(BridgeTraits::TypeOf(r_element), BridgeTraits::IdOf(r_element),
                                          BridgeTraits::ConnectivityOf(r_element),
                                          BridgeTraits::PropertiesIdOf(r_element));
        }
    }
    if (WriteConditions) {
        for (const auto& r_condition : rSource.Conditions()) {
            rDestination.CreateNewCondition(BridgeTraits::TypeOf(r_condition), BridgeTraits::IdOf(r_condition),
                                            BridgeTraits::ConnectivityOf(r_condition),
                                            BridgeTraits::PropertiesIdOf(r_condition));
        }
    }

    CopySubModelParts(rSource, rDestination, WriteElements, WriteConditions);
}

/**
 * @brief Collects the listed variables of a container as flat data arrays.
 * @details Supported types (mirroring VtkOutput): double, int, bool (scalar),
 * array_1d<double, 3/4/6/9> and Vector (multi-component, Vector size taken
 * from the first entity). Unknown or unsupported variables are skipped with a
 * warning. rGetValue(entity, variable) provides the value, rValidate(variable)
 * runs an optional per-variable check (e.g. historical availability).
 */
template <class TContainer, class TGetter, class TValidator>
void CollectVariableDataArrays(
    const std::vector<std::string>& rVariableNames,
    const TContainer& rEntities,
    const std::size_t NumberOfEntities,
    const TGetter& rGetValue,
    const TValidator& rValidate,
    std::vector<DataArray>& rOutput
    )
{
    for (const auto& r_variable_name : rVariableNames) {
        DataArray data;
        data.Name = r_variable_name;

        auto collect_scalar = [&](const auto& rVariable) {
            rValidate(rVariable);
            data.NumberOfComponents = 1;
            data.Values.reserve(NumberOfEntities);
            for (const auto& r_entity : rEntities) {
                data.Values.push_back(static_cast<double>(rGetValue(r_entity, rVariable)));
            }
        };
        auto collect_vector = [&](const auto& rVariable, const std::size_t NumberOfComponents) {
            rValidate(rVariable);
            data.NumberOfComponents = NumberOfComponents;
            data.Values.reserve(NumberOfEntities * NumberOfComponents);
            for (const auto& r_entity : rEntities) {
                const auto& r_value = rGetValue(r_entity, rVariable);
                for (std::size_t i = 0; i < NumberOfComponents; ++i) {
                    data.Values.push_back(i < r_value.size() ? static_cast<double>(r_value[i]) : 0.0);
                }
            }
        };

        if (KratosComponents<Variable<double>>::Has(r_variable_name)) {
            collect_scalar(KratosComponents<Variable<double>>::Get(r_variable_name));
        } else if (KratosComponents<Variable<int>>::Has(r_variable_name)) {
            collect_scalar(KratosComponents<Variable<int>>::Get(r_variable_name));
        } else if (KratosComponents<Variable<bool>>::Has(r_variable_name)) {
            collect_scalar(KratosComponents<Variable<bool>>::Get(r_variable_name));
        } else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(r_variable_name)) {
            collect_vector(KratosComponents<Variable<array_1d<double, 3>>>::Get(r_variable_name), 3);
        } else if (KratosComponents<Variable<array_1d<double, 4>>>::Has(r_variable_name)) {
            collect_vector(KratosComponents<Variable<array_1d<double, 4>>>::Get(r_variable_name), 4);
        } else if (KratosComponents<Variable<array_1d<double, 6>>>::Has(r_variable_name)) {
            collect_vector(KratosComponents<Variable<array_1d<double, 6>>>::Get(r_variable_name), 6);
        } else if (KratosComponents<Variable<array_1d<double, 9>>>::Has(r_variable_name)) {
            collect_vector(KratosComponents<Variable<array_1d<double, 9>>>::Get(r_variable_name), 9);
        } else if (KratosComponents<Variable<Vector>>::Has(r_variable_name)) {
            const auto& r_variable = KratosComponents<Variable<Vector>>::Get(r_variable_name);
            const std::size_t number_of_components =
                NumberOfEntities > 0 ? rGetValue(*rEntities.begin(), r_variable).size() : 0;
            if (number_of_components == 0) {
                KRATOS_WARNING_ONCE("MeshioPlusPlusIO") << "Vector variable \"" << r_variable_name
                    << "\" has no components on the first entity - skipping it" << std::endl;
                continue;
            }
            collect_vector(r_variable, number_of_components);
        } else {
            KRATOS_WARNING_ONCE("MeshioPlusPlusIO") << "Variable \"" << r_variable_name
                << "\" is not registered with a type suitable for MeshioPlusPlusIO - skipping it" << std::endl;
            continue;
        }

        rOutput.push_back(std::move(data));
    }
}

/**
 * @brief Collects the listed flags of a container as 1/0/-1 scalar arrays
 * (VtkOutput convention: -1 when the flag is not defined on the entity).
 */
template <class TContainer>
void CollectFlagDataArrays(
    const std::vector<std::string>& rFlagNames,
    const TContainer& rEntities,
    const std::size_t NumberOfEntities,
    std::vector<DataArray>& rOutput
    )
{
    for (const auto& r_flag_name : rFlagNames) {
        if (!KratosComponents<Flags>::Has(r_flag_name)) {
            KRATOS_WARNING_ONCE("MeshioPlusPlusIO") << "Flag \"" << r_flag_name
                << "\" is not registered - skipping it" << std::endl;
            continue;
        }
        const Flags& r_flag = KratosComponents<Flags>::Get(r_flag_name);

        DataArray data;
        data.Name = r_flag_name;
        data.NumberOfComponents = 1;
        data.Values.reserve(NumberOfEntities);
        for (const auto& r_entity : rEntities) {
            data.Values.push_back(r_entity.IsDefined(r_flag) ? (r_entity.Is(r_flag) ? 1.0 : 0.0) : -1.0);
        }
        rOutput.push_back(std::move(data));
    }
}

/**
 * @brief Collects entity ids as a scalar array named rName.
 */
template <class TContainer>
DataArray CollectIdsArray(
    const TContainer& rEntities,
    const std::size_t NumberOfEntities,
    const std::string& rName
    )
{
    DataArray data;
    data.Name = rName;
    data.NumberOfComponents = 1;
    data.Values.reserve(NumberOfEntities);
    for (const auto& r_entity : rEntities) {
        data.Values.push_back(static_cast<double>(r_entity.Id()));
    }
    return data;
}

/**
 * @brief Collects the entities' properties ids as the PROPERTIES_ID array.
 */
template <class TContainer>
DataArray CollectPropertiesIdsArray(
    const TContainer& rEntities,
    const std::size_t NumberOfEntities
    )
{
    DataArray data;
    data.Name = "PROPERTIES_ID";
    data.NumberOfComponents = 1;
    data.Values.reserve(NumberOfEntities);
    for (const auto& r_entity : rEntities) {
        data.Values.push_back(static_cast<double>(r_entity.GetProperties().Id()));
    }
    return data;
}

/**
 * @brief Collects gauss point results averaged over the integration points
 * (VtkOutput convention). Variables whose CalculateOnIntegrationPoints returns
 * nothing (e.g. the generic core entities) are skipped with a warning.
 * @note CalculateOnIntegrationPoints is non-const, hence the non-const container.
 */
template <class TContainer>
void CollectGaussPointDataArrays(
    const std::vector<std::string>& rVariableNames,
    TContainer& rEntities,
    const std::size_t NumberOfEntities,
    const ProcessInfo& rProcessInfo,
    std::vector<DataArray>& rOutput
    )
{
    for (const auto& r_variable_name : rVariableNames) {
        DataArray data;
        data.Name = r_variable_name;

        // Averages the integration point values of one entity into rAppendTo
        auto collect_scalar = [&]<class TDataType>(const Variable<TDataType>& rVariable) -> bool {
            std::vector<TDataType> gp_values;
            if (NumberOfEntities == 0) {
                return false;
            }
            rEntities.begin()->CalculateOnIntegrationPoints(rVariable, gp_values, rProcessInfo);
            if (gp_values.empty()) {
                KRATOS_WARNING_ONCE("MeshioPlusPlusIO") << "Gauss point variable \"" << r_variable_name
                    << "\" returns no integration point values - skipping it" << std::endl;
                return false;
            }
            data.NumberOfComponents = 1;
            data.Values.reserve(NumberOfEntities);
            for (auto& r_entity : rEntities) {
                r_entity.CalculateOnIntegrationPoints(rVariable, gp_values, rProcessInfo);
                double average = 0.0;
                for (const auto& r_value : gp_values) {
                    average += static_cast<double>(r_value);
                }
                data.Values.push_back(gp_values.empty() ? 0.0 : average / gp_values.size());
            }
            return true;
        };
        auto collect_vector = [&]<class TDataType>(const Variable<TDataType>& rVariable) -> bool {
            std::vector<TDataType> gp_values;
            if (NumberOfEntities == 0) {
                return false;
            }
            rEntities.begin()->CalculateOnIntegrationPoints(rVariable, gp_values, rProcessInfo);
            if (gp_values.empty() || gp_values[0].size() == 0) {
                KRATOS_WARNING_ONCE("MeshioPlusPlusIO") << "Gauss point variable \"" << r_variable_name
                    << "\" returns no integration point values - skipping it" << std::endl;
                return false;
            }
            const std::size_t number_of_components = gp_values[0].size();
            data.NumberOfComponents = number_of_components;
            data.Values.reserve(NumberOfEntities * number_of_components);
            std::vector<double> average(number_of_components);
            for (auto& r_entity : rEntities) {
                r_entity.CalculateOnIntegrationPoints(rVariable, gp_values, rProcessInfo);
                std::fill(average.begin(), average.end(), 0.0);
                for (const auto& r_value : gp_values) {
                    for (std::size_t i = 0; i < number_of_components && i < r_value.size(); ++i) {
                        average[i] += r_value[i];
                    }
                }
                for (std::size_t i = 0; i < number_of_components; ++i) {
                    data.Values.push_back(gp_values.empty() ? 0.0 : average[i] / gp_values.size());
                }
            }
            return true;
        };

        bool collected = false;
        if (KratosComponents<Variable<double>>::Has(r_variable_name)) {
            collected = collect_scalar(KratosComponents<Variable<double>>::Get(r_variable_name));
        } else if (KratosComponents<Variable<int>>::Has(r_variable_name)) {
            collected = collect_scalar(KratosComponents<Variable<int>>::Get(r_variable_name));
        } else if (KratosComponents<Variable<bool>>::Has(r_variable_name)) {
            collected = collect_scalar(KratosComponents<Variable<bool>>::Get(r_variable_name));
        } else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(r_variable_name)) {
            collected = collect_vector(KratosComponents<Variable<array_1d<double, 3>>>::Get(r_variable_name));
        } else if (KratosComponents<Variable<array_1d<double, 6>>>::Has(r_variable_name)) {
            collected = collect_vector(KratosComponents<Variable<array_1d<double, 6>>>::Get(r_variable_name));
        } else if (KratosComponents<Variable<Vector>>::Has(r_variable_name)) {
            collected = collect_vector(KratosComponents<Variable<Vector>>::Get(r_variable_name));
        } else {
            KRATOS_WARNING_ONCE("MeshioPlusPlusIO") << "Gauss point variable \"" << r_variable_name
                << "\" is not registered with a type suitable for MeshioPlusPlusIO - skipping it" << std::endl;
        }

        if (collected) {
            rOutput.push_back(std::move(data));
        }
    }
}

/**
 * @brief Merges per-kind cell arrays into combined ones covering the full cell
 * range of the written mesh (element rows first, then condition rows),
 * zero-filling the rows of the entity kind an array does not apply to.
 */
std::vector<DataArray> MergeCellDataParts(
    std::vector<DataArray>&& rElementPart,
    const std::size_t NumberOfElementRows,
    std::vector<DataArray>&& rConditionPart,
    const std::size_t NumberOfConditionRows
    )
{
    std::vector<DataArray> merged;
    merged.reserve(rElementPart.size() + rConditionPart.size());

    for (auto& r_element_data : rElementPart) {
        const auto it_condition = std::find_if(rConditionPart.begin(), rConditionPart.end(),
            [&r_element_data](const DataArray& rData) { return rData.Name == r_element_data.Name; });
        if (it_condition != rConditionPart.end() &&
            it_condition->NumberOfComponents == r_element_data.NumberOfComponents) {
            r_element_data.Values.insert(r_element_data.Values.end(),
                                         it_condition->Values.begin(), it_condition->Values.end());
            rConditionPart.erase(it_condition);
        } else {
            r_element_data.Values.resize(
                r_element_data.Values.size() + NumberOfConditionRows * r_element_data.NumberOfComponents, 0.0);
        }
        merged.push_back(std::move(r_element_data));
    }

    for (auto& r_condition_data : rConditionPart) {
        r_condition_data.Values.insert(r_condition_data.Values.begin(),
                                       NumberOfElementRows * r_condition_data.NumberOfComponents, 0.0);
        merged.push_back(std::move(r_condition_data));
    }

    return merged;
}

/**
 * @brief Builds a meshio++ NDArray ({n} or {n, components} Float64) from a DataArray.
 */
mio::NDArray ToNDArray(const DataArray& rData)
{
    const std::size_t number_of_rows =
        rData.NumberOfComponents > 0 ? rData.Values.size() / rData.NumberOfComponents : 0;
    mio::NDArray array = mio::NDArray::Uninit(
        mio::DType::Float64,
        rData.NumberOfComponents == 1 ? std::vector<std::size_t>{number_of_rows}
                                      : std::vector<std::size_t>{number_of_rows, rData.NumberOfComponents});
    std::copy(rData.Values.begin(), rData.Values.end(), array.As<double>());
    return array;
}

// Format enum <-> canonical meshio++ format name. Keep in sync with the
// Format enum in meshioplusplus_io.h and meshio++'s registry.cpp.
#define KRATOS_MESHIOPLUSPLUS_FORMATS(X) \
    X(ABAQUS, "abaqus")                  \
    X(ANSYS, "ansys")                    \
    X(ANSYSINP, "ansysinp")              \
    X(AVSUCD, "avsucd")                  \
    X(CGNS, "cgns")                      \
    X(DEX, "dex")                        \
    X(DOLFIN, "dolfin")                  \
    X(ENSIGHT, "ensight")                \
    X(EXODUS, "exodus")                  \
    X(FLAC3D, "flac3d")                  \
    X(FLUX, "flux")                      \
    X(FREEFEM, "freefem")                \
    X(GMSH, "gmsh")                      \
    X(H5M, "h5m")                        \
    X(HMF, "hmf")                        \
    X(IP, "ip")                          \
    X(MED, "med")                        \
    X(MEDIT, "medit")                    \
    X(MFF, "mff")                        \
    X(MFM, "mfm")                        \
    X(MPHTXT, "mphtxt")                  \
    X(NASTRAN, "nastran")                \
    X(NETGEN, "netgen")                  \
    X(OBJ, "obj")                        \
    X(OFF, "off")                        \
    X(OPENFOAM, "openfoam")              \
    X(PERMAS, "permas")                  \
    X(PLY, "ply")                        \
    X(STL, "stl")                        \
    X(SU2, "su2")                        \
    X(SVG, "svg")                        \
    X(TECPLOT, "tecplot")                \
    X(TETGEN, "tetgen")                  \
    X(TIKZ, "tikz")                      \
    X(TRIANGLE, "triangle")              \
    X(UGRID, "ugrid")                    \
    X(UNV, "unv")                        \
    X(VTK, "vtk")                        \
    X(VTP, "vtp")                        \
    X(VTU, "vtu")                        \
    X(WKT, "wkt")                        \
    X(XDMF, "xdmf")

const std::unordered_map<std::string, MeshioPlusPlusIO::Format>& GetFormatNameMap()
{
    static const std::unordered_map<std::string, MeshioPlusPlusIO::Format> format_map = {
#define KRATOS_MESHIOPLUSPLUS_FORMAT_LOOKUP(EnumName, Name) {Name, MeshioPlusPlusIO::Format::EnumName},
        KRATOS_MESHIOPLUSPLUS_FORMATS(KRATOS_MESHIOPLUSPLUS_FORMAT_LOOKUP)
#undef KRATOS_MESHIOPLUSPLUS_FORMAT_LOOKUP
    };
    return format_map;
}

/**
 * @brief Writes rMesh honoring an ascii/binary "file_format" override for the
 * formats exposing such a flag, bypassing the meshio++ registry defaults.
 * @param Skin Whether stl/ply extract and write the boundary skin of volume
 *        meshes (the "skin" setting); ignored by the other formats.
 * @return false when the format has no ascii/binary variant (caller falls
 *         back to the registry writer).
 */
bool WriteWithFileFormatOverride(
    const std::string& rFormatName,
    const bool Binary,
    const bool Skin,
    const std::filesystem::path& rPath,
    const mio::Mesh& rMesh
    )
{
    const std::string path = rPath.string();
    if (rFormatName == "vtu") {
        mio::write_vtu(path, rMesh, Binary, /*zlib=*/Binary);
        return true;
    }
    if (rFormatName == "vtp") {
        mio::write_vtp(path, rMesh, Binary, /*zlib=*/Binary);
        return true;
    }
    if (rFormatName == "vtk") {
        mio::write_vtk(path, rMesh, Binary, /*v51=*/true);
        return true;
    }
    if (rFormatName == "gmsh") {
        mio::write_gmsh41(path, rMesh, Binary);
        return true;
    }
    if (rFormatName == "stl") {
        mio::write_stl(path, rMesh, Binary, Skin);
        return true;
    }
    if (rFormatName == "ply") {
        mio::write_ply(path, rMesh, Binary, Skin);
        return true;
    }
    if (rFormatName == "ansys") {
        mio::write_ansys(path, rMesh, Binary);
        return true;
    }
    if (rFormatName == "flac3d") {
        mio::write_flac3d(path, rMesh, ".16e", Binary);
        return true;
    }
    if (rFormatName == "ensight") {
        mio::write_ensight(path, rMesh, Binary);
        return true;
    }
    return false;
}

/**
 * @brief XDMF mixed-topology cell type index (XDMF v3 specification).
 */
int MeshioToXdmfIndex(const std::string& rMeshioType)
{
    static const std::unordered_map<std::string, int> index_map = {
        {"vertex", 0x1},        {"line", 0x2},          {"triangle", 0x4},     {"quad", 0x5},
        {"tetra", 0x6},         {"pyramid", 0x7},       {"wedge", 0x8},        {"hexahedron", 0x9},
        {"line3", 0x22},        {"quad9", 0x23},        {"triangle6", 0x24},   {"quad8", 0x25},
        {"tetra10", 0x26},      {"pyramid13", 0x27},    {"wedge15", 0x28},     {"wedge18", 0x29},
        {"hexahedron20", 0x30}, {"hexahedron24", 0x31}, {"hexahedron27", 0x32}};
    const auto it = index_map.find(rMeshioType);
    KRATOS_ERROR_IF(it == index_map.end())
        << "XDMF: cell type '" << rMeshioType << "' cannot be written into a mixed topology" << std::endl;
    return it->second;
}

std::pair<const char*, const char*> XdmfDTypeOf(const mio::DType DataType)
{
    switch (DataType) {
        case mio::DType::Int8:    return {"Int", "1"};
        case mio::DType::Int16:   return {"Int", "2"};
        case mio::DType::Int32:   return {"Int", "4"};
        case mio::DType::Int64:   return {"Int", "8"};
        case mio::DType::UInt8:   return {"UInt", "1"};
        case mio::DType::UInt16:  return {"UInt", "2"};
        case mio::DType::UInt32:  return {"UInt", "4"};
        case mio::DType::UInt64:  return {"UInt", "8"};
        case mio::DType::Float32: return {"Float", "4"};
        case mio::DType::Float64: return {"Float", "8"};
    }
    return {"Float", "8"};
}

std::string XdmfAttributeTypeOf(const std::vector<std::size_t>& rShape)
{
    if (rShape.size() == 1 || (rShape.size() == 2 && rShape[1] == 1)) return "Scalar";
    if (rShape.size() == 2 && (rShape[1] == 2 || rShape[1] == 3)) return "Vector";
    if (rShape.size() == 2 && rShape[1] == 9) return "Tensor";
    if (rShape.size() == 2 && rShape[1] == 6) return "Tensor6";
    return "Matrix";
}

/**
 * @brief Counts the <DataItem> nodes of a document (recursive, no XPath).
 */
int CountDataItems(const pugi::xml_node Node)
{
    int count = 0;
    for (const pugi::xml_node child : Node.children()) {
        if (std::string(child.name()) == "DataItem") {
            ++count;
        }
        count += CountDataItems(child);
    }
    return count;
}

} // namespace

// Forwards the HDF5 file handle to AddDataItem only when HDF5 is compiled in.
#ifdef MESHIOPLUSPLUS_HAS_HDF5
#define KRATOS_MESHIOPLUSPLUS_H5_ARG(h5_file) , h5_file
#else
#define KRATOS_MESHIOPLUSPLUS_H5_ARG(h5_file)
#endif

/***********************************************************************************/
/***********************************************************************************/

class XdmfTimeSeriesWriter::Impl
{
public:
    Impl(
        const std::filesystem::path& rFileName,
        const std::string& rDataFormat,
        const std::size_t Precision
        ) : mFileName(rFileName),
            mDataFormat(rDataFormat),
            mPrecision(Precision)
    {
        KRATOS_ERROR_IF(mDataFormat != "XML" && mDataFormat != "Binary" && mDataFormat != "HDF")
            << "Unknown XDMF data format \"" << mDataFormat << "\" (use \"XML\", \"Binary\" or \"HDF\")" << std::endl;
#ifndef MESHIOPLUSPLUS_HAS_HDF5
        KRATOS_ERROR_IF(mDataFormat == "HDF")
            << "The \"HDF\" XDMF data format requires an HDF5-enabled build" << std::endl;
#endif
        mBase = mFileName.string();
        const std::size_t dot = mBase.find_last_of('.');
        if (dot != std::string::npos) {
            mBase = mBase.substr(0, dot);
        }
    }

    bool MeshWritten() const
    {
        if (!std::filesystem::exists(mFileName)) {
            return false;
        }
        pugi::xml_document document;
        if (!document.load_file(mFileName.string().c_str())) {
            return false;
        }
        const pugi::xml_node domain = document.child("Xdmf").child("Domain");
        return FindMeshGrid(domain) && FindCollectionGrid(domain);
    }

    void WriteMesh(
        const ModelPart& rModelPart,
        const bool WriteElements,
        const bool WriteConditions,
        const bool WriteDeformedConfiguration
        )
    {
        KRATOS_TRY

        // Materialize the mesh through the meshio++ Kratos backend: filling the
        // ModelPart view and invalidating the blocks makes the staging accessors
        // (Points()/CellRange()) rebuild point/cell blocks from the entities.
        mio::Mesh mesh;
        FillMeshioModelPart(rModelPart, mesh.GetModelPart(), WriteElements, WriteConditions,
                            WriteDeformedConfiguration);
        mesh.InvalidateBlocks();

        pugi::xml_document document;
        pugi::xml_node xdmf = document.append_child("Xdmf");
        xdmf.append_attribute("Version") = "3.0";
        xdmf.append_attribute("xmlns:xi") = "http://www.w3.org/2001/XInclude";
        pugi::xml_node domain = xdmf.append_child("Domain");

        pugi::xml_node grid = domain.append_child("Grid");
        grid.append_attribute("Name") = "mesh";
        grid.append_attribute("GridType") = "Uniform";

        int data_counter = 0;
#ifdef MESHIOPLUSPLUS_HAS_HDF5
        mio::h5::Hid h5_file;
        if (mDataFormat == "HDF") {
            h5_file = mio::h5::create_file(mBase + ".h5"); // truncates a stale companion
        }
#endif

        // Geometry (points are always materialized as n x 3 Float64)
        pugi::xml_node geometry = grid.append_child("Geometry");
        geometry.append_attribute("GeometryType") = "XYZ";
        AddDataItem(geometry, mesh.Points(), data_counter KRATOS_MESHIOPLUSPLUS_H5_ARG(h5_file));

        // Topology: plain for a single cell block, Mixed otherwise
        const std::size_t number_of_blocks = mesh.NumCellBlocks();
        if (number_of_blocks == 1) {
            const auto cell_block = mesh.Cells(0);
            pugi::xml_node topology = grid.append_child("Topology");
            topology.append_attribute("TopologyType") = mio::xdmfcommon::meshio_to_xdmf(cell_block.Type());
            topology.append_attribute("NumberOfElements") = std::to_string(cell_block.NumCells()).c_str();
            topology.append_attribute("NodesPerElement") = std::to_string(cell_block.NodesPerCell()).c_str();
            AddDataItem(topology, cell_block.Conn(), data_counter KRATOS_MESHIOPLUSPLUS_H5_ARG(h5_file));
        } else if (number_of_blocks > 1) {
            std::size_t total_cells = 0;
            std::size_t total_length = 0;
            for (const auto cell_block : mesh.CellRange()) {
                // vertex/line entries carry a second index token (their count)
                const std::size_t prefix = (cell_block.Type() == "vertex" || cell_block.Type() == "line") ? 2 : 1;
                total_cells += cell_block.NumCells();
                total_length += cell_block.NumCells() * (prefix + cell_block.NodesPerCell());
            }
            mio::NDArray mixed(mio::DType::Int64, {total_length});
            std::int64_t* p_mixed = mixed.As<std::int64_t>();
            std::size_t position = 0;
            for (const auto cell_block : mesh.CellRange()) {
                const mio::NDArray& r_connectivities = cell_block.Conn();
                const std::int64_t* p_connectivities = r_connectivities.As<std::int64_t>();
                const std::size_t nodes_per_cell = cell_block.NodesPerCell();
                const int type_index = MeshioToXdmfIndex(cell_block.Type());
                const std::size_t prefix = (cell_block.Type() == "vertex" || cell_block.Type() == "line") ? 2 : 1;
                for (std::size_t i_cell = 0; i_cell < cell_block.NumCells(); ++i_cell) {
                    for (std::size_t i_prefix = 0; i_prefix < prefix; ++i_prefix) {
                        p_mixed[position++] = type_index;
                    }
                    for (std::size_t i_node = 0; i_node < nodes_per_cell; ++i_node) {
                        p_mixed[position++] = p_connectivities[i_cell * nodes_per_cell + i_node];
                    }
                }
            }
            pugi::xml_node topology = grid.append_child("Topology");
            topology.append_attribute("TopologyType") = "Mixed";
            topology.append_attribute("NumberOfElements") = std::to_string(total_cells).c_str();
            AddDataItem(topology, mixed, data_counter KRATOS_MESHIOPLUSPLUS_H5_ARG(h5_file));
        }

        // The (initially empty) temporal collection the steps are appended to
        pugi::xml_node collection = domain.append_child("Grid");
        collection.append_attribute("Name") = "TimeSeries_Kratos";
        collection.append_attribute("GridType") = "Collection";
        collection.append_attribute("CollectionType") = "Temporal";

        KRATOS_ERROR_IF_NOT(document.save_file(mFileName.string().c_str(), "  "))
            << "Could not write XDMF file " << mFileName << std::endl;

        KRATOS_CATCH("")
    }

    void AppendStep(
        const double TimeValue,
        const std::vector<DataArray>& rPointData,
        const std::vector<DataArray>& rCellData
        )
    {
        KRATOS_TRY

        pugi::xml_document document;
        KRATOS_ERROR_IF_NOT(document.load_file(mFileName.string().c_str()))
            << "Could not parse existing XDMF file " << mFileName
            << " - write the mesh first (WriteMesh)" << std::endl;

        pugi::xml_node domain = document.child("Xdmf").child("Domain");
        pugi::xml_node collection = FindCollectionGrid(domain);
        KRATOS_ERROR_IF_NOT(collection)
            << "XDMF file " << mFileName << " has no temporal collection grid to extend" << std::endl;

        // Continue the DataItem numbering so Binary/HDF appends never clash,
        // even when a new writer extends a file from a previous run.
        int data_counter = CountDataItems(document);

#ifdef MESHIOPLUSPLUS_HAS_HDF5
        mio::h5::Hid h5_file;
        if (mDataFormat == "HDF") {
            const std::string h5_path = mBase + ".h5";
            h5_file = mio::h5::Hid(H5Fopen(h5_path.c_str(), H5F_ACC_RDWR, H5P_DEFAULT), H5Fclose);
            KRATOS_ERROR_IF_NOT(h5_file.Valid())
                << "Could not open HDF5 companion file " << h5_path << " for appending" << std::endl;
        }
#endif

        pugi::xml_node grid = collection.append_child("Grid");
        pugi::xml_node include = grid.append_child("xi:include");
        include.append_attribute("xpointer") =
            "xpointer(//Grid[@Name=\"mesh\"]/*[self::Topology or self::Geometry])";
        pugi::xml_node time = grid.append_child("Time");
        time.append_attribute("Value") = TimeValue;

        auto add_attribute = [&](const DataArray& rData, const char* pCenter) {
            if (rData.Values.empty()) {
                return; // nothing to write (e.g. no entities of this kind)
            }
            const mio::NDArray array = ToNDArray(rData);
            pugi::xml_node attribute = grid.append_child("Attribute");
            attribute.append_attribute("Name") = rData.Name.c_str();
            attribute.append_attribute("AttributeType") = XdmfAttributeTypeOf(array.Shape()).c_str();
            attribute.append_attribute("Center") = pCenter;
            AddDataItem(attribute, array, data_counter KRATOS_MESHIOPLUSPLUS_H5_ARG(h5_file));
        };

        for (const auto& r_data : rPointData) {
            add_attribute(r_data, "Node");
        }
        for (const auto& r_data : rCellData) {
            add_attribute(r_data, "Cell");
        }

        KRATOS_ERROR_IF_NOT(document.save_file(mFileName.string().c_str(), "  "))
            << "Could not write XDMF file " << mFileName << std::endl;

        KRATOS_CATCH("")
    }

private:
    std::filesystem::path mFileName;
    std::string mDataFormat;
    std::size_t mPrecision;
    std::string mBase; /// The file path without extension (for .bin/.h5 companions)

    static pugi::xml_node FindMeshGrid(const pugi::xml_node Domain)
    {
        for (const pugi::xml_node grid : Domain.children("Grid")) {
            if (std::string(grid.attribute("GridType").as_string()) == "Uniform") {
                return grid;
            }
        }
        return pugi::xml_node();
    }

    static pugi::xml_node FindCollectionGrid(const pugi::xml_node Domain)
    {
        for (const pugi::xml_node grid : Domain.children("Grid")) {
            if (std::string(grid.attribute("GridType").as_string()) == "Collection" &&
                std::string(grid.attribute("CollectionType").as_string()) == "Temporal") {
                return grid;
            }
        }
        return pugi::xml_node();
    }

    /**
     * @brief Appends a <DataItem> carrying rArray in the configured data format.
     */
    void AddDataItem(
        pugi::xml_node Parent,
        const mio::NDArray& rArray,
        int& rDataCounter
#ifdef MESHIOPLUSPLUS_HAS_HDF5
        , mio::h5::Hid& rH5File
#endif
        ) const
    {
        const auto [dtype_name, precision] = XdmfDTypeOf(rArray.Dtype());
        std::string dimensions;
        for (std::size_t i = 0; i < rArray.Shape().size(); ++i) {
            if (i > 0) dimensions += " ";
            dimensions += std::to_string(rArray.Shape()[i]);
        }

        pugi::xml_node data_item = Parent.append_child("DataItem");
        data_item.append_attribute("DataType") = dtype_name;
        data_item.append_attribute("Dimensions") = dimensions.c_str();
        data_item.append_attribute("Format") = mDataFormat.c_str();
        data_item.append_attribute("Precision") = precision;

        if (mDataFormat == "Binary") {
            const std::string bin_name = mBase + std::to_string(rDataCounter++) + ".bin";
            std::ofstream bin_file(bin_name, std::ios::binary);
            KRATOS_ERROR_IF_NOT(bin_file.good()) << "Could not open " << bin_name << " for writing" << std::endl;
            bin_file.write(reinterpret_cast<const char*>(rArray.Data()),
                           static_cast<std::streamsize>(rArray.Nbytes()));
            data_item.text().set(bin_name.c_str());
            return;
        }
#ifdef MESHIOPLUSPLUS_HAS_HDF5
        if (mDataFormat == "HDF") {
            const std::string dataset_name = "data" + std::to_string(rDataCounter++);
            mio::h5::write_dataset(rH5File, dataset_name, rArray);
            data_item.text().set(
                (std::filesystem::path(mBase + ".h5").filename().string() + ":/" + dataset_name).c_str());
            return;
        }
#endif
        // XML inline text
        const std::size_t number_of_rows = rArray.Shape().empty() ? 0 : rArray.Shape()[0];
        const std::size_t number_of_columns = number_of_rows > 0 ? rArray.Size() / number_of_rows : 0;
        const bool is_float = rArray.Dtype() == mio::DType::Float32 || rArray.Dtype() == mio::DType::Float64;
        std::string text = "\n";
        char buffer[64];
        for (std::size_t i_row = 0; i_row < number_of_rows; ++i_row) {
            for (std::size_t i_column = 0; i_column < number_of_columns; ++i_column) {
                const std::size_t index = i_row * number_of_columns + i_column;
                if (is_float) {
                    const double value = rArray.Dtype() == mio::DType::Float32
                        ? static_cast<double>(rArray.As<float>()[index])
                        : rArray.As<double>()[index];
                    std::snprintf(buffer, sizeof(buffer), "%.*e", static_cast<int>(mPrecision), value);
                } else {
                    std::snprintf(buffer, sizeof(buffer), "%lld",
                                  static_cast<long long>(rArray.As<std::int64_t>()[index]));
                }
                if (i_column > 0) text += " ";
                text += buffer;
            }
            text += "\n";
        }
        data_item.text().set(text.c_str());
    }
};

#undef KRATOS_MESHIOPLUSPLUS_H5_ARG

/***********************************************************************************/
/***********************************************************************************/

XdmfTimeSeriesWriter::XdmfTimeSeriesWriter(
    const std::filesystem::path& rFileName,
    const std::string& rDataFormat,
    const std::size_t Precision
    ) : mpImpl(std::make_unique<Impl>(rFileName, rDataFormat, Precision))
{
}

XdmfTimeSeriesWriter::~XdmfTimeSeriesWriter() = default;

bool XdmfTimeSeriesWriter::MeshWritten() const
{
    return mpImpl->MeshWritten();
}

void XdmfTimeSeriesWriter::WriteMesh(
    const ModelPart& rModelPart,
    const bool WriteElements,
    const bool WriteConditions,
    const bool WriteDeformedConfiguration
    )
{
    mpImpl->WriteMesh(rModelPart, WriteElements, WriteConditions, WriteDeformedConfiguration);
}

void XdmfTimeSeriesWriter::AppendStep(
    const double TimeValue,
    const std::vector<DataArray>& rPointData,
    const std::vector<DataArray>& rCellData
    )
{
    mpImpl->AppendStep(TimeValue, rPointData, rCellData);
}

/***********************************************************************************/
/***********************************************************************************/

MeshioPlusPlusIO::MeshioPlusPlusIO(
    const std::filesystem::path& rFileName,
    Parameters ThisParameters
    ) : mFileName(rFileName),
        mParameters(ThisParameters)
{
    KRATOS_TRY

    mParameters.ValidateAndAssignDefaults(GetDefaultParameters());

    // Eagerly validate the enumerated settings so misconfiguration fails at
    // construction and not in the middle of a simulation.
    const std::string time_series = mParameters["time_series"].GetString();
    KRATOS_ERROR_IF(time_series != "automatic" && time_series != "file_series" && time_series != "single_file")
        << "Unknown \"time_series\" setting \"" << time_series
        << "\" (use \"automatic\", \"file_series\" or \"single_file\")" << std::endl;

    const std::string output_control = mParameters["output_control_type"].GetString();
    KRATOS_ERROR_IF(output_control != "step" && output_control != "time")
        << "Unknown \"output_control_type\" setting \"" << output_control
        << "\" (use \"step\" or \"time\")" << std::endl;

    const std::string xdmf_data_format = mParameters["xdmf_data_format"].GetString();
    KRATOS_ERROR_IF(xdmf_data_format != "auto" && xdmf_data_format != "XML" &&
                    xdmf_data_format != "Binary" && xdmf_data_format != "HDF")
        << "Unknown \"xdmf_data_format\" setting \"" << xdmf_data_format
        << "\" (use \"auto\", \"XML\", \"Binary\" or \"HDF\")" << std::endl;

    const std::string entity_type = mParameters["entity_type"].GetString();
    KRATOS_ERROR_IF(entity_type != "automatic" && entity_type != "element" && entity_type != "condition")
        << "Unknown \"entity_type\" setting \"" << entity_type
        << "\" (use \"automatic\", \"element\" or \"condition\")" << std::endl;

    const std::string file_format = mParameters["file_format"].GetString();
    KRATOS_ERROR_IF(file_format != "default" && file_format != "ascii" && file_format != "binary")
        << "Unknown \"file_format\" setting \"" << file_format
        << "\" (use \"default\", \"ascii\" or \"binary\")" << std::endl;

    // Resolves the format name, throwing a descriptive error for unknown names
    // and formats compiled out of this build.
    FormatFromString(mParameters["format"].GetString());

    // Gauss point results extrapolated to the nodes are written as
    // non-historical nodal data (VtkOutput convention); the extrapolation
    // process itself is created lazily at the first write (it needs the
    // model part, which the IO does not hold at construction).
    for (const auto& r_name : mParameters["gauss_point_variables_extrapolated_to_nodes"].GetStringArray()) {
        mParameters["nodal_data_value_variables"].Append(r_name);
    }

    KRATOS_CATCH("")
}

MeshioPlusPlusIO::~MeshioPlusPlusIO() = default;

/***********************************************************************************/
/***********************************************************************************/

Parameters MeshioPlusPlusIO::GetDefaultParameters()
{
    return Parameters(R"({
        "format"                                      : "auto",
        "file_format"                                 : "default",
        "skin"                                        : true,
        "time_series"                                 : "automatic",
        "output_control_type"                         : "step",
        "output_precision"                            : 7,
        "label_precision"                             : 4,
        "custom_name_prefix"                          : "",
        "custom_name_postfix"                         : "",
        "entity_type"                                 : "automatic",
        "output_sub_model_parts"                      : false,
        "write_deformed_configuration"                : false,
        "write_ids"                                   : false,
        "xdmf_data_format"                            : "auto",
        "nodal_solution_step_data_variables"          : [],
        "nodal_data_value_variables"                  : [],
        "nodal_flags"                                 : [],
        "element_data_value_variables"                : [],
        "element_flags"                               : [],
        "condition_data_value_variables"              : [],
        "condition_flags"                             : [],
        "gauss_point_variables_extrapolated_to_nodes" : [],
        "gauss_point_variables_in_elements"           : []
    })");
}

/***********************************************************************************/
/***********************************************************************************/

std::vector<std::string> MeshioPlusPlusIO::GetSupportedFormats()
{
    std::set<std::string> names;
    for (const auto& r_entry : mio::registry_readers()) {
        names.insert(r_entry.first);
    }
    for (const auto& r_entry : mio::registry_writers()) {
        names.insert(r_entry.first);
    }
    return std::vector<std::string>(names.begin(), names.end());
}

std::vector<std::string> MeshioPlusPlusIO::GetSupportedReadFormats()
{
    std::vector<std::string> names;
    names.reserve(mio::registry_readers().size());
    for (const auto& r_entry : mio::registry_readers()) {
        names.push_back(r_entry.first);
    }
    return names;
}

std::vector<std::string> MeshioPlusPlusIO::GetSupportedWriteFormats()
{
    std::vector<std::string> names;
    names.reserve(mio::registry_writers().size());
    for (const auto& r_entry : mio::registry_writers()) {
        names.push_back(r_entry.first);
    }
    return names;
}

/***********************************************************************************/
/***********************************************************************************/

MeshioPlusPlusIO::Format MeshioPlusPlusIO::FormatFromString(const std::string& rFormatName)
{
    KRATOS_TRY

    std::string lower_name = rFormatName;
    std::transform(lower_name.begin(), lower_name.end(), lower_name.begin(),
                   [](unsigned char Character) { return std::tolower(Character); });

    if (lower_name.empty() || lower_name == "auto" || lower_name == "automatic") {
        return Format::AUTOMATIC;
    }

    const auto& r_format_map = GetFormatNameMap();
    const auto it = r_format_map.find(lower_name);
    if (it == r_format_map.end()) {
        std::ostringstream supported;
        for (const auto& r_name : GetSupportedFormats()) {
            supported << r_name << " ";
        }
        KRATOS_ERROR << "Unknown format \"" << rFormatName
                     << "\". Supported formats in this build: " << supported.str() << std::endl;
    }

    // Known format compiled out of this build: fail with the missing dependency.
    const char* p_missing_dependency = mio::registry_compiled_out(lower_name);
    KRATOS_ERROR_IF(p_missing_dependency != nullptr)
        << "Format \"" << lower_name << "\" is not available in this build (requires "
        << p_missing_dependency << ")" << std::endl;

    return it->second;

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

std::string MeshioPlusPlusIO::FormatName(const Format Value)
{
    switch (Value) {
        case Format::AUTOMATIC: return "auto";
#define KRATOS_MESHIOPLUSPLUS_FORMAT_NAME(EnumName, Name) \
        case Format::EnumName: return Name;
        KRATOS_MESHIOPLUSPLUS_FORMATS(KRATOS_MESHIOPLUSPLUS_FORMAT_NAME)
#undef KRATOS_MESHIOPLUSPLUS_FORMAT_NAME
    }
    KRATOS_ERROR << "Unknown Format enum value" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

MeshioPlusPlusIO::Format MeshioPlusPlusIO::ResolveFormat(const std::filesystem::path& rPath)
{
    KRATOS_TRY

    std::string format_name;
    try {
        format_name = mio::resolve_format(rPath.string(), "");
    } catch (const std::exception& r_exception) {
        KRATOS_ERROR << "Cannot resolve a format from the extension of " << rPath
                     << ": " << r_exception.what() << std::endl;
    }

    const auto& r_format_map = GetFormatNameMap();
    const auto it = r_format_map.find(format_name);
    KRATOS_ERROR_IF(it == r_format_map.end())
        << "The extension of " << rPath << " resolves to unknown format \"" << format_name << "\"" << std::endl;
    return it->second;

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

bool MeshioPlusPlusIO::IsFormatAvailable(const Format Value)
{
    if (Value == Format::AUTOMATIC) {
        return true;
    }
    const std::string name = FormatName(Value);
    if (mio::registry_compiled_out(name) != nullptr) {
        return false;
    }
    return mio::registry_readers().count(name) > 0 || mio::registry_writers().count(name) > 0;
}

/***********************************************************************************/
/***********************************************************************************/

std::string MeshioPlusPlusIO::ResolveEffectiveFormat(const bool CheckWritable) const
{
    KRATOS_TRY

    std::string format_name = mParameters["format"].GetString();
    std::transform(format_name.begin(), format_name.end(), format_name.begin(),
                   [](unsigned char Character) { return std::tolower(Character); });

    if (format_name.empty() || format_name == "auto" || format_name == "automatic") {
        try {
            format_name = mio::resolve_format(mFileName.string(), "");
        } catch (const std::exception& r_exception) {
            KRATOS_ERROR << "Cannot resolve a format from the extension of " << mFileName
                         << ": " << r_exception.what() << ". Set the \"format\" setting explicitly." << std::endl;
        }
    }

    const char* p_missing_dependency = mio::registry_compiled_out(format_name);
    KRATOS_ERROR_IF(p_missing_dependency != nullptr)
        << "Format \"" << format_name << "\" is not available in this build (requires "
        << p_missing_dependency << ")" << std::endl;

    if (CheckWritable) {
        KRATOS_ERROR_IF(mio::registry_writers().count(format_name) == 0)
            << "Format \"" << format_name << "\" is read-only or unknown - it cannot be written" << std::endl;
    } else {
        KRATOS_ERROR_IF(mio::registry_readers().count(format_name) == 0)
            << "Format \"" << format_name << "\" is write-only or unknown - it cannot be read" << std::endl;
    }

    return format_name;

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void MeshioPlusPlusIO::ReadModelPart(ModelPart& rThisModelPart)
{
    KRATOS_TRY

    KRATOS_ERROR_IF(rThisModelPart.IsDistributed())
        << "MeshioPlusPlusIO does not support distributed model parts" << std::endl;
    KRATOS_ERROR_IF(rThisModelPart.NumberOfNodes() > 0)
        << "The target model part \"" << rThisModelPart.FullName() << "\" must be empty" << std::endl;
    KRATOS_ERROR_IF_NOT(std::filesystem::exists(mFileName))
        << "The file " << mFileName << " does not exist" << std::endl;

    const std::string format_name = ResolveEffectiveFormat(false);

    // Read into the meshio++ Kratos backend: the materialized model part view
    // splits max-dimension cell blocks into elements, lower-dimension ones into
    // conditions, and turns integer tag arrays (gmsh physical groups, ...) into
    // sub model parts.
    mio::Mesh mesh = mio::registry_readers().at(format_name)(mFileName.string());
    mio::ModelPart& r_source = mesh.GetModelPart();

    // Bridge into the real Kratos model part (one bulk O(n) creation pass).
    mio::to_model_part(r_source, rThisModelPart, [&rThisModelPart](std::size_t PropertiesId) {
        return rThisModelPart.HasProperties(PropertiesId)
            ? rThisModelPart.pGetProperties(PropertiesId)
            : rThisModelPart.CreateNewProperties(PropertiesId);
    });

    KRATOS_INFO("MeshioPlusPlusIO") << "Read " << rThisModelPart.NumberOfNodes() << " nodes, "
        << rThisModelPart.NumberOfElements() << " elements and " << rThisModelPart.NumberOfConditions()
        << " conditions from " << mFileName << " (format \"" << format_name << "\")" << std::endl;

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void MeshioPlusPlusIO::WriteModelPart(const ModelPart& rThisModelPart)
{
    KRATOS_TRY

    KRATOS_ERROR_IF(rThisModelPart.IsDistributed())
        << "MeshioPlusPlusIO does not support distributed model parts" << std::endl;

    // Extrapolate the requested gauss point results to nodal data first (the
    // process mutates nodal values, hence the contained const_cast - identical
    // effective behavior to VtkOutput, which holds a non-const reference).
    if (mParameters["gauss_point_variables_extrapolated_to_nodes"].size() > 0) {
        if (mpGaussToNodesProcess == nullptr) {
            Parameters gauss_parameters(R"({
                "echo_level"                 : 0,
                "area_average"               : true,
                "average_variable"           : "NODAL_AREA",
                "list_of_variables"          : [],
                "extrapolate_non_historical" : true
            })");
            gauss_parameters["list_of_variables"].SetStringArray(
                mParameters["gauss_point_variables_extrapolated_to_nodes"].GetStringArray());
            mpGaussToNodesProcess = std::make_unique<IntegrationValuesExtrapolationToNodesProcess>(
                const_cast<ModelPart&>(rThisModelPart), gauss_parameters);
        }
        mpGaussToNodesProcess->Execute();
    }

    WriteTarget(rThisModelPart, "");

    if (mParameters["output_sub_model_parts"].GetBool()) {
        for (const auto& r_sub_model_part : rThisModelPart.SubModelParts()) {
            WriteTarget(r_sub_model_part, "_" + rThisModelPart.Name() + "_" + r_sub_model_part.Name());
        }
    }

    ++mOutputStep;

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void MeshioPlusPlusIO::WriteTarget(
    const ModelPart& rModelPart,
    const std::string& rTargetSuffix
    )
{
    KRATOS_TRY

    const std::string format_name = ResolveEffectiveFormat(true);
    const std::string time_series = mParameters["time_series"].GetString();

    try {
        if (format_name == "xdmf" && time_series == "automatic") {
            // Transient XDMF: extend the target file's temporal collection in place.
            WriteXdmfStep(rModelPart, rTargetSuffix);
        } else if (time_series == "single_file") {
            WriteStatic(ComposeOutputPath(rTargetSuffix, ""), format_name, rModelPart);
        } else { // "file_series", or "automatic" for formats without in-file appending
            WriteStatic(ComposeOutputPath(rTargetSuffix, GetOutputLabel(rModelPart)), format_name, rModelPart);
        }
    } catch (const std::out_of_range& r_exception) {
        // meshio++ throws std::out_of_range when an entity references a node
        // that is not part of the written model part (possible for
        // inconsistent sub model parts when "output_sub_model_parts" is on).
        KRATOS_ERROR << "Could not write model part \"" << rModelPart.FullName()
                     << "\": an entity references a node the model part does not contain ("
                     << r_exception.what() << ")" << std::endl;
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

bool MeshioPlusPlusIO::WritesElements() const
{
    return mParameters["entity_type"].GetString() != "condition";
}

bool MeshioPlusPlusIO::WritesConditions() const
{
    return mParameters["entity_type"].GetString() != "element";
}

/***********************************************************************************/
/***********************************************************************************/

std::filesystem::path MeshioPlusPlusIO::ComposeOutputPath(
    const std::string& rTargetSuffix,
    const std::string& rLabel
    ) const
{
    const std::filesystem::path directory = mFileName.parent_path();
    std::string name = mParameters["custom_name_prefix"].GetString()
        + mFileName.stem().string()
        + rTargetSuffix
        + mParameters["custom_name_postfix"].GetString()
        + (rLabel.empty() ? "" : "_" + rLabel)
        + mFileName.extension().string();
    return directory / name;
}

/***********************************************************************************/
/***********************************************************************************/

std::string MeshioPlusPlusIO::GetOutputLabel(const ModelPart& rModelPart) const
{
    const auto& r_process_info = rModelPart.GetProcessInfo();
    if (mParameters["output_control_type"].GetString() == "time") {
        const double time = r_process_info.Has(TIME) ? r_process_info[TIME] : static_cast<double>(mOutputStep);
        std::ostringstream label;
        label << std::fixed << std::setprecision(mParameters["label_precision"].GetInt()) << time;
        return label.str();
    }
    const int step = r_process_info.Has(STEP) ? r_process_info[STEP] : static_cast<int>(mOutputStep);
    return std::to_string(step);
}

/***********************************************************************************/
/***********************************************************************************/

double MeshioPlusPlusIO::GetOutputTimeValue(const ModelPart& rModelPart) const
{
    const auto& r_process_info = rModelPart.GetProcessInfo();
    if (mParameters["output_control_type"].GetString() == "time") {
        return r_process_info.Has(TIME) ? r_process_info[TIME] : static_cast<double>(mOutputStep);
    }
    return static_cast<double>(r_process_info.Has(STEP) ? r_process_info[STEP] : static_cast<int>(mOutputStep));
}

/***********************************************************************************/
/***********************************************************************************/

std::vector<XdmfTimeSeriesWriter::DataArray> MeshioPlusPlusIO::CollectPointData(const ModelPart& rThisModelPart) const
{
    KRATOS_TRY

    std::vector<DataArray> point_data;
    const auto& r_nodes = rThisModelPart.Nodes();
    const std::size_t number_of_nodes = rThisModelPart.NumberOfNodes();

    // Historical nodal variables (validated against the solution step data)
    CollectVariableDataArrays(
        mParameters["nodal_solution_step_data_variables"].GetStringArray(), r_nodes, number_of_nodes,
        [](const auto& rNode, const auto& rVariable) -> decltype(auto) {
            return rNode.FastGetSolutionStepValue(rVariable);
        },
        [&rThisModelPart](const auto& rVariable) {
            KRATOS_ERROR_IF(!rThisModelPart.HasNodalSolutionStepVariable(rVariable))
                << "Variable " << rVariable.Name() << " is not a nodal solution step variable of model part \""
                << rThisModelPart.FullName() << "\"" << std::endl;
        },
        point_data);

    // Non-historical nodal variables (includes the extrapolated gauss point results)
    CollectVariableDataArrays(
        mParameters["nodal_data_value_variables"].GetStringArray(), r_nodes, number_of_nodes,
        [](const auto& rNode, const auto& rVariable) -> decltype(auto) { return rNode.GetValue(rVariable); },
        [](const auto&) {}, point_data);

    // Nodal flags (1/0, -1 when undefined - VtkOutput convention)
    CollectFlagDataArrays(mParameters["nodal_flags"].GetStringArray(), r_nodes, number_of_nodes, point_data);

    if (mParameters["write_ids"].GetBool()) {
        point_data.push_back(CollectIdsArray(r_nodes, number_of_nodes, "KRATOS_NODE_ID"));
    }

    return point_data;

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

std::vector<XdmfTimeSeriesWriter::DataArray> MeshioPlusPlusIO::CollectCellData(const ModelPart& rThisModelPart) const
{
    KRATOS_TRY

    const bool write_elements = WritesElements();
    const bool write_conditions = WritesConditions();
    const std::size_t number_of_elements = write_elements ? rThisModelPart.NumberOfElements() : 0;
    const std::size_t number_of_conditions = write_conditions ? rThisModelPart.NumberOfConditions() : 0;
    const bool write_ids = mParameters["write_ids"].GetBool();
    const auto gauss_point_variables = mParameters["gauss_point_variables_in_elements"].GetStringArray();

    auto non_historical_getter = [](const auto& rEntity, const auto& rVariable) -> decltype(auto) {
        return rEntity.GetValue(rVariable);
    };
    auto no_validation = [](const auto&) {};

    std::vector<DataArray> element_part;
    if (write_elements && number_of_elements > 0) {
        const auto& r_elements = rThisModelPart.Elements();
        CollectVariableDataArrays(mParameters["element_data_value_variables"].GetStringArray(), r_elements,
                                  number_of_elements, non_historical_getter, no_validation, element_part);
        CollectFlagDataArrays(mParameters["element_flags"].GetStringArray(), r_elements, number_of_elements,
                              element_part);
        if (write_ids) {
            element_part.push_back(CollectIdsArray(r_elements, number_of_elements, "KRATOS_ELEMENT_ID"));
            element_part.push_back(CollectPropertiesIdsArray(r_elements, number_of_elements));
        }
        if (!gauss_point_variables.empty()) {
            // CalculateOnIntegrationPoints is non-const, hence the contained
            // const_cast (identical effective behavior to VtkOutput)
            auto& r_mutable_elements = const_cast<ModelPart&>(rThisModelPart).Elements();
            CollectGaussPointDataArrays(gauss_point_variables, r_mutable_elements, number_of_elements,
                                        rThisModelPart.GetProcessInfo(), element_part);
        }
    }

    std::vector<DataArray> condition_part;
    if (write_conditions && number_of_conditions > 0) {
        const auto& r_conditions = rThisModelPart.Conditions();
        CollectVariableDataArrays(mParameters["condition_data_value_variables"].GetStringArray(), r_conditions,
                                  number_of_conditions, non_historical_getter, no_validation, condition_part);
        CollectFlagDataArrays(mParameters["condition_flags"].GetStringArray(), r_conditions, number_of_conditions,
                              condition_part);
        if (write_ids) {
            condition_part.push_back(CollectIdsArray(r_conditions, number_of_conditions, "KRATOS_CONDITION_ID"));
            condition_part.push_back(CollectPropertiesIdsArray(r_conditions, number_of_conditions));
        }
        if (!gauss_point_variables.empty()) {
            auto& r_mutable_conditions = const_cast<ModelPart&>(rThisModelPart).Conditions();
            CollectGaussPointDataArrays(gauss_point_variables, r_mutable_conditions, number_of_conditions,
                                        rThisModelPart.GetProcessInfo(), condition_part);
        }
    }

    return MergeCellDataParts(std::move(element_part), number_of_elements,
                              std::move(condition_part), number_of_conditions);

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void MeshioPlusPlusIO::WriteStatic(
    const std::filesystem::path& rPath,
    const std::string& rFormatName,
    const ModelPart& rThisModelPart
    ) const
{
    KRATOS_TRY

    const bool write_elements = WritesElements();
    const bool write_conditions = WritesConditions();

    mio::Mesh mesh;
    FillMeshioModelPart(rThisModelPart, mesh.GetModelPart(), write_elements, write_conditions,
                        mParameters["write_deformed_configuration"].GetBool());

    // Nodal data set on the meshio++ model part becomes point data when the
    // staging is rebuilt from the model part view (same node container order).
    for (const auto& r_data : CollectPointData(rThisModelPart)) {
        mesh.GetModelPart().SetNodalData(r_data.Name, ToNDArray(r_data));
    }

    // Cell data: the combined arrays cover element rows first, then condition
    // rows; split them into the per-kind containers meshio++ restores per block.
    const std::size_t number_of_elements = write_elements ? rThisModelPart.NumberOfElements() : 0;
    const std::size_t number_of_conditions = write_conditions ? rThisModelPart.NumberOfConditions() : 0;
    for (const auto& r_data : CollectCellData(rThisModelPart)) {
        if (number_of_elements > 0) {
            DataArray element_slice;
            element_slice.Name = r_data.Name;
            element_slice.NumberOfComponents = r_data.NumberOfComponents;
            element_slice.Values.assign(r_data.Values.begin(),
                                        r_data.Values.begin() + number_of_elements * r_data.NumberOfComponents);
            mesh.GetModelPart().SetElementalData(r_data.Name, ToNDArray(element_slice));
        }
        if (number_of_conditions > 0) {
            DataArray condition_slice;
            condition_slice.Name = r_data.Name;
            condition_slice.NumberOfComponents = r_data.NumberOfComponents;
            condition_slice.Values.assign(r_data.Values.begin() + number_of_elements * r_data.NumberOfComponents,
                                          r_data.Values.end());
            mesh.GetModelPart().SetConditionalData(r_data.Name, ToNDArray(condition_slice));
        }
    }

    // The model part view was mutated directly: rebuild the point/cell staging
    // the format writers read from.
    mesh.InvalidateBlocks();

    // Honor an ascii/binary override where the format supports it
    const std::string file_format = mParameters["file_format"].GetString();
    const bool skin = mParameters["skin"].GetBool();
    if (file_format != "default") {
        if (WriteWithFileFormatOverride(rFormatName, file_format == "binary", skin, rPath, mesh)) {
            return;
        }
        KRATOS_WARNING_ONCE("MeshioPlusPlusIO") << "The \"file_format\" setting is not supported for format \""
            << rFormatName << "\" - using its default" << std::endl;
    }

    // The registry writers extract the boundary skin of volume meshes for the
    // surface-only formats; honor a "skin" opt-out by calling those writers
    // directly with their registry ascii/binary defaults.
    if (!skin) {
        if (rFormatName == "stl") {
            mio::write_stl(rPath.string(), mesh, /*binary=*/false, /*skin=*/false);
            return;
        }
        if (rFormatName == "ply") {
            mio::write_ply(rPath.string(), mesh, /*binary=*/true, /*skin=*/false);
            return;
        }
    }

    mio::registry_writers().at(rFormatName)(rPath.string(), mesh);

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void MeshioPlusPlusIO::WriteXdmfStep(
    const ModelPart& rThisModelPart,
    const std::string& rTargetSuffix
    )
{
    KRATOS_TRY

    auto& p_writer = mXdmfWriters[rTargetSuffix];
    if (p_writer == nullptr) {
        std::string data_format = mParameters["xdmf_data_format"].GetString();
        if (data_format == "auto") {
#ifdef MESHIOPLUSPLUS_HAS_HDF5
            data_format = "HDF";
#else
            data_format = "Binary";
#endif
        }
        p_writer = std::make_unique<XdmfTimeSeriesWriter>(
            ComposeOutputPath(rTargetSuffix, ""), data_format,
            static_cast<std::size_t>(mParameters["output_precision"].GetInt()));
    }

    // First write of this IO: extend the file if it already holds a valid time
    // series ("the buffer already exists"), otherwise create it with the mesh.
    if (mOutputStep == 0 && !p_writer->MeshWritten()) {
        p_writer->WriteMesh(rThisModelPart, WritesElements(), WritesConditions(),
                            mParameters["write_deformed_configuration"].GetBool());
    }

    p_writer->AppendStep(GetOutputTimeValue(rThisModelPart), CollectPointData(rThisModelPart),
                         CollectCellData(rThisModelPart));

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

std::string MeshioPlusPlusIO::Info() const
{
    return "MeshioPlusPlusIO";
}

void MeshioPlusPlusIO::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "MeshioPlusPlusIO";
}

void MeshioPlusPlusIO::PrintData(std::ostream& rOStream) const
{
    rOStream << "File name: " << mFileName << "\n"
             << "Settings: " << mParameters.PrettyPrintJsonString() << std::endl;
}

} // namespace Kratos
