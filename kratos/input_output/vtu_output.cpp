//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher (based on vtk_output.h)
//
//

// System includes

// External includes
#ifdef KRATOS_USE_VTU11
#include "inc/vtu11.hpp"
#else
namespace vtu11 { // dummy namespace that mimics vtu11 such that we don't need conditionals in the code

using DataSet = std::tuple<std::string, size_t, std::vector<double>>;

using VtkCellType = std::uint8_t;

template<typename MeshGenerator, typename Writer>
void write( const std::string& filename,
            MeshGenerator& mesh,
            const std::vector<DataSet>& pointData,
            const std::vector<DataSet>& cellData,
            Writer writer = Writer( ) );

struct Vtu11UnstructuredMesh
{
  std::vector<double>& points_;
  std::vector<size_t>& connectivity_;
  std::vector<size_t>& offsets_;
  std::vector<VtkCellType>& types_;
};

}
#endif

// Project includes
#include "vtu_output.h"

namespace Kratos {

namespace { // helpers namespace

void GetNodalCoordinates(const ModelPart& rModelPart, std::vector<double>& rCoordinates, const bool WriteDeformedConfiguration)
{
    // NOTE: also in MPI all nodes (local and ghost) have to be written, because
    // they might be needed by the elements/conditions due to the connectivity

    if (rCoordinates.size() != rModelPart.NumberOfNodes()*3) {
        rCoordinates.resize(rModelPart.NumberOfNodes()*3);
    }

    std::size_t index(0);

    // Write nodes TODO loops should be OMP
    if (WriteDeformedConfiguration) {
        for (const auto& r_node : rModelPart.Nodes()) {
            rCoordinates[index++] = r_node.X();
            rCoordinates[index++] = r_node.Y();
            rCoordinates[index++] = r_node.Z();
        }
    } else {
        for (const auto& r_node : rModelPart.Nodes()) {
            rCoordinates[index++] = r_node.X0();
            rCoordinates[index++] = r_node.Y0();
            rCoordinates[index++] = r_node.Z0();
        }
    }
}

template <typename TContainerType>
void GetCellInformation(const TContainerType& rContainer, const std::unordered_map<int, int>& rKratosIdToVtuId, std::vector<std::size_t>& rConnectivities, std::vector<std::size_t>& rOffsets, std::vector<vtu11::VtkCellType>& rTypes)
{
        // IMPORTANT: The map geo_type_vtk_cell_type_map is to be extended to support new geometries
    // NOTE: See https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf
    const std::map<GeometryData::KratosGeometryType, int> geo_type_vtk_cell_type_map = {
        { GeometryData::KratosGeometryType::Kratos_Point2D,          1 },
        { GeometryData::KratosGeometryType::Kratos_Point3D,          1 },
        { GeometryData::KratosGeometryType::Kratos_Line2D2,          3 },
        { GeometryData::KratosGeometryType::Kratos_Line3D2,          3 },
        { GeometryData::KratosGeometryType::Kratos_Triangle2D3,      5 },
        { GeometryData::KratosGeometryType::Kratos_Triangle3D3,      5 },
        { GeometryData::KratosGeometryType::Kratos_Quadrilateral2D4, 9 },
        { GeometryData::KratosGeometryType::Kratos_Quadrilateral3D4, 9 },
        { GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4,    10 },
        { GeometryData::KratosGeometryType::Kratos_Hexahedra3D8,     12 },
        { GeometryData::KratosGeometryType::Kratos_Prism3D6,         13 },
        { GeometryData::KratosGeometryType::Kratos_Line2D3,          21 },
        { GeometryData::KratosGeometryType::Kratos_Line3D3,          21 },
        { GeometryData::KratosGeometryType::Kratos_Triangle2D6,      22 },
        { GeometryData::KratosGeometryType::Kratos_Triangle3D6,      22 },
        { GeometryData::KratosGeometryType::Kratos_Quadrilateral2D8, 23 },
        { GeometryData::KratosGeometryType::Kratos_Quadrilateral3D8, 23 },
        { GeometryData::KratosGeometryType::Kratos_Tetrahedra3D10,   24 }
//         { GeometryData::KratosGeometryType::Kratos_Hexahedra3D20,    25 } // NOTE: Quadratic hexahedra (20) requires a conversor, order does not coincide with VTK
    };

    const std::size_t constainer_size(rContainer.size());
    rConnectivities.resize(0);
    if (constainer_size > 0) { // this can happen when a domain in MPI has no local Elements/Consitions
        // reserving a guess for the necessary entries
        rConnectivities.reserve(constainer_size * rContainer.begin()->GetGeometry().PointsNumber());
    }

    if (rOffsets.size() != constainer_size) {
        rOffsets.resize(constainer_size);
    }

    if (rTypes.size() != constainer_size) {
        rTypes.resize(constainer_size);
    }

    std::size_t container_index(0);
    for (const auto& r_entity : rContainer) {
        const auto& r_geom = r_entity.GetGeometry();
        rOffsets[container_index] = r_geom.PointsNumber();
        const auto& r_kratos_cell = r_geom.GetGeometryType();
        if (geo_type_vtk_cell_type_map.count(r_kratos_cell) > 0) {
            rTypes[container_index++] = geo_type_vtk_cell_type_map.at(r_kratos_cell);
        } else {
            KRATOS_ERROR << "Modelpart contains elements or conditions with "
             << "geometries for which no VTK-output is implemented!" << std::endl
             << "Cell type: " << static_cast<int>(r_kratos_cell) << std::endl;
        }
        for (const auto& r_node : r_geom) {
            rConnectivities.push_back(rKratosIdToVtuId.at(r_node.Id()));
        }
    }
}

void CreateMapFromKratosIdToVtuId(const ModelPart& rModelPart, std::unordered_map<int, int>& rKratosIdToVtuId)
{
    int vtk_id(0);
    for(const auto& r_node : rModelPart.Nodes()) {
        rKratosIdToVtuId[r_node.Id()] = vtk_id++;
    }
}

} // helpers namespace



VtuOutput::VtuOutput(
    ModelPart& rModelPart,
    Parameters ThisParameters
    ) : mrModelPart(rModelPart),
        mOutputSettings(ThisParameters)
{
#ifndef KRATOS_USE_VTU11
    KRATOS_ERROR << "The vtu output process has to be enabled at compile time by adding \"VTU11_DIR\" in the configure-script" << std::endl;
#endif

    // The default parameters
    Parameters default_parameters = GetDefaultParameters();
    mOutputSettings.ValidateAndAssignDefaults(default_parameters);

    // Initialize other variables
    const std::string file_format = mOutputSettings["file_format"].GetString();
    if (file_format == "ascii") {
        mFileFormat = VtuOutput::FileFormat::VTU_ASCII;
    } else if (file_format == "binary_raw") {
        mFileFormat = VtuOutput::FileFormat::VTU_BINARY_RAW;
    } else if (file_format == "binary_raw_compressed") {
        mFileFormat = VtuOutput::FileFormat::VTU_BINARY_RAW_COMPRESSED;
    } else if (file_format == "binary_base64") {
        mFileFormat = VtuOutput::FileFormat::VTU_BINARY_BASE64;
    } else if (file_format == "binary_base64_appended") {
        mFileFormat = VtuOutput::FileFormat::VTU_BINARY_BASE64_APPENDED;
    } else {
        KRATOS_ERROR << "Option for \"file_format\": " << file_format << " not recognised!\n Possible output formats options are: \"ascii\", \"binary_raw\", \"binary_raw_compressed\", \"binary_base64\", \"binary_base64_appended\"" << std::endl;
    }

    const auto& r_local_mesh = rModelPart.GetCommunicator().LocalMesh();
    const auto& r_data_comm = rModelPart.GetCommunicator().GetDataCommunicator();

    const int num_elements = r_data_comm.SumAll(static_cast<int>(r_local_mesh.NumberOfElements()));
    const int num_conditions = r_data_comm.SumAll(static_cast<int>(r_local_mesh.NumberOfConditions()));

    KRATOS_WARNING_IF("VtuOutput", num_elements > 0 && num_conditions > 0) << r_data_comm << "Modelpart \"" << rModelPart.Name() << "\" has both elements and conditions.\nGiving precedence to elements and writing only elements!" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

void VtuOutput::PrintOutput()
{
    KRATOS_TRY;
    std::vector<double> coordinates;
    KRATOS_WATCH("1")
    GetNodalCoordinates(mrModelPart, coordinates, mOutputSettings["write_deformed_configuration"].GetBool());
    KRATOS_WATCH("2")

    CreateMapFromKratosIdToVtuId(mrModelPart, mKratosIdToVtkId);
    KRATOS_WATCH("3")

    std::vector<std::size_t> connectivities;
    std::vector<std::size_t> offsets;
    std::vector<vtu11::VtkCellType> types;

    const auto& r_local_mesh = mrModelPart.GetCommunicator().LocalMesh();

    const int num_elements = mrModelPart.GetCommunicator().GetDataCommunicator().SumAll(static_cast<int>(r_local_mesh.NumberOfElements()));
    const int num_conditions = mrModelPart.GetCommunicator().GetDataCommunicator().SumAll(static_cast<int>(r_local_mesh.NumberOfConditions()));

    if (num_elements > 0) {
        GetCellInformation(r_local_mesh.Elements(), mKratosIdToVtkId, connectivities, offsets, types);
    } else if (num_conditions > 0) {
        GetCellInformation(r_local_mesh.Conditions(), mKratosIdToVtkId, connectivities, offsets, types);
    }


    vtu11::Vtu11UnstructuredMesh vtu_mesh{ coordinates, connectivities, offsets, types };
    std::vector<vtu11::DataSet> data_dummy;


    vtu11::write("vtu_file_test.vtu", vtu_mesh, data_dummy, data_dummy);
    KRATOS_CATCH("VTU PrintOutput");
}

/***********************************************************************************/
/***********************************************************************************/

std::string VtuOutput::GetOutputFileName(const ModelPart& rModelPart) const
{
    const int rank = rModelPart.GetCommunicator().MyPID();
    std::string model_part_name(rModelPart.Name());

    std::string label;
    std::stringstream ss;
    const std::string output_control = mOutputSettings["output_control_type"].GetString();
    if (output_control == "step") {
        ss << std::fixed << std::setfill('0')
           << rModelPart.GetProcessInfo()[STEP];
        label = ss.str();
    } else if(output_control == "time") {
        ss << std::fixed << std::setfill('0')
           << rModelPart.GetProcessInfo()[TIME];
        label = ss.str();
    } else {
        KRATOS_ERROR << "Option for \"output_control_type\": " << output_control
            <<" not recognised!\nPossible output_control_type options "
            << "are: \"step\", \"time\"" << std::endl;
    }

    // Putting everything together
    std::string output_file_name;
    if (mOutputSettings["save_output_files_in_folder"].GetBool()) {
        output_file_name += mOutputSettings["folder_name"].GetString() + "/";
    }
    const std::string& custom_name_prefix = mOutputSettings["custom_name_prefix"].GetString();
    output_file_name += custom_name_prefix + model_part_name + "_" + std::to_string(rank) + "_" + label + ".vtk";

    return output_file_name;
}

Parameters VtuOutput::GetDefaultParameters()
{
    // IMPORTANT: when "output_control_type" is "time", then paraview will not be able to group them
    Parameters default_parameters = Parameters(R"(
    {
        "model_part_name"                    : "PLEASE_SPECIFY_MODEL_PART_NAME",
        "file_format"                        : "ascii",
        "output_control_type"                : "step",
        "output_frequency"                   : 1.0,
        "folder_name"                        : "VTU_Output",
        "custom_name_prefix"                 : "",
        "save_output_files_in_folder"        : true,
        "write_deformed_configuration"       : false,
        "nodal_solution_step_data_variables" : [],
        "nodal_data_value_variables"         : [],
        "nodal_flags"                        : [],
        "element_data_value_variables"       : [],
        "element_flags"                      : [],
        "condition_data_value_variables"     : [],
        "condition_flags"                    : []
    })" );

    return default_parameters;
}

} // namespace Kratos
