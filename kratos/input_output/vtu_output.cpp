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
#include <numeric>

// External includes
#include "vtu11.hpp"

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
void GetCellInformation(const TContainerType& rContainer, const std::unordered_map<int, int>& rKratosIdToVtuId, std::vector<vtu11::VtkIndexType>& rConnectivities, std::vector<vtu11::VtkIndexType>& rOffsets, std::vector<vtu11::VtkCellType>& rTypes)
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

    std::partial_sum(rOffsets.begin(), rOffsets.end(), rOffsets.begin());
}

void CreateMapFromKratosIdToVtuId(const ModelPart& rModelPart, std::unordered_map<int, int>& rKratosIdToVtuId)
{
    int vtk_id(0);
    for(const auto& r_node : rModelPart.Nodes()) {
        rKratosIdToVtuId[r_node.Id()] = vtk_id++;
    }
}

void GetData(
    const ModelPart& rModelPart,
    std::vector<vtu11::DataSet>& rDataSetContainer,
    const std::string& rVariableName)
{
    KRATOS_TRY;

    // TODO sync vars in MPI
    if (KratosComponents<Variable<double>>::Has(rVariableName)){
        const auto& r_var_to_write = KratosComponents<Variable<double>>::Get(rVariableName);
        std::vector<double> results(rModelPart.NumberOfNodes());
        std::size_t counter = 0;
        for (const auto& r_node : rModelPart.Nodes()) {
            results[counter++] = r_node.FastGetSolutionStepValue(r_var_to_write);
        }

        rDataSetContainer.push_back(vtu11::DataSet(rVariableName, 1, results));

    } else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(rVariableName)){
        const auto& r_var_to_write = KratosComponents<Variable<array_1d<double, 3>>>::Get(rVariableName);
        std::vector<double> results(rModelPart.NumberOfNodes()*3);
        std::size_t counter = 0;
        for (const auto& r_node : rModelPart.Nodes()) {
            const auto& r_vals = r_node.FastGetSolutionStepValue(r_var_to_write);
            results[counter++] = r_vals[0];
            results[counter++] = r_vals[1];
            results[counter++] = r_vals[2];
        }

        rDataSetContainer.push_back(vtu11::DataSet(rVariableName, 3, results));


    } else {
        KRATOS_WARNING_ONCE(rVariableName) << rModelPart.GetCommunicator().GetDataCommunicator() << "Variable \"" << rVariableName << "\" is "
            << "not suitable for VtkOutput, skipping it" << std::endl;
    }

    KRATOS_CATCH("VTU GetData");
}

} // helpers namespace


VtuOutput::VtuOutput(
    ModelPart& rModelPart,
    Parameters ThisParameters
    ) : mrModelPart(rModelPart),
        mOutputSettings(ThisParameters)
{
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

void VtuOutput::PrintOutput(const std::string& rOutputFilename)
{
    KRATOS_TRY;

    std::vector<double> coordinates;

    GetNodalCoordinates(mrModelPart, coordinates, mOutputSettings["write_deformed_configuration"].GetBool());

    CreateMapFromKratosIdToVtuId(mrModelPart, mKratosIdToVtkId);

    std::vector<vtu11::VtkIndexType> connectivities;
    std::vector<vtu11::VtkIndexType> offsets;
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
    std::vector<vtu11::DataSet> point_data;
    std::vector<vtu11::DataSet> cell_data; // unused for now

    for (const std::string var_name :  mOutputSettings["nodal_solution_step_data_variables"].GetStringArray()) {
        GetData(mrModelPart, point_data, var_name);
    }

    const std::string output_file_name = GetOutputFileName(mrModelPart, false, rOutputFilename);

    const int my_pid = mrModelPart.GetCommunicator().MyPID();
    const int total_processes = mrModelPart.GetCommunicator().TotalProcesses();

    if (mFileFormat == VtuOutput::FileFormat::VTU_ASCII) {
        auto writer = vtu11::AsciiWriter();
        vtu11::write(output_file_name, vtu_mesh, point_data, cell_data, writer);
    } else if (mFileFormat == VtuOutput::FileFormat::VTU_BINARY_RAW) {
        auto writer = vtu11::RawBinaryAppendedWriter();
        vtu11::write(output_file_name, vtu_mesh, point_data, cell_data, writer);
    } else if (mFileFormat == VtuOutput::FileFormat::VTU_BINARY_RAW_COMPRESSED) {
        auto writer = vtu11::CompressedRawBinaryAppendedWriter();
        vtu11::write(output_file_name, vtu_mesh, point_data, cell_data, writer);
    } else if (mFileFormat == VtuOutput::FileFormat::VTU_BINARY_BASE64) {
        auto writer = vtu11::Base64BinaryWriter();
        vtu11::write(output_file_name, vtu_mesh, point_data, cell_data, writer);
    } else if (mFileFormat == VtuOutput::FileFormat::VTU_BINARY_BASE64_APPENDED) {
        auto writer = vtu11::Base64BinaryAppendedWriter();
        vtu11::parallel_write(output_file_name, vtu_mesh, point_data, cell_data, writer, my_pid, total_processes);
    }

    KRATOS_CATCH("VTU PrintOutput");
}

/***********************************************************************************/
/***********************************************************************************/

std::string VtuOutput::GetOutputFileName(const ModelPart& rModelPart, const bool IsSubModelPart, const std::string& rOutputFilename)
{
    // Putting everything together
    std::string output_file_name = "";
    if (mOutputSettings["save_output_files_in_folder"].GetBool()) {
        output_file_name = mOutputSettings["folder_name"].GetString() + "/";
    }

    if (rOutputFilename != "")
    {
        output_file_name += rOutputFilename + ".vtu";
    }
    else
    {
        const int rank = rModelPart.GetCommunicator().MyPID();
        std::string model_part_name;

        if (IsSubModelPart) {
            model_part_name = rModelPart.GetParentModelPart()->Name() + "_" + rModelPart.Name();
        } else {
            model_part_name = rModelPart.Name();
        }

        std::string label;
        std::stringstream ss;
        const std::string output_control = mOutputSettings["output_control_type"].GetString();
        if (output_control == "step") {
            ss << std::fixed << std::setprecision(mDefaultPrecision)<< std::setfill('0')
            << rModelPart.GetProcessInfo()[STEP];
            label = ss.str();
        } else if(output_control == "time") {
            ss << std::fixed << std::setprecision(mDefaultPrecision) << std::setfill('0')
            << rModelPart.GetProcessInfo()[TIME];
            label = ss.str();
        } else {
            KRATOS_ERROR << "Option for \"output_control_type\": " << output_control
                <<" not recognised!\nPossible output_control_type options "
                << "are: \"step\", \"time\"" << std::endl;
        }


        const std::string& r_custom_name_prefix = mOutputSettings["custom_name_prefix"].GetString();
        const std::string& r_custom_name_postfix = mOutputSettings["custom_name_postfix"].GetString();
        output_file_name += r_custom_name_prefix + model_part_name + r_custom_name_postfix + "_" + std::to_string(rank) + "_" + label + ".vtu";
    }

    return output_file_name;
}

Parameters VtuOutput::GetDefaultParameters()
{
    // IMPORTANT: when "output_control_type" is "time", then paraview will not be able to group them
    Parameters default_parameters = Parameters(R"(
    {
        "model_part_name"                             : "PLEASE_SPECIFY_MODEL_PART_NAME",
        "file_format"                                 : "ascii",
        "output_precision"                            : 7,
        "output_control_type"                         : "step",
        "output_frequency"                            : 1.0,
        "output_sub_model_parts"                      : false,
        "folder_name"                                 : "VTU_Output",
        "custom_name_prefix"                          : "",
        "custom_name_postfix"                         : "",
        "save_output_files_in_folder"                 : true,
        "write_deformed_configuration"                : false,
        "write_ids"                                   : false,
        "nodal_solution_step_data_variables"          : [],
        "nodal_data_value_variables"                  : [],
        "nodal_flags"                                 : [],
        "element_data_value_variables"                : [],
        "element_flags"                               : [],
        "condition_data_value_variables"              : [],
        "condition_flags"                             : [],
        "gauss_point_variables_extrapolated_to_nodes" : [],
        "gauss_point_variables_in_elements"           : []
    })" );

    return default_parameters;
}

} // namespace Kratos
