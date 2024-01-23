//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Nicolò Crescenzio
//

// System includes

// External includes

// Project includes
#include "mpm_application_variables.h"
#include "mpm_vtk_output.h"
#include "includes/kratos_filesystem.h"

namespace Kratos
{

Parameters MPMVtkOutput::GetDefaultParameters()
{
    // IMPORTANT: when "output_control_type" is "time", then paraview will not be able to group them
    Parameters default_parameters = Parameters(R"(
    {
        "model_part_name"                             : "PLEASE_SPECIFY_MODEL_PART_NAME",
        "file_format"                                 : "binary",
        "output_precision"                            : 7,
        "output_control_type"                         : "step",
        "output_interval"                             : 1.0,
        "output_sub_model_parts"                      : false,
        "output_path"                                 : "MPM_VTK_Output",
        "custom_name_prefix"                          : "",
        "custom_name_postfix"                         : "",
        "save_output_files_in_folder"                 : true,
        "entity_type"                                 : "automatic",
        "write_ids"                                   : false,
        "element_flags"                               : [],
        "condition_flags"                             : [],
        "gauss_point_variables_in_elements"           : []
    })" );

    return default_parameters;
}

///***********************************************************************************/
///***********************************************************************************/

MPMVtkOutput::MPMVtkOutput(
    ModelPart& rModelPart,
    Parameters ThisParameters
    ) : VtkOutput(rModelPart, ThisParameters)
{
}

///***********************************************************************************/
///***********************************************************************************/

void MPMVtkOutput::WriteNodesToFile(
    const ModelPart& rModelPart,
    std::ofstream& rFileStream
    ) const
{
    const auto& r_local_mesh = rModelPart.GetCommunicator().LocalMesh();

    if (GetEntityType(rModelPart) == EntityType::ELEMENT) {
        rFileStream << "POINTS " << r_local_mesh.NumberOfElements() << " float\n";
        for (auto itr_element = r_local_mesh.ElementsBegin(); itr_element != r_local_mesh.ElementsEnd(); itr_element++) {
            std::vector<array_1d<double, 3>> mp_coord = { ZeroVector(3) };
            itr_element->CalculateOnIntegrationPoints(MP_COORD, mp_coord, rModelPart.GetProcessInfo());
            WriteVectorDataToFile(mp_coord[0], rFileStream);
            if (mFileFormat == MPMVtkOutput::FileFormat::VTK_ASCII) rFileStream << "\n";
        }
    } else if (GetEntityType(rModelPart) == EntityType::CONDITION) {
        rFileStream << "POINTS " << r_local_mesh.NumberOfConditions() << " float\n";
        for (auto itr_condition = r_local_mesh.ConditionsBegin(); itr_condition != r_local_mesh.ConditionsEnd(); itr_condition++) {
            std::vector<array_1d<double, 3>> mpc_coord = { ZeroVector(3) };
            itr_condition->CalculateOnIntegrationPoints(MPC_COORD, mpc_coord, rModelPart.GetProcessInfo());
            WriteVectorDataToFile(mpc_coord[0], rFileStream);
            if (mFileFormat == MPMVtkOutput::FileFormat::VTK_ASCII) rFileStream << "\n";
        }
    } else if (GetEntityType(rModelPart) == EntityType::NONE) {
        rFileStream << "POINTS 0 float\n";
    }
}

///***********************************************************************************/
///***********************************************************************************/

void MPMVtkOutput::WriteConditionsAndElementsToFile(const ModelPart& rModelPart, std::ofstream& rFileStream) const
{
    const auto& r_local_mesh = rModelPart.GetCommunicator().LocalMesh();
    const auto entity_type = GetEntityType(rModelPart);

    if (entity_type == EntityType::ELEMENT) {
        // write cells header
        rFileStream << "\nCELLS " << r_local_mesh.NumberOfElements() << " "
            << 2*r_local_mesh.NumberOfElements() << "\n";
        WriteConnectivity(r_local_mesh.Elements(), rFileStream);
        // write cell types header
        rFileStream << "\nCELL_TYPES " << r_local_mesh.NumberOfElements() << "\n";
        WriteCellType(r_local_mesh.Elements(), rFileStream);
    } else if (entity_type == EntityType::CONDITION) {
        // write cells header
        rFileStream << "\nCELLS " << r_local_mesh.NumberOfConditions() << " "
            << 2*r_local_mesh.NumberOfConditions() << "\n";
        WriteConnectivity(r_local_mesh.Conditions(), rFileStream);
        // write cell types header
        rFileStream << "\nCELL_TYPES " << r_local_mesh.NumberOfConditions() << "\n";
        WriteCellType(r_local_mesh.Conditions(), rFileStream);
    }
}

///***********************************************************************************/
///***********************************************************************************/

template <typename TContainerType>
void MPMVtkOutput::WriteCellType(
    const TContainerType& rContainer,
    std::ofstream& rFileStream
    ) const
{
    // Write entity types
    for (IndexType itr_entity = 0; itr_entity < rContainer.size(); ++itr_entity) {
        WriteScalarDataToFile((unsigned int)1, rFileStream);
        if (mFileFormat == MPMVtkOutput::FileFormat::VTK_ASCII) rFileStream << "\n";
    }
}

///***********************************************************************************/
///***********************************************************************************/

template <typename TContainerType>
void MPMVtkOutput::WriteConnectivity(
    const TContainerType& rContainer,
    std::ofstream& rFileStream
    ) const
{
    for (IndexType itr_entity = 0; itr_entity < rContainer.size(); ++itr_entity) {
        WriteScalarDataToFile((unsigned int)1, rFileStream);
        if (mFileFormat == MPMVtkOutput::FileFormat::VTK_ASCII) rFileStream << " ";
        WriteScalarDataToFile((int)itr_entity, rFileStream);
        if (mFileFormat == MPMVtkOutput::FileFormat::VTK_ASCII) rFileStream << "\n";
    }
}

///***********************************************************************************/
///***********************************************************************************/

void MPMVtkOutput::WriteNodalResultsToFile(
    const ModelPart& rModelPart,
    std::ofstream& rFileStream
    )
{
}

} // namespace Kratos
