// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Anne van de Graaf
//
#include "dgeosettlement.h"
#include "input_output/logger.h"
#include "custom_utilities/input_utilities.h"
#include "utilities/variable_utils.h"
#include "includes/model_part_io.h"


namespace Kratos
{

KratosGeoSettlement::KratosGeoSettlement()
{
    KRATOS_INFO("KratosGeoSettlement") << "Setting up Kratos" << std::endl;

    if (!mKernel.IsImported("GeoMechanicsApplication"))
    {
        KRATOS_INFO("KratosGeoSettlement") << "Importing GeoMechanicsApplication" << std::endl;
        mpGeoApp = Kratos::make_shared<KratosGeoMechanicsApplication>();
        mKernel.ImportApplication(mpGeoApp);
    }
}

int KratosGeoSettlement::RunStage(const std::string&                      rWorkingDirectory,
                                  const std::string&                      rProjectParametersFileName,
                                  const std::function<void(const char*)>& ,
                                  const std::function<void(double)>&      ,
                                  const std::function<void(const char*)>& ,
                                  const std::function<bool()>&            )
{
    KRATOS_INFO("KratosGeoSettlement") << "About to run a stage..." << std::endl;

    const auto project_parameters_file_path = rWorkingDirectory + "/" + rProjectParametersFileName;
    const auto project_parameters = InputUtilities::ProjectParametersFrom(project_parameters_file_path);
    KRATOS_INFO("KratosGeoSettlement") << "Parsed project parameters file " << project_parameters_file_path << std::endl;

    const auto model_part_name = project_parameters["solver_settings"]["model_part_name"].GetString();
    ModelPart& model_part = mModel.CreateModelPart(model_part_name);
    model_part.SetBufferSize(2);
    KRATOS_INFO("KratosGeoSettlement") << "Created a model part" << std::endl;

    AddNodalSolutionStepVariablesTo(model_part);
    KRATOS_INFO("KratosGeoSettlement") << "Added nodal solution step variables" << std::endl;

    AddDegreesOfFreedomTo(model_part);
    KRATOS_INFO("KratosGeoSettlement") << "Added degrees of freedom" << std::endl;

    // Don't include the file extension of the mesh file name, since that is automatically appended by the
    // constructor of class ModelPartIO
    const auto mesh_file_name = project_parameters["solver_settings"]["model_import_settings"]["input_filename"].GetString();
    const auto mesh_file_path = rWorkingDirectory + "/" + mesh_file_name;
    ModelPartIO reader{mesh_file_path};
    reader.ReadModelPart(model_part);
    KRATOS_INFO("KratosGeoSettlement") << "Read the mesh data from " << mesh_file_path << std::endl;

    const auto material_file_name = project_parameters["solver_settings"]["material_import_settings"]["materials_filename"].GetString();
    const auto material_file_path = rWorkingDirectory + "/" + material_file_name;
    InputUtilities::AddMaterialsFrom(material_file_path, mModel);
    KRATOS_INFO("KratosGeoSettlement") << "Read the materials from " << material_file_path << std::endl;

    return 0;
}

void KratosGeoSettlement::AddNodalSolutionStepVariablesTo(ModelPart& rModelPart)
{
    rModelPart.AddNodalSolutionStepVariable(VELOCITY);
    rModelPart.AddNodalSolutionStepVariable(ACCELERATION);

    // Displacement
    rModelPart.AddNodalSolutionStepVariable(DISPLACEMENT);
    rModelPart.AddNodalSolutionStepVariable(TOTAL_DISPLACEMENT);
    rModelPart.AddNodalSolutionStepVariable(REACTION);
    rModelPart.AddNodalSolutionStepVariable(POINT_LOAD);
    rModelPart.AddNodalSolutionStepVariable(LINE_LOAD);
    rModelPart.AddNodalSolutionStepVariable(SURFACE_LOAD);
    rModelPart.AddNodalSolutionStepVariable(VOLUME_ACCELERATION);
    rModelPart.AddNodalSolutionStepVariable(NORMAL_CONTACT_STRESS);
    rModelPart.AddNodalSolutionStepVariable(TANGENTIAL_CONTACT_STRESS);

    // Water
    rModelPart.AddNodalSolutionStepVariable(WATER_PRESSURE);
    rModelPart.AddNodalSolutionStepVariable(REACTION_WATER_PRESSURE);
    rModelPart.AddNodalSolutionStepVariable(DT_WATER_PRESSURE);
    rModelPart.AddNodalSolutionStepVariable(NORMAL_FLUID_FLUX);
    rModelPart.AddNodalSolutionStepVariable(HYDRAULIC_DISCHARGE);
}

void KratosGeoSettlement::AddDegreesOfFreedomTo(Kratos::ModelPart &rModelPart)
{
    VariableUtils().AddDofWithReaction(DISPLACEMENT_X, REACTION_X, rModelPart);
    VariableUtils().AddDofWithReaction(DISPLACEMENT_Y, REACTION_Y, rModelPart);
    VariableUtils().AddDofWithReaction(DISPLACEMENT_Z, REACTION_Z, rModelPart);

    VariableUtils().AddDofWithReaction(WATER_PRESSURE, REACTION_WATER_PRESSURE, rModelPart);
    VariableUtils().AddDof(VOLUME_ACCELERATION_X, rModelPart);
    VariableUtils().AddDof(VOLUME_ACCELERATION_Y, rModelPart);
    VariableUtils().AddDof(VOLUME_ACCELERATION_Z, rModelPart);
}

}