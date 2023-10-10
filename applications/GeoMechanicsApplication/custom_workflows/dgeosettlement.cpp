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
#include "time_loop_executor.h"

#include "utilities/variable_utils.h"

#include "custom_processes/apply_scalar_constraints_table_process.h"
#include "custom_processes/apply_vector_constraints_table_process.hpp"
#include "custom_processes/set_parameter_field_process.hpp"
#include "custom_processes/apply_k0_procedure_process.hpp"
#include "custom_processes/apply_excavation_process.hpp"

#include "custom_utilities/input_utility.h"
#include "custom_utilities/process_factory.hpp"
#include "custom_utilities/process_info_parser.h"

namespace Kratos
{

KratosGeoSettlement::KratosGeoSettlement(std::unique_ptr<InputUtility> pInputUtility,
                                         std::unique_ptr<ProcessInfoParser> pProcessInfoParser,
                                         std::unique_ptr<TimeLoopExecutor> pTimeLoopExecutor) :
    mpInputUtility{std::move(pInputUtility)},
    mpProcessInfoParser{std::move(pProcessInfoParser)},
    mpTimeLoopExecutor{std::move(pTimeLoopExecutor)}
{
    KRATOS_INFO("KratosGeoSettlement") << "Setting up Kratos" << std::endl;
    KRATOS_ERROR_IF_NOT(mpInputUtility) << "Invalid Input Utility";

    if (!mKernel.IsImported("GeoMechanicsApplication"))
    {
        KRATOS_INFO("KratosGeoSettlement") << "Importing GeoMechanicsApplication" << std::endl;
        mpGeoApp = Kratos::make_shared<KratosGeoMechanicsApplication>();
        mKernel.ImportApplication(mpGeoApp);
    }

    InitializeProcessFactory();
}

void KratosGeoSettlement::InitializeProcessFactory()
{
    mProcessFactory->AddCreator("ApplyScalarConstraintsTableProcess",
                                [this](const Parameters& rParameters)
                                {
                                    return std::make_unique<ApplyScalarConstraintsTableProcess>(mModel.GetModelPart(mModelPartName),
                                                                                                rParameters);
                                });

    mProcessFactory->AddCreator("ApplyVectorConstraintsTableProcess",
                                [this](const Parameters& rParameters)
                                {
                                    return std::make_unique<ApplyVectorConstraintsTableProcess>(mModel.GetModelPart(mModelPartName),
                                                                                                rParameters);
                                });

    mProcessFactory->AddCreator("SetParameterFieldProcess",
                                [this](const Parameters& rParameters)
                                {
                                    return std::make_unique<SetParameterFieldProcess>(mModel.GetModelPart(mModelPartName),
                                                                                      rParameters);
                                });

    mProcessFactory->AddCreator("ApplyExcavationProcess",
                                [this](const Parameters& rParameters)
                                {
                                    return std::make_unique<ApplyExcavationProcess>(mModel.GetModelPart(mModelPartName),
                                                                                    rParameters);
                                });

    mProcessFactory->AddCreator("ApplyK0ProcedureProcess",
                                [this](const Parameters& rParameters)
                                {
                                    return std::make_unique<ApplyK0ProcedureProcess>(mModel.GetModelPart(mModelPartName),
                                                                                     rParameters);
                                });

    mProcessFactory->SetCallBackWhenProcessIsUnknown([](const std::string& rProcessName)
    {
        KRATOS_ERROR << "Unexpected process (" << rProcessName << "), calculation is aborted";
    });
}

int KratosGeoSettlement::RunStage(const std::filesystem::path&            rWorkingDirectory,
                                  const std::filesystem::path&            rProjectParametersFile,
                                  const std::function<void(const char*)>& ,
                                  const std::function<void(double)>&      ,
                                  const std::function<void(const char*)>& ,
                                  const std::function<bool()>&            )
{
    KRATOS_INFO("KratosGeoSettlement") << "About to run a stage..." << std::endl;

    const auto project_parameters_file_path = rWorkingDirectory / rProjectParametersFile;
    const auto project_parameters = mpInputUtility->ProjectParametersFromFile(
            project_parameters_file_path.generic_string());
    KRATOS_INFO("KratosGeoSettlement") << "Parsed project parameters file " << project_parameters_file_path << std::endl;

    mModelPartName = project_parameters["solver_settings"]["model_part_name"].GetString();
    if (const auto model_part_name = project_parameters["solver_settings"]["model_part_name"].GetString();
        !mModel.HasModelPart(model_part_name)) {
        auto& model_part = AddNewModelPart(model_part_name);
        const auto mesh_file_name = project_parameters["solver_settings"]["model_import_settings"]["input_filename"].GetString();
        mpInputUtility->ReadModelFromFile(rWorkingDirectory / mesh_file_name, model_part);
    }

    if (project_parameters["solver_settings"].Has("material_import_settings")) {
        const auto material_file_name = project_parameters["solver_settings"]["material_import_settings"]["materials_filename"].GetString();
        const auto material_file_path = rWorkingDirectory / material_file_name;
        mpInputUtility->AddMaterialsFromFile(material_file_path.generic_string(), mModel);
        KRATOS_INFO("KratosGeoSettlement") << "Read the materials from " << material_file_path << std::endl;
    }

    std::vector<std::shared_ptr<Process>> processes = GetProcesses(project_parameters);
    std::vector<std::weak_ptr<Process>> process_observables(processes.begin(), processes.end());

    if (mpTimeLoopExecutor)
    {
        mpTimeLoopExecutor->SetProcessReferences(process_observables);
    }

    return 0;
}

ModelPart& KratosGeoSettlement::AddNewModelPart(const std::string& rModelPartName)
{
    auto& result = mModel.CreateModelPart(rModelPartName);
    KRATOS_INFO("KratosGeoSettlement") << "Created a new model part named '" << rModelPartName << "'" << std::endl;

    result.SetBufferSize(2);

    AddNodalSolutionStepVariablesTo(result);
    KRATOS_INFO("KratosGeoSettlement") << "Added nodal solution step variables" << std::endl;

    AddDegreesOfFreedomTo(result);
    KRATOS_INFO("KratosGeoSettlement") << "Added degrees of freedom" << std::endl;

    return result;
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

const InputUtility* KratosGeoSettlement::GetInterfaceInputUtility() const
{
    return mpInputUtility.get();
}

std::vector<std::shared_ptr<Process>> KratosGeoSettlement::GetProcesses(const Parameters& project_parameters) const
{
    std::vector<std::shared_ptr<Process>> result;
    if (project_parameters.Has("processes")) {
        const auto processes = mpProcessInfoParser->GetProcessList(project_parameters["processes"]);
        for (const auto &process: processes) {
            result.emplace_back(mProcessFactory->Create(process.name, process.parameters));
        }

    }

    return result;
}

// This default destructor is added in the cpp to be able to forward member variables in a unique_ptr
KratosGeoSettlement::~KratosGeoSettlement() = default;

}