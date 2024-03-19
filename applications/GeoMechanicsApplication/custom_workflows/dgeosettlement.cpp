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
#include "time_loop_executor_interface.h"

#include "custom_processes/apply_excavation_process.h"
#include "custom_processes/apply_k0_procedure_process.h"
#include "custom_processes/apply_normal_load_table_process.h"
#include "custom_processes/apply_scalar_constraint_table_process.h"
#include "custom_processes/apply_vector_constraint_table_process.h"
#include "custom_processes/set_parameter_field_process.hpp"

#include "adaptive_time_incrementor.h"
#include "custom_processes/deactivate_conditions_on_inactive_elements_process.hpp"
#include "custom_processes/find_neighbour_elements_of_conditions_process.hpp"
#include "custom_utilities/input_utility.h"
#include "custom_utilities/process_info_parser.h"
#include "custom_utilities/solving_strategy_factory.hpp"
#include "solving_strategy_wrapper.hpp"
#include "spaces/ublas_space.h"

namespace
{

using namespace Kratos;

using SparseSpaceType  = UblasSpace<double, CompressedMatrix, Vector>;
using DenseSpaceType   = UblasSpace<double, Matrix, Vector>;
using LinearSolverType = LinearSolver<SparseSpaceType, DenseSpaceType>;
using SolvingStrategyFactoryType = SolvingStrategyFactory<SparseSpaceType, DenseSpaceType, LinearSolverType>;

double GetStartTimeFrom(const Parameters& rProjectParameters)
{
    return rProjectParameters["problem_data"]["start_time"].GetDouble();
}

double GetEndTimeFrom(const Parameters& rProjectParameters)
{
    return rProjectParameters["problem_data"]["end_time"].GetDouble();
}

double GetTimeIncrementFrom(const Parameters& rProjectParameters)
{
    return rProjectParameters["solver_settings"]["time_stepping"]["time_step"].GetDouble();
}

std::size_t GetMaxNumberOfCyclesFrom(const Parameters& rProjectParameters)
{
    return static_cast<std::size_t>(rProjectParameters["solver_settings"]["number_cycles"].GetInt());
}

double GetReductionFactorFrom(const Parameters& rProjectParameters)
{
    return rProjectParameters["solver_settings"]["reduction_factor"].GetDouble();
}

double GetIncreaseFactorFrom(const Parameters& rProjectParameters)
{
    return rProjectParameters["solver_settings"]["increase_factor"].GetDouble();
}

double GetMaxDeltaTimeFactorFrom(const Parameters& rProjectParameters)
{
    return rProjectParameters["solver_settings"]["time_stepping"]["max_delta_time_factor"].GetDouble();
}

std::size_t GetMinNumberOfIterationsFrom(const Parameters& rProjectParameters)
{
    return static_cast<std::size_t>(rProjectParameters["solver_settings"]["min_iterations"].GetInt());
}

std::size_t GetMaxNumberOfIterationsFrom(const Parameters& rProjectParameters)
{
    return static_cast<std::size_t>(rProjectParameters["solver_settings"]["max_iterations"].GetInt());
}

bool GetResetDisplacementsFrom(const Parameters& rProjectParameters)
{
    return rProjectParameters["solver_settings"]["reset_displacements"].GetBool();
}

} // namespace

namespace Kratos
{

KratosGeoSettlement::KratosGeoSettlement(std::unique_ptr<InputUtility>      pInputUtility,
                                         std::unique_ptr<ProcessInfoParser> pProcessInfoParser,
                                         std::unique_ptr<TimeLoopExecutorInterface> pTimeLoopExecutorInterface)
    : mpInputUtility{std::move(pInputUtility)},
      mpProcessInfoParser{std::move(pProcessInfoParser)},
      mpTimeLoopExecutor{std::move(pTimeLoopExecutorInterface)}
{
    Kernel::GetApplicationsList().clear();
    mModel.Reset();
    KRATOS_INFO("KratosGeoSettlement") << "Setting up Kratos" << std::endl;
    KRATOS_ERROR_IF_NOT(mpInputUtility) << "Invalid Input Utility";

    if (!mKernel.IsImported("GeoMechanicsApplication")) {
        KRATOS_INFO("KratosGeoSettlement") << "Importing GeoMechanicsApplication" << std::endl;
        mpGeoApp = Kratos::make_shared<KratosGeoMechanicsApplication>();
        mKernel.ImportApplication(mpGeoApp);
    }
    if (!mKernel.IsImported("LinearSolversApplication")) {
        KRATOS_INFO("KratosGeoSettlement") << "Importing LinearSolversApplication" << std::endl;
        mpLinearSolversApp = Kratos::make_shared<KratosLinearSolversApplication>();
        mKernel.ImportApplication(mpLinearSolversApp);
    }
    if (!mKernel.IsImported("StructuralMechanicsApplication")) {
        KRATOS_INFO("KratosGeoSettlement") << "Importing StructuralMechanicsApplication" << std::endl;
        mpStructuralMechanicsApp = Kratos::make_shared<KratosStructuralMechanicsApplication>();
        mKernel.ImportApplication(mpStructuralMechanicsApp);
    }

    ParallelUtilities::SetNumThreads(1);

    InitializeProcessFactory();
}

void KratosGeoSettlement::InitializeProcessFactory()
{
    mProcessFactory->AddCreator("ApplyScalarConstraintTableProcess",
                                MakeCreatorFor<ApplyScalarConstraintTableProcess>());
    mProcessFactory->AddCreator("ApplyNormalLoadTableProcess", MakeCreatorFor<ApplyNormalLoadTableProcess>());
    mProcessFactory->AddCreator("ApplyVectorConstraintTableProcess",
                                MakeCreatorFor<ApplyVectorConstraintTableProcess>());
    mProcessFactory->AddCreator("SetParameterFieldProcess", MakeCreatorFor<SetParameterFieldProcess>());
    mProcessFactory->AddCreator("ApplyExcavationProcess", MakeCreatorFor<ApplyExcavationProcess>());
    mProcessFactory->AddCreator("ApplyK0ProcedureProcess", MakeCreatorFor<ApplyK0ProcedureProcess>());

    mProcessFactory->SetCallBackWhenProcessIsUnknown([](const std::string& rProcessName) {
        KRATOS_ERROR << "Unexpected process (" << rProcessName << "), calculation is aborted";
    });
}

int KratosGeoSettlement::RunStage(const std::filesystem::path&            rWorkingDirectory,
                                  const std::filesystem::path&            rProjectParametersFile,
                                  const std::function<void(const char*)>& rLogCallback,
                                  const std::function<void(double)>&      rProgressDelegate,
                                  const std::function<void(const char*)>&,
                                  const std::function<bool()>& rShouldCancel)
{
    std::stringstream     kratos_log_buffer;
    LoggerOutput::Pointer logger_output = CreateLoggingOutput(kratos_log_buffer);

    try {
        KRATOS_INFO("KratosGeoSettlement") << "About to run a stage..." << std::endl;

        const auto project_parameters_file_path = rWorkingDirectory / rProjectParametersFile;
        const auto project_parameters =
            mpInputUtility->ProjectParametersFromFile(project_parameters_file_path.generic_string());
        KRATOS_INFO("KratosGeoSettlement")
            << "Parsed project parameters file " << project_parameters_file_path << std::endl;

        mModelPartName = project_parameters["solver_settings"]["model_part_name"].GetString();
        if (const auto model_part_name = project_parameters["solver_settings"]["model_part_name"].GetString();
            !mModel.HasModelPart(model_part_name)) {
            auto&      model_part = AddNewModelPart(model_part_name);
            const auto mesh_file_name =
                project_parameters["solver_settings"]["model_import_settings"]["input_filename"].GetString();
            mpInputUtility->ReadModelFromFile(rWorkingDirectory / mesh_file_name, model_part);
            AddDegreesOfFreedomTo(model_part);
            KRATOS_INFO("KratosGeoSettlement") << "Added degrees of freedom" << std::endl;
        }

        PrepareModelPart(project_parameters["solver_settings"]);

        if (project_parameters["solver_settings"].Has("material_import_settings")) {
            const auto material_file_name =
                project_parameters["solver_settings"]["material_import_settings"]
                                  ["materials_filename"]
                                      .GetString();
            const auto material_file_path = rWorkingDirectory / material_file_name;
            mpInputUtility->AddMaterialsFromFile(material_file_path.generic_string(), mModel);
            KRATOS_INFO("KratosGeoSettlement") << "Read the materials from " << material_file_path << std::endl;
        }

        std::vector<std::shared_ptr<Process>> processes = GetProcesses(project_parameters);
        std::vector<std::weak_ptr<Process>> process_observables(processes.begin(), processes.end());

        if (mpTimeLoopExecutor) {
            mpTimeLoopExecutor->SetCancelDelegate(rShouldCancel);
            mpTimeLoopExecutor->SetProgressDelegate(rProgressDelegate);
            mpTimeLoopExecutor->SetProcessObservables(process_observables);
            mpTimeLoopExecutor->SetTimeIncrementor(MakeTimeIncrementor(project_parameters));
            mpTimeLoopExecutor->SetSolverStrategyWrapper(MakeStrategyWrapper(project_parameters, rWorkingDirectory));
        }

        for (const auto& process : processes) {
            process->ExecuteInitialize();
        }

        for (const auto& process : processes) {
            process->ExecuteBeforeSolutionLoop();
        }

        if (mpTimeLoopExecutor) {
            // For now, pass a dummy state. THIS PROBABLY NEEDS TO BE REFINED AT SOME POINT!
            TimeStepEndState start_of_loop_state;
            start_of_loop_state.convergence_state = TimeStepEndState::ConvergenceState::converged;
            start_of_loop_state.time              = GetStartTimeFrom(project_parameters);
            mpTimeLoopExecutor->Run(start_of_loop_state);
        }

        for (const auto& process : processes) {
            process->ExecuteFinalize();
        }

        FlushLoggingOutput(rLogCallback, logger_output, kratos_log_buffer);

        return 0;
    } catch (const std::exception& exc) {
        KRATOS_INFO("KratosGeoSettlement") << exc.what();

        FlushLoggingOutput(rLogCallback, logger_output, kratos_log_buffer);

        return 1;
    }
}

ModelPart& KratosGeoSettlement::AddNewModelPart(const std::string& rModelPartName)
{
    auto& result = mModel.CreateModelPart(rModelPartName);
    KRATOS_INFO("KratosGeoSettlement")
        << "Created a new model part named '" << rModelPartName << "'" << std::endl;

    result.SetBufferSize(2);

    AddNodalSolutionStepVariablesTo(result);
    KRATOS_INFO("KratosGeoSettlement") << "Added nodal solution step variables" << std::endl;

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
    rModelPart.AddNodalSolutionStepVariable(LINE_LOAD_Y);
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

void KratosGeoSettlement::AddDegreesOfFreedomTo(Kratos::ModelPart& rModelPart)
{
    VariableUtils().AddDofWithReaction(DISPLACEMENT_X, REACTION_X, rModelPart);
    VariableUtils().AddDofWithReaction(DISPLACEMENT_Y, REACTION_Y, rModelPart);
    VariableUtils().AddDofWithReaction(DISPLACEMENT_Z, REACTION_Z, rModelPart);

    VariableUtils().AddDofWithReaction(WATER_PRESSURE, REACTION_WATER_PRESSURE, rModelPart);
    VariableUtils().AddDof(VOLUME_ACCELERATION_X, rModelPart);
    VariableUtils().AddDof(VOLUME_ACCELERATION_Y, rModelPart);
    VariableUtils().AddDof(VOLUME_ACCELERATION_Z, rModelPart);
}

LoggerOutput::Pointer KratosGeoSettlement::CreateLoggingOutput(std::stringstream& rKratosLogBuffer) const
{
    auto logger_output = std::make_shared<LoggerOutput>(rKratosLogBuffer);
    Logger::AddOutput(logger_output);
    return logger_output;
}

void KratosGeoSettlement::FlushLoggingOutput(const std::function<void(const char*)>& rLogCallback,
                                             LoggerOutput::Pointer                   pLoggerOutput,
                                             const std::stringstream& rKratosLogBuffer) const
{
    rLogCallback(rKratosLogBuffer.str().c_str());
    Logger::RemoveOutput(pLoggerOutput);
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
        for (const auto& process : processes) {
            result.emplace_back(mProcessFactory->Create(process.name, process.parameters));
        }
    }

    return result;
}

std::unique_ptr<TimeIncrementor> KratosGeoSettlement::MakeTimeIncrementor(const Parameters& rProjectParameters)
{
    // For now, we can create adaptive time incrementors only
    return std::make_unique<AdaptiveTimeIncrementor>(
        GetStartTimeFrom(rProjectParameters), GetEndTimeFrom(rProjectParameters),
        GetTimeIncrementFrom(rProjectParameters), GetMaxNumberOfCyclesFrom(rProjectParameters),
        GetReductionFactorFrom(rProjectParameters), GetIncreaseFactorFrom(rProjectParameters),
        GetMaxDeltaTimeFactorFrom(rProjectParameters), GetMinNumberOfIterationsFrom(rProjectParameters),
        GetMaxNumberOfIterationsFrom(rProjectParameters));
}

std::shared_ptr<StrategyWrapper> KratosGeoSettlement::MakeStrategyWrapper(const Parameters& rProjectParameters,
                                                                          const std::filesystem::path& rWorkingDirectory)
{
    auto solving_strategy = SolvingStrategyFactoryType::Create(
        rProjectParameters["solver_settings"], GetComputationalModelPart());
    KRATOS_ERROR_IF_NOT(solving_strategy) << "No solving strategy was created!" << std::endl;

    solving_strategy->Initialize();
    GetMainModelPart().CloneTimeStep();

    if (rProjectParameters["solver_settings"]["reset_displacements"].GetBool()) {
        constexpr auto source_index      = std::size_t{0};
        constexpr auto destination_index = std::size_t{1};
        RestoreValuesOfNodalVariable(DISPLACEMENT, source_index, destination_index);
        RestoreValuesOfNodalVariable(ROTATION, source_index, destination_index);

        VariableUtils{}.UpdateCurrentToInitialConfiguration(GetComputationalModelPart().Nodes());
    }

    FindNeighbourElementsOfConditionsProcess{GetComputationalModelPart()}.Execute();
    DeactivateConditionsOnInactiveElements{GetComputationalModelPart()}.Execute();

    // For now, we can create solving strategy wrappers only
    using SolvingStrategyWrapperType = SolvingStrategyWrapper<SparseSpaceType, DenseSpaceType>;
    return std::make_shared<SolvingStrategyWrapperType>(std::move(solving_strategy),
                                                        GetResetDisplacementsFrom(rProjectParameters),
                                                        rWorkingDirectory, rProjectParameters);
}

void KratosGeoSettlement::PrepareModelPart(const Parameters& rSolverSettings)
{
    auto& main_model_part = GetMainModelPart();
    if (!main_model_part.HasSubModelPart(mComputationalSubModelPartName)) {
        main_model_part.CreateSubModelPart(mComputationalSubModelPartName);
    }

    if (rSolverSettings.Has("nodal_smoothing")) {
        main_model_part.GetProcessInfo().SetValue(NODAL_SMOOTHING,
                                                  rSolverSettings["nodal_smoothing"].GetBool());
    }

    // Note that the computing part and the main model part _share_ their process info and properties
    GetComputationalModelPart().SetProcessInfo(main_model_part.GetProcessInfo());
    for (auto i = ModelPart::SizeType{0}; i < main_model_part.NumberOfMeshes(); ++i) {
        auto& mesh = main_model_part.GetMesh(i);
        for (const auto& property : mesh.Properties()) {
            GetComputationalModelPart().AddProperties(mesh.pGetProperties(property.GetId()));
        }
    }

    GetComputationalModelPart().Set(ACTIVE);

    const auto problem_domain_sub_model_part_list =
        rSolverSettings["problem_domain_sub_model_part_list"];
    std::vector<std::string> domain_part_names;
    for (const auto& sub_model_part : problem_domain_sub_model_part_list) {
        domain_part_names.emplace_back(sub_model_part.GetString());
    }

    // Add nodes to computing model part
    std::set<Node::IndexType> node_id_set;
    for (const auto& name : domain_part_names) {
        auto& domain_part = main_model_part.GetSubModelPart(name);
        for (const auto& node : domain_part.Nodes()) {
            node_id_set.insert(node.Id());
        }
    }
    GetComputationalModelPart().AddNodes(
        std::vector<Node::IndexType>{node_id_set.begin(), node_id_set.end()});

    std::set<IndexedObject::IndexType> element_id_set;
    for (const auto& name : domain_part_names) {
        auto& domain_part = main_model_part.GetSubModelPart(name);
        for (const auto& element : domain_part.Elements()) {
            element_id_set.insert(element.Id());
        }
    }
    GetComputationalModelPart().AddElements(
        std::vector<IndexedObject::IndexType>{element_id_set.begin(), element_id_set.end()});

    GetComputationalModelPart().Conditions().clear();
    const auto processes_sub_model_part_list = rSolverSettings["processes_sub_model_part_list"];
    std::vector<std::string> domain_condition_names;
    for (const auto& sub_model_part : processes_sub_model_part_list) {
        domain_condition_names.emplace_back(sub_model_part.GetString());
    }

    std::set<IndexedObject::IndexType> condition_id_set;
    for (const auto& name : domain_condition_names) {
        auto& domain_part = main_model_part.GetSubModelPart(name);
        for (const auto& condition : domain_part.Conditions()) {
            condition_id_set.insert(condition.Id());
        }
    }
    GetComputationalModelPart().AddConditions(
        std::vector<IndexedObject::IndexType>{condition_id_set.begin(), condition_id_set.end()});

    // NOTE TO SELF: here we don't yet "Adding Computing Sub Sub Model Parts" (see check_and_prepare_model_process_geo.py)
    // We are not sure yet whether that piece of code is relevant or not.
}

ModelPart& KratosGeoSettlement::GetMainModelPart() { return mModel.GetModelPart(mModelPartName); }

ModelPart& KratosGeoSettlement::GetComputationalModelPart()
{
    return GetMainModelPart().GetSubModelPart(mComputationalSubModelPartName);
}

// This default destructor is added in the cpp to be able to forward member variables in a unique_ptr
KratosGeoSettlement::~KratosGeoSettlement() = default;

} // namespace Kratos