// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Jonathan Nuttall
//

#include "dgeoflow.h"
#include "custom_utilities/file_input_utility.h"
#include "geo_output_writer.h"
#include "includes/model_part_io.h"
#include "input_output/logger.h"
#include "input_output/logger_output.h"
#include "input_output/logger_table_output.h"
#include "processes/apply_constant_scalarvalue_process.h"
#include <iomanip>
#include <sstream>

class GeoFlowApplyConstantScalarValueProcess : public Kratos::ApplyConstantScalarValueProcess
{
public:
    GeoFlowApplyConstantScalarValueProcess(Kratos::ModelPart&              rModelPart,
                                           const Kratos::Variable<double>& rVariable,
                                           double                          DoubleValue,
                                           std::size_t                     MeshId,
                                           const Flags&                    rOptions)
        : Kratos::ApplyConstantScalarValueProcess(rModelPart, rVariable, DoubleValue, MeshId, rOptions)
    {
    }

    bool hasWaterPressure() const { return mvariable_name == "WATER_PRESSURE"; }

    Kratos::ModelPart& GetModelPart() { return mr_model_part; }
};

class GeoFlowApplyConstantHydrostaticPressureProcess : public Kratos::ApplyConstantHydrostaticPressureProcess
{
public:
    GeoFlowApplyConstantHydrostaticPressureProcess(Kratos::ModelPart& rModelPart, const Kratos::Parameters& rSettings)
        : Kratos::ApplyConstantHydrostaticPressureProcess(rModelPart, rSettings)
    {
    }

    Kratos::ModelPart& GetModelPart() { return mrModelPart; }

    double GetReferenceCoord() const { return mReferenceCoordinate; }

    void SetReferenceCoord(double value) { mReferenceCoordinate = value; }

    bool hasWaterPressure() const { return mVariableName == "WATER_PRESSURE"; }
};

namespace Kratos
{
KratosExecute::KratosExecute()
{
    KRATOS_INFO("KratosExecute") << "Setting Up Kratos" << std::endl;

    if (!mKernel.IsImported("GeoMechanicsApplication")) {
        KRATOS_INFO("KratosExecute") << "Importing GeoMechanicsApplication" << std::endl;
        mpGeoApp = Kratos::make_shared<KratosGeoMechanicsApplication>();
        mKernel.ImportApplication(mpGeoApp);
    }

    ParallelUtilities::SetNumThreads(1);
    if (this->GetEchoLevel() > 0) {
        Kratos::OpenMPUtils::PrintOMPInfo();
    }

    this->SetEchoLevel(0);
}

int KratosExecute::GetEchoLevel() const { return mEchoLevel; }

void KratosExecute::SetEchoLevel(int level) { mEchoLevel = level; }

void KratosExecute::ResetModelParts()
{
    KRATOS_INFO("Resetting Model") << "Setting Up Execution" << std::endl;
    mCurrentModel.Reset();
}

KratosExecute::ConvergenceCriteriaType::Pointer KratosExecute::setup_criteria_dgeoflow()
{
    const double  rel_tol      = 1.0e-4;
    const double  abs_tol      = 1.0e-9;
    VariableData* p_water_pres = &WATER_PRESSURE;
    KratosExecute::ConvergenceVariableListType convergence_settings{std::make_tuple(p_water_pres, rel_tol, abs_tol)};
    return std::make_shared<KratosExecute::MixedGenericCriteriaType>(convergence_settings);
}

KratosExecute::LinearSolverType::Pointer KratosExecute::setup_solver_dgeoflow()
{
    return Kratos::make_shared<SkylineLUFactorizationSolverType>();
}

KratosExecute::GeoMechanicsNewtonRaphsonErosionProcessStrategyType::Pointer KratosExecute::setup_strategy_dgeoflow(ModelPart& rModelPart)
{
    // Create the linear strategy
    auto p_solver = setup_solver_dgeoflow();

    Scheme<SparseSpaceType, LocalSpaceType>::Pointer p_scheme =
        Kratos::make_shared<BackwardEulerQuasistaticPwScheme<SparseSpaceType, LocalSpaceType>>();

    auto p_builder_and_solver =
        Kratos::make_shared<ResidualBasedBlockBuilderAndSolver<SparseSpaceType, LocalSpaceType, KratosExecute::LinearSolverType>>(
            p_solver);
    p_builder_and_solver->SetEchoLevel(0);

    auto p_criteria = setup_criteria_dgeoflow();
    p_criteria->SetEchoLevel(0);

    Parameters p_parameters(R"(
    {
        "min_iteration":    6,
        "number_cycles":    100,
        "increase_factor":  2.0,
        "reduction_factor": 0.5,
		"max_piping_iterations": 500,
        "desired_iterations": 4,
        "max_radius_factor": 10.0,
        "min_radius_factor": 0.1,
        "search_neighbours_step": false,
        "body_domain_sub_model_part_list": [],
        "loads_sub_model_part_list": [],
        "loads_variable_list" : []
    }  )");

    int  MaxIterations          = 15;
    bool CalculateReactions     = true;
    bool ReformDofSetAtEachStep = false;
    bool MoveMeshFlag           = false;

    auto pSolvingStrategy =
        Kratos::make_unique<GeoMechanicsNewtonRaphsonErosionProcessStrategy<SparseSpaceType, LocalSpaceType, KratosExecute::LinearSolverType>>(
            rModelPart, p_scheme, p_solver, p_criteria, p_builder_and_solver, p_parameters,
            MaxIterations, CalculateReactions, ReformDofSetAtEachStep, MoveMeshFlag);

    pSolvingStrategy->Check();
    return pSolvingStrategy;
}

void KratosExecute::ParseProcesses(ModelPart& rModelPart, Parameters projFile)
{
    // Currently: In DGeoflow only fixed hydrostatic head has been , also need load of gravity.

    auto constraints_processes = projFile["processes"]["constraints_process_list"];
    for (Parameters process : constraints_processes) {
        // we only support fixed hydrostatic head
        auto name          = process["Parameters"]["model_part_name"].GetString();
        auto pressure_type = process["Parameters"]["fluid_pressure_type"].GetString();

        std::size_t found   = name.find_last_of('.');
        std::string subname = name.substr(found + 1);

        ModelPart& part = rModelPart.GetSubModelPart(subname);

        if (pressure_type == "Uniform") {
            auto value = process["Parameters"]["value"].GetDouble();
            mProcesses.push_back(make_shared<GeoFlowApplyConstantScalarValueProcess>(
                part, WATER_PRESSURE, value, 0, ApplyConstantScalarValueProcess::VARIABLE_IS_FIXED));
        } else if (pressure_type == "Hydrostatic") {
            auto cProcesses = process.Clone();
            cProcesses["Parameters"].RemoveValue("fluid_pressure_type");
            mProcesses.push_back(make_shared<GeoFlowApplyConstantHydrostaticPressureProcess>(
                part, cProcesses["Parameters"]));
        } else {
            KRATOS_ERROR << "Reading Processing - Not Implemented - Pressure_type" << std::endl;
        }
    }

    auto loads_processes = projFile["processes"]["loads_process_list"];
    // Should only have one.
    auto        name = loads_processes.GetArrayItem(0)["Parameters"]["model_part_name"].GetString();
    std::size_t found   = name.find_last_of('.');
    std::string subname = name.substr(found + 1);
    ModelPart&  part    = rModelPart.GetSubModelPart(subname);
    mProcesses.push_back(make_shared<ApplyConstantScalarValueProcess>(
        part, VOLUME_ACCELERATION_X, 0.0, 0, ApplyConstantScalarValueProcess::VARIABLE_IS_FIXED));

    mProcesses.push_back(make_shared<ApplyConstantScalarValueProcess>(
        part, VOLUME_ACCELERATION_Y, -9.81, 0, ApplyConstantScalarValueProcess::VARIABLE_IS_FIXED));

    mProcesses.push_back(make_shared<ApplyConstantScalarValueProcess>(
        part, VOLUME_ACCELERATION_Z, 0.0, 0, ApplyConstantScalarValueProcess::VARIABLE_IS_FIXED));
}

int KratosExecute::MainExecution(ModelPart& rModelPart,
                                 const GeoMechanicsNewtonRaphsonErosionProcessStrategyType::Pointer rpSolvingStrategy,
                                 double       Time,
                                 double       DeltaTime,
                                 unsigned int NumberOfIterations) const
{
    // Initialize
    for (const auto& process : mProcesses) {
        process->ExecuteInitialize();
    }

    for (const auto& process : mProcesses) {
        process->ExecuteBeforeSolutionLoop();
    }

    for (unsigned int iter = 0; iter < NumberOfIterations; ++iter) {
        Time += DeltaTime;
        rModelPart.CloneTimeStep(Time);
        rpSolvingStrategy->Initialize();
        rpSolvingStrategy->InitializeSolutionStep();

        for (const auto& process : mProcesses) {
            process->ExecuteInitializeSolutionStep();
        }

        rpSolvingStrategy->Predict();
        rpSolvingStrategy->SolveSolutionStep();

        for (const auto& process : mProcesses) {
            process->ExecuteFinalizeSolutionStep();
        }

        rpSolvingStrategy->FinalizeSolutionStep();
    }

    for (const auto& process : mProcesses) {
        process->ExecuteFinalize();
    }

    return 0;
}

int KratosExecute::ExecuteFlowAnalysis(std::string_view         WorkingDirectory,
                                       const std::string&       rProjectParamsFileName,
                                       const CriticalHeadInfo&  rCriticalHeadInfo,
                                       std::string_view         CriticalHeadBoundaryModelPartName,
                                       const CallBackFunctions& rCallBackFunctions)
{
    mWorkingDirectory                  = WorkingDirectory;
    mCriticalHeadBoundaryModelPartName = CriticalHeadBoundaryModelPartName;

    this->SetEchoLevel(1);

    std::stringstream kratos_log_buffer;
    auto              pOutput = std::make_shared<LoggerOutput>(kratos_log_buffer);
    Logger::AddOutput(pOutput);

    try {
        rCallBackFunctions.ReportProgress(0.0);

        std::string            projectpath = mWorkingDirectory + "/" + rProjectParamsFileName;
        const FileInputUtility input_utility;
        auto                   projectfile = input_utility.ProjectParametersFromFile(projectpath);

        auto materialname =
            projectfile["solver_settings"]["material_import_settings"]["materials_filename"].GetString();
        std::string materialpath = mWorkingDirectory + "/" + materialname;

        auto modelName = projectfile["solver_settings"]["model_part_name"].GetString();

        ModelPart& rModelPart = mCurrentModel.CreateModelPart(modelName);
        rModelPart.SetBufferSize(2);

        KRATOS_INFO_IF("GeoFlowKernel", this->GetEchoLevel() > 0)
            << "Working Directory: " << mWorkingDirectory << std::endl;
        KRATOS_INFO_IF("GeoFlowKernel", this->GetEchoLevel() > 0)
            << "Project Name: " << rProjectParamsFileName << std::endl;

        const auto pSolvingStrategy = setup_strategy_dgeoflow(rModelPart);
        pSolvingStrategy->SetEchoLevel(0);

        KRATOS_INFO_IF("GeoFlowKernel", this->GetEchoLevel() > 0) << "Setup Solving Strategy" << std::endl;

        AddNodalSolutionStepVariables(rModelPart);

        KRATOS_INFO_IF("GeoFlowKernel", this->GetEchoLevel() > 0)
            << "Nodal Solution Variables Added" << std::endl;

        // Don't include the file extension of the mesh file name, since that is automatically
        // appended by the constructor of class ModelPartIO
        const auto mesh_file_name =
            projectfile["solver_settings"]["model_import_settings"]["input_filename"].GetString();
        const auto  mesh_file_path = mWorkingDirectory + "/" + mesh_file_name;
        ModelPartIO reader{mesh_file_path};
        reader.ReadModelPart(rModelPart);

        KRATOS_INFO_IF("GeoFlowKernel", this->GetEchoLevel() > 0) << "Parsed Mesh" << std::endl;

        input_utility.AddMaterialsFromFile(materialpath, mCurrentModel);

        KRATOS_INFO_IF("GeoFlowKernel", this->GetEchoLevel() > 0) << "Parsed Material" << std::endl;

        // Dofs for Water Pressure
        VariableUtils().AddDofWithReaction(WATER_PRESSURE, REACTION_WATER_PRESSURE, rModelPart);
        VariableUtils().AddDof(VOLUME_ACCELERATION_X, rModelPart);
        VariableUtils().AddDof(VOLUME_ACCELERATION_Y, rModelPart);
        VariableUtils().AddDof(VOLUME_ACCELERATION_Z, rModelPart);

        KRATOS_INFO_IF("GeoFlowKernel", this->GetEchoLevel() > 0) << "Added DoF" << std::endl;

        ParseProcesses(rModelPart, projectfile);

        KRATOS_INFO_IF("GeoFlowKernel", this->GetEchoLevel() > 0) << "Parsed Process Data" << std::endl;

        bool has_piping = rCriticalHeadInfo.stepCriticalHead != 0;

        if (rCallBackFunctions.ShouldCancel()) {
            HandleCleanUp(rCallBackFunctions, pOutput);

            return 0;
        }

        const auto rGidOutputSettings =
            projectfile["output_processes"]["gid_output"][0]["Parameters"];

        if (has_piping) {
            ExecuteWithPiping(rModelPart, rGidOutputSettings, rCriticalHeadInfo, pOutput,
                              rCallBackFunctions, pSolvingStrategy);
        } else {
            ExecuteWithoutPiping(rModelPart, rGidOutputSettings, pSolvingStrategy);
        }

        HandleCleanUp(rCallBackFunctions, pOutput);

        return 0;
    } catch (const std::exception& exc) {
        KRATOS_INFO_IF("GeoFlowKernel", this->GetEchoLevel() > 0) << exc.what();

        HandleCleanUp(rCallBackFunctions, pOutput);

        return 1;
    }
}

void KratosExecute::ExecuteWithoutPiping(ModelPart&                rModelPart,
                                         const Kratos::Parameters& rGidOutputSettings,
                                         const GeoMechanicsNewtonRaphsonErosionProcessStrategyType::Pointer pSolvingStrategy) const
{
    MainExecution(rModelPart, pSolvingStrategy, 0.0, 1.0, 1);

    GeoOutputWriter writer{rGidOutputSettings, mWorkingDirectory, rModelPart};
    writer.WriteGiDOutput(rModelPart, rGidOutputSettings);
}

int KratosExecute::ExecuteWithPiping(ModelPart&                rModelPart,
                                     const Kratos::Parameters& rGidOutputSettings,
                                     const CriticalHeadInfo&   rCriticalHeadInfo,
                                     LoggerOutput::Pointer     pOutput,
                                     const CallBackFunctions&  rCallBackFunctions,
                                     const GeoMechanicsNewtonRaphsonErosionProcessStrategyType::Pointer pSolvingStrategy)
{
    KRATOS_INFO_IF("GeoFlowKernel", this->GetEchoLevel() > 0) << "Critical head search started." << std::endl;
    KRATOS_INFO_IF("GeoFlowKernel", this->GetEchoLevel() > 0)
        << "Critical head min head: " << rCriticalHeadInfo.minCriticalHead << std::endl;
    KRATOS_INFO_IF("GeoFlowKernel", this->GetEchoLevel() > 0)
        << "Critical head max head: " << rCriticalHeadInfo.maxCriticalHead << std::endl;
    KRATOS_INFO_IF("GeoFlowKernel", this->GetEchoLevel() > 0)
        << "Critical head step size: " << rCriticalHeadInfo.stepCriticalHead << std::endl;

    shared_ptr<Process> p_river_boundary;
    if (mCriticalHeadBoundaryModelPartName.empty()) {
        p_river_boundary = FindRiverBoundaryAutomatically(pSolvingStrategy);
    } else {
        p_river_boundary = FindRiverBoundaryByName(mCriticalHeadBoundaryModelPartName);
    }

    if (!p_river_boundary) {
        KRATOS_ERROR << "No river boundary found.";
    }

    FindCriticalHead(rModelPart, rGidOutputSettings, rCriticalHeadInfo, pOutput, p_river_boundary,
                     pSolvingStrategy, rCallBackFunctions);

    WriteCriticalHeadResultToFile();

    KRATOS_INFO_IF("GeoFlowKernel", this->GetEchoLevel() > 0) << "Finished writing result" << std::endl;
    return 0;
}

void KratosExecute::WriteCriticalHeadResultToFile() const
{
    KRATOS_INFO_IF("GeoFlowKernel", this->GetEchoLevel() > 0)
        << "Writing result to: " << mWorkingDirectory << "\\criticalHead.json" << std::endl;

    // output critical head_json
    std::ofstream critical_head_file(mWorkingDirectory + "\\criticalHead.json");

    critical_head_file << "{\n";
    critical_head_file << "\t \"PipeData\":\t{\n";
    if (mPipingSuccess) {
        critical_head_file << "\t\t \"Success\": \"True\",\n";
        critical_head_file << "\t\t \"CriticalHead\": \"" + std::to_string(mCriticalHead) + "\"\n";
    } else {
        critical_head_file << "\t\t \"Success\": \"False\"\n";
    }
    critical_head_file << "\t }\n";
    critical_head_file << "}\n";

    critical_head_file.close();
}

void KratosExecute::AddNodalSolutionStepVariables(ModelPart& rModelPart) const
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

    // Smoothing
    rModelPart.AddNodalSolutionStepVariable(NODAL_AREA);
    rModelPart.AddNodalSolutionStepVariable(NODAL_CAUCHY_STRESS_TENSOR);
    rModelPart.AddNodalSolutionStepVariable(NODAL_DAMAGE_VARIABLE);
    rModelPart.AddNodalSolutionStepVariable(NODAL_JOINT_AREA);
    rModelPart.AddNodalSolutionStepVariable(NODAL_JOINT_WIDTH);
    rModelPart.AddNodalSolutionStepVariable(NODAL_JOINT_DAMAGE);
}

int KratosExecute::FindCriticalHead(ModelPart&                 rModelPart,
                                    const Kratos::Parameters&  rGidOutputSettings,
                                    const CriticalHeadInfo&    rCriticalHeadInfo,
                                    LoggerOutput::Pointer      pOutput,
                                    const shared_ptr<Process>& pRiverBoundary,
                                    const GeoMechanicsNewtonRaphsonErosionProcessStrategyType::Pointer pSolvingStrategy,
                                    const CallBackFunctions& rCallBackFunctions)
{
    auto current_process =
        std::static_pointer_cast<GeoFlowApplyConstantHydrostaticPressureProcess>(pRiverBoundary);
    KRATOS_INFO_IF("GeoFlowKernel", this->GetEchoLevel() > 0)
        << "River boundary name: " << current_process->GetName() << std::endl;

    current_process->SetReferenceCoord(rCriticalHeadInfo.minCriticalHead);
    mCurrentHead  = rCriticalHeadInfo.minCriticalHead;
    mCriticalHead = mCurrentHead;

    std::vector<Element*> pipe_elements;
    pipe_elements             = pSolvingStrategy->GetPipingElements();
    const auto noPipeElements = pipe_elements.size();

    int        step = 1;
    const auto max_steps =
        static_cast<int>(std::ceil((rCriticalHeadInfo.maxCriticalHead - rCriticalHeadInfo.minCriticalHead) /
                                   rCriticalHeadInfo.stepCriticalHead)) +
        2;

    while (!AreExceedingMaxCriticalHead(mCurrentHead, rCriticalHeadInfo.maxCriticalHead)) {
        KRATOS_INFO_IF("GeoFlowKernel", this->GetEchoLevel() > 0)
            << "Searching at head: " << mCurrentHead << std::endl;

        std::ostringstream current_head_stream;
        current_head_stream << std::setprecision(8) << std::noshowpoint << mCurrentHead;
        std::string current_head_string = current_head_stream.str();

        std::string progress = "Calculating head level " + current_head_string + "m (" +
                               std::to_string(step) + "/" + std::to_string(max_steps) + ")";
        rCallBackFunctions.ReportTextualProgress(progress.data());
        rCallBackFunctions.ReportProgress(((double)step) / ((double)max_steps));

        MainExecution(rModelPart, pSolvingStrategy, 0.0, 1.0, 1);

        auto count = std::size_t{0};
        for (Element* element : pipe_elements) {
            if (element->GetValue(PIPE_ACTIVE)) count += 1;
        }

        KRATOS_INFO_IF("GeoFlowKernel", this->GetEchoLevel() > 0)
            << "Open pipe elements: " << count << std::endl;

        if (count == noPipeElements) {
            HandleCriticalHeadFound(rCriticalHeadInfo);
            break;
        }

        GeoOutputWriter writer{rGidOutputSettings, mWorkingDirectory, rModelPart};
        writer.WriteGiDOutput(rModelPart, rGidOutputSettings);

        // Update boundary conditions for next search head.
        if (pRiverBoundary->Info() == "ApplyConstantScalarValueProcess") {
            ResetModelParts();
            KRATOS_ERROR << "ApplyConstantScalarValueProcess process search is not implemented.";
        }

        if (pRiverBoundary->Info() == "ApplyConstantHydrostaticPressureProcess") {
            mCriticalHead = current_process->GetReferenceCoord();
            mCurrentHead  = mCriticalHead + rCriticalHeadInfo.stepCriticalHead;
            current_process->SetReferenceCoord(mCurrentHead);
            step++;
        }

        if (rCallBackFunctions.ShouldCancel()) {
            HandleCleanUp(rCallBackFunctions, pOutput);

            return 0;
        }
    }
    return 0;
}

void KratosExecute::HandleCriticalHeadFound(const CriticalHeadInfo& rCriticalHeadInfo)
{
    if (std::abs(mCurrentHead - rCriticalHeadInfo.minCriticalHead) < 1e-9) {
        KRATOS_INFO_IF("GeoFlowKernel", this->GetEchoLevel() > 0)
            << "Critical head undetermined: All pipe elements open at initial search value :"
            << rCriticalHeadInfo.minCriticalHead << std::endl;
    } else {
        mPipingSuccess = true;
        KRATOS_INFO_IF("GeoFlowKernel", this->GetEchoLevel() > 0)
            << "Critical head found: " << mCriticalHead << std::endl;
    }
}

void KratosExecute::HandleCleanUp(const CallBackFunctions& rCallBackFunctions, LoggerOutput::Pointer pOutput)
{
    std::stringstream kratos_log_buffer;

    rCallBackFunctions.LogCallback(kratos_log_buffer.str().c_str());
    Logger::RemoveOutput(pOutput);
    ResetModelParts();
}

bool KratosExecute::AreExceedingMaxCriticalHead(double CurrentHead, double MaxCriticalHead) const
{
    const auto result = (CurrentHead > MaxCriticalHead + 1e-9);
    KRATOS_INFO_IF("GeoFlowKernel", result && (this->GetEchoLevel() > 0))
        << "Critical head undetermined at " << CurrentHead
        << ", max search head reached: " << MaxCriticalHead << std::endl;
    return result;
}

shared_ptr<Process> KratosExecute::FindRiverBoundaryByName(const std::string& CriticalHeadBoundaryModelPartName) const
{
    shared_ptr<Process> p_river_boundary;

    for (const auto& process : mProcesses) {
        if (process->Info() == "ApplyConstantHydrostaticPressureProcess") {
            auto current_process =
                std::static_pointer_cast<GeoFlowApplyConstantHydrostaticPressureProcess>(process);
            if (current_process->hasWaterPressure() &&
                (current_process->GetName() == CriticalHeadBoundaryModelPartName)) {
                p_river_boundary = current_process;
            }
        }
    }

    KRATOS_ERROR_IF_NOT(p_river_boundary) << "No boundary found with the model part name "
                                          << CriticalHeadBoundaryModelPartName << "." << std::endl;

    return p_river_boundary;
}

shared_ptr<Process> KratosExecute::FindRiverBoundaryAutomatically(
    const KratosExecute::GeoMechanicsNewtonRaphsonErosionProcessStrategyType::Pointer rpSolvingStrategy) const
{
    shared_ptr<Process> p_river_boundary;

    std::vector<Element*> pipe_elements;
    pipe_elements = rpSolvingStrategy->GetPipingElements();

    double firstNode_A = pipe_elements.front()->GetGeometry().GetPoint(0).X0();
    double firstNode_B = pipe_elements.front()->GetGeometry().GetPoint(1).X0();
    double lastNode_A  = pipe_elements.back()->GetGeometry().GetPoint(0).X0();

    IndexType RiverNode;

    if ((firstNode_A < lastNode_A) && (firstNode_A < firstNode_B)) {
        RiverNode = pipe_elements.back()->GetGeometry().GetPoint(1).Id();
    } else {
        RiverNode = pipe_elements.back()->GetGeometry().GetPoint(0).Id();
    }

    // Get Find boundary in Processes
    for (const auto& process : mProcesses) {
        ModelPart* currentModelPart = nullptr;

        if (process->Info() == "ApplyConstantScalarValueProcess") {
            auto current_process = std::static_pointer_cast<GeoFlowApplyConstantScalarValueProcess>(process);
            currentModelPart = &current_process->GetModelPart();
            if (current_process->hasWaterPressure()) {
                currentModelPart->GetNode(RiverNode);
                p_river_boundary = current_process;
            }
        } else if (process->Info() == "ApplyConstantHydrostaticPressureProcess") {
            auto current_process =
                std::static_pointer_cast<GeoFlowApplyConstantHydrostaticPressureProcess>(process);
            currentModelPart = &current_process->GetModelPart();
            if (current_process->hasWaterPressure()) {
                currentModelPart->GetNode(RiverNode);
                p_river_boundary = current_process;
            }
        }
    }

    KRATOS_ERROR_IF_NOT(p_river_boundary)
        << "No boundary found on the river side at node " << RiverNode << "." << std::endl;

    return p_river_boundary;
}

} // namespace Kratos
