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

#pragma once

#include <sstream>
#include <iomanip>
#include "dgeoflow.h"
#include "processes/apply_constant_scalarvalue_process.h"
#include "utilities/read_materials_utility.h"
#include "input_output/logger.h"
#include "input_output/logger_output.h"
#include "input_output/logger_table_output.h"
#include "includes/model_part_io.h"
#include "write_output.h"


class GeoFlowApplyConstantScalarValueProcess : public Kratos::ApplyConstantScalarValueProcess
{
public:
    using ApplyConstantScalarValueProcess ::ApplyConstantScalarValueProcess;

    bool hasWaterPressure()
    {
        return mvariable_name == "WATER_PRESSURE";
    }

    Kratos::ModelPart &GetModelPart()
    {
        return mr_model_part;
    }
};

class GeoFlowApplyConstantHydrostaticPressureProcess : public Kratos::ApplyConstantHydrostaticPressureProcess
{
    using ApplyConstantHydrostaticPressureProcess::ApplyConstantHydrostaticPressureProcess;

public:
    Kratos::ModelPart &GetModelPart()
    {
        return mrModelPart;
    }

    double GetReferenceCoord()
    {
        return mReferenceCoordinate;
    }

    void SetReferenceCoord(double value)
    {
        mReferenceCoordinate = value;
    }

    bool hasWaterPressure()
    {
        return mVariableName == "WATER_PRESSURE";
    }
};


namespace Kratos
{
    KratosExecute::KratosExecute()
    {
        KRATOS_INFO("KratosExecute") << "Setting Up Kratos" << std::endl;

    	if (!kernel.IsImported("GeoMechanicsApplication"))
        {
            KRATOS_INFO("KratosExecute") << "Importing GeoMechanicsApplication" << std::endl;
    		geoApp = Kratos::make_shared<KratosGeoMechanicsApplication>();
            kernel.ImportApplication(geoApp);
        }

        Kratos::OpenMPUtils::SetNumThreads(1);
        if (this->GetEchoLevel() > 0)
        {
            Kratos::OpenMPUtils::PrintOMPInfo();
        }

        this->SetEchoLevel(0);
    }

    int KratosExecute::GetEchoLevel() const
    {
        return echoLevel;
    }

    void KratosExecute::SetEchoLevel(int level)
    {
        echoLevel = level;
    }

    void KratosExecute::ResetModelParts()
    {
        KRATOS_INFO("Resetting Model") << "Setting Up Execution" << std::endl;
        current_model.Reset();
    }

    KratosExecute::ConvergenceCriteriaType::Pointer KratosExecute::setup_criteria_dgeoflow()
    {
        const double rel_tol = 1.0e-4;
        const double abs_tol = 1.0e-9;
        VariableData *p_water_pres = &WATER_PRESSURE;
        KratosExecute::ConvergenceVariableListType convergence_settings{std::make_tuple(p_water_pres, rel_tol, abs_tol)};
        return KratosExecute::ConvergenceCriteriaType::Pointer(new KratosExecute::MixedGenericCriteriaType(convergence_settings));
    }

    KratosExecute::LinearSolverType::Pointer KratosExecute::setup_solver_dgeoflow()
    {
        return Kratos::make_shared<SkylineLUFactorizationSolverType>();
    }

    KratosExecute::GeoMechanicsNewtonRaphsonErosionProcessStrategyType::Pointer KratosExecute::setup_strategy_dgeoflow(ModelPart &model_part)
    {
        // Create the linear strategy
        auto p_solver = setup_solver_dgeoflow();

        Scheme<SparseSpaceType, LocalSpaceType>::Pointer p_scheme = Kratos::make_shared<BackwardEulerQuasistaticPwScheme<SparseSpaceType, LocalSpaceType>>();

        auto p_builder_and_solver = Kratos::make_shared<ResidualBasedBlockBuilderAndSolver<SparseSpaceType, LocalSpaceType, KratosExecute::LinearSolverType>>(p_solver);
        p_builder_and_solver->SetEchoLevel(0);

        auto p_criteria = setup_criteria_dgeoflow();
        p_criteria->SetEchoLevel(0);

        Parameters p_parameters(R"(
    {
        "min_iteration":    6,
        "number_cycles":    100,
        "increase_factor":  2.0,
        "reduction_factor": 0.5,
        "end_time": 1.0,
		"max_piping_iterations": 500,
        "desired_iterations": 4,
        "max_radius_factor": 10.0,
        "min_radius_factor": 0.1,
        "search_neighbours_step": false,
        "body_domain_sub_model_part_list": [],
        "loads_sub_model_part_list": [],
        "loads_variable_list" : []
    }  )");

        int MaxIterations = 15;
        bool CalculateReactions = true;
        bool ReformDofSetAtEachStep = false;
        bool MoveMeshFlag = false;

        auto p_solving_strategy = Kratos::make_unique<GeoMechanicsNewtonRaphsonErosionProcessStrategy<SparseSpaceType, LocalSpaceType, KratosExecute::LinearSolverType>>(
            model_part,
            p_scheme,
            p_solver,
            p_criteria,
            p_builder_and_solver,
            p_parameters,
            MaxIterations, CalculateReactions, ReformDofSetAtEachStep, MoveMeshFlag);

        p_solving_strategy->Check();
        return p_solving_strategy;
    }

    void KratosExecute::parseMaterial(Model &model, const std::string& rMaterialFilePath)
    {
        std::string parameters = R"({ "Parameters" : { "materials_filename" :")" + rMaterialFilePath + R"("}})";
        Parameters material_file{parameters};
        ReadMaterialsUtility(material_file, model);
    }

    Parameters KratosExecute::openProjectParamsFile(const std::string& rProjectParamsFilePath)
    {
        std::ifstream t(rProjectParamsFilePath);
        std::stringstream buffer;
        buffer << t.rdbuf();
        Parameters projFile{buffer.str()};
        return projFile;
    }

    std::vector<std::shared_ptr<Process>> KratosExecute::parseProcess(ModelPart &model_part, Parameters projFile)
    {
        // Currently: In DGeoflow only fixed hydrostatic head has been , also need load of gravity.

        std::vector<std::shared_ptr<Process>> processes;

        auto constraints_processes = projFile["processes"]["constraints_process_list"];
        for (Parameters process : constraints_processes)
        {
            // we only support fixed hydrostatic head
            auto name = process["Parameters"]["model_part_name"].GetString();
            auto pressure_type = process["Parameters"]["fluid_pressure_type"].GetString();

            std::size_t found = name.find_last_of('.');
            std::string subname = name.substr(found + 1);

            ModelPart &part = model_part.GetSubModelPart(subname);

            if (pressure_type == "Uniform")
            {
                auto value = process["Parameters"]["value"].GetDouble();
                processes.push_back(make_shared<GeoFlowApplyConstantScalarValueProcess>(GeoFlowApplyConstantScalarValueProcess(part, WATER_PRESSURE,
                                                                                                                               value, 0, GeoFlowApplyConstantScalarValueProcess::VARIABLE_IS_FIXED)));
            }
            else if (pressure_type == "Hydrostatic")
            {
                auto cProcesses = process.Clone();
                cProcesses["Parameters"].RemoveValue("fluid_pressure_type");
                processes.push_back(make_shared<GeoFlowApplyConstantHydrostaticPressureProcess>(GeoFlowApplyConstantHydrostaticPressureProcess(part, cProcesses["Parameters"])));
            }
            else
            {
                KRATOS_ERROR << "Reading Processing - Not Implemented - Pressure_type" << std::endl;
            }
        }

        auto loads_processes = projFile["processes"]["loads_process_list"];
        // Should only have one.
        auto name = loads_processes.GetArrayItem(0)["Parameters"]["model_part_name"].GetString();
        std::size_t found = name.find_last_of('.');
        std::string subname = name.substr(found + 1);
        ModelPart &part = model_part.GetSubModelPart(subname);
        processes.push_back(make_shared<ApplyConstantScalarValueProcess>(ApplyConstantScalarValueProcess(part, VOLUME_ACCELERATION_X,
                                                                                                         0.0, 0, ApplyConstantScalarValueProcess::VARIABLE_IS_FIXED)));

        processes.push_back(make_shared<ApplyConstantScalarValueProcess>(ApplyConstantScalarValueProcess(part, VOLUME_ACCELERATION_Y, -9.81,
                                                                                                         0, ApplyConstantScalarValueProcess::VARIABLE_IS_FIXED)));

        processes.push_back(make_shared<Process>(ApplyConstantScalarValueProcess(part, VOLUME_ACCELERATION_Z, 0.0,
                                                                                 0, ApplyConstantScalarValueProcess::VARIABLE_IS_FIXED)));

        return processes;
    }

    int KratosExecute::MainExecution(ModelPart&                                                          rModelPart,
                                     const std::vector<std::shared_ptr<Process>>&                        rProcesses,
                                     const GeoMechanicsNewtonRaphsonErosionProcessStrategyType::Pointer& rpSolvingStrategy,
                                     double                                                              Time,
                                     double                                                              DeltaTime,
                                     unsigned int                                                        NumberOfIterations)
    {

    	// Initialize
        for (const auto& process : rProcesses)
        {
            process->ExecuteInitialize();
        }

        for (const auto& process : rProcesses)
        {
            process->ExecuteBeforeSolutionLoop();
        }

        for (unsigned int iter = 0; iter < NumberOfIterations; ++iter)
        {
            Time += DeltaTime;
            rModelPart.CloneTimeStep(Time);
            rpSolvingStrategy->Initialize();
            rpSolvingStrategy->InitializeSolutionStep();

            for (const auto& process : rProcesses)
            {
                process->ExecuteInitializeSolutionStep();
            }

            rpSolvingStrategy->Predict();
            rpSolvingStrategy->SolveSolutionStep();

            for (const auto& process : rProcesses)
            {
                process->ExecuteFinalizeSolutionStep();
            }

            rpSolvingStrategy->FinalizeSolutionStep();
        }

        for (const auto& process : rProcesses)
        {
            process->ExecuteFinalize();
        }

        return 0;
    }

    int KratosExecute::ExecuteFlowAnalysis(const std::string&                       rWorkingDirectory,
                                           const std::string&                       rProjectParamsFileName,
                                           double                                   minCriticalHead,
                                           double                                   maxCriticalHead,
                                           double                                   stepCriticalHead,
                                           const std::string&                       rCriticalHeadBoundaryModelPartName,
                                           const std::function<void(const char*)>&  rLogCallback,
                                           const std::function<void(double)>&       rReportProgress,
                                           const std::function<void(const char*)>&  rReportTextualProgress,
                                           const std::function<bool()>&             rShouldCancel)
    {
        this->SetEchoLevel(1);

        std::stringstream kratosLogBuffer;
        LoggerOutput::Pointer p_output(new LoggerOutput(kratosLogBuffer));
        Logger::AddOutput(p_output);
        
        try
        {
            rReportProgress(0.0);

            std::string projectpath = rWorkingDirectory + "/" + rProjectParamsFileName;
            auto projectfile = openProjectParamsFile(projectpath);

            auto materialname = projectfile["solver_settings"]["material_import_settings"]["materials_filename"].GetString();
            std::string materialpath = rWorkingDirectory + "/" + materialname;

            auto modelName = projectfile["solver_settings"]["model_part_name"].GetString();

            ModelPart &model_part = current_model.CreateModelPart(modelName);
            model_part.SetBufferSize(2);

            KRATOS_INFO_IF("GeoFlowKernel", this->GetEchoLevel() > 0) << "Working Directory: " << rWorkingDirectory << std::endl;
            KRATOS_INFO_IF("GeoFlowKernel", this->GetEchoLevel() > 0) << "Project Name: " << rProjectParamsFileName << std::endl;

            const auto p_solving_strategy = setup_strategy_dgeoflow(model_part);
            p_solving_strategy->SetEchoLevel(0);

            KRATOS_INFO_IF("GeoFlowKernel", this->GetEchoLevel() > 0) << "Setup Solving Strategy" << std::endl;

            model_part.AddNodalSolutionStepVariable(VELOCITY);
            model_part.AddNodalSolutionStepVariable(ACCELERATION);

            // Displacement
            model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
            model_part.AddNodalSolutionStepVariable(TOTAL_DISPLACEMENT);
            model_part.AddNodalSolutionStepVariable(REACTION);
            model_part.AddNodalSolutionStepVariable(POINT_LOAD);
            model_part.AddNodalSolutionStepVariable(LINE_LOAD);
            model_part.AddNodalSolutionStepVariable(SURFACE_LOAD);
            model_part.AddNodalSolutionStepVariable(VOLUME_ACCELERATION);
            model_part.AddNodalSolutionStepVariable(NORMAL_CONTACT_STRESS);
            model_part.AddNodalSolutionStepVariable(TANGENTIAL_CONTACT_STRESS);

            // Water
            model_part.AddNodalSolutionStepVariable(WATER_PRESSURE);
            model_part.AddNodalSolutionStepVariable(REACTION_WATER_PRESSURE);
            model_part.AddNodalSolutionStepVariable(DT_WATER_PRESSURE);
            model_part.AddNodalSolutionStepVariable(NORMAL_FLUID_FLUX);
            model_part.AddNodalSolutionStepVariable(HYDRAULIC_DISCHARGE);

            // Smoothing
            model_part.AddNodalSolutionStepVariable(NODAL_AREA);
            model_part.AddNodalSolutionStepVariable(NODAL_CAUCHY_STRESS_TENSOR);
            model_part.AddNodalSolutionStepVariable(NODAL_DAMAGE_VARIABLE);
            model_part.AddNodalSolutionStepVariable(NODAL_JOINT_AREA);
            model_part.AddNodalSolutionStepVariable(NODAL_JOINT_WIDTH);
            model_part.AddNodalSolutionStepVariable(NODAL_JOINT_DAMAGE);

            KRATOS_INFO_IF("GeoFlowKernel", this->GetEchoLevel() > 0) << "Nodal Solution Variables Added" << std::endl;

            // Don't include the file extension of the mesh file name, since that is automatically appended by the
            // constructor of class ModelPartIO
            const auto mesh_file_name = projectfile["solver_settings"]["model_import_settings"]["input_filename"].GetString();
            const auto mesh_file_path = rWorkingDirectory + "/" + mesh_file_name;
            ModelPartIO reader{mesh_file_path};
            reader.ReadModelPart(model_part);

            KRATOS_INFO_IF("GeoFlowKernel", this->GetEchoLevel() > 0) << "Parsed Mesh" << std::endl;

            parseMaterial(current_model, materialpath);

            KRATOS_INFO_IF("GeoFlowKernel", this->GetEchoLevel() > 0) << "Parsed Material" << std::endl;

            // Dofs for Water Pressure
            VariableUtils().AddDofWithReaction(WATER_PRESSURE, REACTION_WATER_PRESSURE, model_part);
            VariableUtils().AddDof(VOLUME_ACCELERATION_X, model_part);
            VariableUtils().AddDof(VOLUME_ACCELERATION_Y, model_part);
            VariableUtils().AddDof(VOLUME_ACCELERATION_Z, model_part);

            KRATOS_INFO_IF("GeoFlowKernel", this->GetEchoLevel() > 0) << "Added DoF" << std::endl;

            std::vector<std::shared_ptr<Process>> processes = parseProcess(model_part, projectfile);

            KRATOS_INFO_IF("GeoFlowKernel", this->GetEchoLevel() > 0) << "Parsed Process Data" << std::endl;

            bool hasPiping = stepCriticalHead != 0;

            if (rShouldCancel())
            {
                rLogCallback(kratosLogBuffer.str().c_str());
                Logger::RemoveOutput(p_output);
                ResetModelParts();
                return 0;
            }

            if (!hasPiping)
            {
                MainExecution(model_part, processes, p_solving_strategy, 0.0, 1.0, 1);
                GeoOutputWriter::WriteGiDOutput(model_part, projectfile, rWorkingDirectory);
            }
            else
            {
                KRATOS_INFO_IF("GeoFlowKernel", this->GetEchoLevel() > 0) << "Critical head search started." << std::endl;
                KRATOS_INFO_IF("GeoFlowKernel", this->GetEchoLevel() > 0) << "Critical head min head: " << minCriticalHead << std::endl;
                KRATOS_INFO_IF("GeoFlowKernel", this->GetEchoLevel() > 0) << "Critical head max head: " << maxCriticalHead << std::endl;
                KRATOS_INFO_IF("GeoFlowKernel", this->GetEchoLevel() > 0) << "Critical head step size: " << stepCriticalHead << std::endl;

                shared_ptr<Process> RiverBoundary;
                if (rCriticalHeadBoundaryModelPartName.empty())
                {
                    RiverBoundary = FindRiverBoundaryAutomatically(p_solving_strategy, processes);
                }
                else
                {
                    RiverBoundary = FindRiverBoundaryByName(rCriticalHeadBoundaryModelPartName, processes);
                }

                if (!RiverBoundary)
                {
                    throw std::logic_error("No river boundary found.");
                }

                double criticalHead;
                double currentHead;
                bool pipingSuccess = false;

                auto currentProcess = std::static_pointer_cast<GeoFlowApplyConstantHydrostaticPressureProcess>(RiverBoundary);
                KRATOS_INFO_IF("GeoFlowKernel", this->GetEchoLevel() > 0) << "River boundary name: " << currentProcess->GetName() << std::endl;

                currentProcess->SetReferenceCoord(minCriticalHead);
                currentHead = minCriticalHead;
                criticalHead = currentHead;

                std::vector<Element *> pipeElements;
                pipeElements = p_solving_strategy->GetPipingElements();
                const auto noPipeElements = pipeElements.size();

                int step = 1;
                const auto maxSteps = static_cast<int>(std::ceil((maxCriticalHead - minCriticalHead) / stepCriticalHead)) + 2;

                while (true)
                {
                    if (maxCriticalHead - criticalHead < -1e-9)
                    {
                        KRATOS_INFO_IF("GeoFlowKernel", this->GetEchoLevel() > 0) << "Critical head undetermined at " << criticalHead << ", max search head reached: " << maxCriticalHead << std::endl;
                        break;
                    }

                    KRATOS_INFO_IF("GeoFlowKernel", this->GetEchoLevel() > 0) << "Searching at head: " << currentHead << std::endl;

                    std::ostringstream currentHeadStream;
                    currentHeadStream << std::setprecision(8) << std::noshowpoint << currentHead;
                    std::string currentHeadString = currentHeadStream.str();

                    std::string progress = "Calculating head level " + currentHeadString + "m (" + std::to_string(step) + "/" + std::to_string(maxSteps) + ")";
                    rReportTextualProgress(progress.data());
                    rReportProgress(((double)step) / ((double)maxSteps));

                    MainExecution(model_part, processes, p_solving_strategy, 0.0, 1.0, 1);

                    auto count = std::size_t{0};
                    for (Element *element : pipeElements)
                    {
                        if (element->GetValue(PIPE_ACTIVE))
                            count += 1;
                    }

                    KRATOS_INFO_IF("GeoFlowKernel", this->GetEchoLevel() > 0) << "Open pipe elements: " << count << std::endl;

                    if (count == noPipeElements)
                    {
                        if (abs(currentHead - minCriticalHead) < 1e-9)
                        {
                            KRATOS_INFO_IF("GeoFlowKernel", this->GetEchoLevel() > 0) << "Critical head undetermined: All pipe elements open at initial search value :" << minCriticalHead << std::endl;
                        }
                        else
                        {
                            pipingSuccess = true;
                            KRATOS_INFO_IF("GeoFlowKernel", this->GetEchoLevel() > 0) << "Critical head found: " << criticalHead << std::endl;
                        }
                        break;
                    }

                    GeoOutputWriter::WriteGiDOutput(model_part, projectfile, rWorkingDirectory);

                    // Update boundary conditions for next search head.
                    if (RiverBoundary->Info() == "ApplyConstantScalarValueProcess")
                    {
                        ResetModelParts();
                        throw std::logic_error("ApplyConstantScalarValueProcess process search is not Implemented");
                    }

                    if (RiverBoundary->Info() == "ApplyConstantHydrostaticPressureProcess")
                    {
                        criticalHead = currentProcess->GetReferenceCoord();
                        currentHead = criticalHead + stepCriticalHead;
                        currentProcess->SetReferenceCoord(currentHead);
                        step++;
                    }

                    if (rShouldCancel())
                    {
                        rLogCallback(kratosLogBuffer.str().c_str());
                        Logger::RemoveOutput(p_output);
                        ResetModelParts();
                        return 0;
                    }
                }

                KRATOS_INFO_IF("GeoFlowKernel", this->GetEchoLevel() > 0) << "Writing result to: " << rWorkingDirectory << "\\criticalHead.json" << std::endl;

                // output critical head_json
                std::ofstream CriticalHeadFile(rWorkingDirectory + "\\criticalHead.json");

                CriticalHeadFile << "{\n";
                CriticalHeadFile << "\t \"PipeData\":\t{\n";
                if (pipingSuccess)
                {
                    CriticalHeadFile << "\t\t \"Success\": \"True\",\n";
                    CriticalHeadFile << "\t\t \"CriticalHead\": \"" + std::to_string(criticalHead) + "\"\n";
                }
                else
                {
                    CriticalHeadFile << "\t\t \"Success\": \"False\"\n";
                }
                CriticalHeadFile << "\t }\n";
                CriticalHeadFile << "}\n";

                // Close the file
                CriticalHeadFile.close();

                KRATOS_INFO_IF("GeoFlowKernel", this->GetEchoLevel() > 0) << "Finished writing result" << std::endl;
            }

            rLogCallback(kratosLogBuffer.str().c_str());
            Logger::RemoveOutput(p_output);

            ResetModelParts();
            return 0;
        }
        catch (const std::exception &exc)
        {
            KRATOS_INFO_IF("GeoFlowKernel", this->GetEchoLevel() > 0) << exc.what();

            rLogCallback(kratosLogBuffer.str().c_str());
            Logger::RemoveOutput(p_output);

            ResetModelParts();
            return 1;
        }
    }

    shared_ptr<Process> KratosExecute::FindRiverBoundaryByName(const std::string& rCriticalHeadBoundaryModelPartName,
                                                               const std::vector<std::shared_ptr<Process>>& rProcesses)
    {
        shared_ptr<Process> RiverBoundary;

        for (const auto& process : rProcesses)
        {
            if (process->Info() == "ApplyConstantHydrostaticPressureProcess")
            {
                auto current_process = std::static_pointer_cast<GeoFlowApplyConstantHydrostaticPressureProcess>(process);
                if (current_process->hasWaterPressure() &&
                    (current_process->GetName() == rCriticalHeadBoundaryModelPartName))
                {
                    RiverBoundary = current_process;
                }
            }
        }

        if (!RiverBoundary)
        {
            std::cerr << "No boundary found with the model part name " << rCriticalHeadBoundaryModelPartName << "." << std::endl;
            return nullptr;
        }

        return RiverBoundary;
    }

    shared_ptr<Process> KratosExecute::FindRiverBoundaryAutomatically(const KratosExecute::GeoMechanicsNewtonRaphsonErosionProcessStrategyType::Pointer& rpSolvingStrategy,
                                                                      const std::vector<std::shared_ptr<Process>>& rProcesses)
    {
        shared_ptr<Process> RiverBoundary;

        std::vector<Element *> pipeElements;
        pipeElements = rpSolvingStrategy->GetPipingElements();

        double firstNode_A = pipeElements.front()->GetGeometry().GetPoint(0).X0();
        double firstNode_B = pipeElements.front()->GetGeometry().GetPoint(1).X0();
        double lastNode_A = pipeElements.back()->GetGeometry().GetPoint(0).X0();

        IndexType RiverNode;

        if ((firstNode_A < lastNode_A) && (firstNode_A < firstNode_B))
        {
            RiverNode = pipeElements.back()->GetGeometry().GetPoint(1).Id();
        }
        else
        {
            RiverNode = pipeElements.back()->GetGeometry().GetPoint(0).Id();
        }

        // Get Find boundary in Processes
        for (const auto& process : rProcesses)
        {
            ModelPart *currentModelPart;

            if (process->Info() == "ApplyConstantScalarValueProcess")
            {
                auto current_process = std::static_pointer_cast<GeoFlowApplyConstantScalarValueProcess>(process);
                if (current_process->hasWaterPressure())
                {
                    currentModelPart = &current_process->GetModelPart();
                    try
                    {
                        currentModelPart->GetNode(RiverNode);
                        RiverBoundary = current_process;
                    }
                    catch (...)
                    {
                    }
                }
            }
            else if (process->Info() == "ApplyConstantHydrostaticPressureProcess")
            {
                auto current_process = std::static_pointer_cast<GeoFlowApplyConstantHydrostaticPressureProcess>(process);
                currentModelPart = &current_process->GetModelPart();
                if (current_process->hasWaterPressure())
                {
                    try
                    {
                        currentModelPart->GetNode(RiverNode);
                        RiverBoundary = current_process;
                    }
                    catch (...)
                    {
                    }
                }
            }
        }

        if (!RiverBoundary)
        {
            std::cerr << "No boundary found on the river side at node " << RiverNode << "." << std::endl;
            return nullptr;
        }

        return RiverBoundary;
    }
}
