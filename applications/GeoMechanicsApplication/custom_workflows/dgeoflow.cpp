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


#include "dgeoflow.h"
#include "processes/apply_constant_scalarvalue_process.h"
#include "utilities/read_materials_utility.h"
#include "input_output/logger.h"
#include "input_output/logger_output.h"
#include "input_output/logger_table_output.h"
#include "dgeooutput.h"
#include "dgeoparser.h"

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

    double GetProcessDoubleValue()
    {
        return mdouble_value;
    }

    void SetProcessDoubleValue(double value)
    {
        mdouble_value = value;
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
    KratosGeoFlow::KratosGeoFlow()
	    : KratosGeoApplication()
    {
    }

    KratosGeoFlow::ConvergenceCriteriaType::Pointer KratosGeoFlow::convergence_criteria()
    {
        const double rel_tol = 1.0e-4;
        const double abs_tol = 1.0e-9;
        VariableData *p_water_pres = &WATER_PRESSURE;
        ConvergenceVariableListType convergence_settings;
        convergence_settings.push_back(std::make_tuple(p_water_pres, rel_tol, abs_tol));
        return ConvergenceCriteriaType::Pointer(new MixedGenericCriteriaType(convergence_settings));
    }

    KratosGeoFlow::LinearSolverType::Pointer KratosGeoFlow::solver_settings()
    {
        // Parameters linear_solver_settings(R"({"solver_type": "sparse_lu"})");
        // return linear_solver_factory.Create(linear_solver_settings);
        // LinearSolverType::Pointer p_solver = LinearSolverFactoryType().Create(linear_solver_settings);
        LinearSolverType::Pointer p_solver = Kratos::make_shared<SkylineLUFactorizationSolverType>();
        // LinearSolverType::Pointer p_solver = Kratos::make_shared<EigenSparseLUSolverType>();
        return p_solver;
    }

    KratosGeoFlow::ImplicitSolvingStrategyType::Pointer KratosGeoFlow::strategy_settings(ModelPart &model_part)
    {
        // Create the linear strategy
        auto p_solver = solver_settings();

        Scheme<SparseSpaceType, LocalSpaceType>::Pointer p_scheme = Kratos::make_shared<BackwardEulerQuasistaticPwScheme<SparseSpaceType, LocalSpaceType>>();

        auto p_builder_and_solver = Kratos::make_shared<ResidualBasedBlockBuilderAndSolver<SparseSpaceType, LocalSpaceType, KratosGeoFlow::LinearSolverType>>(p_solver);
        p_builder_and_solver->SetEchoLevel(0);

        auto p_criteria = convergence_criteria();
        p_criteria->SetEchoLevel(0);

        Parameters p_parameters(R"(
    {
        "min_iteration":    6,
        "number_cycles":    100,
        "increase_factor":  2.0,
        "reduction_factor": 0.5,
        "end_time": 1.0,
        "realised_factor": 1.0,
		
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

        auto p_solving_strategy = Kratos::make_unique<GeoMechanicsNewtonRaphsonErosionProcessStrategy<SparseSpaceType, LocalSpaceType, KratosGeoFlow::LinearSolverType>>(
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

    int KratosGeoFlow::execute_application(std::string workingDirectory, std::string projectName,
                                           double minCriticalHead, double maxCriticalHead, double stepCriticalHead,
                                           std::string criticalHeadBoundaryModelPartName,
                                           std::function<void(char *)> logCallback,
                                           std::function<void(double)> reportProgress,
                                           std::function<void(char *)> reportTextualProgress,
                                           std::function<bool()> shouldCancel)
    {
        this->SetEchoLevel(1);

        std::stringstream kratosLogBuffer;
        LoggerOutput::Pointer p_output(new LoggerOutput(kratosLogBuffer));
        Logger::AddOutput(p_output);
        
        try
        {
            reportProgress(0.0);

            std::string projectpath = workingDirectory + "/" + projectName;
            auto projectfile = KratosGeoParser::openProjectParamsFile(projectpath);

            auto materialname = projectfile["solver_settings"]["material_import_settings"]["materials_filename"].GetString();
            auto meshname = projectfile["solver_settings"]["model_import_settings"]["input_filename"].GetString() + "." +
                            projectfile["solver_settings"]["model_import_settings"]["input_type"].GetString();

            std::string meshpath = workingDirectory + "/" + meshname;
            std::string materialpath = workingDirectory + "/" + materialname;

            auto modelName = projectfile["solver_settings"]["model_part_name"].GetString();

            Model& current_model = GetModelPointer();
        	ModelPart &model_part = current_model.CreateModelPart(modelName);
            model_part.SetBufferSize(2);

            KRATOS_INFO_IF("GeoFlowKernel", this->GetEchoLevel() > 0) << "Working Directory: " << workingDirectory << std::endl;
            KRATOS_INFO_IF("GeoFlowKernel", this->GetEchoLevel() > 0) << "Project Name: " << projectName << std::endl;

            const auto p_solving_strategy = strategy_settings(model_part);
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

            KratosGeoParser::parseMesh(model_part, meshpath);

            KRATOS_INFO_IF("GeoFlowKernel", this->GetEchoLevel() > 0) << "Parsed Mesh" << std::endl;

            KratosGeoParser::parseMaterial(current_model, materialpath);

            KRATOS_INFO_IF("GeoFlowKernel", this->GetEchoLevel() > 0) << "Parsed Material" << std::endl;

            // Dofs for Water Pressure
            VariableUtils().AddDofWithReaction(WATER_PRESSURE, REACTION_WATER_PRESSURE, model_part);
            VariableUtils().AddDof(VOLUME_ACCELERATION_X, model_part);
            VariableUtils().AddDof(VOLUME_ACCELERATION_Y, model_part);
            VariableUtils().AddDof(VOLUME_ACCELERATION_Z, model_part);

            KRATOS_INFO_IF("GeoFlowKernel", this->GetEchoLevel() > 0) << "Added DoF" << std::endl;

            std::vector<std::shared_ptr<Process>> processes = KratosGeoParser::parseProcess(model_part, projectfile);

            KRATOS_INFO_IF("GeoFlowKernel", this->GetEchoLevel() > 0) << "Parsed Process Data" << std::endl;

            bool hasPiping = stepCriticalHead != 0;

            if (shouldCancel())
            {
                logCallback(strdup(kratosLogBuffer.str().c_str()));
                Logger::RemoveOutput(p_output);
                ResetModelParts();
                return 0;
            }

            if (!hasPiping)
            {
                execute_kratos_calculation(model_part, processes, p_solving_strategy, 0.0, 1.0, 1);
                KratosGeoOutput::outputGiD(current_model, model_part, projectfile, workingDirectory);
            }
            else
            {
                KRATOS_INFO_IF("GeoFlowKernel", this->GetEchoLevel() > 0) << "Critical head search started." << std::endl;
                KRATOS_INFO_IF("GeoFlowKernel", this->GetEchoLevel() > 0) << "Critical head min head: " << minCriticalHead << std::endl;
                KRATOS_INFO_IF("GeoFlowKernel", this->GetEchoLevel() > 0) << "Critical head max head: " << maxCriticalHead << std::endl;
                KRATOS_INFO_IF("GeoFlowKernel", this->GetEchoLevel() > 0) << "Critical head step size: " << stepCriticalHead << std::endl;

                shared_ptr<Process> RiverBoundary;
            	RiverBoundary = FindRiverBoundaryByName(criticalHeadBoundaryModelPartName, processes);

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
                pipeElements = std::static_pointer_cast<GeoMechanicsNewtonRaphsonErosionProcessStrategyType>(p_solving_strategy)->GetPipingElements();
                int noPipeElements = pipeElements.size();

                int step = 1;
                int maxSteps = std::ceil((maxCriticalHead - minCriticalHead) / stepCriticalHead) + 2;

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
                    reportTextualProgress(progress.data());
                    reportProgress(((double)step) / ((double)maxSteps));

                    execute_kratos_calculation(model_part, processes, p_solving_strategy, 0.0, 1.0, 1);

                    int count = 0;
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

                    KratosGeoOutput::outputGiD(current_model, model_part, projectfile, workingDirectory);

                    // Update boundary conditions for next search head.
                    if (RiverBoundary->Info() == "ApplyConstantScalarValueProcess")
                    {
                        ResetModelParts();
                        throw std::logic_error("ApplyConstantScalarValueProcess process search is not Implemented");
                    }

                    if (RiverBoundary->Info() == "ApplyConstantHydrostaticPressureProcess")
                    {
                        auto currentProcess = std::static_pointer_cast<GeoFlowApplyConstantHydrostaticPressureProcess>(RiverBoundary);
                        criticalHead = currentProcess->GetReferenceCoord();
                        currentHead = criticalHead + stepCriticalHead;
                        currentProcess->SetReferenceCoord(currentHead);
                        step++;
                    }

                    if (shouldCancel())
                    {
                        logCallback(strdup(kratosLogBuffer.str().c_str()));
                        Logger::RemoveOutput(p_output);
                        ResetModelParts();
                        return 0;
                    }
                }

                KRATOS_INFO_IF("GeoFlowKernel", this->GetEchoLevel() > 0) << "Writing result to: " << workingDirectory << "\\criticalHead.json" << std::endl;

                // output critical head_json
                std::ofstream CriticalHeadFile(workingDirectory + "\\criticalHead.json");

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

            logCallback(strdup(kratosLogBuffer.str().c_str()));
            Logger::RemoveOutput(p_output);

            ResetModelParts();
            return 0;
        }
        catch (const std::exception &exc)
        {
            KRATOS_INFO_IF("GeoFlowKernel", this->GetEchoLevel() > 0) << exc.what();

            logCallback(strdup(kratosLogBuffer.str().c_str()));
            Logger::RemoveOutput(p_output);

            ResetModelParts();
            return 1;
        }
    };

    shared_ptr<Process> KratosGeoFlow::FindRiverBoundaryByName(std::string criticalHeadBoundaryModelPartName,
                                                               std::vector<std::shared_ptr<Process>> processes)
    {
        shared_ptr<Process> RiverBoundary;

        for (shared_ptr<Process> process : processes)
        {
            if (process->Info() == "ApplyConstantHydrostaticPressureProcess")
            {
                auto current_process = std::static_pointer_cast<GeoFlowApplyConstantHydrostaticPressureProcess>(process);
                if (current_process->hasWaterPressure())
                {
                    if (current_process->GetName() == criticalHeadBoundaryModelPartName)
                    {
                        RiverBoundary = current_process;
                    }
                }
            }
        }

        if (!RiverBoundary)
        {
            std::cerr << "No boundary found with the model part name " << criticalHeadBoundaryModelPartName << "." << std::endl;
            return NULL;
        }

        return RiverBoundary;
    }

}
