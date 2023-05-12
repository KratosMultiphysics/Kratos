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

void NodeOperation::write(Kratos::GidIO<> &gid_io, Kratos::ModelPart &model_part){};

void NodeDISPLACEMENT::write(Kratos::GidIO<> &gid_io, Kratos::ModelPart &model_part) { gid_io.WriteNodalResults(Kratos::DISPLACEMENT, model_part.Nodes(), 0, 0); }

void NodeTOTAL_DISPLACEMENT::write(Kratos::GidIO<> &gid_io, Kratos::ModelPart &model_part) { gid_io.WriteNodalResults(Kratos::TOTAL_DISPLACEMENT, model_part.Nodes(), 0, 0); }

void NodeWATER_PRESSURE::write(Kratos::GidIO<> &gid_io, Kratos::ModelPart &model_part) { gid_io.WriteNodalResults(Kratos::WATER_PRESSURE, model_part.Nodes(), 0, 0); }

void NodeNORMAL_FLUID_FLUX::write(Kratos::GidIO<> &gid_io, Kratos::ModelPart &model_part) { gid_io.WriteNodalResults(Kratos::NORMAL_FLUID_FLUX, model_part.Nodes(), 0, 0); }

void NodeVOLUME_ACCELERATION::write(Kratos::GidIO<> &gid_io, Kratos::ModelPart &model_part) { gid_io.WriteNodalResults(Kratos::VOLUME_ACCELERATION, model_part.Nodes(), 0, 0); }

void NodeHYDRAULIC_DISCHARGE::write(Kratos::GidIO<> &gid_io, Kratos::ModelPart &model_part) { gid_io.WriteNodalResults(Kratos::HYDRAULIC_DISCHARGE, model_part.Nodes(), 0, 0); }

void NodeHYDRAULIC_HEAD::write(Kratos::GidIO<> &gid_io, Kratos::ModelPart &model_part) { gid_io.WriteNodalResults(Kratos::HYDRAULIC_HEAD, model_part.Nodes(), 0, 0); }

void GaussOperation::write(Kratos::GidIO<> &gid_io, Kratos::ModelPart &model_part){};

void GaussFLUID_FLUX_VECTOR::write(Kratos::GidIO<> &gid_io, Kratos::ModelPart &model_part) { gid_io.PrintOnGaussPoints(Kratos::FLUID_FLUX_VECTOR, model_part, 0, 0); }

void GaussHYDRAULIC_HEAD::write(Kratos::GidIO<> &gid_io, Kratos::ModelPart &model_part) { gid_io.PrintOnGaussPoints(Kratos::HYDRAULIC_HEAD, model_part, 0, 0); }

void GaussLOCAL_FLUID_FLUX_VECTOR::write(Kratos::GidIO<> &gid_io, Kratos::ModelPart &model_part) { gid_io.PrintOnGaussPoints(Kratos::LOCAL_FLUID_FLUX_VECTOR, model_part, 0, 0); }

void GaussLOCAL_PERMEABILITY_MATRIX::write(Kratos::GidIO<> &gid_io, Kratos::ModelPart &model_part) { gid_io.PrintOnGaussPoints(Kratos::LOCAL_PERMEABILITY_MATRIX, model_part, 0, 0); }

void GaussPERMEABILITY_MATRIX::write(Kratos::GidIO<> &gid_io, Kratos::ModelPart &model_part) { gid_io.PrintOnGaussPoints(Kratos::PERMEABILITY_MATRIX, model_part, 0, 0); }

void GaussDEGREE_OF_SATURATION::write(Kratos::GidIO<> &gid_io, Kratos::ModelPart &model_part) { gid_io.PrintOnGaussPoints(Kratos::DEGREE_OF_SATURATION, model_part, 0, 0); }

void GaussDERIVATIVE_OF_SATURATION::write(Kratos::GidIO<> &gid_io, Kratos::ModelPart &model_part) { gid_io.PrintOnGaussPoints(Kratos::DERIVATIVE_OF_SATURATION, model_part, 0, 0); }

void GaussRELATIVE_PERMEABILITY::write(Kratos::GidIO<> &gid_io, Kratos::ModelPart &model_part) { gid_io.PrintOnGaussPoints(Kratos::RELATIVE_PERMEABILITY, model_part, 0, 0); }

void GaussPIPE_ACTIVE::write(Kratos::GidIO<> &gid_io, Kratos::ModelPart &model_part) { gid_io.PrintOnGaussPoints(Kratos::PIPE_ACTIVE, model_part, 0, 0); }

void GaussPIPE_HEIGHT::write(Kratos::GidIO<> &gid_io, Kratos::ModelPart &model_part) { gid_io.PrintOnGaussPoints(Kratos::PIPE_HEIGHT, model_part, 0, 0); }

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

    int KratosExecute::GetEchoLevel()
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
        KratosExecute::ConvergenceVariableListType convergence_settings;
        convergence_settings.push_back(std::make_tuple(p_water_pres, rel_tol, abs_tol));
        return KratosExecute::ConvergenceCriteriaType::Pointer(new KratosExecute::MixedGenericCriteriaType(convergence_settings));
    }

    KratosExecute::LinearSolverType::Pointer KratosExecute::setup_solver_dgeoflow()
    {
        // Parameters linear_solver_settings(R"({"solver_type": "sparse_lu"})");
        // return linear_solver_factory.Create(linear_solver_settings);
        // LinearSolverType::Pointer p_solver = LinearSolverFactoryType().Create(linear_solver_settings);
        LinearSolverType::Pointer p_solver = Kratos::make_shared<SkylineLUFactorizationSolverType>();
        // LinearSolverType::Pointer p_solver = Kratos::make_shared<EigenSparseLUSolverType>();
        return p_solver;
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

    void KratosExecute::parseMaterial(Model &model, std::string filepath)
    {
        std::string parameters = "{ \"Parameters\" : { \"materials_filename\" :\"" + filepath + "\"}}";
        Parameters material_file{parameters};
        ReadMaterialsUtility(material_file, model);
    }

    Parameters KratosExecute::openProjectParamsFile(std::string filepath)
    {
        std::ifstream t(filepath);
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

            std::size_t found = name.find_last_of(".");
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
        std::size_t found = name.find_last_of(".");
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

    void KratosExecute::outputGiD(Model &model, ModelPart &model_part, Parameters parameters, std::string workingDirectory)
    {
        std::map<std::string, GiD_PostMode> PostMode;
        PostMode["GiD_PostAscii"] = GiD_PostAscii;
        PostMode["GiD_PostAsciiZipped"] = GiD_PostAsciiZipped;
        PostMode["GiD_PostBinary"] = GiD_PostBinary;
        PostMode["GiD_PostHDF5"] = GiD_PostHDF5;

        std::map<std::string, MultiFileFlag> MultiFiles;
        MultiFiles["SingleFile"] = SingleFile;
        MultiFiles["MultipleFiles"] = MultipleFiles;

        std::map<std::string, WriteDeformedMeshFlag> DeformedFlag;
        DeformedFlag["WriteDeformed"] = WriteDeformed;
        DeformedFlag["WriteUndeformed"] = WriteUndeformed;

        std::map<std::string, WriteConditionsFlag> ConditionFlag;
        ConditionFlag["WriteConditions"] = WriteConditions;
        ConditionFlag["WriteElementsOnly"] = WriteElementsOnly;
        ConditionFlag["WriteConditionsOnly"] = WriteConditionsOnly;

        Parameters gid_out = parameters["output_processes"]["gid_output"].GetArrayItem(0);
        Parameters outputParameters = gid_out["Parameters"];
        std::string filename = outputParameters["output_name"].GetString();
        GiD_PostMode gid_output_type = PostMode[outputParameters["postprocess_parameters"]["result_file_configuration"]["gidpost_flags"]["GiDPostMode"].GetString()];
        MultiFileFlag multifiles_output = MultiFiles[outputParameters["postprocess_parameters"]["result_file_configuration"]["gidpost_flags"]["MultiFileFlag"].GetString()];
        WriteDeformedMeshFlag deformed_output = DeformedFlag[outputParameters["postprocess_parameters"]["result_file_configuration"]["gidpost_flags"]["WriteDeformedMeshFlag"].GetString()];
        WriteConditionsFlag condition_output = ConditionFlag[outputParameters["postprocess_parameters"]["result_file_configuration"]["gidpost_flags"]["WriteConditionsFlag"].GetString()];

        filename = workingDirectory + "/" + filename;
        GidIO<> gid_io(filename, gid_output_type, multifiles_output, deformed_output, condition_output);

        gid_io.InitializeMesh(0.0);
        gid_io.WriteMesh(model_part.GetMesh());
        gid_io.FinalizeMesh();
        gid_io.InitializeResults(0, model_part.GetMesh());

        std::unordered_map<std::string, unique_ptr<NodeOperation>> NodeOutput;
        NodeOutput["DISPLACEMENT"] = Kratos::make_unique<NodeDISPLACEMENT>();
        NodeOutput["TOTAL_DISPLACEMENT"] = Kratos::make_unique<NodeTOTAL_DISPLACEMENT>();
        NodeOutput["WATER_PRESSURE"] = Kratos::make_unique<NodeWATER_PRESSURE>();
        NodeOutput["NORMAL_FLUID_FLUX"] = Kratos::make_unique<NodeNORMAL_FLUID_FLUX>();
        NodeOutput["VOLUME_ACCELERATION"] = Kratos::make_unique<NodeVOLUME_ACCELERATION>();
        NodeOutput["HYDRAULIC_DISCHARGE"] = Kratos::make_unique<NodeHYDRAULIC_DISCHARGE>();
        NodeOutput["HYDRAULIC_HEAD"] = Kratos::make_unique<NodeHYDRAULIC_HEAD>();

        // Calculate hydraulic head on the nodes
        auto gauss_outputs = outputParameters["postprocess_parameters"]["result_file_configuration"]["gauss_point_results"].GetStringArray();
        auto nodal_outputs = outputParameters["postprocess_parameters"]["result_file_configuration"]["nodal_results"].GetStringArray();

        if (std::find(gauss_outputs.begin(), gauss_outputs.end(), "HYDRAULIC_HEAD") != gauss_outputs.end())
        {
            calculateNodalHydraulicHead(gid_io, model_part);
        }

        for (std::string var : nodal_outputs)
        {
            NodeOutput[var]->write(gid_io, model_part);
        }

        std::unordered_map<std::string, std::unique_ptr<GaussOperation>> GaussOutput;
        GaussOutput["FLUID_FLUX_VECTOR"] = Kratos::make_unique<GaussFLUID_FLUX_VECTOR>();
        GaussOutput["HYDRAULIC_HEAD"] = Kratos::make_unique<GaussHYDRAULIC_HEAD>();
        GaussOutput["LOCAL_FLUID_FLUX_VECTOR"] = Kratos::make_unique<GaussLOCAL_FLUID_FLUX_VECTOR>();
        GaussOutput["LOCAL_PERMEABILITY_MATRIX"] = Kratos::make_unique<GaussLOCAL_PERMEABILITY_MATRIX>();
        GaussOutput["PERMEABILITY_MATRIX"] = Kratos::make_unique<GaussPERMEABILITY_MATRIX>();
        GaussOutput["DEGREE_OF_SATURATION"] = Kratos::make_unique<GaussDEGREE_OF_SATURATION>();
        GaussOutput["DERIVATIVE_OF_SATURATION"] = Kratos::make_unique<GaussDERIVATIVE_OF_SATURATION>();
        GaussOutput["RELATIVE_PERMEABILITY"] = Kratos::make_unique<GaussRELATIVE_PERMEABILITY>();
        GaussOutput["PIPE_ACTIVE"] = Kratos::make_unique<GaussPIPE_ACTIVE>();
        GaussOutput["PIPE_HEIGHT"] = Kratos::make_unique<GaussPIPE_HEIGHT>();

        // Now Output Gauss Point Results on Gauss Points
        for (std::string var : gauss_outputs)
        {
            GaussOutput[var]->write(gid_io, model_part);
        }

        gid_io.FinalizeResults();
    }

    void KratosExecute::calculateNodalHydraulicHead(GidIO<> &gid_io, ModelPart &model_part) {
            const auto& element_var = KratosComponents<Variable<double>>::Get("HYDRAULIC_HEAD");

            for (Element element : model_part.Elements())
            {
                auto& rGeom = element.GetGeometry();
                const auto& rProp = element.GetProperties();
                
                const auto NodalHydraulicHead = GeoElementUtilities::CalculateNodalHydraulicHeadFromWaterPressures<3>(rGeom, rProp);

            	for (unsigned int node = 0; node < 3; ++node)
                {
                    rGeom[node].SetValue(element_var, NodalHydraulicHead[node]);
                }
            }
            gid_io.WriteNodalResultsNonHistorical(element_var, model_part.Nodes(), 0);
    }

    int KratosExecute::mainExecution(ModelPart &model_part,
                                     std::vector<std::shared_ptr<Process>> processes,
                                     GeoMechanicsNewtonRaphsonErosionProcessStrategyType::Pointer p_solving_strategy,
                                     double time, double delta_time, double number_iterations)
    {

    	// Initialize
        for (auto process : processes)
        {
            process->ExecuteInitialize();
        }

        for (auto process : processes)
        {
            process->ExecuteBeforeSolutionLoop();
        }

        for (unsigned int iter = 0; iter < number_iterations; ++iter)
        {
            time += delta_time;
            model_part.CloneTimeStep(time);
            p_solving_strategy->Initialize();
            p_solving_strategy->InitializeSolutionStep();

            for (auto process : processes)
            {
                process->ExecuteInitializeSolutionStep();
            }

            p_solving_strategy->Predict();
            p_solving_strategy->SolveSolutionStep();

            for (auto process : processes)
            {
                process->ExecuteFinalizeSolutionStep();
            }

            p_solving_strategy->FinalizeSolutionStep();
        }

        for (auto process : processes)
        {
            process->ExecuteFinalize();
        }

        return 0;
    }

    int KratosExecute::execute_flow_analysis(std::string workingDirectory, std::string projectName,
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
            auto projectfile = openProjectParamsFile(projectpath);

            auto materialname = projectfile["solver_settings"]["material_import_settings"]["materials_filename"].GetString();
            std::string materialpath = workingDirectory + "/" + materialname;

            auto modelName = projectfile["solver_settings"]["model_part_name"].GetString();

            ModelPart &model_part = current_model.CreateModelPart(modelName);
            model_part.SetBufferSize(2);

            KRATOS_INFO_IF("GeoFlowKernel", this->GetEchoLevel() > 0) << "Working Directory: " << workingDirectory << std::endl;
            KRATOS_INFO_IF("GeoFlowKernel", this->GetEchoLevel() > 0) << "Project Name: " << projectName << std::endl;

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
            const auto mesh_file_path = workingDirectory + "/" + mesh_file_name;
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

            if (shouldCancel())
            {
                logCallback(strdup(kratosLogBuffer.str().c_str()));
                Logger::RemoveOutput(p_output);
                ResetModelParts();
                return 0;
            }

            if (!hasPiping)
            {
                mainExecution(model_part, processes, p_solving_strategy, 0.0, 1.0, 1);
                outputGiD(current_model, model_part, projectfile, workingDirectory);
            }
            else
            {
                KRATOS_INFO_IF("GeoFlowKernel", this->GetEchoLevel() > 0) << "Critical head search started." << std::endl;
                KRATOS_INFO_IF("GeoFlowKernel", this->GetEchoLevel() > 0) << "Critical head min head: " << minCriticalHead << std::endl;
                KRATOS_INFO_IF("GeoFlowKernel", this->GetEchoLevel() > 0) << "Critical head max head: " << maxCriticalHead << std::endl;
                KRATOS_INFO_IF("GeoFlowKernel", this->GetEchoLevel() > 0) << "Critical head step size: " << stepCriticalHead << std::endl;

                shared_ptr<Process> RiverBoundary;
                if (criticalHeadBoundaryModelPartName.empty())
                {
                    RiverBoundary = FindRiverBoundaryAutomatically(p_solving_strategy, processes);
                }
                else
                {
                    RiverBoundary = FindRiverBoundaryByName(criticalHeadBoundaryModelPartName, processes);
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

                    mainExecution(model_part, processes, p_solving_strategy, 0.0, 1.0, 1);

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

                    outputGiD(current_model, model_part, projectfile, workingDirectory);

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

    shared_ptr<Process> KratosExecute::FindRiverBoundaryByName(std::string criticalHeadBoundaryModelPartName,
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

    shared_ptr<Process> KratosExecute::FindRiverBoundaryAutomatically(KratosExecute::GeoMechanicsNewtonRaphsonErosionProcessStrategyType::Pointer p_solving_strategy,
                                                                      std::vector<std::shared_ptr<Process>> processes)
    {
        shared_ptr<Process> RiverBoundary;

        std::vector<Element *> pipeElements;
        pipeElements = p_solving_strategy->GetPipingElements();

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
        for (shared_ptr<Process> process : processes)
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
            return NULL;
        }

        return RiverBoundary;
    }
}
