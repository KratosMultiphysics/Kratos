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

// System includes
#include <limits>
#include <map>

/* External includes */

/* Utility includes */
#include "includes/define.h"
#include "containers/model.h"
#include "includes/model_part.h"
#include "spaces/ublas_space.h"

// Linear solvers
#include "linear_solvers/linear_solver.h"

// The most basic scheme (static)
#include "custom_strategies/schemes/backward_euler_quasistatic_Pw_scheme.hpp"

// The most extended criteria (residual criteria)
#include "solving_strategies/convergencecriterias/residual_criteria.h"

// The most builder and solver (the block builder and solver)
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"

// The strategies to test
#include <custom_processes/apply_component_table_process.hpp>
#include <custom_processes/apply_constant_hydrostatic_pressure_process.hpp>
#include <factories/linear_solver_factory.h>
#include <ghc/filesystem.hpp>
#include <includes/gid_io.h>

#include <processes/apply_constant_scalarvalue_process.h>
#include <solving_strategies/convergencecriterias/mixed_generic_criteria.h>
#include <solving_strategies/strategies/implicit_solving_strategy.h>
#include <solving_strategies/strategies/residualbased_newton_raphson_strategy.h>

#include "custom_strategies/strategies/geo_mechanics_newton_raphson_erosion_process_strategy.hpp"

#include "utilities/variable_utils.h"
#include "utilities/read_materials_utility.h"

using namespace std;

#pragma region NodeVariables

class NodeOperation
{
public:
    virtual void write(Kratos::GidIO<>& gid_io, Kratos::ModelPart& model_part) = 0;
};
class NodeDISPLACEMENT : public NodeOperation {
public:
    void write(Kratos::GidIO<>& gid_io, Kratos::ModelPart& model_part) { gid_io.WriteNodalResults(Kratos::DISPLACEMENT, model_part.Nodes(), 0, 0); }
};
class NodeTOTAL_DISPLACEMENT : public NodeOperation {
public:
    void write(Kratos::GidIO<>& gid_io, Kratos::ModelPart& model_part) { gid_io.WriteNodalResults(Kratos::TOTAL_DISPLACEMENT, model_part.Nodes(), 0, 0); }
};
class NodeWATER_PRESSURE : public NodeOperation {
public:
    void write(Kratos::GidIO<>& gid_io, Kratos::ModelPart& model_part) { gid_io.WriteNodalResults(Kratos::WATER_PRESSURE, model_part.Nodes(), 0, 0); }
};
class NodeNORMAL_FLUID_FLUX : public NodeOperation {
public:
    void write(Kratos::GidIO<>& gid_io, Kratos::ModelPart& model_part) { gid_io.WriteNodalResults(Kratos::NORMAL_FLUID_FLUX, model_part.Nodes(), 0, 0); }
};
class NodeVOLUME_ACCELERATION : public NodeOperation {
public:
    void write(Kratos::GidIO<>& gid_io, Kratos::ModelPart& model_part) { gid_io.WriteNodalResults(Kratos::VOLUME_ACCELERATION, model_part.Nodes(), 0, 0); }
};
class NodeHYDRAULIC_DISCHARGE : public NodeOperation {
public:
    void write(Kratos::GidIO<>& gid_io, Kratos::ModelPart& model_part) { gid_io.WriteNodalResults(Kratos::HYDRAULIC_DISCHARGE, model_part.Nodes(), 0, 0); }
};
#pragma endregion NodeVariables

#pragma region GaussVariables
class GaussOperation
{
public:
    virtual void write(Kratos::GidIO<>& gid_io, Kratos::ModelPart& model_part) = 0;
};
class GaussFLUID_FLUX_VECTOR : public GaussOperation {
public:
    void write(Kratos::GidIO<>& gid_io, Kratos::ModelPart& model_part) { gid_io.PrintOnGaussPoints(Kratos::FLUID_FLUX_VECTOR, model_part, 0, 0); }
};
class GaussHYDRAULIC_HEAD : public GaussOperation {
public:
    void write(Kratos::GidIO<>& gid_io, Kratos::ModelPart& model_part) { gid_io.PrintOnGaussPoints(Kratos::HYDRAULIC_HEAD, model_part, 0, 0); }
};
class GaussLOCAL_FLUID_FLUX_VECTOR : public GaussOperation {
public:
    void write(Kratos::GidIO<>& gid_io, Kratos::ModelPart& model_part) { gid_io.PrintOnGaussPoints(Kratos::LOCAL_FLUID_FLUX_VECTOR, model_part, 0, 0); }
};
class GaussLOCAL_PERMEABILITY_MATRIX : public GaussOperation {
public:
    void write(Kratos::GidIO<>& gid_io, Kratos::ModelPart& model_part) { gid_io.PrintOnGaussPoints(Kratos::LOCAL_PERMEABILITY_MATRIX, model_part, 0, 0); }
};
#pragma endregion GaussVariables

namespace Kratos
{

    typedef Node<3> NodeType;
    typedef Geometry<NodeType> GeometryType;
    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

    // The direct solver
    typedef LinearSolver<SparseSpaceType, LocalSpaceType> LinearSolverType;

    // The convergence criteria type
    typedef ConvergenceCriteria< SparseSpaceType, LocalSpaceType > ConvergenceCriteriaType;
    typedef MixedGenericCriteria<SparseSpaceType, LocalSpaceType> MixedGenericCriteriaType;
    typedef typename MixedGenericCriteriaType::ConvergenceVariableListType ConvergenceVariableListType;

    // The builder ans solver type
    typedef BuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > BuilderAndSolverType;
    typedef ResidualBasedBlockBuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > ResidualBasedBlockBuilderAndSolverType;

    // The time scheme
    typedef Scheme< SparseSpaceType, LocalSpaceType >  SchemeType;
    typedef BackwardEulerQuasistaticPwScheme< SparseSpaceType, LocalSpaceType > BackwardEulerQuasistaticPwSchemeType;

    // The strategies
    typedef ImplicitSolvingStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >
        ImplicitSolvingStrategyType;

    typedef ResidualBasedNewtonRaphsonStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >
        ResidualBasedNewtonRaphsonStrategyType;

    typedef GeoMechanicsNewtonRaphsonErosionProcessStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >
        GeoMechanicsNewtonRaphsonErosionProcessStrategyType;

    // Dof arrays
    typedef PointerVectorSet<Dof<double>, SetIdentityFunction<Dof<double>>, std::less<SetIdentityFunction<Dof<double>>::result_type>,
        std::equal_to<SetIdentityFunction<Dof<double>>::result_type>, Dof<double>* > DofsArrayType;

	ConvergenceCriteriaType::Pointer setup_criteria_dgeoflow()
    {
        const double rel_tol = 1.0e-4;
        const double abs_tol = 1.0e-9;
        VariableData* p_water_pres = &WATER_PRESSURE;
        ConvergenceVariableListType convergence_settings;
        convergence_settings.push_back(std::make_tuple(p_water_pres, rel_tol, abs_tol));
        auto mixed_generic_criteria = MixedGenericCriteriaType(convergence_settings);
        return ConvergenceCriteriaType::Pointer(new MixedGenericCriteriaType(convergence_settings));
    }

    LinearSolverType::Pointer setup_solver_dgeoflow()
    {
        Parameters linear_solver_settings(R"({"solver_type": "skyline_lu_factorization"})");
        LinearSolverFactory<SparseSpaceType, LocalSpaceType> linear_solver_factory;
        return linear_solver_factory.Create(linear_solver_settings);
    }

    GeoMechanicsNewtonRaphsonErosionProcessStrategyType::Pointer setup_strategy_dgeoflow(ModelPart& model_part)
    {
        // Create the linear strategy
        auto p_solver = setup_solver_dgeoflow();
        Scheme<SparseSpaceType, LocalSpaceType>::Pointer p_scheme = Kratos::make_shared<BackwardEulerQuasistaticPwScheme<SparseSpaceType, LocalSpaceType>>();
        auto p_builder_and_solver = Kratos::make_shared<ResidualBasedBlockBuilderAndSolver<SparseSpaceType, LocalSpaceType, LinearSolverType>>(p_solver);

        Parameters p_parameters(R"(
            {
                "min_iteration":    6,
                "number_cycles":    100,
                "increase_factor":  2.0,
                "reduction_factor": 0.5,
                "end_time": 1.0,
                "realised_factor": 1.0,
				
				"max_piping_iterations": 50,

                "desired_iterations": 4,
                "max_radius_factor": 20.0,
                "min_radius_factor": 0.5,
                "search_neighbours_step": false,
                "body_domain_sub_model_part_list": [],
                "loads_sub_model_part_list": [],
                "loads_variable_list" : []
            }  )");

        int MaxIterations = 15;
        bool CalculateReactions = false;
        bool ReformDofSetAtEachStep = true;
        bool MoveMeshFlag = false;

        auto p_solving_strategy = Kratos::make_unique<GeoMechanicsNewtonRaphsonErosionProcessStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType>>(
            model_part,
            p_scheme,
            p_solver,
            setup_criteria_dgeoflow(),
            p_builder_and_solver,
            p_parameters,
            MaxIterations, CalculateReactions, ReformDofSetAtEachStep, MoveMeshFlag);

        p_solving_strategy->Check();
        return p_solving_strategy;
    }

    void parseMesh(ModelPart& model_part, string filepath)
    {
        // Parses MDPA file into model_part
        std::ifstream input(filepath);
        bool read_properties = false;
        bool read_nodes = false;
        bool read_elements = false;

        bool read_subparts = false;
        bool read_subparts_table = false;
        bool read_subparts_nodes = false;
        bool read_subparts_elements = false;
        bool read_subparts_conditions = false;

        string element_type;
        string part_name;
        string nodeStr;

        for (std::string line; getline(input, line); )
        {

            //===================== Properties =========================
            if (line.substr(0, 16) == "Begin Properties")
            {
                read_properties = true;
                std::size_t found = line.find_last_of(" ");
                int property_id = stoi(line.substr(found + 1));
                model_part.CreateNewProperties(property_id);
                continue;
            }
            if (line == "End Properties")
            {
                read_properties = false;
                continue;
            }
            if (read_properties)
            {
                KRATOS_ERROR << "Reading Properties - Not Implemented " << std::endl;
            }
            //=====================   Nodes   =========================
            if (line == "Begin Nodes")
            {
                read_nodes = true;
                continue;
            }
            if (line == "End Nodes")
            {
                read_nodes = false;
                continue;
            }
            if (read_nodes)
            {
                std::istringstream iss(line);
                int nodeId;
                double x, y, z;
                iss >> nodeId >> x >> y >> z;
                model_part.CreateNewNode(nodeId, x, y, z);
            }
            //====================   Element   ==========================
            if (line.substr(0, 14) == "Begin Elements")
            {
                read_elements = true;
                std::size_t found = line.find_last_of(" ");
                element_type = line.substr(found + 1);

                std::size_t posD = element_type.find_last_of("D");
                nodeStr = element_type.substr(posD + 1);

                continue;
            }
            if (line == "End Elements")
            {
                read_elements = false;
                continue;
            }

            if (read_elements)
            {
                unsigned long elementId, propertyId, node1, node2, node3;
                std::vector<ModelPart::IndexType> element_nodes;
                std::istringstream iss(line);

                if (nodeStr == "4N")
                {
                    unsigned long node4;
                    iss >> elementId >> propertyId >> node1 >> node2 >> node3 >> node4;
                    element_nodes = { node1, node2, node3, node4 };
                }
                else if (nodeStr == "3N")
                {
                    iss >> elementId >> propertyId >> node1 >> node2 >> node3;
                    element_nodes = { node1, node2, node3 };
                }
                else
                {
                    KRATOS_ERROR << "Element Type Unknown / Not Implemented " << std::endl;
                }
                auto p_elem_prop = model_part.pGetProperties(propertyId);
                model_part.CreateNewElement(element_type, elementId, element_nodes, p_elem_prop);
            }

            //===================== Properties =========================

            if (line.substr(0, 18) == "Begin SubModelPart")
            {
                read_subparts = true;
                std::size_t found = line.find_last_of(" ");
                part_name = line.substr(found + 1);
                model_part.CreateSubModelPart(part_name);
                continue;
            }
            if (line == "End SubModelPart")
            {
                read_subparts = false;
                continue;
            }
            if (read_subparts)
            {
                auto subpart = model_part.pGetSubModelPart(part_name);
                //===========  Sub-Tables  ===============
                if (line == "  Begin SubModelPartTables")
                {
                    read_subparts_table = true;
                    continue;
                }
                if (line == "  End SubModelPartTables")
                {
                    read_subparts_table = false;
                    continue;
                }
                if (read_subparts_table)
                {
                    KRATOS_ERROR << "Subpart Tables - Not Implemented " << std::endl;
                }

                //===========  Sub-Nodes  ===============

                if (line == "  Begin SubModelPartNodes")
                {
                    read_subparts_nodes = true;
                    continue;
                }
                if (line == "  End SubModelPartNodes")
                {
                    read_subparts_nodes = false;
                    continue;
                }
                if (read_subparts_nodes)
                {
                    auto node = model_part.pGetNode(stoi(line));
                    subpart->AddNode(node);
                }

                //===========  Sub-Elements  ===============

                if (line == "  Begin SubModelPartElements")
                {
                    read_subparts_elements = true;
                    continue;
                }
                if (line == "  End SubModelPartElements")
                {
                    read_subparts_elements = false;
                    continue;
                }
                if (read_subparts_elements)
                {
                    auto element = model_part.pGetElement(stoi(line));
                    subpart->AddElement(element);
                }

                //===========  Sub-Elements  ===============

                if (line == "  Begin SubModelPartConditions")
                {
                    read_subparts_conditions = true;
                    continue;
                }
                if (line == "  End SubModelPartConditions")
                {
                    read_subparts_conditions = false;
                    continue;
                }
                if (read_subparts_conditions)
                {
                    KRATOS_ERROR << "Subpart Conditions - Not Implemented " << std::endl;
                }
            }
        }
    }

    void parseMaterial(Model& model, string filepath)
    {

        string parameters = "{ \"Parameters\" : { \"materials_filename\" :\"" + filepath + "\"}}";
        Parameters material_file{ parameters };
        ReadMaterialsUtility(material_file, model);

    }

    Parameters openProjectParamsFile(string filepath)
    {
        std::ifstream t(filepath);
        std::stringstream buffer;
        buffer << t.rdbuf();
        Parameters projFile{ buffer.str() };
        return projFile;
    }

    std::vector<std::shared_ptr<Process>> parseProcess(ModelPart& model_part, Parameters projFile)
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
            string subname = name.substr(found + 1);

            ModelPart& part = model_part.GetSubModelPart(subname);

            if (pressure_type == "Uniform")
            {
                auto value = process["Parameters"]["value"].GetDouble();
                processes.push_back(make_shared<ApplyConstantScalarValueProcess>(ApplyConstantScalarValueProcess(part, WATER_PRESSURE,
                    value, 0, ApplyConstantScalarValueProcess::VARIABLE_IS_FIXED)));
            }
            else if (pressure_type == "Hydrostatic")
            {
                auto cProcesses = process.Clone();
                cProcesses["Parameters"].RemoveValue("fluid_pressure_type");
                processes.push_back(make_shared<ApplyConstantHydrostaticPressureProcess>(ApplyConstantHydrostaticPressureProcess(part, cProcesses["Parameters"])));
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
        string subname = name.substr(found + 1);
        ModelPart& part = model_part.GetSubModelPart(subname);
        processes.push_back(make_shared<ApplyConstantScalarValueProcess>(ApplyConstantScalarValueProcess(part, VOLUME_ACCELERATION_X,
            0.0, 0, ApplyConstantScalarValueProcess::VARIABLE_IS_FIXED)));

        processes.push_back(make_shared<ApplyConstantScalarValueProcess>(ApplyConstantScalarValueProcess(part, VOLUME_ACCELERATION_Y, 10.0,
            0, ApplyConstantScalarValueProcess::VARIABLE_IS_FIXED)));

        processes.push_back(make_shared<Process>(ApplyConstantScalarValueProcess(part, VOLUME_ACCELERATION_Z, 0.0,
            0, ApplyConstantScalarValueProcess::VARIABLE_IS_FIXED)));

        return processes;

    }

    void outputGiD(Model& model, ModelPart& model_part, Parameters parameters)
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
        string filename = outputParameters["output_name"].GetString();
        GiD_PostMode gid_output_type = PostMode[outputParameters["postprocess_parameters"]["result_file_configuration"]["gidpost_flags"]["GiDPostMode"].GetString()];
        MultiFileFlag multifiles_output = MultiFiles[outputParameters["postprocess_parameters"]["result_file_configuration"]["gidpost_flags"]["MultiFileFlag"].GetString()];
        WriteDeformedMeshFlag deformed_output = DeformedFlag[outputParameters["postprocess_parameters"]["result_file_configuration"]["gidpost_flags"]["WriteDeformedMeshFlag"].GetString()];
        WriteConditionsFlag condition_output = ConditionFlag[outputParameters["postprocess_parameters"]["result_file_configuration"]["gidpost_flags"]["WriteConditionsFlag"].GetString()];

        GidIO<> gid_io(filename, gid_output_type, multifiles_output, deformed_output, condition_output);

        gid_io.InitializeMesh(0.0);
        gid_io.WriteMesh(model_part.GetMesh());
        gid_io.FinalizeMesh();
        gid_io.InitializeResults(0, model_part.GetMesh());

        std::unordered_map<string, unique_ptr<NodeOperation>> NodeOutput;
        NodeOutput["DISPLACEMENT"] = std::make_unique<NodeDISPLACEMENT>();
        NodeOutput["TOTAL_DISPLACEMENT"] = std::make_unique<NodeTOTAL_DISPLACEMENT>();
        NodeOutput["WATER_PRESSURE"] = std::make_unique<NodeWATER_PRESSURE>();
        NodeOutput["NORMAL_FLUID_FLUX"] = std::make_unique<NodeNORMAL_FLUID_FLUX>();
        NodeOutput["VOLUME_ACCELERATION"] = std::make_unique<NodeVOLUME_ACCELERATION>();
        NodeOutput["HYDRAULIC_DISCHARGE"] = std::make_unique<NodeHYDRAULIC_DISCHARGE>();

        auto nodal_outputs = outputParameters["postprocess_parameters"]["result_file_configuration"]["nodal_results"].GetStringArray();
        for (string var : nodal_outputs)
        {
            NodeOutput[var]->write(gid_io, model_part);
        }

        std::unordered_map<string, unique_ptr<GaussOperation>> GaussOutput;
        GaussOutput["FLUID_FLUX_VECTOR"] = std::make_unique<GaussFLUID_FLUX_VECTOR>();
        GaussOutput["HYDRAULIC_HEAD"] = std::make_unique<GaussHYDRAULIC_HEAD>();
        GaussOutput["LOCAL_FLUID_FLUX_VECTOR"] = std::make_unique<GaussLOCAL_FLUID_FLUX_VECTOR>();
        GaussOutput["LOCAL_PERMEABILITY_MATRIX"] = std::make_unique<GaussLOCAL_PERMEABILITY_MATRIX>();
        auto gauss_outputs = outputParameters["postprocess_parameters"]["result_file_configuration"]["gauss_point_results"].GetStringArray();
        for (string var : gauss_outputs)
        {
            GaussOutput[var]->write(gid_io, model_part);
        }

        gid_io.FinalizeResults();
    }

    int cpp_geomechanics(string meshpath, string projectpath, string materialpath)
    {
        auto projectfile = openProjectParamsFile(projectpath);
        auto modelName = projectfile["solver_settings"]["model_part_name"].GetString();


        // Initial Setup
        Model current_model;
        constexpr double tolerance = 1e-6;
        ModelPart& model_part = current_model.CreateModelPart(modelName);
        model_part.SetBufferSize(2);

        // ReSharper disable once CppMsExtBindingRValueToLvalueReference
        auto& p_solving_strategy = setup_strategy_dgeoflow(model_part);
        p_solving_strategy->SetEchoLevel(0);


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

        parseMesh(model_part, meshpath);

        parseMaterial(current_model, materialpath);

        // Dofs for Water Pressure
        VariableUtils().AddDofWithReaction(WATER_PRESSURE, REACTION_WATER_PRESSURE, model_part);
        VariableUtils().AddDof(VOLUME_ACCELERATION_X, model_part);
        VariableUtils().AddDof(VOLUME_ACCELERATION_Y, model_part);
        VariableUtils().AddDof(VOLUME_ACCELERATION_Z, model_part);

        std::vector<std::shared_ptr<Process>> processes = parseProcess(model_part, projectfile);

        // Initialize
        for (auto process : processes) {
            process->ExecuteInitialize();
        }

        for (auto process : processes) {
            process->ExecuteBeforeSolutionLoop();
        }

        double time = 0.0;
        const double delta_time = 1.0;
        const unsigned int number_iterations = 1;
        for (unsigned int iter = 0; iter < number_iterations; ++iter) {
            time += delta_time;
            model_part.CloneTimeStep(time);
            p_solving_strategy->Initialize();
            p_solving_strategy->InitializeSolutionStep();

            for (auto process : processes) {
                process->ExecuteInitializeSolutionStep();
            }

            p_solving_strategy->Predict();
            p_solving_strategy->SolveSolutionStep();

            for (auto process : processes) {
                process->ExecuteFinalizeSolutionStep();
            }

            p_solving_strategy->FinalizeSolutionStep();
        }
        for (auto process : processes) {
            process->ExecuteFinalize();
        }

        outputGiD(current_model, model_part, projectfile);
        return 0;
    }

}
