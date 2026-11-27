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

// System includes
#include <map>

/* Project includes */
#include "containers/model.h"
#include "custom_conditions/geo_seepage_condition.h"
#include "custom_strategies/schemes/backward_euler_quasistatic_Pw_scheme.hpp"
#include "custom_workflows/dgeoflow.h"
#include "geo_mechanics_application_variables.h"
#include "includes/kratos_parameters.h"
#include "tests/cpp_tests/flow_stubs.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"

namespace Kratos::Testing
{

using namespace Kratos;

namespace
{

class SeepageTestElement : public Element
{
public:
    using Element::Element;

    void GetDofList(DofsVectorType& rDofs, const ProcessInfo&) const override
    {
        rDofs = {GetGeometry()[0].pGetDof(WATER_PRESSURE),
                 GetGeometry()[1].pGetDof(WATER_PRESSURE)};
    }

    void EquationIdVector(EquationIdVectorType& rIds, const ProcessInfo&) const override
    {
        rIds = {GetGeometry()[0].GetDof(WATER_PRESSURE).EquationId(),
                GetGeometry()[1].GetDof(WATER_PRESSURE).EquationId()};
    }

    void CalculateRightHandSide(VectorType& rRhs, const ProcessInfo&) override
    {
        rRhs.resize(2, false);
        for (auto i = std::size_t{0}; i < 2; ++i) {
            rRhs[i] = GetValue(PRESSURE) + GetGeometry()[i].FastGetSolutionStepValue(WATER_PRESSURE);
        }
    }

    void CalculateLocalSystem(MatrixType& rLhs, VectorType& rRhs, const ProcessInfo& rProcessInfo) override
    {
        rLhs = -IdentityMatrix(2);
        CalculateRightHandSide(rRhs, rProcessInfo);
    }
};

class TestableDGeoFlowStrategy : public KratosExecute::GeoMechanicsNewtonRaphsonErosionProcessStrategyType
{
public:
    using BaseType = KratosExecute::GeoMechanicsNewtonRaphsonErosionProcessStrategyType;
    using BaseType::BaseType;
    using BaseType::Recalculate;
};

} // namespace

class ParametrizedDGeoFlowTests : public ::testing::TestWithParam<std::string>
{
};

TEST_P(ParametrizedDGeoFlowTests, ErosionProcessStrategy)
{
    // Arrange
    auto working_directory =
        "./applications/GeoMechanicsApplication/tests/test_compare_sellmeijer/" + GetParam();
    auto project_file = "ProjectParameters.json";

    auto execute = KratosExecute();

    std::vector<double> progress_report_values;
    auto                reportProgress = [&progress_report_values](double Progress) {
        progress_report_values.push_back(Progress);
    };

    std::vector<std::string> reported_textual_progress;
    auto report_textual_progress = [&reported_textual_progress](const char* pMessage) {
        reported_textual_progress.emplace_back(pMessage);
    };

    const KratosExecute::CriticalHeadInfo  critical_head_info(3, 4, 0.1);
    const KratosExecute::CallBackFunctions call_back_functions(
        &flow_stubs::emptyLog, reportProgress, report_textual_progress, &flow_stubs::emptyCancel);

    // Act
    const auto status = execute.ExecuteFlowAnalysis(working_directory, project_file, critical_head_info,
                                                    "PorousDomain.Left_head", call_back_functions);

    // Assert
    KRATOS_EXPECT_EQ(status, 0);
    KRATOS_EXPECT_EQ(reported_textual_progress.front(), "Calculating head level 3m (1/12)");
    KRATOS_EXPECT_EQ(reported_textual_progress.back(), "Calculating head level 3.8m (9/12)");
    KRATOS_EXPECT_EQ(reported_textual_progress.size(), 9);
    KRATOS_EXPECT_DOUBLE_EQ(progress_report_values.front(), 0.0);
    KRATOS_EXPECT_DOUBLE_EQ(progress_report_values.back(), 0.75);
    KRATOS_EXPECT_EQ(progress_report_values.size(), 10);
}

KRATOS_TEST_CASE_IN_SUITE(DGeoFlowParsesSeepageAndNormalFluxLoads, KratosGeoMechanicsFastSuite)
{
    Model model;
    auto& r_model_part = model.CreateModelPart("PorousDomain");
    r_model_part.AddNodalSolutionStepVariable(NORMAL_FLUID_FLUX);
    r_model_part.AddNodalSolutionStepVariable(VOLUME_ACCELERATION_X);
    r_model_part.AddNodalSolutionStepVariable(VOLUME_ACCELERATION_Y);
    r_model_part.AddNodalSolutionStepVariable(VOLUME_ACCELERATION_Z);
    auto& r_boundaries = r_model_part.CreateSubModelPart("boundaries");
    r_boundaries.CreateSubModelPart("flux_boundary");
    r_model_part.CreateSubModelPart("seepage_boundary");
    r_model_part.CreateSubModelPart("gravity");

    Parameters project_parameters(R"(
    {
        "processes": {
            "constraints_process_list": [],
            "loads_process_list": [
                {
                    "process_name": "ApplyScalarConstraintTableProcess",
                    "Parameters": {
                        "model_part_name": "PorousDomain.boundaries.flux_boundary",
                        "variable_name": "NORMAL_FLUID_FLUX",
                        "value": -5.0,
                        "table": 1
                    }
                },
                {
                    "process_name": "ApplySeepageBoundaryProcess",
                    "Parameters": {
                        "model_part_name": "PorousDomain.seepage_boundary"
                    }
                },
                {
                    "process_name": "ApplyVectorConstraintTableProcess",
                    "Parameters": {
                        "model_part_name": "PorousDomain.gravity",
                        "variable_name": "VOLUME_ACCELERATION"
                    }
                }
            ]
        }
    })");

    KratosExecute execute;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(execute.ParseProcesses(r_model_part, project_parameters),
                                      "DGeoFlow only supports constant NORMAL_FLUID_FLUX values.");

    project_parameters["processes"]["loads_process_list"][0]["Parameters"]["table"].SetInt(0);
    execute.ParseProcesses(r_model_part, project_parameters);
}

KRATOS_TEST_CASE_IN_SUITE(DGeoFlowParsesSeepageConstraints, KratosGeoMechanicsFastSuite)
{
    Model model;
    auto& r_model_part = model.CreateModelPart("PorousDomain");
    r_model_part.AddNodalSolutionStepVariable(WATER_PRESSURE);
    r_model_part.AddNodalSolutionStepVariable(VOLUME_ACCELERATION_X);
    r_model_part.AddNodalSolutionStepVariable(VOLUME_ACCELERATION_Y);
    r_model_part.AddNodalSolutionStepVariable(VOLUME_ACCELERATION_Z);
    r_model_part.CreateSubModelPart("fixed_head");
    r_model_part.CreateSubModelPart("gravity");

    Parameters project_parameters(R"(
    {
        "processes": {
            "constraints_process_list": [
                {
                    "Parameters": {
                        "model_part_name": "PorousDomain.fixed_head",
                        "fluid_pressure_type": "Uniform",
                        "value": 0.0
                    }
                },
                {
                    "process_name": "ApplySeepageBoundaryProcess",
                    "Parameters": {
                        "model_part_name": "PorousDomain.seepage_boundary"
                    }
                }
            ],
            "loads_process_list": [
                {
                    "process_name": "ApplyVectorConstraintTableProcess",
                    "Parameters": {
                        "model_part_name": "PorousDomain.gravity",
                        "variable_name": "VOLUME_ACCELERATION"
                    }
                }
            ]
        }
    })");

    KratosExecute execute;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        execute.ParseProcesses(r_model_part, project_parameters),
        "ApplySeepageBoundaryProcess: \"PorousDomain.seepage_boundary\" does not correspond to any model part");

    r_model_part.CreateSubModelPart("seepage_boundary");
    KratosExecute valid_execute;
    valid_execute.ParseProcesses(r_model_part, project_parameters);
}

KRATOS_TEST_CASE_IN_SUITE(DGeoFlowSwitchesSeepageDuringSolveAndErosionRecalculation, KratosGeoMechanicsFastSuite)
{
    Model model;
    auto& r_model_part = model.CreateModelPart("PorousDomain", 2);
    KratosExecute::AddNodalSolutionStepVariables(r_model_part);
    auto p_first_node = r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    auto p_second_node = r_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    for (auto& r_node : r_model_part.Nodes()) {
        r_node.AddDof(WATER_PRESSURE, REACTION_WATER_PRESSURE);
    }
    auto p_geometry = make_shared<Line2D2<Node>>(p_first_node, p_second_node);
    auto p_properties = r_model_part.CreateNewProperties(0);
    auto p_element = make_intrusive<SeepageTestElement>(1, p_geometry, p_properties);
    p_element->SetValue(PRESSURE, 1.0);
    r_model_part.AddElement(p_element);
    r_model_part.AddCondition(make_intrusive<GeoSeepageCondition>(1, p_geometry, p_properties));
    r_model_part.CloneTimeStep(1.0);

    using SparseSpaceType = KratosExecute::SparseSpaceType;
    using LocalSpaceType = KratosExecute::LocalSpaceType;
    using LinearSolverType = KratosExecute::LinearSolverType;
    auto settings = Parameters{R"({"number_cycles":1,"max_piping_iterations":2})"};
    TestableDGeoFlowStrategy strategy{
        r_model_part,
        make_shared<BackwardEulerQuasistaticPwScheme<SparseSpaceType, LocalSpaceType>>(),
        KratosExecute::setup_criteria_dgeoflow(),
        make_shared<ResidualBasedBlockBuilderAndSolver<SparseSpaceType, LocalSpaceType, LinearSolverType>>(
            KratosExecute::setup_solver_dgeoflow()),
        settings, 15, true, false, false};
    strategy.SetEchoLevel(0);
    strategy.Initialize();
    strategy.InitializeSolutionStep();
    strategy.Predict();

    KRATOS_EXPECT_TRUE(strategy.SolveSolutionStep());
    for (const auto& r_node : r_model_part.Nodes()) {
        KRATOS_EXPECT_TRUE(r_node.IsFixed(WATER_PRESSURE));
        KRATOS_EXPECT_NEAR(r_node.FastGetSolutionStepValue(WATER_PRESSURE), 0.0, 1e-12);
    }

    // Reversing the flow must release the seepage nodes inside the erosion recalculation.
    p_element->SetValue(PRESSURE, -1.0);
    KRATOS_EXPECT_TRUE(strategy.Recalculate());
    for (const auto& r_node : r_model_part.Nodes()) {
        KRATOS_EXPECT_FALSE(r_node.IsFixed(WATER_PRESSURE));
        KRATOS_EXPECT_NEAR(r_node.FastGetSolutionStepValue(WATER_PRESSURE), 1.0, 1e-12);
    }
    strategy.FinalizeSolutionStep();
}

INSTANTIATE_TEST_SUITE_P(KratosGeoMechanicsIntegrationSuite,
                         ParametrizedDGeoFlowTests,
                         ::testing::Values("HeightAquiferD10L30line"));

} // namespace Kratos::Testing
