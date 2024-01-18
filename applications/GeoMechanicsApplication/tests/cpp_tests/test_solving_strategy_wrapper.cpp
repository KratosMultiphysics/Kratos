// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Richard Faasse
//

#include "testing/testing.h"
#include "containers/model.h"
#include "spaces/ublas_space.h"
#include "custom_workflows/solving_strategy_wrapper.hpp"
#include "linear_solvers/linear_solver.h"
#include "custom_utilities/solving_strategy_factory.hpp"
#include "geo_mechanics_application_variables.h"

namespace
{

const std::string testParameters = R"(
    {
        "solver_type":                        "U_Pw",
        "model_part_name":                    "PorousDomain",
        "domain_size":                        2,
        "model_import_settings":              {
            "input_type":       "mdpa",
            "input_filename":   "test_model"
        },
        "material_import_settings":              {
            "materials_filename":       "MaterialParameters.json"
        },
        "time_stepping":              {
            "time_step":                1.0,
            "max_delta_time_factor":    1000
        },
        "buffer_size":                        2,
        "echo_level":                         1,
        "clear_storage":                      false,
        "compute_reactions":                  true,
        "move_mesh_flag":                     false,
        "reform_dofs_at_each_step":           false,
        "nodal_smoothing":                    false,
        "block_builder":                      true,
        "solution_type":                      "Quasi-Static",
        "scheme_type":                        "Backward_Euler",
        "reset_displacements":                false,
        "newmark_beta":                       0.25,
        "newmark_gamma":                      0.5,
        "newmark_theta":                      0.5,
        "rayleigh_m":                         0.0,
        "rayleigh_k":                         0.0,
        "strategy_type":                      "newton_raphson",
        "convergence_criterion":              "displacement_criterion",
        "displacement_relative_tolerance":    1.0E-4,
        "displacement_absolute_tolerance":    1.0E-9,
        "residual_relative_tolerance":        1.0E-4,
        "residual_absolute_tolerance":        1.0E-9,
        "water_pressure_relative_tolerance":  1.0E-4,
        "water_pressure_absolute_tolerance":  1.0E-9,
        "min_iterations":                     6,
        "max_iterations":                     15,
        "number_cycles":                      100,
        "reduction_factor":                   0.5,
        "increase_factor":                    2.0,
        "desired_iterations":                 4,
        "max_radius_factor":                  10.0,
        "min_radius_factor":                  0.1,
        "calculate_reactions":                true,
        "max_line_search_iterations":         5,
        "first_alpha_value":                  0.5,
        "second_alpha_value":                 1.0,
        "min_alpha":                          0.1,
        "max_alpha":                          2.0,
        "line_search_tolerance":              0.5,
        "rotation_dofs":                      true,
        "linear_solver_settings":             {
            "solver_type":     "sparse_lu",
            "scaling":         true
        },
        "problem_domain_sub_model_part_list": ["Clay_after_excavation","Excavated_clay"],
        "processes_sub_model_part_list":      ["Side_sliders","Bottom_fixed","Head_line","Bottom_pressure","Gravitational_load"],
        "body_domain_sub_model_part_list":    ["Clay_after_excavation","Excavated_clay"]
    }
)";

}

using namespace Kratos;

namespace Kratos::Testing
{

using SparseSpaceType = UblasSpace<double, CompressedMatrix, Vector>;
using DenseSpaceType = UblasSpace<double, Matrix, Vector>;
using LinearSolverType = LinearSolver<SparseSpaceType, DenseSpaceType>;
using SolvingStrategyFactoryType = SolvingStrategyFactory<SparseSpaceType, DenseSpaceType, LinearSolverType>;
using SolvingStrategyWrapperType = SolvingStrategyWrapper<SparseSpaceType, DenseSpaceType>;


SolvingStrategyWrapperType CreateWrapperWithDefaultProcessInfoEntries(ModelPart& rModelPart)
{
    auto info = std::make_shared<ProcessInfo>();
    (*info)[NL_ITERATION_NUMBER] = 5;
    (*info)[TIME] = 17.0;
    (*info)[DELTA_TIME] = 3.4;
    (*info)[STEP] = 3;
    rModelPart.SetProcessInfo(info);
    return SolvingStrategyWrapperType{SolvingStrategyFactoryType::Create(Parameters{testParameters}, rModelPart)};
}

ModelPart& CreateDummyModelPart(Model& rModel)
{
    const auto buffer_size = ModelPart::IndexType{2};
    auto& r_model_part = rModel.CreateModelPart("dummy", buffer_size);
    r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    return r_model_part;
}

SolvingStrategyWrapperType CreateWrapperWithEmptyProcessInfo(ModelPart& rModelPart)
{
    const auto reset_displacements = true;
    return SolvingStrategyWrapperType{SolvingStrategyFactoryType::Create(Parameters{testParameters}, rModelPart), reset_displacements};
}

KRATOS_TEST_CASE_IN_SUITE(GetNumberOfIterationsFromStrategyWrapper_ReturnsCorrectNumber, KratosGeoMechanicsFastSuite)
{
    Model model;
    auto& model_part = CreateDummyModelPart(model);
    auto  wrapper    = CreateWrapperWithDefaultProcessInfoEntries(model_part);

    KRATOS_EXPECT_EQ(wrapper.GetNumberOfIterations(), 5);
}

KRATOS_TEST_CASE_IN_SUITE(GetEndTimeFromStrategyWrapper_ReturnsCorrectNumber, KratosGeoMechanicsFastSuite)
{
    Model model;
    auto& model_part = CreateDummyModelPart(model);
    auto  wrapper    = CreateWrapperWithDefaultProcessInfoEntries(model_part);

    KRATOS_EXPECT_DOUBLE_EQ(wrapper.GetEndTime(), 17.0);
}

KRATOS_TEST_CASE_IN_SUITE(SolveSolutionStepFromStrategyWrapper_ReturnsCorrectState, KratosGeoMechanicsFastSuite)
{
    Model model;
    auto& model_part = CreateDummyModelPart(model);
    auto  wrapper    = CreateWrapperWithDefaultProcessInfoEntries(model_part);

    KRATOS_EXPECT_EQ(wrapper.SolveSolutionStep(), TimeStepEndState::ConvergenceState::converged);
}

KRATOS_TEST_CASE_IN_SUITE(SetEndTimeFromStrategyWrapper, KratosGeoMechanicsFastSuite)
{
    Model      model;
    auto&      model_part = CreateDummyModelPart(model);
    auto       wrapper    = CreateWrapperWithEmptyProcessInfo(model_part);
    const auto end_time   = 0.4;

    wrapper.SetEndTime(end_time);

    KRATOS_EXPECT_DOUBLE_EQ(model_part.GetProcessInfo()[TIME], end_time);
}

KRATOS_TEST_CASE_IN_SUITE(GetTimeIncrementFromStrategyWrapper_ReturnsCorrectNumber, KratosGeoMechanicsFastSuite)
{
    Model model;
    auto& model_part = CreateDummyModelPart(model);
    auto  wrapper    = CreateWrapperWithDefaultProcessInfoEntries(model_part);

    KRATOS_EXPECT_DOUBLE_EQ(wrapper.GetTimeIncrement(), 3.4);
}

KRATOS_TEST_CASE_IN_SUITE(SetTimeIncrementFromStrategyWrapper, KratosGeoMechanicsFastSuite)
{
    Model      model;
    auto&      model_part     = CreateDummyModelPart(model);
    auto       wrapper        = CreateWrapperWithEmptyProcessInfo(model_part);
    const auto time_increment = 0.8;

    wrapper.SetTimeIncrement(time_increment);

    KRATOS_EXPECT_EQ(model_part.GetProcessInfo()[DELTA_TIME], time_increment);
}

KRATOS_TEST_CASE_IN_SUITE(GetStepNumberFromStrategyWrapper_ReturnsCorrectNumber, KratosGeoMechanicsFastSuite)
{
    Model model;
    auto& model_part = CreateDummyModelPart(model);
    auto  wrapper    = CreateWrapperWithDefaultProcessInfoEntries(model_part);

    KRATOS_EXPECT_EQ(wrapper.GetStepNumber(), 3);
}

KRATOS_TEST_CASE_IN_SUITE(IncrementStepNumberFromStrategyWrapper, KratosGeoMechanicsFastSuite)
{
    Model model;
    auto& model_part = CreateDummyModelPart(model);
    auto  wrapper    = CreateWrapperWithDefaultProcessInfoEntries(model_part);

    wrapper.IncrementStepNumber();
    KRATOS_EXPECT_EQ(wrapper.GetStepNumber(), 4);
    wrapper.IncrementStepNumber();
    KRATOS_EXPECT_EQ(wrapper.GetStepNumber(), 5);
}

KRATOS_TEST_CASE_IN_SUITE(SaveAndAccumulateTotalDisplacementField, KratosGeoMechanicsFastSuite)
{
    Model model;
    auto& model_part       = CreateDummyModelPart(model);
    auto  strategy_wrapper = CreateWrapperWithEmptyProcessInfo(model_part);
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    model_part.AddNodalSolutionStepVariable(TOTAL_DISPLACEMENT);

    auto p_node = model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    const auto original_total_displacement = array_1d<double, 3>{1.0, 2.0, 3.0};
    p_node->GetSolutionStepValue(TOTAL_DISPLACEMENT) = original_total_displacement;

    // Saving the total displacement field twice should still result in the same expected result.
    strategy_wrapper.SaveTotalDisplacementFieldAtStartOfTimeLoop();
    strategy_wrapper.SaveTotalDisplacementFieldAtStartOfTimeLoop();

    const auto displacement_in_time_step = array_1d<double, 3>{3.0, 2.0, 1.0};
    p_node->GetSolutionStepValue(DISPLACEMENT) = displacement_in_time_step;
    strategy_wrapper.AccumulateTotalDisplacementField();

    const auto expected_total_displacement = array_1d<double, 3>{original_total_displacement + displacement_in_time_step};

    KRATOS_EXPECT_EQ(p_node->GetSolutionStepValue(TOTAL_DISPLACEMENT), expected_total_displacement);
}

KRATOS_TEST_CASE_IN_SUITE(RestorePositionsAndDOFVectorToStartOfStep_UpdatesPosition, KratosGeoMechanicsFastSuite)
{
    Model model;
    auto& model_part       = CreateDummyModelPart(model);
    auto  strategy_wrapper = CreateWrapperWithEmptyProcessInfo(model_part);
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT);

    const auto initial_position = array_1d<double, 3>{1.0, 2.0, 3.0};
    auto p_node = model_part.CreateNewNode(1, initial_position[0], initial_position[1], initial_position[2]);

    const auto displacement_in_time_step = array_1d<double, 3>{3.0, 2.0, 1.0};
    p_node->GetSolutionStepValue(DISPLACEMENT, 1) = displacement_in_time_step;
    strategy_wrapper.RestorePositionsAndDOFVectorToStartOfStep();

    const auto expected_position_after_displacement = array_1d<double, 3>{initial_position + displacement_in_time_step};

    KRATOS_EXPECT_VECTOR_EQ(p_node->Coordinates(), expected_position_after_displacement)
}

KRATOS_TEST_CASE_IN_SUITE(RestoreNodalDisplacementsAndWaterPressuresOnRequest, KratosGeoMechanicsFastSuite)
{
    Model model;
    auto& model_part       = CreateDummyModelPart(model);
    auto  strategy_wrapper = CreateWrapperWithEmptyProcessInfo(model_part);
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    model_part.AddNodalSolutionStepVariable(WATER_PRESSURE);

    // Define the old Degrees of Freedom
    auto p_node = model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    const auto old_displacement = array_1d<double, 3>{1.0, 2.0, 3.0};
    p_node->GetSolutionStepValue(DISPLACEMENT, 1) = old_displacement;
    const auto old_water_pressure = 4.0;
    p_node->GetSolutionStepValue(WATER_PRESSURE, 1) = old_water_pressure;

    // Set some new Degrees of Freedom (emulating a time step calculation)
    const auto new_displacement = array_1d<double, 3>{10.0, 20.0, 30.0};
    p_node->GetSolutionStepValue(DISPLACEMENT, 0) = new_displacement;
    const auto new_water_pressure = 40.0;
    p_node->GetSolutionStepValue(WATER_PRESSURE, 0) = new_water_pressure;

    strategy_wrapper.RestorePositionsAndDOFVectorToStartOfStep();

    KRATOS_EXPECT_VECTOR_EQ(p_node->GetSolutionStepValue(DISPLACEMENT, 0), old_displacement)
    KRATOS_EXPECT_EQ(p_node->GetSolutionStepValue(WATER_PRESSURE, 0), old_water_pressure);
}

KRATOS_TEST_CASE_IN_SUITE(RestoreNodalDisplacementsAndRotationsOnRequest, KratosGeoMechanicsFastSuite)
{
    Model model;
    auto& model_part       = CreateDummyModelPart(model);
    auto  strategy_wrapper = CreateWrapperWithEmptyProcessInfo(model_part);
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    model_part.AddNodalSolutionStepVariable(ROTATION);

    // Define the old Degrees of Freedom
    auto p_node = model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    const auto old_displacement = array_1d<double, 3>{1.0, 2.0, 3.0};
    p_node->GetSolutionStepValue(DISPLACEMENT, 1) = old_displacement;
    const auto old_rotation = array_1d<double, 3>{5.0, 6.0, 7.0};
    p_node->GetSolutionStepValue(ROTATION, 1) = old_rotation;

    // Set some new Degrees of Freedom (emulating a time step calculation)
    const auto new_displacement = array_1d<double, 3>{10.0, 20.0, 30.0};
    p_node->GetSolutionStepValue(DISPLACEMENT, 0) = new_displacement;
    const auto new_rotation = array_1d<double, 3>{50.0, 60.0, 70.0};
    p_node->GetSolutionStepValue(ROTATION, 0) = new_rotation;

    strategy_wrapper.RestorePositionsAndDOFVectorToStartOfStep();

    KRATOS_EXPECT_VECTOR_EQ(p_node->GetSolutionStepValue(DISPLACEMENT, 0), old_displacement)
    KRATOS_EXPECT_VECTOR_EQ(p_node->GetSolutionStepValue(ROTATION, 0), old_rotation)
}

}
