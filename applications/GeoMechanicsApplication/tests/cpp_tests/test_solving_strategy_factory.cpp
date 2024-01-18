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
#include "custom_utilities/solving_strategy_factory.hpp"
#include "containers/model.h"

using namespace Kratos;

namespace Kratos::Testing
{

using SparseSpaceType = UblasSpace<double, CompressedMatrix, Vector>;
using DenseSpaceType = UblasSpace<double, Matrix, Vector>;
using LinearSolverType = LinearSolver<SparseSpaceType, DenseSpaceType>;
using SolvingStrategyFactoryType = SolvingStrategyFactory<SparseSpaceType, DenseSpaceType, LinearSolverType>;

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
        "move_mesh_flag":                     true,
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


KRATOS_TEST_CASE_IN_SUITE(CreateSolvingStrategy_Throws_WhenNoStrategyTypeIsDefined, KratosGeoMechanicsFastSuite)
{
    Model model;
    auto& dummy_model_part = model.CreateModelPart("dummy");

    const Parameters parameters{};

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(auto test = SolvingStrategyFactoryType::Create(
            parameters, dummy_model_part), "The parameter strategy_type is undefined, aborting.")
}

KRATOS_TEST_CASE_IN_SUITE(Create_ReturnsSolvingStrategy_ForNewtonRaphsonStrategy, KratosGeoMechanicsFastSuite)
{
    Model model;
    const int buffer_size = 2;
    auto& dummy_model_part = model.CreateModelPart("dummy", buffer_size);
    const Parameters parameters{testParameters};

    const auto created_strategy = SolvingStrategyFactoryType::Create(
            parameters, dummy_model_part);

    KRATOS_EXPECT_TRUE(created_strategy->GetMoveMeshFlag()) // since it is set to true in the parameters
    KRATOS_EXPECT_NE(created_strategy, nullptr);
    KRATOS_EXPECT_EQ(created_strategy->Check(), 0);
}

KRATOS_TEST_CASE_IN_SUITE(Create_ReturnsSolvingStrategy_ForLineSearchStrategy, KratosGeoMechanicsFastSuite)
{
    Model model;
    const int buffer_size = 2;
    auto& dummy_model_part = model.CreateModelPart("dummy", buffer_size);
    Parameters parameters{testParameters};
    parameters["strategy_type"].SetString("line_search");

    const auto created_strategy = SolvingStrategyFactoryType::Create(
            parameters, dummy_model_part);

    KRATOS_EXPECT_TRUE(created_strategy->GetMoveMeshFlag()) // since it is set to true in the parameters
    KRATOS_EXPECT_NE(created_strategy, nullptr);
    KRATOS_EXPECT_EQ(created_strategy->Check(), 0);
}

}
