// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Gennady Markelov
//

// Project includes
#include "custom_elements/Pw_element.hpp"
#include "custom_elements/geo_steady_state_Pw_piping_element.h"
#include "custom_strategies/strategies/geo_mechanics_newton_raphson_erosion_process_strategy.hpp"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "tests/cpp_tests/test_utilities.h"

#include <geo_mechanics_application.h>
#include <linear_solvers/skyline_lu_factorization_solver.h>
#include <solving_strategies/convergencecriterias/mixed_generic_criteria.h>

namespace Kratos
{
template <class TSparseSpace, class TDenseSpace, class TLinearSolver>
class MockGeoMechanicsNewtonRaphsonErosionProcessStrategy
    : public GeoMechanicsNewtonRaphsonErosionProcessStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:
    using SolvingStrategyType = ImplicitSolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>;
    using TConvergenceCriteriaType = ConvergenceCriteria<TSparseSpace, TDenseSpace>;
    using TBuilderAndSolverType    = typename SolvingStrategyType::TBuilderAndSolverType;
    using TSchemeType              = typename SolvingStrategyType::TSchemeType;

    MockGeoMechanicsNewtonRaphsonErosionProcessStrategy(ModelPart&                    rModelPart,
                                                        typename TSchemeType::Pointer pScheme,
                                                        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
                                                        typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver,
                                                        Parameters& rParameters)
        : GeoMechanicsNewtonRaphsonErosionProcessStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(
              rModelPart, std::move(pScheme), std::move(pNewConvergenceCriteria), std::move(pNewBuilderAndSolver), rParameters),
          mrModelPart(rModelPart)
    {
    }

private:
    bool Recalculate() override
    {
        // The function provides a pressure change over iterations that leads to a non-linear growth of an equilibrium pipe height.
        // As a result the height curve crosses twice a pipe height growth curve in the search algorithm implemented
        // in  GeoMechanicsNewtonRaphsonErosionProcessStrategy. The search stops at the second point.
        for (auto& node : mrModelPart.Nodes()) {
            node.FastGetSolutionStepValue(WATER_PRESSURE) =
                std::pow(node.FastGetSolutionStepValue(WATER_PRESSURE), 0.9725);
        }
        return true;
    }

    void BaseClassFinalizeSolutionStep() override
    {
        // to avoid calling the BaseClassFinalizeSolutionStep in tested class.
    }

    ModelPart& mrModelPart;
};

Node::Pointer CreateNodeWithSolutionStepValues(ModelPart& rModelPart, int Id, double CoordinateX, double CoordinateY, double WaterPressure)
{
    constexpr double coordinate_z = 0.0;
    auto             p_node = rModelPart.CreateNewNode(Id, CoordinateX, CoordinateY, coordinate_z);

    p_node->FastGetSolutionStepValue(WATER_PRESSURE) = WaterPressure;
    const array_1d<double, 3> GravitationalAcceleration{0, -9.81, 0};
    p_node->FastGetSolutionStepValue(VOLUME_ACCELERATION) = GravitationalAcceleration;

    return p_node;
}

Geometry<Node>::Pointer CreateQuadrilateral2D4N(ModelPart&              rModelPart,
                                                const std::vector<int>& rNodeIds,
                                                double                  Xmin,
                                                double                  Xmax,
                                                double                  WaterPressureLeft,
                                                double                  WaterPressureRight)
{
    constexpr double y_min = -0.1;
    constexpr double y_max = 0.1;

    return make_shared<Quadrilateral2D4<Node>>(
        CreateNodeWithSolutionStepValues(rModelPart, rNodeIds[0], Xmin, y_min, WaterPressureLeft),
        CreateNodeWithSolutionStepValues(rModelPart, rNodeIds[1], Xmax, y_min, WaterPressureRight),
        CreateNodeWithSolutionStepValues(rModelPart, rNodeIds[2], Xmax, y_max, WaterPressureRight),
        CreateNodeWithSolutionStepValues(rModelPart, rNodeIds[3], Xmin, y_max, WaterPressureLeft));
}

Geometry<Node>::Pointer CreateLine2D2N(ModelPart&              rModelPart,
                                       const std::vector<int>& rNodeIds,
                                       double                  Xmin,
                                       double                  Xmax,
                                       double                  WaterPressureLeft,
                                       double                  WaterPressureRight)
{
    constexpr double y = 0.0;

    return make_shared<Line2D2<Node>>(
        CreateNodeWithSolutionStepValues(rModelPart, rNodeIds[0], Xmin, y, WaterPressureLeft),
        CreateNodeWithSolutionStepValues(rModelPart, rNodeIds[1], Xmax, y, WaterPressureRight));
}

auto SetupPipingStrategy(Model& rModel)
{
    using SparseSpaceType             = UblasSpace<double, CompressedMatrix, Vector>;
    using LocalSpaceType              = UblasSpace<double, Matrix, Vector>;
    using LinearSolverType            = LinearSolver<SparseSpaceType, LocalSpaceType>;
    using ConvergenceCriteriaType     = ConvergenceCriteria<SparseSpaceType, LocalSpaceType>;
    using MixedGenericCriteriaType    = MixedGenericCriteria<SparseSpaceType, LocalSpaceType>;
    using ConvergenceVariableListType = MixedGenericCriteriaType::ConvergenceVariableListType;

    auto& r_model_part = rModel.CreateModelPart("ModelPart", 1);
    r_model_part.AddNodalSolutionStepVariable(WATER_PRESSURE);
    r_model_part.AddNodalSolutionStepVariable(VOLUME_ACCELERATION);
    const auto& r_process_info = r_model_part.GetProcessInfo();

    auto p_element_props = r_model_part.CreateNewProperties(0);
    p_element_props->SetValue(CONSTITUTIVE_LAW, LinearElastic2DInterfaceLaw().Clone());

    auto contributions = {CalculationContribution::Permeability, CalculationContribution::FluidBodyFlow};
    auto p_element = make_intrusive<PwElement<2, 4>>(
        0, CreateQuadrilateral2D4N(r_model_part, std::vector<int>{13, 14, 15, 16}, 3.0, 4.0, 2000.0, 2000.0),
        p_element_props, contributions, nullptr);
    p_element->Initialize(r_process_info);
    r_model_part.AddElement(p_element);

    // Set the pipe element properties
    auto p_piping_element_props = r_model_part.CreateNewProperties(1);
    p_piping_element_props->SetValue(PIPE_D_70, 2e-4);
    p_piping_element_props->SetValue(PIPE_ETA, 0.25);
    p_piping_element_props->SetValue(PIPE_THETA, 30);
    p_piping_element_props->SetValue(DENSITY_SOLID, 2650);
    p_piping_element_props->SetValue(DENSITY_WATER, 1000);
    p_piping_element_props->SetValue(PIPE_ELEMENT_LENGTH, 1.0);
    p_piping_element_props->SetValue(PIPE_MODIFIED_D, false);
    p_piping_element_props->SetValue(PIPE_MODEL_FACTOR, 1);
    p_piping_element_props->SetValue(PIPE_START_ELEMENT, 1);

    // Create the start piping element
    auto p_piping_element = make_intrusive<GeoSteadyStatePwPipingElement<2, 2>>(
        1, CreateLine2D2N(r_model_part, std::vector{1, 2}, 0.0, 1.0, 0.0, 500.0), p_piping_element_props);
    p_piping_element->Initialize(r_process_info);
    r_model_part.AddElement(p_piping_element);

    // Create other piping elements
    p_piping_element = make_intrusive<GeoSteadyStatePwPipingElement<2, 2>>(
        2, CreateLine2D2N(r_model_part, std::vector{3, 4}, 2.0, 3.0, 500.0, 1000.0), p_piping_element_props);
    p_piping_element->Initialize(r_process_info);
    r_model_part.AddElement(p_piping_element);

    p_piping_element = make_intrusive<GeoSteadyStatePwPipingElement<2, 2>>(
        3, CreateLine2D2N(r_model_part, std::vector{5, 6}, 1.0, 2.0, 1000.0, 1500.0), p_piping_element_props);
    p_piping_element->Initialize(r_process_info);
    r_model_part.AddElement(p_piping_element);

    Parameters p_parameters(R"(
    {
 		"max_piping_iterations": 25
    }  )");

    const auto p_builder_and_solver =
        Kratos::make_shared<ResidualBasedBlockBuilderAndSolver<SparseSpaceType, LocalSpaceType, LinearSolverType>>(nullptr);
    const ConvergenceVariableListType      convergence_settings{};
    const ConvergenceCriteriaType::Pointer p_criteria =
        std::make_shared<MixedGenericCriteriaType>(convergence_settings);

    return std::make_unique<MockGeoMechanicsNewtonRaphsonErosionProcessStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType>>(
        r_model_part, nullptr, p_criteria, p_builder_and_solver, p_parameters);
}
} // namespace Kratos

namespace Kratos::Testing
{
KRATOS_TEST_CASE_IN_SUITE(CheckGetPipingElements, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    Model model;
    auto  p_solving_strategy = SetupPipingStrategy(model);

    // Act
    const auto piping_elements = p_solving_strategy->GetPipingElements();

    // Assert
    const auto    elements           = model.GetModelPart("ModelPart").Elements();
    constexpr int number_of_elements = 4;
    KRATOS_EXPECT_EQ(elements.size(), number_of_elements);

    constexpr int number_of_piping_elements = 3;
    KRATOS_EXPECT_EQ(piping_elements.size(), number_of_piping_elements);
    KRATOS_EXPECT_EQ(piping_elements[0]->GetId(), 1);
    KRATOS_EXPECT_EQ(piping_elements[1]->GetId(), 3);
    KRATOS_EXPECT_EQ(piping_elements[2]->GetId(), 2);
}

KRATOS_TEST_CASE_IN_SUITE(CheckFinalizeSolutionStep, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    Model model;
    auto  p_solving_strategy = SetupPipingStrategy(model);

    // Act
    p_solving_strategy->FinalizeSolutionStep();

    // Assert
    auto elements = model.GetModelPart("ModelPart").Elements();

    constexpr int active = 1;
    KRATOS_EXPECT_EQ(elements[1].GetValue(PIPE_ACTIVE), active);

    constexpr double expected_active_pipe_height = 0.0183333;
    KRATOS_EXPECT_NEAR(elements[1].GetValue(PIPE_HEIGHT), expected_active_pipe_height, Defaults::relative_tolerance);

    constexpr int inactive = 0;
    KRATOS_EXPECT_EQ(elements[2].GetValue(PIPE_ACTIVE), inactive);
    KRATOS_EXPECT_EQ(elements[3].GetValue(PIPE_ACTIVE), inactive);

    constexpr double expected_inactive_pipe_height = 0.0;
    KRATOS_EXPECT_NEAR(elements[2].GetValue(PIPE_HEIGHT), expected_inactive_pipe_height, Defaults::relative_tolerance);
    KRATOS_EXPECT_NEAR(elements[3].GetValue(PIPE_HEIGHT), expected_inactive_pipe_height, Defaults::relative_tolerance);
}
} // namespace Kratos::Testing
