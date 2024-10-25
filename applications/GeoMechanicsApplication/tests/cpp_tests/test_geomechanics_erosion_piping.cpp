// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Aron Noordam
//

#include <string>

// Project includes
#include "containers/model.h"
#include "custom_constitutive/linear_elastic_2D_interface_law.h"
#include "custom_elements/U_Pw_small_strain_element.hpp"
#include "custom_elements/plane_strain_stress_state.h"
#include "custom_elements/steady_state_Pw_piping_element.hpp"
#include "custom_strategies/strategies/geo_mechanics_newton_raphson_erosion_process_strategy.hpp"
#include "custom_workflows/dgeoflow.h"
#include "geo_mechanics_fast_suite.h"

namespace Kratos
{
template <class TSparseSpace, class TDenseSpace, class TLinearSolver>
class MockGeoMechanicsNewtonRaphsonErosionProcessStrategy
    : public GeoMechanicsNewtonRaphsonErosionProcessStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:
    using BaseType = ImplicitSolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>;
    using TConvergenceCriteriaType = ConvergenceCriteria<TSparseSpace, TDenseSpace>;
    using TBuilderAndSolverType    = typename BaseType::TBuilderAndSolverType;
    using TSchemeType              = typename BaseType::TSchemeType;

    MockGeoMechanicsNewtonRaphsonErosionProcessStrategy(ModelPart&                    model_part,
                                                        typename TSchemeType::Pointer pScheme,
                                                        typename TLinearSolver::Pointer pNewLinearSolver,
                                                        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
                                                        typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver,
                                                        Parameters& rParameters)
        : GeoMechanicsNewtonRaphsonErosionProcessStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(
              model_part, pScheme, pNewLinearSolver, pNewConvergenceCriteria, pNewBuilderAndSolver, rParameters),
          mModelPart(model_part){};

private:
    bool Recalculate() override
    {
        // The function provides a pressure change over iterations that leads to a non-linear growth of an equilibrium pipe height.
        // As a result the height curve crosses twice a pipe height growth curve in the search algorithm implemented
        // in  GeoMechanicsNewtonRaphsonErosionProcessStrategy. The search stops at the second point.
        for (auto& node : mModelPart.Nodes()) {
            node.FastGetSolutionStepValue(WATER_PRESSURE) =
                pow(node.FastGetSolutionStepValue(WATER_PRESSURE), 0.9725);
        }
        return true;
    }

    void       BaseClassFinalizeSolutionStep() override{};
    ModelPart& mModelPart;
};

const array_1d<double, 3> GravityArray{0, 9.81, 0};

Node::Pointer CreateNodeWithFastSolutionStepValues(ModelPart& rModelPart, int Id, double X, double Y, double WaterPressure)
{
    constexpr double z      = 0.0;
    auto             p_node = rModelPart.CreateNewNode(Id, X, Y, z);

    p_node->SetLock();

    p_node->FastGetSolutionStepValue(WATER_PRESSURE)      = WaterPressure;
    p_node->FastGetSolutionStepValue(VOLUME_ACCELERATION) = GravityArray;

    p_node->UnSetLock();

    return p_node;
}

Geometry<Node>::Pointer CreateQuadrilateral2D4N(
    ModelPart& rModelPart, int Id1, int Id2, int Id3, int Id4, double Xmin, double Xmax, double WaterPressureLeft, double WaterPressureRight)
{
    constexpr double y_min = -0.1;
    constexpr double y_max = 0.1;

    return Kratos::make_shared<Quadrilateral2D4<Node>>(
        CreateNodeWithFastSolutionStepValues(rModelPart, Id1, Xmin, y_min, WaterPressureLeft),
        CreateNodeWithFastSolutionStepValues(rModelPart, Id2, Xmax, y_min, WaterPressureRight),
        CreateNodeWithFastSolutionStepValues(rModelPart, Id3, Xmax, y_max, WaterPressureRight),
        CreateNodeWithFastSolutionStepValues(rModelPart, Id4, Xmin, y_max, WaterPressureLeft));
}

KratosExecute::GeoMechanicsNewtonRaphsonErosionProcessStrategyType::Pointer SetupPipingStrategy(Model& rModel)
{
    using SparseSpaceType = UblasSpace<double, CompressedMatrix, Vector>;
    using LocalSpaceType  = UblasSpace<double, Matrix, Vector>;

    auto& r_model_part = rModel.CreateModelPart("ModelPart", 1);
    r_model_part.AddNodalSolutionStepVariable(WATER_PRESSURE);
    r_model_part.AddNodalSolutionStepVariable(VOLUME_ACCELERATION);
    const auto& r_process_info = r_model_part.GetProcessInfo();

    auto p_element_props = r_model_part.CreateNewProperties(0);
    p_element_props->SetValue(CONSTITUTIVE_LAW, LinearElastic2DInterfaceLaw().Clone());

    auto p_element = make_intrusive<UPwSmallStrainElement<2, 4>>(
        0, CreateQuadrilateral2D4N(r_model_part, 13, 14, 15, 16, 3.0, 4.0, 2000.0, 2000.0),
        p_element_props, std::make_unique<PlaneStrainStressState>());
    p_element->Initialize(r_process_info);
    r_model_part.AddElement(p_element);

    // Set the pipe element properties
    auto p_piping_element_prop = r_model_part.CreateNewProperties(1);
    p_piping_element_prop->SetValue(PIPE_D_70, 2e-4);
    p_piping_element_prop->SetValue(PIPE_ETA, 0.25);
    p_piping_element_prop->SetValue(PIPE_THETA, 30);
    p_piping_element_prop->SetValue(DENSITY_SOLID, 2650);
    p_piping_element_prop->SetValue(DENSITY_WATER, 1000);
    p_piping_element_prop->SetValue(PIPE_ELEMENT_LENGTH, 1.0);
    p_piping_element_prop->SetValue(PIPE_MODIFIED_D, false);
    p_piping_element_prop->SetValue(PIPE_MODEL_FACTOR, 1);
    p_piping_element_prop->SetValue(PIPE_START_ELEMENT, 1);
    p_piping_element_prop->SetValue(CONSTITUTIVE_LAW, LinearElastic2DInterfaceLaw().Clone());

    // Create the start piping element
    auto p_piping_element = make_intrusive<SteadyStatePwPipingElement<2, 4>>(
        1, CreateQuadrilateral2D4N(r_model_part, 1, 2, 3, 4, 0.0, 1.0, 0.0, 500.0),
        p_piping_element_prop, std::make_unique<PlaneStrainStressState>());
    p_piping_element->Initialize(r_process_info);
    r_model_part.AddElement(p_piping_element);

    // Create other piping elements
    p_piping_element = make_intrusive<SteadyStatePwPipingElement<2, 4>>(
        2, CreateQuadrilateral2D4N(r_model_part, 5, 6, 7, 8, 2.0, 3.0, 500.0, 1000.0),
        p_piping_element_prop, std::make_unique<PlaneStrainStressState>());
    p_piping_element->Initialize(r_process_info);
    r_model_part.AddElement(p_piping_element);

    p_piping_element = make_intrusive<SteadyStatePwPipingElement<2, 4>>(
        3, CreateQuadrilateral2D4N(r_model_part, 9, 10, 11, 12, 1.0, 2.0, 1000.0, 1500.0),
        p_piping_element_prop, std::make_unique<PlaneStrainStressState>());
    p_piping_element->Initialize(r_process_info);
    r_model_part.AddElement(p_piping_element);

    Parameters p_parameters(R"(
    {
 		"max_piping_iterations": 25
    }  )");

    const auto p_solver =
        Kratos::make_shared<SkylineLUFactorizationSolver<SparseSpaceType, LocalSpaceType>>();
    const auto p_scheme =
        Kratos::make_shared<BackwardEulerQuasistaticPwScheme<SparseSpaceType, LocalSpaceType>>();
    const auto p_builder_and_solver =
        Kratos::make_shared<ResidualBasedBlockBuilderAndSolver<SparseSpaceType, LocalSpaceType, KratosExecute::LinearSolverType>>(
            p_solver);
    const KratosExecute::ConvergenceVariableListType convergence_settings{};
    const auto p_criteria = std::make_shared<KratosExecute::MixedGenericCriteriaType>(convergence_settings);

    auto p_solving_strategy =
        Kratos::make_unique<MockGeoMechanicsNewtonRaphsonErosionProcessStrategy<SparseSpaceType, LocalSpaceType, KratosExecute::LinearSolverType>>(
            r_model_part, p_scheme, p_solver, p_criteria, p_builder_and_solver, p_parameters);

    return p_solving_strategy;
}
} // namespace Kratos

namespace Kratos::Testing
{
KRATOS_TEST_CASE_IN_SUITE(CheckGetPipingElements, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Model model;
    auto  p_solving_strategy = SetupPipingStrategy(model);

    const auto piping_elements = p_solving_strategy->GetPipingElements();

    constexpr int number_of_piping_elements = 3;

    KRATOS_EXPECT_EQ(piping_elements.size(), number_of_piping_elements);
    KRATOS_EXPECT_EQ(piping_elements[0]->GetId(), 1);
    KRATOS_EXPECT_EQ(piping_elements[1]->GetId(), 3);
    KRATOS_EXPECT_EQ(piping_elements[2]->GetId(), 2);
}

KRATOS_TEST_CASE_IN_SUITE(CheckFinalizeSolutionStep, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Model model;
    auto  p_solving_strategy = SetupPipingStrategy(model);

    p_solving_strategy->FinalizeSolutionStep();

    auto elements = model.GetModelPart("ModelPart").Elements();

    constexpr int active = 1;
    KRATOS_EXPECT_EQ(elements[1].GetValue(PIPE_ACTIVE), active);

    constexpr double expected_active_pipe_height = 0.0183333;
    KRATOS_EXPECT_NEAR(elements[1].GetValue(PIPE_HEIGHT), expected_active_pipe_height, 1.0e-6);

    constexpr int inactive = 0;
    KRATOS_EXPECT_EQ(elements[2].GetValue(PIPE_ACTIVE), inactive);
    KRATOS_EXPECT_EQ(elements[3].GetValue(PIPE_ACTIVE), inactive);

    constexpr double expected_inactive_pipe_height = 0.0;
    KRATOS_EXPECT_NEAR(elements[2].GetValue(PIPE_HEIGHT), expected_inactive_pipe_height, 1.0e-6);
    KRATOS_EXPECT_NEAR(elements[3].GetValue(PIPE_HEIGHT), expected_inactive_pipe_height, 1.0e-6);
}
} // namespace Kratos::Testing