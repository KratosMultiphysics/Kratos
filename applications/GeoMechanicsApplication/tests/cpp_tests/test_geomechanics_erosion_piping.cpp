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
#include "custom_elements/plane_strain_stress_state.h"
#include "custom_elements/steady_state_Pw_piping_element.hpp"
#include "geo_mechanics_fast_suite.h"
#include "custom_elements/U_Pw_small_strain_element.hpp"
#include "custom_strategies/strategies/geo_mechanics_newton_raphson_erosion_process_strategy.hpp"
#include "custom_workflows/dgeoflow.h"

namespace Kratos
{
template <class TSparseSpace, class TDenseSpace, class TLinearSolver>
class MockGeoMechanicsNewtonRaphsonErosionProcessStrategy
    : public GeoMechanicsNewtonRaphsonErosionProcessStrategy<
        TSparseSpace, TDenseSpace, TLinearSolver>
{
public:
    using BaseType = ImplicitSolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>;
    using TConvergenceCriteriaType = ConvergenceCriteria<TSparseSpace, TDenseSpace>;
    using TBuilderAndSolverType = typename BaseType::TBuilderAndSolverType;
    using TSchemeType = typename BaseType::TSchemeType;
    using DofsArrayType = typename BaseType::DofsArrayType;
    using TSystemMatrixType = typename BaseType::TSystemMatrixType;
    using TSystemVectorType = typename BaseType::TSystemVectorType;
    using filtered_elements =
    typename boost::range_detail::filtered_range<
        std::function<bool(Element*)>, std::vector<Element*>>;
    using PropertiesType = Properties;
    using NodeType       = Node;
    using GeometryType   = Geometry<NodeType>;

    MockGeoMechanicsNewtonRaphsonErosionProcessStrategy(ModelPart&                    model_part,
                                                        typename TSchemeType::Pointer pScheme,
                                                        typename TLinearSolver::Pointer
                                                        pNewLinearSolver,
                                                        typename TConvergenceCriteriaType::Pointer
                                                        pNewConvergenceCriteria,
                                                        typename TBuilderAndSolverType::Pointer
                                                        pNewBuilderAndSolver,
                                                        Parameters& rParameters,
                                                        int         MaxIterations          = 30,
                                                        bool        CalculateReactions     = false,
                                                        bool        ReformDofSetAtEachStep = false,
                                                        bool        MoveMeshFlag           = false)
        : GeoMechanicsNewtonRaphsonErosionProcessStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(
            model_part,
            pScheme,
            pNewLinearSolver,
            pNewConvergenceCriteria,
            pNewBuilderAndSolver,
            rParameters,
            MaxIterations,
            CalculateReactions,
            ReformDofSetAtEachStep,
            MoveMeshFlag)
    {
    };

private:
    bool Recalculate() override
    {
        std::cout << " Recalculate test" << std::endl;
        return true;
    }

    void BaseClassFinalizeSolutionStep() override
    {
        std::cout << " BaseClassFinalizeSolutionStep test" << std::endl;
    };
};

KratosExecute::GeoMechanicsNewtonRaphsonErosionProcessStrategyType::Pointer SetupPipingStrategy(
    Model& pModel)
{
    using SparseSpaceType = UblasSpace<double, CompressedMatrix, Vector>;
    using LocalSpaceType  = UblasSpace<double, Matrix, Vector>;

    // The direct solver
    using LinearSolverType = LinearSolver<SparseSpaceType, LocalSpaceType>;

    std::cout << " begins" << std::endl;

    // initialize modelpart
    auto& r_model_part = pModel.CreateModelPart("ModelPart", 1);
    r_model_part.GetProcessInfo().SetValue(DOMAIN_SIZE, 5);
    r_model_part.AddNodalSolutionStepVariable(WATER_PRESSURE);
    r_model_part.AddNodalSolutionStepVariable(VOLUME_ACCELERATION);
    const auto& r_process_info = r_model_part.GetProcessInfo();

    auto p_non_piping_elem_prop = r_model_part.CreateNewProperties(0);
    p_non_piping_elem_prop->SetValue(CONSTITUTIVE_LAW, LinearElastic2DInterfaceLaw().Clone());

    // Create the test non-piping element
    auto p_node_13 = r_model_part.CreateNewNode(13, 3.0, -0.1, 0.0);
    auto p_node_14 = r_model_part.CreateNewNode(14, 4.0, -0.1, 0.0);
    auto p_node_15 = r_model_part.CreateNewNode(15, 4.0, 0.1, 0.0);
    auto p_node_16 = r_model_part.CreateNewNode(16, 3.0, 0.1, 0.0);

    auto p_non_piping_element = make_intrusive<UPwSmallStrainElement<2, 4>>(
        0, Kratos::make_shared<Quadrilateral2D4<Node>>(p_node_13, p_node_14, p_node_15, p_node_16),
        p_non_piping_elem_prop, std::make_unique<PlaneStrainStressState>());
    p_non_piping_element->Initialize(r_process_info);
    r_model_part.AddElement(p_non_piping_element);

    // Set the pipe element properties
    auto p_elem_prop = r_model_part.CreateNewProperties(1);
    p_elem_prop->SetValue(PIPE_D_70, 2e-4);
    p_elem_prop->SetValue(PIPE_ETA, 0.25);
    p_elem_prop->SetValue(PIPE_THETA, 30);
    p_elem_prop->SetValue(DENSITY_SOLID, 2650);
    p_elem_prop->SetValue(DENSITY_WATER, 1000);

    p_elem_prop->SetValue(PIPE_ELEMENT_LENGTH, 0.5);
    p_elem_prop->SetValue(PIPE_MODIFIED_D, false);
    p_elem_prop->SetValue(PIPE_MODEL_FACTOR, 1);
    p_elem_prop->SetValue(PIPE_START_ELEMENT, 1);
    p_elem_prop->SetValue(CONSTITUTIVE_LAW, LinearElastic2DInterfaceLaw().Clone());

    auto p_node_1 = r_model_part.CreateNewNode(1, 0.0, -0.1, 0.0);
    auto p_node_2 = r_model_part.CreateNewNode(2, 1.0, -0.1, 0.0);
    auto p_node_3 = r_model_part.CreateNewNode(3, 1.0, 0.1, 0.0);
    auto p_node_4 = r_model_part.CreateNewNode(4, 0.0, 0.1, 0.0);

    // Create the start piping element
    auto p_element = make_intrusive<SteadyStatePwPipingElement<2, 4>>(
        1, Kratos::make_shared<Quadrilateral2D4<Node>>(p_node_1, p_node_2, p_node_3, p_node_4),
        p_elem_prop, std::make_unique<PlaneStrainStressState>());
    p_element->Initialize(r_process_info);
    r_model_part.AddElement(p_element);

    auto p_node_5 = r_model_part.CreateNewNode(5, 2.0, -0.1, 0.0);
    auto p_node_6 = r_model_part.CreateNewNode(6, 3.0, -0.1, 0.0);
    auto p_node_7 = r_model_part.CreateNewNode(7, 3.0, 0.1, 0.0);
    auto p_node_8 = r_model_part.CreateNewNode(8, 2.0, 0.1, 0.0);

    // Create the test piping element
    p_element = make_intrusive<SteadyStatePwPipingElement<2, 4>>(
        2, Kratos::make_shared<Quadrilateral2D4<Node>>(p_node_5, p_node_6, p_node_7, p_node_8),
        p_elem_prop, std::make_unique<PlaneStrainStressState>());
    p_element->Initialize(r_process_info);
    r_model_part.AddElement(p_element);

    auto p_node_9  = r_model_part.CreateNewNode(9, 1.0, -0.1, 0.0);
    auto p_node_10 = r_model_part.CreateNewNode(10, 2.0, -0.1, 0.0);
    auto p_node_11 = r_model_part.CreateNewNode(11, 2.0, 0.1, 0.0);
    auto p_node_12 = r_model_part.CreateNewNode(12, 1.0, 0.1, 0.0);

    // Create the test piping element
    p_element = make_intrusive<SteadyStatePwPipingElement<2, 4>>(
        3, Kratos::make_shared<Quadrilateral2D4<Node>>(p_node_9, p_node_10, p_node_11, p_node_12),
        p_elem_prop, std::make_unique<PlaneStrainStressState>());
    p_element->Initialize(r_process_info);
    r_model_part.AddElement(p_element);

    Parameters p_parameters(R"(
    {
        "min_iteration":    6,
        "number_cycles":    100,
        "increase_factor":  2.0,
        "reduction_factor": 0.5,
		"max_piping_iterations": 50,
        "desired_iterations": 4,
        "max_radius_factor": 10.0,
        "min_radius_factor": 0.1,
        "search_neighbours_step": false,
        "body_domain_sub_model_part_list": [],
        "loads_sub_model_part_list": [],
        "loads_variable_list" : []
    }  )");

    int  max_iterations            = 15;
    bool calculate_reactions       = true;
    bool reform_dof_set_atEachStep = false;
    bool move_mesh_flag            = false;

    auto p_solver = Kratos::make_shared<SkylineLUFactorizationSolver<
        SparseSpaceType, LocalSpaceType>>();
    Scheme<SparseSpaceType, LocalSpaceType>::Pointer p_scheme =
        Kratos::make_shared<BackwardEulerQuasistaticPwScheme<SparseSpaceType, LocalSpaceType>>();
    auto p_builder_and_solver =
        Kratos::make_shared<ResidualBasedBlockBuilderAndSolver<
            SparseSpaceType, LocalSpaceType, KratosExecute::LinearSolverType>>(
            p_solver);
    p_builder_and_solver->SetEchoLevel(0);
    const double                               rel_tol      = 1.0e-4;
    const double                               abs_tol      = 1.0e-9;
    VariableData*                              p_water_pres = &WATER_PRESSURE;
    KratosExecute::ConvergenceVariableListType convergence_settings{
        std::make_tuple(p_water_pres, rel_tol, abs_tol)};
    auto p_criteria = std::make_shared<KratosExecute::MixedGenericCriteriaType>(
        convergence_settings);
    p_criteria->SetEchoLevel(0);

    auto p_solving_strategy =
        Kratos::make_unique<MockGeoMechanicsNewtonRaphsonErosionProcessStrategy<
            SparseSpaceType, LocalSpaceType, KratosExecute::LinearSolverType>>(
            pModel.GetModelPart("ModelPart"), p_scheme, p_solver, p_criteria, p_builder_and_solver,
            p_parameters,
            max_iterations, calculate_reactions, reform_dof_set_atEachStep, move_mesh_flag);

    return p_solving_strategy;
}
}

namespace Kratos::Testing
{
KRATOS_TEST_CASE_IN_SUITE(GetPipingElements, KratosGeoMechanicsFastSuite)
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

KRATOS_TEST_CASE_IN_SUITE(FinalizeSolutionStep, KratosGeoMechanicsFastSuite)
{
    Model model;
    auto  p_solving_strategy = SetupPipingStrategy(model);

    const auto piping_elements = p_solving_strategy->GetPipingElements();

    p_solving_strategy->FinalizeSolutionStep();

    constexpr double expected_pipe_height = 3;

    KRATOS_EXPECT_NEAR(piping_elements[0]->GetProperties().GetValue(PIPE_HEIGHT),
                       expected_pipe_height, 1.0e-6);
    KRATOS_EXPECT_NEAR(piping_elements[1]->GetProperties().GetValue(PIPE_HEIGHT),
                       expected_pipe_height, 1.0e-6);
    KRATOS_EXPECT_NEAR(piping_elements[2]->GetProperties().GetValue(PIPE_HEIGHT),
                       expected_pipe_height, 1.0e-6);
}
} // namespace Kratos::Testing