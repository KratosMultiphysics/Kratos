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

#include "custom_strategies/builder_and_solvers/residualbased_block_builder_and_solver_linear_elastic_dynamic.h"
#include "custom_strategies/schemes/incremental_newmark_linear_elastic_U_scheme.hpp"
#include "custom_strategies/schemes/newmark_dynamic_U_Pw_scheme.hpp"
#include "custom_strategies/strategies/residualbased_newton_raphson_strategy_linear_elastic_dynamic.hpp"
#include "tests/cpp_tests/test_utilities/geo_mock_condition.h"
#include "tests/cpp_tests/test_utilities/geo_mock_element.h"

#include "geo_mechanics_fast_suite.h"
#include "spaces/ublas_space.h"

#include "factories/linear_solver_factory.h"
#include "linear_solvers/linear_solver.h"

#include "solving_strategies/convergencecriterias/displacement_criteria.h"

namespace Kratos::Testing
{

using namespace Kratos;
using SparseSpaceType  = UblasSpace<double, CompressedMatrix, Vector>;
using LocalSpaceType   = UblasSpace<double, Matrix, Vector>;
using LinearSolverType = LinearSolver<SparseSpaceType, LocalSpaceType>;

/// <summary>
/// Test class for the GeoMechanicNewtonRaphsonStrategyLinearElasticDynamic class
/// The tests are based on "Finite Element Procedures, Bathe K.J., 1996, page 782"
/// </summary>
class NewtonRaphsonStrategyLinearElasticDynamicTester
{
public:
    NewtonRaphsonStrategyLinearElasticDynamicTester()
    {
        CreateValidModelPartWithTwoDegreesOfFreedom();
    }

    void CreateValidModelPartWithTwoDegreesOfFreedom()
    {
        // add solution step variables
        auto& result = mModel.CreateModelPart("model_part", 2);
        result.AddNodalSolutionStepVariable(DISPLACEMENT);
        result.AddNodalSolutionStepVariable(VELOCITY);
        result.AddNodalSolutionStepVariable(ACCELERATION);
		result.AddNodalSolutionStepVariable(TOTAL_DISPLACEMENT);

        // add nodes, elements and conditions
        auto       p_node   = result.CreateNewNode(0, 0.0, 0.0, 0.0);
        const auto node_ids = std::vector<ModelPart::IndexType>{0};

        auto p_geometry = result.CreateNewGeometry("Point2D", node_ids);
        auto p_element  = Kratos::make_intrusive<GeoMockElement>(0, p_geometry);
        result.AddElement(p_element, 0);

        auto p_condition = Kratos::make_intrusive<GeoMockCondition>(0, p_geometry);
        result.AddCondition(p_condition, 0);

        // add dofs
        p_node->AddDof(DISPLACEMENT_X);
   //     p_node->AddDof(TOTAL_DISPLACEMENT_X);
        p_node->pGetDof(DISPLACEMENT_X)->SetEquationId(0);

        p_node->AddDof(DISPLACEMENT_Y);
//		p_node->AddDof(TOTAL_DISPLACEMENT_Y);
        p_node->pGetDof(DISPLACEMENT_Y)->SetEquationId(1);

        // initialize nodal values
        p_node->FastGetSolutionStepValue(VELOCITY, 0) = Kratos::array_1d<double, 3>{0.0, 0.0, 0.0};
        p_node->FastGetSolutionStepValue(ACCELERATION, 0) = Kratos::array_1d<double, 3>{0.0, 0.0, 0.0};
        p_node->FastGetSolutionStepValue(DISPLACEMENT, 0) = Kratos::array_1d<double, 3>{0.0, 0.0, 0.0};

        p_node->FastGetSolutionStepValue(VELOCITY, 1) = Kratos::array_1d<double, 3>{0.0, 0.0, 0.0};
        p_node->FastGetSolutionStepValue(ACCELERATION, 1) = Kratos::array_1d<double, 3>{0.0, 0.0, 0.0};
        p_node->FastGetSolutionStepValue(DISPLACEMENT, 1) = Kratos::array_1d<double, 3>{0.0, 0.0, 0.0};
    }

    ModelPart& GetModelPart() { return mModel.GetModelPart("model_part"); }

    /// <summary>
    /// Creates a valid strategy with either a direct linear solve or an iterative linear solver
    /// </summary>
    /// <param name="rModelPart">model part to be solved</param>
    /// <param name="RelativeTollerance">relative tolerance for the displacement criteria</param>
    /// <param name="AbsoluteTollerance">absolute tolerance for the displacement criteria</param>
    /// <param name="CalculateInitialSecondDerivative">When true, the initial value for the second derivative is calculated which
    /// fullfills equilibrium in the very first timestep</param>
    /// <param name="UseDirectSolver"> When true, the sparse lu linear solver is used, when false, the cg linear solver is used</param>
    /// <param name="SingleIteration"> When true, the solver is set to only perform a single iteration</param>
    /// <returns></returns>
    static GeoMechanicNewtonRaphsonStrategyLinearElasticDynamic<SparseSpaceType, LocalSpaceType, LinearSolverType> CreateValidStrategy(
        ModelPart&   rModelPart,
        const double RelativeTolerance,
        const double AbsoluteTolerance,
        const bool   CalculateInitialSecondDerivative,
        const bool   UseDirectSolver,
        const bool   SingleIteration)
    {
        double beta  = 0.25;
        double gamma = 0.5;
        // create strategy
        auto pScheme =
            std::make_shared<IncrementalNewmarkLinearElasticUScheme<SparseSpaceType, LocalSpaceType>>(beta, gamma);

        auto factory = LinearSolverFactory<SparseSpaceType, LocalSpaceType>{};

        std::string solver_type_parameters;
        if (UseDirectSolver) {
            solver_type_parameters = R"({ "solver_type": "sparse_lu" })";
        } else {
            solver_type_parameters = R"({ "solver_type": "cg" })";
        }

        auto pLinearSolver = factory.Create(Parameters(solver_type_parameters));

        auto pBuilderAndSolver =
            std::make_shared<ResidualBasedBlockBuilderAndSolverLinearElasticDynamic<SparseSpaceType, LocalSpaceType, LinearSolverType>>(
                pLinearSolver, 0.25, 0.5, CalculateInitialSecondDerivative);
        auto pConvergenceCriteria = std::make_shared<DisplacementCriteria<SparseSpaceType, LocalSpaceType>>(
            RelativeTolerance, AbsoluteTolerance);

        int max_iterations = 30;
        if (SingleIteration) {
            max_iterations = 1;
        }

        return GeoMechanicNewtonRaphsonStrategyLinearElasticDynamic<SparseSpaceType, LocalSpaceType, LinearSolverType>(
            rModelPart, pScheme, pConvergenceCriteria, pBuilderAndSolver, max_iterations, false, false);
    }

private:
    Model mModel;
};

/// <summary>
/// Runs the tests for NewtonRaphsonLinearElasticDynamic
/// </summary>
/// <param name="DeltaTime"> time step size</param>
/// <param name="UseIterations"> When true, low convergence tolerance is used such that multiple non linear iterations are perfomed</param>
/// <param name="CalculateInitialAcceleration">When true, the initial value for the second derivative is calculated which
/// fullfills equilibrium in the very first timestep</param>
/// <param name="UseDirectSolver"> When true, the sparse lu linear solver is used, when false, the cg linear solver is used </param>
/// <param name="ConvergeOutside"> When true, non-linear convergence is checked outside the solution step </param>
/// <param name="rExpectedDisplacementX"> expected displacements in x direction</param>
/// <param name="rExpectedDisplacementY"> expected dusplacements in y direction</param>
void TestNewtonRaphsonLinearElasticDynamic(const double               DeltaTime,
                                           const bool                 UseIterations,
                                           const bool                 CalculateInitialAcceleration,
                                           const bool                 UseDirectSolver,
                                           const bool                 ConvergeOutside,
                                           const std::vector<double>& rExpectedDisplacementX,
                                           const std::vector<double>& rExpectedDisplacementY)
{
    double RelativeConvergenceTolerance = 100;
    double AbsoluteConvergenceTolerance = 100;

    if (UseIterations) {
        RelativeConvergenceTolerance = 1e-12;
        AbsoluteConvergenceTolerance = 1e-12;
    }

    NewtonRaphsonStrategyLinearElasticDynamicTester tester;

    // set up model and solver
    ModelPart& model_part = tester.GetModelPart();

    auto& r_current_process_info = model_part.GetProcessInfo();

    r_current_process_info[TIME]       = 0.0;
    r_current_process_info[DELTA_TIME] = DeltaTime;

    auto stiffness_matrix  = Matrix(2, 2);
    stiffness_matrix(0, 0) = 6.0;
    stiffness_matrix(0, 1) = -2.0;
    stiffness_matrix(1, 0) = -2.0;
    stiffness_matrix(1, 1) = 4.0;

    auto mass_matrix  = Matrix(2, 2);
    mass_matrix(0, 0) = 2.0;
    mass_matrix(0, 1) = 0.0;
    mass_matrix(1, 0) = 0.0;
    mass_matrix(1, 1) = 1.0;

    auto damping_matrix = Matrix(ZeroMatrix(2, 2));

    auto rhs = Vector(2);
    rhs(0)   = 0.0;
    rhs(1)   = 10.0;

    auto& geo_custom_element = dynamic_cast<GeoMockElement&>(model_part.GetElement(0));

    geo_custom_element.SetLeftHandSide(stiffness_matrix);
    geo_custom_element.SetMassMatrix(mass_matrix);
    geo_custom_element.SetDampingMatrix(damping_matrix);

    auto& geo_custom_condition = dynamic_cast<GeoMockCondition&>(model_part.GetCondition(0));

    geo_custom_condition.SetRightHandSide(rhs);

    // create strategy

    bool single_iteration = false;
    if (ConvergeOutside) {
        single_iteration = true;
    }

    auto r_solver = NewtonRaphsonStrategyLinearElasticDynamicTester::CreateValidStrategy(
        model_part, RelativeConvergenceTolerance, AbsoluteConvergenceTolerance,
        CalculateInitialAcceleration, UseDirectSolver, single_iteration);
    // initialize solver
    r_solver.Initialize();

    std::vector<double> calculated_displacement_x;
    std::vector<double> calculated_displacement_y;

    const std::size_t n_steps          = 12;
    const std::size_t max_iterations   = 30;
    std::size_t       iteration_number = 0;
    // run test solution loop
    for (std::size_t i = 0; i < n_steps; i++) {
        double new_time = r_current_process_info[TIME] + r_current_process_info[DELTA_TIME];

        r_current_process_info[STEP] += 1;
        model_part.CloneTimeStep(new_time);
        r_current_process_info[TIME] = new_time;

        r_solver.InitializeSolutionStep();
        r_solver.Predict();

        bool is_converged = false;
        while (!is_converged && iteration_number < max_iterations) {
            is_converged = r_solver.SolveSolutionStep();
            r_solver.SetStiffnessMatrixIsBuilt(true);
            iteration_number++;
        }

        r_solver.FinalizeSolutionStep();

        // store calculated results
        auto p_node = model_part.pGetNode(0);
        calculated_displacement_x.push_back(p_node->FastGetSolutionStepValue(DISPLACEMENT_X, 0));
        calculated_displacement_y.push_back(p_node->FastGetSolutionStepValue(DISPLACEMENT_Y, 0));
    }

    // check if the building functions are called the right amount, as this is largely why this solver is efficient for linear elastic systems
    KRATOS_EXPECT_EQ(geo_custom_element.GetCountCalculateLeftHandSideCalled(), 1);
    KRATOS_EXPECT_EQ(geo_custom_element.GetCountCalculateMassMatrixCalled(), 1);
    KRATOS_EXPECT_EQ(geo_custom_element.GetCalculateDampingMatrixCalled(), 1);
    KRATOS_EXPECT_EQ(geo_custom_element.GetCountCalculateRightHandSideCalled(), 1);

    KRATOS_EXPECT_EQ(geo_custom_condition.GetCountCalculateLeftHandSideCalled(), 1);
    KRATOS_EXPECT_EQ(geo_custom_condition.GetCountCalculateMassMatrixCalled(), 1);
    KRATOS_EXPECT_EQ(geo_custom_condition.GetCalculateDampingMatrixCalled(), 1);

    // rhs for conditions is called each solution step and during initialisation
    if (UseIterations) {
        // Two iterations are performed per solution step, thus rhs for conditions is called twice each solution step step and during initialisation
        KRATOS_EXPECT_EQ(geo_custom_condition.GetCountCalculateRightHandSideCalled(), n_steps * 2 + 1);
    } else {
        KRATOS_EXPECT_EQ(geo_custom_condition.GetCountCalculateRightHandSideCalled(), n_steps + 1);
    }

    // check results
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(calculated_displacement_x, rExpectedDisplacementX, 0.01);
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(calculated_displacement_y, rExpectedDisplacementY, 0.01);
}

/// <summary>
/// Tests NewtonRaphsonLinearElasticDynamic with a calculated initial second derivative and no newton raphson iterations, using a direct solver
/// </summary>
KRATOS_TEST_CASE_IN_SUITE(NewtonRaphsonLinearElasticDynamicCalculatedInitialAccelerationNoIterations,
                          KratosGeoMechanicsFastSuite)
{
    // set up expected results
    const std::vector<double> expected_displacement_x = {
        0.00673, 0.0505, 0.189, 0.485, 0.961, 1.58, 2.23, 2.76, 3.0, 2.85, 2.28, 1.40};
    const std::vector<double> expected_displacement_y = {0.364, 1.35, 2.68, 4.00, 4.95, 5.34,
                                                         5.13,  4.48, 3.64, 2.90, 2.44, 2.31};

    const double delta_time                     = 0.28;
    const bool   use_iterations                 = false;
    const bool   calculate_initial_acceleration = true;
    const bool   use_direct_solver              = true;
    const bool   converge_outside               = false;

    TestNewtonRaphsonLinearElasticDynamic(delta_time, use_iterations, calculate_initial_acceleration,
                                          use_direct_solver, converge_outside,
                                          expected_displacement_x, expected_displacement_y);
}

/// <summary>
/// Tests NewtonRaphsonLinearElasticDynamic with a calculated initial second derivative and no newton raphson iterations, using an iterative solver
/// </summary>
KRATOS_TEST_CASE_IN_SUITE(NewtonRaphsonLinearElasticDynamicCalculatedInitialAccelerationNoIterationsIterativeSolver,
                          KratosGeoMechanicsFastSuite)
{
    // set up expected results
    const std::vector<double> expected_displacement_x = {
        0.00673, 0.0505, 0.189, 0.485, 0.961, 1.58, 2.23, 2.76, 3.0, 2.85, 2.28, 1.40};
    const std::vector<double> expected_displacement_y = {0.364, 1.35, 2.68, 4.00, 4.95, 5.34,
                                                         5.13,  4.48, 3.64, 2.90, 2.44, 2.31};

    const double delta_time                     = 0.28;
    const bool   use_iterations                 = false;
    const bool   calculate_initial_acceleration = true;
    const bool   use_direct_solver              = false;
    const bool   converge_outside               = false;

    TestNewtonRaphsonLinearElasticDynamic(delta_time, use_iterations, calculate_initial_acceleration,
                                          use_direct_solver, converge_outside,
                                          expected_displacement_x, expected_displacement_y);
}

/// <summary>
/// Tests NewtonRaphsonLinearElasticDynamic with a calculated initial second derivative and multiple newton raphson iterations, using an direct solver
/// </summary>
KRATOS_TEST_CASE_IN_SUITE(NewtonRaphsonLinearElasticDynamicCalculatedInitialAccelerationWithIterations,
                          KratosGeoMechanicsFastSuite)
{
    // set up expected results
    const std::vector<double> expected_displacement_x = {
        0.00673, 0.0505, 0.189, 0.485, 0.961, 1.58, 2.23, 2.76, 3.0, 2.85, 2.28, 1.40};
    const std::vector<double> expected_displacement_y = {0.364, 1.35, 2.68, 4.00, 4.95, 5.34,
                                                         5.13,  4.48, 3.64, 2.90, 2.44, 2.31};

    const double delta_time                     = 0.28;
    const bool   use_iterations                 = true;
    const bool   calculate_initial_acceleration = true;
    const bool   use_direct_solver              = true;
    const bool   converge_outside               = false;

    TestNewtonRaphsonLinearElasticDynamic(delta_time, use_iterations, calculate_initial_acceleration,
                                          use_direct_solver, converge_outside,
                                          expected_displacement_x, expected_displacement_y);
}

/// <summary>
/// Tests NewtonRaphsonLinearElasticDynamic with a 0 initial second derivative and no newton raphson iterations, using an direct solver
/// </summary>
KRATOS_TEST_CASE_IN_SUITE(NewtonRaphsonLinearElasticDynamicZeroInitialAccelerationNoIterations, KratosGeoMechanicsFastSuite)
{
    // set up expected results
    const std::vector<double> expected_displacement_x = {0.996, 1.01, 0.982, 1.02, 0.969, 1.04,
                                                         0.957, 1.05, 0.947, 1.06, 0.940, 1.06};
    const std::vector<double> expected_displacement_y = {2.99, 3.02, 2.97, 3.04, 2.95, 3.06,
                                                         2.93, 3.08, 2.91, 3.09, 2.90, 3.11};

    const double delta_time                     = 28;
    const bool   use_iterations                 = false;
    const bool   calculate_initial_acceleration = false;
    const bool   use_direct_solver              = true;
    const bool   converge_outside               = false;

    TestNewtonRaphsonLinearElasticDynamic(delta_time, use_iterations, calculate_initial_acceleration,
                                          use_direct_solver, converge_outside,
                                          expected_displacement_x, expected_displacement_y);
}

/// <summary>
/// Tests NewtonRaphsonLinearElasticDynamic with a calculated initial second derivative and multiple newton raphson iterations, using an direct solver. Convergence
/// is checked outside the solution step. Such that SolveSolutionStep is called multiple times per solution step
/// </summary>
KRATOS_TEST_CASE_IN_SUITE(NewtonRaphsonLinearElasticDynamicCalculatedInitialAccelerationWithIterationsConvergenceOutside,
                          KratosGeoMechanicsFastSuite)
{
    // set up expected results
    const std::vector<double> expected_displacement_x = {
        0.00673, 0.0505, 0.189, 0.485, 0.961, 1.58, 2.23, 2.76, 3.0, 2.85, 2.28, 1.40};
    const std::vector<double> expected_displacement_y = {0.364, 1.35, 2.68, 4.00, 4.95, 5.34,
                                                         5.13,  4.48, 3.64, 2.90, 2.44, 2.31};

    const double delta_time                     = 0.28;
    const bool   use_iterations                 = true;
    const bool   calculate_initial_acceleration = true;
    const bool   use_direct_solver              = true;
    const bool   converge_outside               = true;

    TestNewtonRaphsonLinearElasticDynamic(delta_time, use_iterations, calculate_initial_acceleration,
                                          use_direct_solver, converge_outside,
                                          expected_displacement_x, expected_displacement_y);
}

} // namespace Kratos::Testing
