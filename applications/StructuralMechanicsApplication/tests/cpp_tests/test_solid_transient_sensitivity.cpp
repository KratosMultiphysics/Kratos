// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Michael Andre, https://github.com/msandre
//

// System includes
#include <functional>
#include <stack>
#include <string>
#include <unordered_map>
#include <cmath>
#include <algorithm>

// External includes

// Project includes
#include "containers/model.h"
#include "linear_solvers/skyline_lu_custom_scalar_solver.h"
#include "solving_strategies/convergencecriterias/residual_criteria.h"
#include "solving_strategies/schemes/residual_based_adjoint_bossak_scheme.h"
#include "solving_strategies/schemes/residual_based_bossak_displacement_scheme.hpp"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"
#include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"
#include "spaces/ublas_space.h"
#include "testing/testing.h"
#include "utilities/sensitivity_builder.h"

// Application includes
#include "custom_response_functions/response_utilities/adjoint_nodal_displacement_response_function.h"
#include "tests/cpp_tests/bits/adjoint_test_model_part_factory.h"
#include "tests/cpp_tests/bits/test_model_part_factory.h"

namespace Kratos
{
namespace
{
namespace test_solid_transient_sensitivity_cpp
{ // unity build unity guard
struct PrimalTestSolver
{
    std::function<ModelPart&()> ModelPartFactory;
    unsigned ResponseNodeId;
    double TotalTime;
    double TimeStepSize;

    double CalculateResponseValue();

    double CalculateResponseValue(unsigned NodeToPerturb, char Direction, double Perturbation);
};

struct AdjointTestSolver
{
    std::function<ModelPart&()> ModelPartFactory;
    unsigned ResponseNodeId;
    double TotalTime;
    double TimeStepSize;

    double CalculateSensitivity(unsigned NodeToPerturb, char Direction);
};

} // namespace test_solid_transient_sensitivity_cpp
} // namespace

namespace Testing
{
KRATOS_TEST_CASE_IN_SUITE(TotalLagrangian2D3_SaintVenantPlaneStrain_TransientSensitivity, KratosStructuralMechanicsFastSuite)
{
    if (!KratosComponents<ConstitutiveLaw>::Has("KirchhoffSaintVenantPlaneStrain2DLaw")) {
        // this test can only be run if the ConstitutiveLawsApp is imported
        return;
    }

    using test_solid_transient_sensitivity_cpp::AdjointTestSolver;
    using test_solid_transient_sensitivity_cpp::PrimalTestSolver;
    Model this_model;
    auto model_part_factory = [&this_model]() -> ModelPart& {
        return CreateStructuralMechanicsTestModelPart(
            &this_model, KratosComponents<Element>::Get("TotalLagrangianElement2D3N"),
            KratosComponents<ConstitutiveLaw>::Get(
                "KirchhoffSaintVenantPlaneStrain2DLaw"),
            [](ModelPart* pModelPart) {
                for (unsigned i : {1, 3})
                {
                    pModelPart->GetNode(i).Fix(DISPLACEMENT_X);
                    pModelPart->GetNode(i).Fix(DISPLACEMENT_Y);
                }
            });
    };
    const double time_step_size = 0.004;
    const double total_time = 0.015; // approx. 1/4 of a period.
    const unsigned response_node_id = 2;
    PrimalTestSolver solver{model_part_factory, response_node_id, total_time, time_step_size};
    AdjointTestSolver adjoint_solver{model_part_factory, response_node_id,
                                     total_time, time_step_size};
    const double delta = 1e-5;
    const double response_value0 = solver.CalculateResponseValue();
    for (unsigned i_node : {1, 2, 3})
    {
        for (char dir : {'x', 'y'})
        {
            const double response_value1 =
                solver.CalculateResponseValue(i_node, dir, delta);
            const double finite_diff_sensitivity =
                -(response_value1 - response_value0) / delta;
            const double adjoint_sensitivity =
                adjoint_solver.CalculateSensitivity(i_node, dir);
            const double tol = std::max(0.001 * std::abs(finite_diff_sensitivity), 1e-10);
            KRATOS_EXPECT_NEAR(adjoint_sensitivity, finite_diff_sensitivity, tol);
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(TotalLagrangian3D8_SaintVenant_TransientSensitivity, KratosStructuralMechanicsFastSuite)
{
    if (!KratosComponents<ConstitutiveLaw>::Has("KirchhoffSaintVenant3DLaw")) {
        // this test can only be run if the ConstitutiveLawsApp is imported
        return;
    }

    using test_solid_transient_sensitivity_cpp::AdjointTestSolver;
    using test_solid_transient_sensitivity_cpp::PrimalTestSolver;
    Model this_model;
    auto model_part_factory = [&this_model]() -> ModelPart& {
        return CreateStructuralMechanicsTestModelPart(
            &this_model, KratosComponents<Element>::Get("TotalLagrangianElement3D8N"),
            KratosComponents<ConstitutiveLaw>::Get("KirchhoffSaintVenant3DLaw"),
            [](ModelPart* pModelPart) {
                for (unsigned i : {1, 5, 8, 4})
                {
                    pModelPart->GetNode(i).Fix(DISPLACEMENT_X);
                    pModelPart->GetNode(i).Fix(DISPLACEMENT_Y);
                    pModelPart->GetNode(i).Fix(DISPLACEMENT_Z);
                }
                pModelPart->GetProperties(1)[VOLUME_ACCELERATION](1) = 10.;
            });
    };
    const double time_step_size = 0.05;
    const double total_time = 0.3; // approx. 1/4 of a period.
    const unsigned response_node_id = 6;
    PrimalTestSolver solver{model_part_factory, response_node_id, total_time, time_step_size};
    AdjointTestSolver adjoint_solver{model_part_factory, response_node_id,
                                     total_time, time_step_size};
    const double delta = 1e-5;
    const double response_value0 = solver.CalculateResponseValue();
    for (unsigned i_node : {1, 8})
    {
        for (char dir : {'x', 'y'})
        {
            const double response_value1 =
                solver.CalculateResponseValue(i_node, dir, delta);
            const double finite_diff_sensitivity =
                -(response_value1 - response_value0) / delta;
            const double adjoint_sensitivity =
                adjoint_solver.CalculateSensitivity(i_node, dir);
            const double tol = std::max(0.01 * std::abs(finite_diff_sensitivity), 1e-10);
            KRATOS_EXPECT_NEAR(adjoint_sensitivity, finite_diff_sensitivity, tol);
        }
    }
}
} // namespace Testing

namespace // cpp internals
{
namespace test_solid_transient_sensitivity_cpp
{ // unity build unity guard
typedef TUblasSparseSpace<double> SparseSpaceType;
typedef TUblasDenseSpace<double> LocalSpaceType;
typedef LinearSolver<SparseSpaceType, LocalSpaceType> LinearSolverType;
typedef Scheme<SparseSpaceType, LocalSpaceType> SchemeType;
typedef ConvergenceCriteria<SparseSpaceType, LocalSpaceType> ConvergenceCriteriaType;
typedef ImplicitSolvingStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType> SolvingStrategyType;

AdjointResponseFunction::Pointer ResponseFunctionFactory(ModelPart* pModelPart, unsigned ResponseNodeId)
{
    Parameters params{R"({"traced_dof": "DISPLACEMENT", "direction": [0.0,1.0,0.0], "gradient_mode": "semi_analytic", "step_size": 1e-2})"};
    std::string response_mp_label = "response_mp";
    if (!pModelPart->HasSubModelPart(response_mp_label)) {
        ModelPart& r_sub = pModelPart->CreateSubModelPart(response_mp_label);
        r_sub.AddNode(pModelPart->pGetNode(ResponseNodeId));
    }
    params.AddEmptyValue("response_part_name");
    params["response_part_name"].SetString(response_mp_label);
    return Kratos::make_shared<AdjointNodalDisplacementResponseFunction>(*pModelPart, params);
}

SolvingStrategyType::Pointer CreatePrimalSolvingStrategy(ModelPart* pModelPart)
{
    LinearSolverType::Pointer p_linear_solver =
        Kratos::make_shared<SkylineLUCustomScalarSolver<SparseSpaceType, LocalSpaceType>>();
    SchemeType::Pointer p_scheme =
        Kratos::make_shared<ResidualBasedBossakDisplacementScheme<SparseSpaceType, LocalSpaceType>>(0.0);
    ConvergenceCriteriaType::Pointer p_conv_criteria =
        Kratos::make_shared<ResidualCriteria<SparseSpaceType, LocalSpaceType>>(1e-6, 1e-9);
    return Kratos::make_shared<ResidualBasedNewtonRaphsonStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType>>(
        *pModelPart, p_scheme, p_linear_solver, p_conv_criteria, 30, true, false, true);
}

void SolvePrimal(ModelPart* pModelPart,
                 double DeltaTime,
                 double TotalTime,
                 std::function<void(ModelPart*)> CallerFun)
{
    auto p_solver = CreatePrimalSolvingStrategy(pModelPart);
    p_solver->SetEchoLevel(0);
    p_solver->Initialize();
    const double start_time = 0.0;
    pModelPart->CloneTimeStep(start_time - DeltaTime);
    pModelPart->CloneTimeStep(start_time);
    for (double current_time = start_time; current_time < TotalTime;)
    {
        current_time += DeltaTime;
        pModelPart->CloneTimeStep(current_time);
        p_solver->Solve();
        CallerFun(pModelPart);
    }
}

double CalculateResponseValue(ModelPart* pModelPart, double DeltaTime, double TotalTime, unsigned ResponseNodeId)
{
    auto p_response_function = ResponseFunctionFactory(pModelPart, ResponseNodeId);
    p_response_function->Initialize();
    double response_value = 0.;
    auto sum_response_value = [&](ModelPart* pModelPart) {
        response_value += pModelPart->GetProcessInfo()[DELTA_TIME] *
                          p_response_function->CalculateValue(*pModelPart);
    };
    SolvePrimal(pModelPart, DeltaTime, TotalTime, sum_response_value);
    return response_value;
}

double PrimalTestSolver::CalculateResponseValue()
{
    ModelPart& primal_model_part = ModelPartFactory();
    return test_solid_transient_sensitivity_cpp::CalculateResponseValue(
        &primal_model_part, TimeStepSize, TotalTime, ResponseNodeId);
}

unsigned DirectionIndex(char Direction)
{
    KRATOS_ERROR_IF(Direction != 'x' && Direction != 'y' && Direction != 'z')
        << "invalid direction: '" << Direction << "'";
    if (Direction == 'x')
        return 0;
    else if (Direction == 'y')
        return 1;
    else
        return 2;
}

double PrimalTestSolver::CalculateResponseValue(unsigned NodeToPerturb, char Direction, double Perturbation)
{
    KRATOS_ERROR_IF(Perturbation <= 0.) << "invalid perturbation: " << Perturbation;
    ModelPart& primal_model_part = ModelPartFactory();
    const unsigned i_dir = DirectionIndex(Direction);
    primal_model_part.GetNode(NodeToPerturb).GetInitialPosition()[i_dir] += Perturbation;
    return test_solid_transient_sensitivity_cpp::CalculateResponseValue(
        &primal_model_part, TimeStepSize, TotalTime, ResponseNodeId);
}

SolvingStrategyType::Pointer CreateAdjointSolvingStrategy(ModelPart* pModelPart,
                                                          AdjointResponseFunction::Pointer pResponseFunction)
{
    LinearSolverType::Pointer p_linear_solver =
        Kratos::make_shared<SkylineLUCustomScalarSolver<SparseSpaceType, LocalSpaceType>>();
    Parameters settings(R"({ "alpha_bossak": 0.0 })");
    SchemeType::Pointer p_adjoint_scheme =
        Kratos::make_shared<ResidualBasedAdjointBossakScheme<SparseSpaceType, LocalSpaceType>>(
            settings, pResponseFunction);
    return Kratos::make_shared<ResidualBasedLinearStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType>>(
        *pModelPart, p_adjoint_scheme, p_linear_solver);
}

struct PrimalSolution
{
    using NodeIdType = unsigned;
    double time;
    std::unordered_map<NodeIdType, array_1d<double, 3>> disp;
    std::unordered_map<NodeIdType, array_1d<double, 3>> vel;
    std::unordered_map<NodeIdType, array_1d<double, 3>> acc;
};

void SolveAdjoint(ModelPart* pAdjointModelPart,
                  std::stack<PrimalSolution>* pPrimalSolution,
                  unsigned ResponseNodeId)
{
    auto p_response_function = ResponseFunctionFactory(pAdjointModelPart, ResponseNodeId);
    auto p_adjoint_solver =
        CreateAdjointSolvingStrategy(pAdjointModelPart, p_response_function);
    p_adjoint_solver->SetEchoLevel(0);
    p_adjoint_solver->Initialize();
    SensitivityBuilder sensitivity_builder(
        Parameters(R"(
            {
                "nodal_solution_step_sensitivity_variables": ["SHAPE_SENSITIVITY"],
                "build_mode" : "integrate"
            })"),
        *pAdjointModelPart, p_response_function);
    sensitivity_builder.Initialize();
    while (!pPrimalSolution->empty())
    {
        auto p_sol = &pPrimalSolution->top();
        pAdjointModelPart->CloneTimeStep(p_sol->time);
        for (auto& r_adjoint_node : pAdjointModelPart->Nodes())
        {
            r_adjoint_node.FastGetSolutionStepValue(DISPLACEMENT) =
                p_sol->disp[r_adjoint_node.Id()];
            r_adjoint_node.X() = r_adjoint_node.X0() +
                                 r_adjoint_node.FastGetSolutionStepValue(DISPLACEMENT_X);
            r_adjoint_node.Y() = r_adjoint_node.Y0() +
                                 r_adjoint_node.FastGetSolutionStepValue(DISPLACEMENT_Y);
            r_adjoint_node.FastGetSolutionStepValue(VELOCITY) =
                p_sol->vel[r_adjoint_node.Id()];
            r_adjoint_node.FastGetSolutionStepValue(ACCELERATION) =
                p_sol->acc[r_adjoint_node.Id()];
        }
        p_adjoint_solver->Solve();
        sensitivity_builder.UpdateSensitivities();
        pPrimalSolution->pop();
    }
}

double AdjointTestSolver::CalculateSensitivity(unsigned NodeToPerturb, char Direction)
{
    ModelPart& primal_model_part = ModelPartFactory();
    std::stack<PrimalSolution> sol;
    auto record_primal = [&](ModelPart* pModelPart) {
        sol.push(PrimalSolution{pModelPart->GetProcessInfo()[TIME]});
        auto& current_sol = sol.top();
        for (const auto& r_node : pModelPart->Nodes())
        {
            current_sol.disp[r_node.Id()] = r_node.FastGetSolutionStepValue(DISPLACEMENT);
            current_sol.vel[r_node.Id()] = r_node.FastGetSolutionStepValue(VELOCITY);
            current_sol.acc[r_node.Id()] = r_node.FastGetSolutionStepValue(ACCELERATION);
        };
    };
    SolvePrimal(&primal_model_part, TimeStepSize, TotalTime, record_primal);
    ModelPart& adjoint_model_part =
        CreateStructuralMechanicsAdjointTestModelPart(&primal_model_part);
    SolveAdjoint(&adjoint_model_part, &sol, ResponseNodeId);
    const unsigned i_dir = DirectionIndex(Direction);
    return adjoint_model_part.GetNode(NodeToPerturb).FastGetSolutionStepValue(SHAPE_SENSITIVITY)[i_dir];
}

} // namespace test_solid_transient_sensitivity_cpp
} // namespace
} // namespace Kratos
