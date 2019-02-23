// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:
//

// System includes
#include <string>
#include <vector>

// External includes

// Project includes
#include "containers/model.h"
#include "includes/shared_pointers.h"
#include "linear_solvers/skyline_lu_custom_scalar_solver.h"
#include "solving_strategies/convergencecriterias/residual_criteria.h"
#include "solving_strategies/schemes/residual_based_adjoint_static_scheme.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
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
namespace smtsss
{ // cotire unity guard
using SparseSpaceType = TUblasSparseSpace<double>;
using LocalSpaceType = TUblasDenseSpace<double>;
using LinearSolverType = LinearSolver<SparseSpaceType, LocalSpaceType>;
using SchemeType = Scheme<SparseSpaceType, LocalSpaceType>;
using ConvergenceCriteriaType = ConvergenceCriteria<SparseSpaceType, LocalSpaceType>;
using SolvingStrategyType =
    SolvingStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType>;

class PrimalTestSolver
{
public:
    explicit PrimalTestSolver(ModelPart* pPrimalModelPart, unsigned ResponseNodeId);

    double CalculateResponseValue();

    double CalculateResponseValue(unsigned NodeToPerturb, char Direction, double Perturbation);

private:
    SolvingStrategyType::Pointer CreateSolvingStrategy();

    ModelPart* mpPrimalModelPart;
    unsigned mResponseNodeId;
};

class AdjointTestSolver
{
public:
    explicit AdjointTestSolver(ModelPart* pAdjointModelPart, unsigned ResponseNodeId);

    double CalculateSensitivityValue(unsigned Node, char Direction) const;

private:
    SolvingStrategyType::Pointer CreateAdjointSolvingStrategy(AdjointResponseFunction::Pointer pResponseFunction);

    ModelPart* mpAdjointModelPart;
    unsigned mResponseNodeId;
};
} // namespace smtsss
} // namespace

namespace Testing
{
KRATOS_TEST_CASE_IN_SUITE(TotalLagrangian2D3_StaticSensitivity, KratosStructuralMechanicsFastSuite)
{
    Model this_model;
    ModelPart& primal_model_part = CreateStructuralMechanicsTestModelPart(
        &this_model, KratosComponents<Element>::Get("TotalLagrangianElement2D3N"),
        KratosComponents<ConstitutiveLaw>::Get("LinearElasticPlaneStrain2DLaw"),
        [](ModelPart* pModelPart) {
            pModelPart->GetNode(1).Fix(DISPLACEMENT_X);
            pModelPart->GetNode(1).Fix(DISPLACEMENT_Y);
            pModelPart->GetNode(3).Fix(DISPLACEMENT_X);
            pModelPart->GetNode(3).Fix(DISPLACEMENT_Y);
        });
    const unsigned response_node_id = 2;
    smtsss::PrimalTestSolver solver{&primal_model_part, response_node_id};
    const double delta = 1e-7;
    const double response_value0 = solver.CalculateResponseValue();
    ModelPart& adjoint_model_part = CreateStructuralMechanicsAdjointTestModelPart(&primal_model_part);
    smtsss::AdjointTestSolver adjoint_solver{&adjoint_model_part, response_node_id};
    for (unsigned i_node : {1, 2, 3})
    {
        for (char dir : {'x', 'y'})
        {
            const double response_value1 =
                solver.CalculateResponseValue(i_node, dir, delta);
            const double finite_diff_sensitivity =
                -(response_value1 - response_value0) / delta;
            const double adjoint_sensitivity =
                adjoint_solver.CalculateSensitivityValue(i_node, dir);
            KRATOS_CHECK_NEAR(finite_diff_sensitivity, adjoint_sensitivity, 1e-8);
        }
    }
}
} // namespace Testing

namespace
{
namespace smtsss
{ // cotire unity guard
AdjointResponseFunction::Pointer ResponseFunctionFactory(ModelPart* pModelPart, unsigned ResponseNodeId)
{
    Parameters params{R"({"traced_dof": "DISPLACEMENT_Y", "gradient_mode": "semi_analytic", "step_size": 1e-2})"};
    params.AddEmptyValue("traced_node_id");
    params["traced_node_id"].SetInt(ResponseNodeId);
    return Kratos::make_shared<AdjointNodalDisplacementResponseFunction>(*pModelPart, params);
}

PrimalTestSolver::PrimalTestSolver(ModelPart* pPrimalModelPart, unsigned ResponseNodeId)
    : mpPrimalModelPart(pPrimalModelPart), mResponseNodeId(ResponseNodeId){};

double PrimalTestSolver::CalculateResponseValue()
{
    auto p_solver = CreateSolvingStrategy();
    p_solver->Initialize();
    p_solver->Solve();
    auto p_response_function = ResponseFunctionFactory(mpPrimalModelPart, mResponseNodeId);
    return p_response_function->CalculateValue(*mpPrimalModelPart);
}

double PrimalTestSolver::CalculateResponseValue(const unsigned NodeToPerturb,
                                                const char Direction,
                                                const double Perturbation)
{
    KRATOS_ERROR_IF(Direction != 'x' && Direction != 'y')
        << "invalid direction: '" << Direction << "'";
    KRATOS_ERROR_IF(Perturbation <= 0.) << "invalid perturbation: " << Perturbation;
    const unsigned i_dir = (Direction == 'x') ? 0 : 1;
    auto& r_node = mpPrimalModelPart->GetNode(NodeToPerturb);
    r_node.GetInitialPosition()[i_dir] += Perturbation;
    const double response_value = CalculateResponseValue();
    r_node.GetInitialPosition()[i_dir] -= Perturbation;
    return response_value;
}

SolvingStrategyType::Pointer PrimalTestSolver::CreateSolvingStrategy()
{
    LinearSolverType::Pointer p_linear_solver =
        Kratos::make_shared<SkylineLUCustomScalarSolver<SparseSpaceType, LocalSpaceType>>();
    SchemeType::Pointer p_scheme =
        Kratos::make_shared<ResidualBasedIncrementalUpdateStaticScheme<SparseSpaceType, LocalSpaceType>>();
    ConvergenceCriteriaType::Pointer p_conv_criteria =
        Kratos::make_shared<ResidualCriteria<SparseSpaceType, LocalSpaceType>>(1e-12, 1e-14);
    return Kratos::make_shared<ResidualBasedNewtonRaphsonStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType>>(
        *mpPrimalModelPart, p_scheme, p_linear_solver, p_conv_criteria, 30,
        true, false, true);
}

AdjointTestSolver::AdjointTestSolver(ModelPart* pAdjointModelPart, unsigned ResponseNodeId)
    : mpAdjointModelPart(pAdjointModelPart), mResponseNodeId(ResponseNodeId)
{
    auto p_adjoint_response_function =
        ResponseFunctionFactory(mpAdjointModelPart, mResponseNodeId);
    auto p_adjoint_solver = CreateAdjointSolvingStrategy(p_adjoint_response_function);
    p_adjoint_solver->Initialize();
    p_adjoint_solver->Solve();
    SensitivityBuilder sensitivity_builder(
        Parameters(R"({"nodal_solution_step_sensitivity_variables": ["SHAPE_SENSITIVITY"], "build_mode": "static"})"),
        *mpAdjointModelPart, p_adjoint_response_function);
    sensitivity_builder.Initialize();
    sensitivity_builder.UpdateSensitivities();
}

double AdjointTestSolver::CalculateSensitivityValue(unsigned Node, char Direction) const
{
    KRATOS_ERROR_IF(Direction != 'x' && Direction != 'y')
        << "invalid direction: '" << Direction << "'";
    const unsigned i_dir = (Direction == 'x') ? 0 : 1;
    return mpAdjointModelPart->GetNode(Node).FastGetSolutionStepValue(SHAPE_SENSITIVITY)[i_dir];
}

SolvingStrategyType::Pointer AdjointTestSolver::CreateAdjointSolvingStrategy(
    AdjointResponseFunction::Pointer pResponseFunction)
{
    LinearSolverType::Pointer p_linear_solver =
        Kratos::make_shared<SkylineLUCustomScalarSolver<SparseSpaceType, LocalSpaceType>>();
    SchemeType::Pointer p_adjoint_scheme =
        Kratos::make_shared<ResidualBasedAdjointStaticScheme<SparseSpaceType, LocalSpaceType>>(
            pResponseFunction);
    return Kratos::make_shared<ResidualBasedLinearStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType>>(
        *mpAdjointModelPart, p_adjoint_scheme, p_linear_solver);
}

} // namespace smtsss
} // namespace
} // namespace Kratos
