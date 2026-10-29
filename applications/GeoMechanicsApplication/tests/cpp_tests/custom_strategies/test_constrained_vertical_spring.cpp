// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//

#include "containers/model.h"
#include "custom_elements/nodal_elements/spring_damper_element.hpp"
#include "custom_utilities/solving_strategy_factory.hpp"
#include "spaces/default_spaces.h"
#include "structural_mechanics_application.h"
#include "structural_mechanics_application_variables.h"
#include "test_setup_utilities/model_setup_utilities.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"

#include <cmath>
#include <string>
#include <vector>

namespace Kratos::Testing
{
namespace
{

using SpringSparseSpace  = TDefaultSparseSpace<double>;
using SpringDenseSpace   = TDefaultDenseSpace<double>;
using SpringLinearSolver = LinearSolver<SpringSparseSpace, SpringDenseSpace>;
using SpringStrategyFactory = SolvingStrategyFactory<SpringSparseSpace, SpringDenseSpace, SpringLinearSolver>;

constexpr double lower_stiffness    = 1000.0;
constexpr double displacement_ratio = 0.5;
constexpr double tolerance          = 1.0e-10;
constexpr double cubic_stiffness    = 1.0e7;

// Test-only extension: f(delta) = k*delta + alpha*delta^3 in the vertical direction.
class CubicVerticalSpringElement : public SpringDamperElement<2>
{
public:
    using SpringDamperElement<2>::SpringDamperElement;

    void CalculateLocalSystem(MatrixType& rLhs, VectorType& rRhs, const ProcessInfo& rProcessInfo) override
    {
        CalculateLeftHandSide(rLhs, rProcessInfo);
        CalculateRightHandSide(rRhs, rProcessInfo);
    }

    void CalculateLeftHandSide(MatrixType& rLhs, const ProcessInfo& rProcessInfo) override
    {
        SpringDamperElement<2>::CalculateLeftHandSide(rLhs, rProcessInfo);
        const auto extension = GetExtension();
        const auto tangent   = 3.0 * cubic_stiffness * extension * extension;
        rLhs(1, 1) += tangent;
        rLhs(4, 4) += tangent;
        rLhs(1, 4) -= tangent;
        rLhs(4, 1) -= tangent;
    }

    void CalculateRightHandSide(VectorType& rRhs, const ProcessInfo& rProcessInfo) override
    {
        SpringDamperElement<2>::CalculateRightHandSide(rRhs, rProcessInfo);
        const auto extension = GetExtension();
        const auto force     = cubic_stiffness * extension * extension * extension;
        rRhs[1] += force;
        rRhs[4] -= force;
    }

private:
    double GetExtension() const
    {
        return GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_Y) -
               GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_Y);
    }
};

ModelPart& CreateConstrainedVerticalSpring(Model& rModel, double UpperStiffness, bool Nonlinear = false)
{
    auto& r_model_part                         = rModel.CreateModelPart("Spring", 2);
    r_model_part.GetProcessInfo()[DOMAIN_SIZE] = 2;
    r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    r_model_part.AddNodalSolutionStepVariable(REACTION);
    r_model_part.AddNodalSolutionStepVariable(ROTATION);
    r_model_part.AddNodalSolutionStepVariable(REACTION_MOMENT);
    r_model_part.AddNodalSolutionStepVariable(POINT_LOAD);

    ModelSetupUtilities::CreateNodes(
        r_model_part, {{1, {0.0, 0.0, 0.0}}, {2, {0.0, 1.0, 0.0}}, {3, {0.0, 2.0, 0.0}}});
    for (auto& r_node : r_model_part.Nodes()) {
        r_node.AddDof(DISPLACEMENT_X, REACTION_X);
        r_node.AddDof(DISPLACEMENT_Y, REACTION_Y);
        r_node.AddDof(DISPLACEMENT_Z, REACTION_Z);
        r_node.AddDof(ROTATION_Z, REACTION_MOMENT_Z);
        r_node.Fix(DISPLACEMENT_X);
        r_node.Fix(DISPLACEMENT_Z);
        r_node.Fix(ROTATION_Z);
    }
    r_model_part.GetNode(1).Fix(DISPLACEMENT_Y);

    const auto p_properties   = r_model_part.CreateNewProperties(1);
    auto       p_lower_spring = r_model_part.CreateNewElement(
        "SpringDamperElement2D", 1, std::vector<ModelPart::IndexType>{1, 2}, p_properties);
    auto p_upper_spring = r_model_part.CreateNewElement(
        "SpringDamperElement2D", 2, std::vector<ModelPart::IndexType>{2, 3}, p_properties);
    p_lower_spring->SetValue(NODAL_DISPLACEMENT_STIFFNESS, array_1d<double, 3>{0.0, lower_stiffness, 0.0});
    p_upper_spring->SetValue(NODAL_DISPLACEMENT_STIFFNESS, array_1d<double, 3>{0.0, UpperStiffness, 0.0});
    p_lower_spring->SetValue(NODAL_ROTATIONAL_STIFFNESS, array_1d<double, 3>{0.0, 0.0, 0.0});
    p_upper_spring->SetValue(NODAL_ROTATIONAL_STIFFNESS, array_1d<double, 3>{0.0, 0.0, 0.0});

    if (Nonlinear) {
        for (const auto& p_spring : {p_lower_spring, p_upper_spring}) {
            auto p_cubic_spring = Kratos::make_intrusive<CubicVerticalSpringElement>(
                p_spring->Id(), p_spring->pGetGeometry(), p_properties);
            p_cubic_spring->GetData() = p_spring->GetData();
            r_model_part.RemoveElement(p_spring->Id());
            r_model_part.AddElement(p_cubic_spring);
        }
    }

    r_model_part.CreateNewCondition("PointLoadCondition2D1N", 1, std::vector<ModelPart::IndexType>{3}, p_properties);
    r_model_part.CreateNewMasterSlaveConstraint(
        "LinearMasterSlaveConstraint", 1, r_model_part.GetNode(3), DISPLACEMENT_Y,
        r_model_part.GetNode(2), DISPLACEMENT_Y, displacement_ratio, 0.0);
    return r_model_part;
}

void SolveConstrainedVerticalSpring(const std::string& rStrategyType,
                                    const std::string& rQuasiNewtonType,
                                    double             UpperStiffness,
                                    bool               Nonlinear = false)
{
    SCOPED_TRACE(rStrategyType + " / " + rQuasiNewtonType + " / upper stiffness = " + std::to_string(UpperStiffness));
    SCOPED_TRACE(Nonlinear ? "cubic springs" : "linear springs");
    Model model;
    auto& r_model_part = CreateConstrainedVerticalSpring(model, UpperStiffness, Nonlinear);
    auto  settings     = Parameters{R"({
        "strategy_type": "newton_raphson",
        "quasi_newton_type": "broyden",
        "solution_type": "static",
        "scheme_type": "static",
        "block_builder": true,
        "linear_solver_settings": {
            "solver_type": "LinearSolversApplication.sparse_lu",
            "scaling": false
        },
        "convergence_criterion": "residual_criterion",
        "residual_relative_tolerance": 1.0e-10,
        "residual_absolute_tolerance": 1.0e-10,
        "max_iterations": 30,
        "compute_reactions": true,
        "reform_dofs_at_each_step": false,
        "move_mesh_flag": false,
        "echo_level": 0
    })"};
    settings["strategy_type"].SetString(rStrategyType);
    settings["quasi_newton_type"].SetString(rQuasiNewtonType);
    auto p_strategy = SpringStrategyFactory::Create(settings, r_model_part);
    ASSERT_NE(p_strategy, nullptr);
    KRATOS_EXPECT_EQ(p_strategy->Check(), 0);
    p_strategy->Initialize();

    // Substituting u_mid = a*u_top gives k_eff = k_lower*a^2 + k_upper*(1-a)^2.
    const auto effective_stiffness = lower_stiffness * displacement_ratio * displacement_ratio +
                                     UpperStiffness * (1.0 - displacement_ratio) * (1.0 - displacement_ratio);
    const auto number_of_steps = Nonlinear ? 1 : 2;
    for (int step = 1; step <= number_of_steps; ++step) {
        r_model_part.CloneTimeStep(static_cast<double>(step));
        r_model_part.GetProcessInfo()[STEP] = step;
        const auto load                     = Nonlinear ? 6.25 : 10.0 * step;
        r_model_part.GetNode(3).FastGetSolutionStepValue(POINT_LOAD_Y) = load;
        p_strategy->InitializeSolutionStep();
        p_strategy->Predict();
        ASSERT_TRUE(p_strategy->SolveSolutionStep());
        p_strategy->FinalizeSolutionStep();

        const auto bottom = r_model_part.GetNode(1).FastGetSolutionStepValue(DISPLACEMENT_Y);
        const auto middle = r_model_part.GetNode(2).FastGetSolutionStepValue(DISPLACEMENT_Y);
        const auto top    = r_model_part.GetNode(3).FastGetSolutionStepValue(DISPLACEMENT_Y);
        KRATOS_EXPECT_NEAR(bottom, 0.0, tolerance);
        // The nonlinear Newton case deliberately characterizes the empty-mask false convergence.
        const auto premature_convergence = Nonlinear && rStrategyType == "newton_raphson";
        const auto expected_top = Nonlinear && !premature_convergence ? 0.01 : load / effective_stiffness;
        KRATOS_EXPECT_NEAR(top, expected_top, tolerance);
        KRATOS_EXPECT_NEAR(middle, displacement_ratio * top, tolerance);
        const auto lower_force =
            lower_stiffness * middle + (Nonlinear ? cubic_stiffness * middle * middle * middle : 0.0);
        KRATOS_EXPECT_NEAR(r_model_part.GetNode(1).FastGetSolutionStepValue(REACTION_Y), -lower_force, tolerance);
        KRATOS_EXPECT_TRUE(r_model_part.GetProcessInfo()[NL_ITERATION_NUMBER] <
                           settings["max_iterations"].GetInt());
        KRATOS_EXPECT_TRUE(std::isfinite(r_model_part.GetProcessInfo()[RESIDUAL_NORM]));
        KRATOS_EXPECT_TRUE(std::isfinite(r_model_part.GetProcessInfo()[CONVERGENCE_RATIO]));

        Vector lower_rhs;
        Vector upper_rhs;
        r_model_part.GetElement(1).CalculateRightHandSide(lower_rhs, r_model_part.GetProcessInfo());
        r_model_part.GetElement(2).CalculateRightHandSide(upper_rhs, r_model_part.GetProcessInfo());
        const auto raw_middle_residual = lower_rhs[4] + upper_rhs[1];
        const auto raw_top_residual    = load + upper_rhs[4];
        const auto reduced_residual = raw_top_residual + displacement_ratio * raw_middle_residual;
        if (premature_convergence) {
            KRATOS_EXPECT_EQ(r_model_part.GetProcessInfo()[NL_ITERATION_NUMBER], 1);
            KRATOS_EXPECT_NEAR(top, 0.0125, tolerance);
            KRATOS_EXPECT_NEAR(reduced_residual, -2.44140625, tolerance);
            KRATOS_EXPECT_NEAR(r_model_part.GetProcessInfo()[RESIDUAL_NORM], 0.0, tolerance);
        } else {
            KRATOS_EXPECT_NEAR(reduced_residual, 0.0, 1.0e-8);
            if (Nonlinear) {
                KRATOS_EXPECT_TRUE(r_model_part.GetProcessInfo()[NL_ITERATION_NUMBER] > 1);
                KRATOS_EXPECT_NEAR(0.5 * lower_stiffness * top + cubic_stiffness * top * top * top / 8.0,
                                   load, 1.0e-8);
            } else {
                KRATOS_EXPECT_NEAR(raw_top_residual,
                                   load * (lower_stiffness - UpperStiffness) / (lower_stiffness + UpperStiffness),
                                   tolerance);
            }
        }
    }
    p_strategy->Clear();
}

} // namespace

KRATOS_TEST_CASE_IN_SUITE(ConstrainedVerticalSpring_NewtonRaphson, KratosGeoMechanicsIntegrationSuite)
{
    auto p_structural_application = std::make_shared<KratosStructuralMechanicsApplication>();
    this->ImportApplicationIntoKernel(p_structural_application);
    for (const auto upper_stiffness : {1000.0, 3000.0}) {
        SolveConstrainedVerticalSpring("newton_raphson", "broyden", upper_stiffness);
    }
    SolveConstrainedVerticalSpring("newton_raphson", "broyden", 1000.0, true);
}

KRATOS_TEST_CASE_IN_SUITE(ConstrainedVerticalSpring_QuasiNewtonBroyden, KratosGeoMechanicsIntegrationSuite)
{
    auto p_structural_application = std::make_shared<KratosStructuralMechanicsApplication>();
    this->ImportApplicationIntoKernel(p_structural_application);
    for (const auto upper_stiffness : {1000.0, 3000.0}) {
        SolveConstrainedVerticalSpring("quasi_newton", "broyden", upper_stiffness);
    }
    SolveConstrainedVerticalSpring("quasi_newton", "broyden", 1000.0, true);
}

KRATOS_TEST_CASE_IN_SUITE(ConstrainedVerticalSpring_QuasiNewtonLbfgs, KratosGeoMechanicsIntegrationSuite)
{
    auto p_structural_application = std::make_shared<KratosStructuralMechanicsApplication>();
    this->ImportApplicationIntoKernel(p_structural_application);
    for (const auto upper_stiffness : {1000.0, 3000.0}) {
        SolveConstrainedVerticalSpring("quasi_newton", "lbfgs", upper_stiffness);
    }
    SolveConstrainedVerticalSpring("quasi_newton", "lbfgs", 1000.0, true);
}

} // namespace Kratos::Testing
