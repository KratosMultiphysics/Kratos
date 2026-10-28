//     |  /           |
//     ' /   __| _` | __|  _ \   __|
//     . \  |   (   | |   (   |\__ `
//    _|\_\_|  \__,_|\__|\___/ ____/
//                    Multi-Physics
//
//   Fluid Dynamics Application
//
//   License:         BSD License
//                    Kratos default license: kratos/license.txt
//
//   Main author:     Nicolò Antonelli
//

// System includes
#include <vector>

// Project includes
#include "tests/cpp_tests/fluid_dynamics_fast_suite.h"
#include "containers/model.h"
#include "fluid_dynamics_application_variables.h"

namespace Kratos::Testing
{
namespace
{
Condition::Pointer create_condition(ModelPart& rModelPart, const double Coefficient,
                                    const bool Quadrilateral = false)
{
    rModelPart.AddNodalSolutionStepVariable(VELOCITY);
    rModelPart.AddNodalSolutionStepVariable(PRESSURE);
    rModelPart.AddNodalSolutionStepVariable(ACCELERATION);
    rModelPart.GetProcessInfo()[BDF_COEFFICIENTS] = ZeroVector(3);
    auto p_material = rModelPart.CreateNewProperties(1);
    p_material->SetValue(DYNAMIC_VISCOSITY, 1.3);
    auto p_interface_properties = rModelPart.CreateNewProperties(2);
    p_interface_properties->SetValue(PENALTY_COEFFICIENT, Coefficient);
    rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
    rModelPart.CreateNewNode(3, 1.0, 1.0, 0.0);
    rModelPart.CreateNewNode(4, 1.0, 0.0, 0.0);
    rModelPart.CreateNewNode(5, 2.0, 0.0, 0.0);
    rModelPart.CreateNewNode(6, 1.0, 1.0, 0.0);
    if (Quadrilateral) {
        // Two non-affine Q1 parents sharing the oriented segment x=1.
        rModelPart.CreateNewNode(7, 0.2, 1.0, 0.0);
        rModelPart.CreateNewNode(8, 2.2, 1.0, 0.0);
    }
    const std::vector<ModelPart::IndexType> first_ids = Quadrilateral
        ? std::vector<ModelPart::IndexType>{1, 2, 3, 7}
        : std::vector<ModelPart::IndexType>{1, 2, 3};
    const std::vector<ModelPart::IndexType> second_ids = Quadrilateral
        ? std::vector<ModelPart::IndexType>{4, 5, 8, 6}
        : std::vector<ModelPart::IndexType>{4, 5, 6};
    const std::string element_name = Quadrilateral ? "Element2D4N" : "Element2D3N";
    auto p_first = rModelPart.CreateNewElement(element_name, 1, first_ids, p_material);
    auto p_second = rModelPart.CreateNewElement(element_name, 2, second_ids, p_material);
    for (auto& r_node : rModelPart.Nodes()) {
        r_node.AddDof(VELOCITY_X).SetEquationId(3 * (r_node.Id() - 1));
        r_node.AddDof(VELOCITY_Y).SetEquationId(3 * (r_node.Id() - 1) + 1);
        r_node.AddDof(PRESSURE).SetEquationId(3 * (r_node.Id() - 1) + 2);
        const bool second = r_node.Id() > 3;
        auto& r_velocity = r_node.FastGetSolutionStepValue(VELOCITY);
        r_velocity[0] = 0.2 + 0.3 * r_node.X() + 0.1 * second;
        r_velocity[1] = -0.4 + 0.2 * r_node.Y() - 0.15 * second;
        r_node.FastGetSolutionStepValue(PRESSURE) = 1.0 + 0.5 * r_node.X() - 0.2 * r_node.Y();
    }
    auto condition_ids = first_ids;
    condition_ids.insert(condition_ids.end(), second_ids.begin(), second_ids.end());
    auto p_condition = rModelPart.CreateNewCondition(
        Quadrilateral ? "FluidPatchCouplingCondition2D8N" : "FluidPatchCouplingCondition2D6N",
        1, condition_ids, p_interface_properties);
    GlobalPointersVector<Element> parents;
    parents.resize(2);
    parents(0) = p_first;
    parents(1) = p_second;
    p_condition->SetValue(NEIGHBOUR_ELEMENTS, parents);
    Matrix coordinates = ZeroMatrix(2, 3);
    coordinates(0, 0) = coordinates(1, 0) = 1.0;
    coordinates(1, 1) = 1.0;
    p_condition->SetValue(COUPLING_INTERFACE_COORDINATES, coordinates);
    return p_condition;
}
} // unnamed namespace

KRATOS_TEST_CASE_IN_SUITE(FluidPatchCouplingQuadrilateral, FluidDynamicsApplicationFastSuite)
{
    for (const double coefficient : {40.0, 0.0}) {
        Model model;
        auto& r_model_part = model.CreateModelPart("Fluid");
        auto p_condition = create_condition(r_model_part, coefficient, true);
        const auto& r_info = r_model_part.GetProcessInfo();
        KRATOS_EXPECT_EQ(p_condition->Check(r_info), 0);
        Matrix lhs;
        Vector rhs;
        p_condition->CalculateLocalSystem(lhs, rhs, r_info);
        KRATOS_EXPECT_EQ(lhs.size1(), 24);
        KRATOS_EXPECT_EQ(rhs.size(), 24);
        const Matrix symmetry_error = coefficient > 0.0
            ? Matrix(lhs - trans(lhs)) : Matrix(lhs + trans(lhs));
        KRATOS_EXPECT_NEAR(norm_frobenius(symmetry_error), 0.0, 1e-12);
        Condition::EquationIdVectorType ids;
        Condition::DofsVectorType dofs;
        p_condition->EquationIdVector(ids, r_info);
        p_condition->GetDofList(dofs, r_info);
        KRATOS_EXPECT_EQ(ids.size(), 24);
        const double epsilon = 1e-6;
        for (std::size_t j = 0; j < 24; ++j) {
            KRATOS_EXPECT_EQ(ids[j], dofs[j]->EquationId());
            auto& r_node = p_condition->GetGeometry()[j / 3];
            double& r_value = j % 3 == 2
                ? r_node.FastGetSolutionStepValue(PRESSURE)
                : r_node.FastGetSolutionStepValue(VELOCITY)[j % 3];
            r_value += epsilon;
            Matrix perturbed_lhs;
            Vector perturbed_rhs;
            p_condition->CalculateLocalSystem(perturbed_lhs, perturbed_rhs, r_info);
            r_value -= epsilon;
            for (std::size_t i = 0; i < 24; ++i) {
                KRATOS_EXPECT_NEAR((perturbed_rhs[i] - rhs[i]) / epsilon, -lhs(i, j), 1e-7);
            }
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(FluidPatchCouplingPenaltyAndPenaltyFree, FluidDynamicsApplicationFastSuite)
{
    Model model;
    auto& r_model_part = model.CreateModelPart("Fluid");
    auto p_condition = create_condition(r_model_part, 40.0);
    const auto& r_info = r_model_part.GetProcessInfo();
    KRATOS_EXPECT_EQ(p_condition->Check(r_info), 0);
    Matrix penalty_lhs;
    Vector rhs;
    p_condition->CalculateLocalSystem(penalty_lhs, rhs, r_info);
    KRATOS_EXPECT_EQ(penalty_lhs.size1(), 18);
    KRATOS_EXPECT_EQ(rhs.size(), 18);
    Matrix symmetric_error = penalty_lhs - trans(penalty_lhs);
    KRATOS_EXPECT_NEAR(norm_frobenius(symmetric_error), 0.0, 1e-12);
    // Pressure is part of traction, including the pressure test equations.
    KRATOS_EXPECT_GT(std::abs(penalty_lhs(3, 5)), 1e-8);
    KRATOS_EXPECT_GT(std::abs(penalty_lhs(5, 3)), 1e-8);
    for (std::size_t i = 2; i < 18; i += 3) {
        for (std::size_t j = 2; j < 18; j += 3) {
            KRATOS_EXPECT_NEAR(penalty_lhs(i, j), 0.0, 1e-14);
        }
    }
    p_condition->GetProperties()[PENALTY_COEFFICIENT] = 0.0;
    Matrix free_lhs;
    p_condition->CalculateLocalSystem(free_lhs, rhs, r_info);
    Matrix skew_error = free_lhs + trans(free_lhs);
    KRATOS_EXPECT_NEAR(norm_frobenius(skew_error), 0.0, 1e-12);
    Matrix difference = penalty_lhs - free_lhs;
    KRATOS_EXPECT_GT(norm_frobenius(difference), 1.0);
    p_condition->GetProperties()[PENALTY_COEFFICIENT] = -1.0;
    Matrix iga_free_lhs;
    p_condition->CalculateLocalSystem(iga_free_lhs, rhs, r_info);
    KRATOS_EXPECT_MATRIX_NEAR(iga_free_lhs, free_lhs, 1e-14);

    Condition::EquationIdVectorType ids;
    Condition::DofsVectorType dofs;
    p_condition->EquationIdVector(ids, r_info);
    p_condition->GetDofList(dofs, r_info);
    KRATOS_EXPECT_EQ(ids.size(), 18);
    for (std::size_t i = 0; i < ids.size(); ++i) {
        KRATOS_EXPECT_EQ(ids[i], i);
        KRATOS_EXPECT_EQ(ids[i], dofs[i]->EquationId());
    }
}

KRATOS_TEST_CASE_IN_SUITE(FluidPatchCouplingResidualLinearization, FluidDynamicsApplicationFastSuite)
{
    for (const double coefficient : {40.0, 0.0}) {
        Model model;
        auto& r_model_part = model.CreateModelPart("Fluid");
        auto p_condition = create_condition(r_model_part, coefficient);
        Matrix lhs;
        Vector rhs;
        p_condition->CalculateLocalSystem(lhs, rhs, r_model_part.GetProcessInfo());
        const double epsilon = 1e-7;
        for (const std::size_t column : {std::size_t(0), std::size_t(5), std::size_t(13), std::size_t(17)}) {
            auto& r_node = p_condition->GetGeometry()[column / 3];
            double& r_value = column % 3 == 2
                ? r_node.FastGetSolutionStepValue(PRESSURE)
                : r_node.FastGetSolutionStepValue(VELOCITY)[column % 3];
            r_value += epsilon;
            Vector perturbed_rhs;
            p_condition->CalculateRightHandSide(perturbed_rhs, r_model_part.GetProcessInfo());
            r_value -= epsilon;
            for (std::size_t row = 0; row < 18; ++row) {
                KRATOS_EXPECT_NEAR(-(perturbed_rhs[row] - rhs[row]) / epsilon, lhs(row, column), 2e-6);
            }
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(FluidPatchCouplingLinearExactFlux, FluidDynamicsApplicationFastSuite)
{
    for (const double coefficient : {40.0, 0.0}) {
        Model model;
        auto& r_model_part = model.CreateModelPart("Fluid");
        auto p_condition = create_condition(r_model_part, coefficient);
        for (auto& r_node : r_model_part.Nodes()) {
            auto& r_velocity = r_node.FastGetSolutionStepValue(VELOCITY);
            r_velocity[0] = 1.0 + 2.0 * r_node.X() + 3.0 * r_node.Y();
            r_velocity[1] = -2.0 + 4.0 * r_node.X() - 2.0 * r_node.Y();
            r_node.FastGetSolutionStepValue(PRESSURE) = 1.0 + 2.0 * r_node.X() - 3.0 * r_node.Y();
        }
        Vector rhs;
        p_condition->CalculateRightHandSide(rhs, r_model_part.GetProcessInfo());
        const double mu = 1.3;
        Vector expected = ZeroVector(18);
        expected[3] = 2.0 * mu - 1.0;
        expected[6] = 2.0 * mu - 0.5;
        expected[4] = expected[7] = 3.5 * mu;
        expected[9] = -expected[3];
        expected[15] = -expected[6];
        expected[10] = expected[16] = -3.5 * mu;
        // The interface residual is the balancing traction, not zero by
        // itself. It cancels each patch's volume residual at exact fields.
        KRATOS_EXPECT_VECTOR_NEAR(rhs, expected, 1e-12);
    }
}

KRATOS_TEST_CASE_IN_SUITE(FluidPatchCouplingTimeAssemblyAndChecks, FluidDynamicsApplicationFastSuite)
{
    Model model;
    auto& r_model_part = model.CreateModelPart("Fluid");
    auto p_condition = create_condition(r_model_part, 0.0);
    auto& r_info = r_model_part.GetProcessInfo();
    Matrix bdf_lhs;
    Vector bdf_rhs;
    p_condition->CalculateLocalSystem(bdf_lhs, bdf_rhs, r_info);
    Matrix damping;
    Vector damping_rhs = bdf_rhs;
    p_condition->CalculateLocalVelocityContribution(damping, damping_rhs, r_info);
    KRATOS_EXPECT_NEAR(norm_frobenius(damping), 0.0, 1e-14);
    KRATOS_EXPECT_VECTOR_NEAR(damping_rhs, bdf_rhs, 1e-14);
    r_info.Erase(BDF_COEFFICIENTS);
    Matrix empty_lhs;
    Vector empty_rhs;
    p_condition->CalculateLocalSystem(empty_lhs, empty_rhs, r_info);
    KRATOS_EXPECT_NEAR(norm_frobenius(empty_lhs), 0.0, 1e-14);
    p_condition->CalculateLocalVelocityContribution(damping, damping_rhs, r_info);
    KRATOS_EXPECT_MATRIX_NEAR(damping, bdf_lhs, 1e-14);
    KRATOS_EXPECT_VECTOR_NEAR(damping_rhs, bdf_rhs, 1e-14);
    p_condition->GetProperties()[PENALTY_COEFFICIENT] = -2.0;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(p_condition->Check(r_info), "PENALTY_COEFFICIENT");
    p_condition->GetProperties()[PENALTY_COEFFICIENT] = 0.0;
    auto& r_coordinates = p_condition->GetValue(COUPLING_INTERFACE_COORDINATES);
    r_coordinates(0, 0) = r_coordinates(1, 0) = 1.1;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(p_condition->Check(r_info), "Coupling endpoints");
}

} // namespace Kratos::Testing
