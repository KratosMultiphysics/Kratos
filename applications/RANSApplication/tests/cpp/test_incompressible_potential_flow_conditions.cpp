//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "testing/testing.h"

// Application includes
#include "custom_utilities/test_utilities.h"
#include "includes/cfd_variables.h"
#include "rans_application_variables.h"

namespace Kratos
{
namespace Testing
{
namespace
{
ModelPart& RansIncompressiblePotentialFlowVelocityInlet2D2N_SetUp(
    Model& rModel)
{
    const auto add_variables_function = [](ModelPart& rModelPart) {
        rModelPart.AddNodalSolutionStepVariable(VELOCITY_POTENTIAL);
        rModelPart.AddNodalSolutionStepVariable(VELOCITY);
    };

    using namespace RansApplicationTestUtilities;

    ModelPart& r_model_part = CreateScalarVariableTestModelPart(
        rModel, "Element2D3N",
        "RansIncompressiblePotentialFlowVelocityInlet2D2N",
        add_variables_function, VELOCITY_POTENTIAL, 1);

    // set nodal historical variables
    RandomFillNodalHistoricalVariable(r_model_part, VELOCITY_POTENTIAL, -10.0, 10.0);
    RandomFillNodalHistoricalVariable(r_model_part, VELOCITY, -5.0, 5.0);

    RandomFillContainerVariable<ModelPart::ConditionsContainerType, array_1d<double, 3>>(
        r_model_part, NORMAL, 1.0, 10.0);

    return r_model_part;
}

ModelPart& RansIncompressiblePotentialFlowPressureBodyForce2D2N_SetUp(
    Model& rModel)
{
    const auto add_variables_function = [](ModelPart& rModelPart) {
        rModelPart.AddNodalSolutionStepVariable(VELOCITY_POTENTIAL);
        rModelPart.AddNodalSolutionStepVariable(PRESSURE_POTENTIAL);
        rModelPart.AddNodalSolutionStepVariable(DENSITY);
        rModelPart.AddNodalSolutionStepVariable(BODY_FORCE);
    };

    using namespace RansApplicationTestUtilities;

    ModelPart& r_model_part = CreateScalarVariableTestModelPart(
        rModel, "Element2D3N",
        "RansIncompressiblePotentialFlowPressureBodyForce2D2N",
        add_variables_function, PRESSURE_POTENTIAL, 1, true, false);

    // set nodal historical variables
    RandomFillNodalHistoricalVariable(r_model_part, DENSITY, 2.0, 2.0);
    RandomFillNodalHistoricalVariable(r_model_part, BODY_FORCE, -2.0, 2.0);
    RandomFillNodalHistoricalVariable(r_model_part, VELOCITY_POTENTIAL, -10.0, 10.0);
    RandomFillNodalHistoricalVariable(r_model_part, PRESSURE_POTENTIAL, -10.0, 10.0);

    RandomFillContainerVariable<ModelPart::ConditionsContainerType, array_1d<double, 3>>(
        r_model_part, NORMAL, 1.0, 10.0);

    return r_model_part;
}

} // namespace

KRATOS_TEST_CASE_IN_SUITE(RansIncompressiblePotentialFlowVelocityInlet2D2N_EquationIdVector,
                          KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansIncompressiblePotentialFlowVelocityInlet2D2N_SetUp(model);

    // Test:
    RansApplicationTestUtilities::TestEquationIdVector<ModelPart::ConditionsContainerType>(r_model_part);
}

KRATOS_TEST_CASE_IN_SUITE(RansIncompressiblePotentialFlowVelocityInlet2D2N_GetDofList, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansIncompressiblePotentialFlowVelocityInlet2D2N_SetUp(model);

    // Test:
    RansApplicationTestUtilities::TestGetDofList<ModelPart::ConditionsContainerType>(
        r_model_part, VELOCITY_POTENTIAL);
}

KRATOS_TEST_CASE_IN_SUITE(RansIncompressiblePotentialFlowVelocityInlet2D2N_CalculateLocalSystem,
                          KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansIncompressiblePotentialFlowVelocityInlet2D2N_SetUp(model);

    // Test:
    Matrix LHS, ref_LHS;
    Vector RHS, ref_RHS;
    auto& r_condition = r_model_part.Conditions().front();

    // checking for non-inlet condition
    r_condition.SetValue(RANS_IS_INLET, 0);
    r_condition.Initialize();
    r_condition.CalculateLocalSystem(
        LHS, RHS, static_cast<const ProcessInfo>(r_model_part.GetProcessInfo()));
    // setting reference values
    ref_RHS = ZeroVector(2);
    ref_LHS = ZeroMatrix(2, 2);

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
    KRATOS_CHECK_MATRIX_NEAR(LHS, ref_LHS, 1e-12);

    // checking for inlet condition
    r_condition.SetValue(RANS_IS_INLET, 1);
    r_condition.Initialize();
    r_condition.CalculateLocalSystem(
        LHS, RHS, static_cast<const ProcessInfo>(r_model_part.GetProcessInfo()));
    // settting reference values
    ref_RHS[0] = -7.4912949377611748e-01;
    ref_RHS[1] = -7.4912949377611748e-01;

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
    KRATOS_CHECK_MATRIX_NEAR(LHS, ref_LHS, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansIncompressiblePotentialFlowVelocityInlet2D2N_CalculateLeftHandSide,
                          KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansIncompressiblePotentialFlowVelocityInlet2D2N_SetUp(model);

    // Test:
    Matrix LHS, ref_LHS(2, 2, 0.0);
    auto& r_condition = r_model_part.Conditions().front();

    r_condition.SetValue(RANS_IS_INLET, 0);
    r_condition.Initialize();
    r_condition.CalculateLeftHandSide(
        LHS, static_cast<const ProcessInfo>(r_model_part.GetProcessInfo()));

    KRATOS_CHECK_MATRIX_NEAR(LHS, ref_LHS, 1e-12);

    r_condition.SetValue(RANS_IS_INLET, 1);
    r_condition.Initialize();
    r_condition.CalculateLeftHandSide(
        LHS, static_cast<const ProcessInfo>(r_model_part.GetProcessInfo()));

    KRATOS_CHECK_MATRIX_NEAR(LHS, ref_LHS, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansIncompressiblePotentialFlowVelocityInlet2D2N_CalculateRightHandSide,
                          KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansIncompressiblePotentialFlowVelocityInlet2D2N_SetUp(model);

    // Test:
    Vector RHS, ref_RHS(2, 0.0);
    auto& r_condition = r_model_part.Conditions().front();

    r_condition.SetValue(RANS_IS_INLET, 0);
    r_condition.Initialize();
    r_condition.CalculateRightHandSide(
        RHS, static_cast<const ProcessInfo>(r_model_part.GetProcessInfo()));

    r_condition.SetValue(RANS_IS_INLET, 1);
    r_condition.Initialize();
    r_condition.CalculateRightHandSide(
        RHS, static_cast<const ProcessInfo>(r_model_part.GetProcessInfo()));
    ref_RHS[0] = -7.4912949377611748e-01;
    ref_RHS[1] = -7.4912949377611748e-01;

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansIncompressiblePotentialFlowPressureBodyForce2D2N_EquationIdVector,
                          KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansIncompressiblePotentialFlowPressureBodyForce2D2N_SetUp(model);

    // Test:
    RansApplicationTestUtilities::TestEquationIdVector<ModelPart::ConditionsContainerType>(r_model_part);
}

KRATOS_TEST_CASE_IN_SUITE(RansIncompressiblePotentialFlowPressureBodyForce2D2N_GetDofList,
                          KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansIncompressiblePotentialFlowPressureBodyForce2D2N_SetUp(model);

    // Test:
    RansApplicationTestUtilities::TestGetDofList<ModelPart::ConditionsContainerType>(
        r_model_part, PRESSURE_POTENTIAL);
}

KRATOS_TEST_CASE_IN_SUITE(RansIncompressiblePotentialFlowPressureBodyForce2D2N_CalculateLocalSystem,
                          KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansIncompressiblePotentialFlowPressureBodyForce2D2N_SetUp(model);

    // Test:
    Matrix LHS, ref_LHS(2, 2, 0.0);
    Vector RHS, ref_RHS(2);
    auto& r_condition = r_model_part.Conditions().front();
    r_condition.Initialize();
    r_condition.CalculateLocalSystem(
        LHS, RHS, static_cast<const ProcessInfo>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_RHS[0] = -9.9559464418787802e-01;
    ref_RHS[1] = -9.9559464418787802e-01;

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
    KRATOS_CHECK_MATRIX_NEAR(LHS, ref_LHS, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansIncompressiblePotentialFlowPressureBodyForce2D2N_CalculateLeftHandSide,
                          KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansIncompressiblePotentialFlowPressureBodyForce2D2N_SetUp(model);

    // Test:
    Matrix LHS, ref_LHS(2, 2, 0.0);
    auto& r_condition = r_model_part.Conditions().front();
    r_condition.Initialize();
    r_condition.CalculateLeftHandSide(
        LHS, static_cast<const ProcessInfo>(r_model_part.GetProcessInfo()));
    KRATOS_CHECK_MATRIX_NEAR(LHS, ref_LHS, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansIncompressiblePotentialFlowPressureBodyForce2D2N_CalculateRightHandSide,
                          KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansIncompressiblePotentialFlowPressureBodyForce2D2N_SetUp(model);

    // Test:
    Vector RHS, ref_RHS(2);
    auto& r_condition = r_model_part.Conditions().front();
    r_condition.Initialize();
    r_condition.CalculateRightHandSide(
        RHS, static_cast<const ProcessInfo>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_RHS[0] = -9.9559464418787802e-01;
    ref_RHS[1] = -9.9559464418787802e-01;

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
}
} // namespace Testing
} // namespace Kratos.
