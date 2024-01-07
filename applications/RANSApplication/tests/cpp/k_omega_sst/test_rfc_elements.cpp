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
#include "../stabilization_method_test_utilities.h"
#include "custom_utilities/fluid_test_utilities.h"
#include "custom_utilities/test_utilities.h"
#include "rans_application_variables.h"
#include "test_utilities.h"

namespace Kratos
{
namespace Testing
{
namespace
{
ModelPart& RansKOmegaSSTKRFC2D3NSetUp(
    Model& rModel)
{
    auto& r_model_part = KOmegaSSTTestUtilities::RansKOmegaSSTK2D3NSetUp(
        rModel, "RansKOmegaSSTKRFC2D3N");

    StabilizationMethodTestUtilities::InitializeResidualBasedFluxCorrectedConstants(
        r_model_part.GetProcessInfo());

    RansApplicationTestUtilities::CheckElementsAndConditions(r_model_part);

    return r_model_part;
}

ModelPart& RansKOmegaSSTOmegaRFC2D3NSetUp(
    Model& rModel)
{
    auto& r_model_part = KOmegaSSTTestUtilities::RansKOmegaSSTOmega2D3NSetUp(
        rModel, "RansKOmegaSSTOmegaRFC2D3N");

    StabilizationMethodTestUtilities::InitializeResidualBasedFluxCorrectedConstants(
        r_model_part.GetProcessInfo());

    RansApplicationTestUtilities::CheckElementsAndConditions(r_model_part);

    return r_model_part;
}

} // namespace

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTKRFC2D3N_EquationIdVector, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTKRFC2D3NSetUp(model);

    // Test:
    FluidTestUtilities::RunEntityEquationIdVectorTest(r_model_part.Elements(), r_model_part.GetProcessInfo(), {&TURBULENT_KINETIC_ENERGY});
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTKRFC2D3N_GetDofList, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTKRFC2D3NSetUp(model);

    // Test:
    FluidTestUtilities::RunEntityGetDofListTest(r_model_part.Elements(), r_model_part.GetProcessInfo(), {&TURBULENT_KINETIC_ENERGY});
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTKRFC2D3N_CalculateLocalSystem, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTKRFC2D3NSetUp(model);

    // Test:
    Matrix LHS, ref_LHS(3, 3);
    Vector RHS, ref_RHS(3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateLocalSystem(
        LHS, RHS, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_RHS[0] =  3.3127898780761695e+00;
    ref_RHS[1] =  5.4533614255485272e+00;
    ref_RHS[2] =  5.6468930884664319e+00;
    ref_LHS = ZeroMatrix(3, 3);

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
    KRATOS_CHECK_MATRIX_NEAR(LHS, ref_LHS, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTKRFC2D3N_CalculateRightHandSide, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTKRFC2D3NSetUp(model);

    // Test:
    Vector RHS, ref_RHS(3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateRightHandSide(
        RHS, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_RHS[0] =  3.3127898780761695e+00;
    ref_RHS[1] =  5.4533614255485272e+00;
    ref_RHS[2] =  5.6468930884664319e+00;

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTKRFC2D3N_CalculateLocalVelocityContribution, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTKRFC2D3NSetUp(model);

    // Test:
    Matrix LHS, ref_LHS(3, 3);
    Vector RHS(3, 0.0), ref_RHS(3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateLocalVelocityContribution(
        LHS, RHS, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_LHS(0, 0) =  3.5469914790686704e+02;
    ref_LHS(0, 1) =  -4.6416943795724485e+01;
    ref_LHS(0, 2) =  -6.5047976071625257e+01;
    ref_LHS(1, 0) =  -5.2413138434578741e+01;
    ref_LHS(1, 1) =  3.0888287548296159e+02;
    ref_LHS(1, 2) =  -6.0931710230200380e+01;
    ref_LHS(2, 0) =  -7.6264727571941393e+01;
    ref_LHS(2, 1) =  -6.8038071892262820e+01;
    ref_LHS(2, 2) =  4.1449416929414042e+02;

    ref_RHS[0] =  1.3545714397014908e+03;
    ref_RHS[1] =  -6.9032705019065743e+03;
    ref_RHS[2] =  -3.5314530765540272e+04;

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
    KRATOS_CHECK_MATRIX_NEAR(LHS, ref_LHS, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTKRFC2D3N_CalculateMassMatrix, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTKRFC2D3NSetUp(model);

    // Test:
    Matrix M, ref_M(3, 3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateMassMatrix(
        M, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_M(0, 0) =  2.4116177005911799e-01;
    ref_M(0, 1) =  3.0349287845057079e-02;
    ref_M(0, 2) =  3.3870458373186078e-02;
    ref_M(1, 0) =  4.0732504067512894e-02;
    ref_M(1, 1) =  2.4642375414775408e-01;
    ref_M(1, 2) =  3.9719726410431998e-02;
    ref_M(2, 0) =  4.8919772188969493e-02;
    ref_M(2, 1) =  5.1830636847301229e-02;
    ref_M(2, 2) =  2.5725633028429495e-01;

    KRATOS_CHECK_MATRIX_NEAR(M, ref_M, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTKRFC2D3N_CalculateDampingMatrix, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTKRFC2D3NSetUp(model);

    // Test:
    Matrix D, ref_D(3, 3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateDampingMatrix(
        D, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_D(0, 0) =  3.5469914790686704e+02;
    ref_D(0, 1) =  -4.6416943795724485e+01;
    ref_D(0, 2) =  -6.5047976071625257e+01;
    ref_D(1, 0) =  -5.2413138434578741e+01;
    ref_D(1, 1) =  3.0888287548296159e+02;
    ref_D(1, 2) =  -6.0931710230200380e+01;
    ref_D(2, 0) =  -7.6264727571941393e+01;
    ref_D(2, 1) =  -6.8038071892262820e+01;
    ref_D(2, 2) =  4.1449416929414042e+02;

    KRATOS_CHECK_MATRIX_NEAR(D, ref_D, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTOmegaRFC2D3N_EquationIdVector, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTOmegaRFC2D3NSetUp(model);

    // Test:
    FluidTestUtilities::RunEntityEquationIdVectorTest(r_model_part.Elements(), r_model_part.GetProcessInfo(), {&TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE});
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTOmegaRFC2D3N_GetDofList, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTOmegaRFC2D3NSetUp(model);

    // Test:
    FluidTestUtilities::RunEntityGetDofListTest(r_model_part.Elements(), r_model_part.GetProcessInfo(), {&TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE});
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTOmegaRFC2D3N_CalculateLocalSystem, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTOmegaRFC2D3NSetUp(model);

    // Test:
    Matrix LHS, ref_LHS(3, 3);
    Vector RHS, ref_RHS(3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateLocalSystem(
        LHS, RHS, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_RHS[0] =  -2.1518154207903831e+03;
    ref_RHS[1] =  -2.2656057944918284e+03;
    ref_RHS[2] =  -2.4311133152853463e+03;

    ref_LHS = ZeroMatrix(3, 3);

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
    KRATOS_CHECK_MATRIX_NEAR(LHS, ref_LHS, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTOmegaRFC2D3N_CalculateRightHandSide, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTOmegaRFC2D3NSetUp(model);

    // Test:
    Vector RHS, ref_RHS(3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateRightHandSide(
        RHS, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_RHS[0] =  -2.1518154207903831e+03;
    ref_RHS[1] =  -2.2656057944918284e+03;
    ref_RHS[2] =  -2.4311133152853463e+03;

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTOmegaRFC2D3N_CalculateLocalVelocityContribution, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTOmegaRFC2D3NSetUp(model);

    // Test:
    Matrix LHS, ref_LHS(3, 3);
    Vector RHS(3, 0.0), ref_RHS(3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateLocalVelocityContribution(
        LHS, RHS, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_LHS(0, 0) =  5.0887437955601195e+02;
    ref_LHS(0, 1) =  -5.9664338962811826e+01;
    ref_LHS(0, 2) =  -8.4002202266447227e+01;
    ref_LHS(1, 0) =  -6.5646377645873713e+01;
    ref_LHS(1, 1) =  3.9177796790362027e+02;
    ref_LHS(1, 2) =  -7.1398449566227896e+01;
    ref_LHS(2, 0) =  -9.5066786607357812e+01;
    ref_LHS(2, 1) =  -7.8314313346957292e+01;
    ref_LHS(2, 2) =  5.5731472652746379e+02;

    ref_RHS[0] =  -1.9779834282070445e+05;
    ref_RHS[1] =  1.6865756847683777e+04;
    ref_RHS[2] =  -1.8426111845645477e+05;

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
    KRATOS_CHECK_MATRIX_NEAR(LHS, ref_LHS, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTOmegaRFC2D3N_CalculateMassMatrix, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTOmegaRFC2D3NSetUp(model);

    // Test:
    Matrix M, ref_M(3, 3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateMassMatrix(
        M, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_M(0, 0) =  2.4384327981160181e-01;
    ref_M(0, 1) =  3.2285095491244052e-02;
    ref_M(0, 2) =  3.6061244277775073e-02;
    ref_M(1, 0) =  4.1041568593648679e-02;
    ref_M(1, 1) =  2.4730258457031049e-01;
    ref_M(1, 2) =  4.0355452571691000e-02;
    ref_M(2, 0) =  4.6992265969397852e-02;
    ref_M(2, 1) =  5.0258623197852181e-02;
    ref_M(2, 2) =  2.5546016764251733e-01;

    KRATOS_CHECK_MATRIX_NEAR(M, ref_M, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTOmegaRFC2D3N_CalculateDampingMatrix, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTOmegaRFC2D3NSetUp(model);

    // Test:
    Matrix D, ref_D(3, 3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateDampingMatrix(
        D, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_D(0, 0) =  5.0887437955601195e+02;
    ref_D(0, 1) =  -5.9664338962811826e+01;
    ref_D(0, 2) =  -8.4002202266447227e+01;
    ref_D(1, 0) =  -6.5646377645873713e+01;
    ref_D(1, 1) =  3.9177796790362027e+02;
    ref_D(1, 2) =  -7.1398449566227896e+01;
    ref_D(2, 0) =  -9.5066786607357812e+01;
    ref_D(2, 1) =  -7.8314313346957292e+01;
    ref_D(2, 2) =  5.5731472652746379e+02;

    KRATOS_CHECK_MATRIX_NEAR(D, ref_D, 1e-12);
}

} // namespace Testing
} // namespace Kratos.
