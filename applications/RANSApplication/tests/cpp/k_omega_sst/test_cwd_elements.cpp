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
#include "custom_utilities/test_utilities.h"
#include "rans_application_variables.h"
#include "test_utilities.h"

namespace Kratos
{
namespace Testing
{
namespace
{
ModelPart& RansKOmegaSSTKCWD2D3NSetUp(
    Model& rModel)
{
    auto& r_model_part = KOmegaSSTTestUtilities::RansKOmegaSSTK2D3NSetUp(
        rModel, "RansKOmegaSSTKCWD2D3N");

    StabilizationMethodTestUtilities::InitializeCrossWindStabilizationConstants(
        r_model_part.GetProcessInfo());

    return r_model_part;
}

ModelPart& RansKOmegaSSTOmegaCWD2D3NSetUp(
    Model& rModel)
{
    auto& r_model_part = KOmegaSSTTestUtilities::RansKOmegaSSTOmega2D3NSetUp(
        rModel, "RansKOmegaSSTOmegaCWD2D3N");

    StabilizationMethodTestUtilities::InitializeCrossWindStabilizationConstants(
        r_model_part.GetProcessInfo());

    return r_model_part;
}

} // namespace

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTKCWD2D3N_EquationIdVector, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTKCWD2D3NSetUp(model);

    // Test:
    RansApplicationTestUtilities::TestEquationIdVector<ModelPart::ElementsContainerType>(r_model_part);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTKCWD2D3N_GetDofList, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTKCWD2D3NSetUp(model);

    // Test:
    RansApplicationTestUtilities::TestGetDofList<ModelPart::ElementsContainerType>(
        r_model_part, TURBULENT_KINETIC_ENERGY);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTKCWD2D3N_CalculateLocalSystem, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTKCWD2D3NSetUp(model);

    // Test:
    Matrix LHS, ref_LHS;
    Vector RHS, ref_RHS(3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateLocalSystem(
        LHS, RHS, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_RHS[0] = 2.4695761774232956e+00;
    ref_RHS[1] = 1.5697917789130236e+00;
    ref_RHS[2] = 1.6364014140218073e+00;
    ref_LHS = ZeroMatrix(3, 3);

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
    KRATOS_CHECK_MATRIX_NEAR(LHS, ref_LHS, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTKCWD2D3N_CalculateRightHandSide, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTKCWD2D3NSetUp(model);

    // Test:
    Vector RHS, ref_RHS(3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateRightHandSide(
        RHS, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_RHS[0] = 2.4695761774232956e+00;
    ref_RHS[1] = 1.5697917789130236e+00;
    ref_RHS[2] = 1.6364014140218073e+00;

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTKCWD2D3N_CalculateLocalVelocityContribution, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTKCWD2D3NSetUp(model);

    // Test:
    Matrix LHS, ref_LHS(3, 3);
    Vector RHS(3, 0.0), ref_RHS(3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateLocalVelocityContribution(
        LHS, RHS, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_LHS(0, 0) = 1.6511106466399615e+03;
    ref_LHS(0, 1) = -1.1991203013598561e+03;
    ref_LHS(0, 2) = 2.5025247807782401e+02;
    ref_LHS(1, 0) = -1.1980827139472685e+03;
    ref_LHS(1, 1) = 3.0794728713090144e+03;
    ref_LHS(1, 2) = -1.1061800119166003e+03;
    ref_LHS(2, 0) = 2.4945420690455288e+02;
    ref_LHS(2, 1) = -1.1068590134706019e+03;
    ref_LHS(2, 2) = 2.1086766565663147e+03;

    ref_RHS[0] = -7.6693427034775901e+04;
    ref_RHS[1] = 1.1565481448237575e+05;
    ref_RHS[2] = -1.5878011887264551e+05;

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
    KRATOS_CHECK_MATRIX_NEAR(LHS, ref_LHS, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTKCWD2D3N_CalculateMassMatrix, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTKCWD2D3NSetUp(model);

    // Test:
    Matrix M, ref_M(3, 3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateMassMatrix(
        M, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_M(0, 0) = 2.5029361543117840e-01;
    ref_M(0, 1) = 4.1823659205420256e-02;
    ref_M(0, 2) = 4.1665035737394232e-02;
    ref_M(1, 0) = 4.1127052781114497e-02;
    ref_M(1, 1) = 2.4972648473291831e-01;
    ref_M(1, 2) = 4.1479037223053218e-02;
    ref_M(2, 0) = 4.1911314404462302e-02;
    ref_M(2, 1) = 4.1782743358597274e-02;
    ref_M(2, 2) = 2.5018882921689373e-01;

    KRATOS_CHECK_MATRIX_NEAR(M, ref_M, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTKCWD2D3N_CalculateDampingMatrix, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTKCWD2D3NSetUp(model);

    // Test:
    Matrix D, ref_D(3, 3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateDampingMatrix(
        D, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_D(0, 0) = 1.6511106466399615e+03;
    ref_D(0, 1) = -1.1991203013598561e+03;
    ref_D(0, 2) = 2.5025247807782401e+02;
    ref_D(1, 0) = -1.1980827139472685e+03;
    ref_D(1, 1) = 3.0794728713090144e+03;
    ref_D(1, 2) = -1.1061800119166003e+03;
    ref_D(2, 0) = 2.4945420690455288e+02;
    ref_D(2, 1) = -1.1068590134706019e+03;
    ref_D(2, 2) = 2.1086766565663147e+03;

    KRATOS_CHECK_MATRIX_NEAR(D, ref_D, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTOmegaCWD2D3N_EquationIdVector, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTOmegaCWD2D3NSetUp(model);

    // Test:
    RansApplicationTestUtilities::TestEquationIdVector<ModelPart::ElementsContainerType>(r_model_part);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTOmegaCWD2D3N_GetDofList, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTOmegaCWD2D3NSetUp(model);

    // Test:
    RansApplicationTestUtilities::TestGetDofList<ModelPart::ElementsContainerType>(
        r_model_part, TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTOmegaCWD2D3N_CalculateLocalSystem, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTOmegaCWD2D3NSetUp(model);

    // Test:
    Matrix LHS, ref_LHS;
    Vector RHS, ref_RHS(3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateLocalSystem(
        LHS, RHS, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_RHS[0] = -4.4653425966448558e+03;
    ref_RHS[1] = -4.4594309479529456e+03;
    ref_RHS[2] = -4.5000424290893034e+03;
    ref_LHS = ZeroMatrix(3, 3);

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
    KRATOS_CHECK_MATRIX_NEAR(LHS, ref_LHS, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTOmegaCWD2D3N_CalculateRightHandSide, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTOmegaCWD2D3NSetUp(model);

    // Test:
    Vector RHS, ref_RHS(3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateRightHandSide(
        RHS, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_RHS[0] = -4.4653425966448558e+03;
    ref_RHS[1] = -4.4594309479529456e+03;
    ref_RHS[2] = -4.5000424290893034e+03;

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTOmegaCWD2D3N_CalculateLocalVelocityContribution, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTOmegaCWD2D3NSetUp(model);

    // Test:
    Matrix LHS, ref_LHS(3, 3);
    Vector RHS(3, 0.0), ref_RHS(3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateLocalVelocityContribution(
        LHS, RHS, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_LHS(0, 0) = 1.1562015794677250e+03;
    ref_LHS(0, 1) = -7.0283984627297491e+02;
    ref_LHS(0, 2) = 1.1057029705825336e+02;
    ref_LHS(1, 0) = -7.0180225886038761e+02;
    ref_LHS(1, 1) = 1.9307105690696326e+03;
    ref_LHS(1, 2) = -7.4095754535836738e+02;
    ref_LHS(2, 0) = 1.0977202588498216e+02;
    ref_LHS(2, 1) = -7.4163654691236877e+02;
    ref_LHS(2, 2) = 9.6115299128284153e+02;

    ref_RHS[0] = -5.2538593169971311e+05;
    ref_RHS[1] = -4.8165689936027513e+05;
    ref_RHS[2] = 2.2773746416511392e+05;

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
    KRATOS_CHECK_MATRIX_NEAR(LHS, ref_LHS, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTOmegaCWD2D3N_CalculateMassMatrix, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTOmegaCWD2D3NSetUp(model);

    // Test:
    Matrix M, ref_M(3, 3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateMassMatrix(
        M, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_M(0, 0) = 2.5004265770405271e-01;
    ref_M(0, 1) = 4.1649022815353244e-02;
    ref_M(0, 2) = 4.0923772294834349e-02;
    ref_M(1, 0) = 4.1197290676099667e-02;
    ref_M(1, 1) = 2.4967444879671088e-01;
    ref_M(1, 2) = 4.1303364091819304e-02;
    ref_M(2, 0) = 4.2091296941878742e-02;
    ref_M(2, 1) = 4.2008202335761254e-02;
    ref_M(2, 2) = 2.5110068448321560e-01;

    KRATOS_CHECK_MATRIX_NEAR(M, ref_M, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTOmegaCWD2D3N_CalculateDampingMatrix, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTOmegaCWD2D3NSetUp(model);

    // Test:
    Matrix D, ref_D(3, 3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateDampingMatrix(
        D, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_D(0, 0) = 1.1562015794677250e+03;
    ref_D(0, 1) = -7.0283984627297491e+02;
    ref_D(0, 2) = 1.1057029705825336e+02;
    ref_D(1, 0) = -7.0180225886038761e+02;
    ref_D(1, 1) = 1.9307105690696326e+03;
    ref_D(1, 2) = -7.4095754535836738e+02;
    ref_D(2, 0) = 1.0977202588498216e+02;
    ref_D(2, 1) = -7.4163654691236877e+02;
    ref_D(2, 2) = 9.6115299128284153e+02;

    KRATOS_CHECK_MATRIX_NEAR(D, ref_D, 1e-12);
}

} // namespace Testing
} // namespace Kratos.
