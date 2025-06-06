//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Andrea Gorgi
//

// Project includes
#include "containers/model.h"
#include "testing/testing.h"
#include "test_creation_utility.h"
#include "custom_conditions/support_solid_iga_condition.h"
#include "iga_structural_mechanics_fast_suite.h"

namespace Kratos::Testing
{
namespace
{
    typedef std::size_t SizeType;

    typename Condition::Pointer GetSupportSolidIGACondition(
        ModelPart& rModelPart, 
        SizeType PolynomialDegree, 
        IntegrationPoint<3> IntegrationPoint)
    {
        // Set the condition properties
        auto p_cond_prop = rModelPart.CreateNewProperties(0);
        p_cond_prop->SetValue(PENALTY_FACTOR, -1.0); // Penalty-free formulation
        p_cond_prop->SetValue(THICKNESS, 1.0);
        p_cond_prop->SetValue(YOUNG_MODULUS, 100);
        p_cond_prop->SetValue(POISSON_RATIO, 0.3);

        // Costitutive law
        const auto& r_clone_cl = KratosComponents<ConstitutiveLaw>::Get("LinearElasticPlaneStrain2DLaw");
        p_cond_prop->SetValue(CONSTITUTIVE_LAW, r_clone_cl.Clone());

        auto p_quadrature_point = TestCreationUtility::GetQuadraturePointGeometryOnCurve(
            rModelPart, PolynomialDegree, IntegrationPoint);

        return Kratos::make_intrusive<SupportSolidIGACondition>(1, p_quadrature_point, p_cond_prop);
    }
}

// Test per SupportSolidIGACondition con p=3
KRATOS_TEST_CASE_IN_SUITE(SupportSolidIGAConditionP3, KratosIgaSMFastSuite)
{
    Model model;
    auto& r_model_part = model.CreateModelPart("ModelPart");

    // Set buffer size
    r_model_part.SetBufferSize(2);

    // Aggiunta variabili
    r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    r_model_part.GetProcessInfo().SetValue(DOMAIN_SIZE, 2);

    const auto& r_process_info = r_model_part.GetProcessInfo();

    IntegrationPoint<3> integration_point(0.333333333333333, 0.05, 0.0, 0.086963711284364);
    auto p_support_condition = GetSupportSolidIGACondition(r_model_part, 3, integration_point);

    Vector u_D = ZeroVector(3);
    u_D[0] = 0.1; // Displacement in x-direction
    u_D[1] = -0.5; // Displacement in y-direction
    p_support_condition->SetValue(DISPLACEMENT, u_D);

    // Aggiunta DOFs
    for (auto& r_node : r_model_part.Nodes()) {
        r_node.AddDof(DISPLACEMENT_X);
        r_node.AddDof(DISPLACEMENT_Y);
    }

    p_support_condition->Initialize(r_process_info);
    p_support_condition->InitializeSolutionStep(r_process_info);
    p_support_condition->FinalizeSolutionStep(r_process_info);

    Matrix left_hand_side_matrix;
    Vector right_hand_side_vector;
    p_support_condition->CalculateLocalSystem(left_hand_side_matrix, right_hand_side_vector, r_process_info);

    // Check LHS e RHS
    const double tolerance = 1.0e-4;

    const std::vector<std::vector<double>> expected_LHS = { 
       {0, 0.596275, 0, 2.68324, 0, 2.23603, 0, 0.521741, 2.93641, 0.0313829, 4.40462, 0.141223, 2.20231, 0.117686, 0.367052, 0.0274601}, 
        {-0.596275, 0, 1.78883, 0, 2.23603, 0, 0.596275, 0, -0.0313829, 10.2774, 0.0941487, 15.4162, 0.117686, 7.70808, 0.0313829, 1.28468}, 
        {0, -1.78883, 0, -1.34055e-15, 0, 1.34162, 1.11022e-16, 0.447207, 4.40462, -0.0941487, 6.60693, -2.57848e-17, 3.30346, 0.0706116, 0.550577, 0.0235372}, 
        {-2.68324, 0, 1.34055e-15, 0, 2.01243, 1.77636e-15, 0.67081, 0, -0.141223, 15.4162, 1.15325e-16, 23.1243, 0.105917, 11.5621, 0.0353058, 1.92702}, 
        {0, -2.23603, 0, -2.01243, 0, -0.335405, 5.55112e-17, 0.0559008, 2.20231, -0.117686, 3.30346, -0.105917, 1.65173, -0.0176529, 0.275289, 0.00294215}, 
        {-2.23603, 0, -1.34162, -1.77636e-15, 0.335405, 0, 0.223603, 0, -0.117686, 7.70808, -0.0706116, 11.5621, 0.0176529, 5.78106, 0.0117686, 0.963511}, 
        {0, -0.596275, -1.11022e-16, -0.67081, -5.55112e-17, -0.223603, 0, -0.0186336, 0.367052, -0.0313829, 0.550577, -0.0353058, 0.275289, -0.0117686, 0.0458815, -0.000980716}, 
        {-0.521741, 0, -0.447207, 0, -0.0559008, 0, 0.0186336, 0, -0.0274601, 1.28468, -0.0235372, 1.92702, -0.00294215, 0.963511, 0.000980716, 0.160585}, 
        {-2.93641, 0.0313829, -4.40462, 0.141223, -2.20231, 0.117686, -0.367052, 0.0274601, 0, 0.00165173, 0, 0.0074328, 0, 0.006194, 0, 0.00144527}, 
        {-0.0313829, -10.2774, 0.0941487, -15.4162, 0.117686, -7.70808, 0.0313829, -1.28468, -0.00165173, 0, 0.0049552, 0, 0.006194, 0, 0.00165173, -1.38778e-17}, 
        {-4.40462, -0.0941487, -6.60693, -1.15325e-16, -3.30346, 0.0706116, -0.550577, 0.0235372, 0, -0.0049552, 0, -3.71343e-18, 0, 0.0037164, 0, 0.0012388}, 
        {-0.141223, -15.4162, 2.57848e-17, -23.1243, 0.105917, -11.5621, 0.0353058, -1.92702, -0.0074328, 0, 3.71343e-18, 0, 0.0055746, 1.11022e-16, 0.0018582, 0}, 
        {-2.20231, -0.117686, -3.30346, -0.105917, -1.65173, -0.0176529, -0.275289, 0.00294215, 0, -0.006194, 0, -0.0055746, 0, -0.000929099, 0, 0.00015485}, 
        {-0.117686, -7.70808, -0.0706116, -11.5621, 0.0176529, -5.78106, 0.0117686, -0.963511, -0.006194, 0, -0.0037164, -1.11022e-16, 0.000929099, 0, 0.0006194, -6.93889e-18}, 
        {-0.367052, -0.0313829, -0.550577, -0.0353058, -0.275289, -0.0117686, -0.0458815, -0.000980716, 0, -0.00165173, 0, -0.0018582, 0, -0.0006194, 0, -5.16166e-05}, 
        {-0.0274601, -1.28468, -0.0235372, -1.92702, -0.00294215, -0.963511, 0.000980716, -0.160585, -0.00144527, 1.38778e-17, -0.0012388, 0, -0.00015485, 6.93889e-18, 5.16166e-05, 0}
    };
    const std::vector<double> expected_RHS = {
        -2.18648,
        -16.9195,
        1.48656,
        -26.0148,
        3.12642,
        -13.3251,
        0.91826,
        -2.27382,
        -1.15828,
        17.3655,
        -1.48656,
        26.0148,
        -0.617851,
        12.9907,
        -0.0820705,
        2.16232
    };

    // std::cout << "LHS matrix:" << std::endl;
    // for (SizeType i = 0; i < left_hand_side_matrix.size1(); ++i) {
    //     std::cout << "{";
    //     for (SizeType j = 0; j < left_hand_side_matrix.size2(); ++j) {
    //         std::cout << left_hand_side_matrix(i,j);
    //         if (j < left_hand_side_matrix.size2() - 1) {
    //             std::cout << ", ";
    //         }
    //     }
    //     std::cout << "}";
    //     if (i < left_hand_side_matrix.size1() - 1) {
    //         std::cout << ", ";
    //     }
    //     std::cout << std::endl;
    // }

    // std::cout << "RHS vector:" << std::endl;
    // for (SizeType i = 0; i < right_hand_side_vector.size(); ++i) {
    //     std::cout << right_hand_side_vector(i) << "," << std::endl;
    // }
    // for (SizeType i = 0; i < left_hand_side_matrix.size1(); ++i) {
    //     for (SizeType j = 0; j < left_hand_side_matrix.size2(); ++j) {
    //         KRATOS_EXPECT_NEAR(left_hand_side_matrix(i,j), expected_LHS[i][j], tolerance);
    //     }
    // }

    for (SizeType i = 0; i < right_hand_side_vector.size(); ++i) {
        KRATOS_EXPECT_NEAR(right_hand_side_vector(i), expected_RHS[i], tolerance);
    }
}
}
