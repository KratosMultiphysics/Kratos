// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                     license: structural_mechanics_application/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//                   Alejandro Cornejo
//

// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "utilities/math_utils.h"
#include "includes/node.h"
#include "containers/model.h"
#include "geometries/triangle_2d_3.h"
#include "custom_utilities/advanced_constitutive_law_utilities.h"
#include "custom_utilities/constitutive_law_utilities.h"

namespace Kratos
{
namespace Testing
{

    // Tolerance
    static constexpr double tolerance = 1.0e-6;

    // NodeType
    typedef Node NodeType;

/**
* Check the correct calculation of the CL utilities principal stress
*/
KRATOS_TEST_CASE_IN_SUITE(ConstitutiveLawUtilitiesPrincipalStresses, KratosConstitutiveLawsFastSuite)
{
    array_1d<double, 6> stress_vector;
    array_1d<double, 3> principal_stress_vector;

    stress_vector[0] = -19.0e6;
    stress_vector[1] = 4.6e6;
    stress_vector[2] = -8.3e6;
    stress_vector[3] = -4.7e6;
    stress_vector[4] = 11.8e6;
    stress_vector[5] = 6.45e6;

    double I1, I2, I3;
    AdvancedConstitutiveLawUtilities<6>::CalculateI1Invariant(stress_vector, I1);
    AdvancedConstitutiveLawUtilities<6>::CalculateI2Invariant(stress_vector, I2);
    AdvancedConstitutiveLawUtilities<6>::CalculateI3Invariant(stress_vector, I3);

    KRATOS_EXPECT_LE((   -2.27e+07 - I1)/I1, tolerance);
    KRATOS_EXPECT_LE((-1.70812e+14 - I2)/I2, tolerance);
    KRATOS_EXPECT_LE(( 2.64752e+21 - I3)/I3, tolerance);

    AdvancedConstitutiveLawUtilities<6>::CalculatePrincipalStresses(principal_stress_vector, stress_vector);

    KRATOS_EXPECT_LE((1.16178e+07 - principal_stress_vector[0])/principal_stress_vector[0], tolerance);
    KRATOS_EXPECT_LE((-2.53163e+07 - principal_stress_vector[1])/principal_stress_vector[1], tolerance);
    KRATOS_EXPECT_LE((-9.00151e+06 - principal_stress_vector[2])/principal_stress_vector[2], tolerance);

    AdvancedConstitutiveLawUtilities<6>::CalculatePrincipalStressesWithCardano(principal_stress_vector, stress_vector);
    KRATOS_EXPECT_LE((1.16178e+07 - principal_stress_vector[0])/principal_stress_vector[0], tolerance);
    KRATOS_EXPECT_LE((-9.00151e+06 - principal_stress_vector[1])/principal_stress_vector[1], tolerance);
    KRATOS_EXPECT_LE((-2.53163e+07 - principal_stress_vector[2])/principal_stress_vector[2], tolerance);
}

/**
* Check the correct calculation of the CL utilities Hencky and Biot strains
*/
KRATOS_TEST_CASE_IN_SUITE(ConstitutiveLawUtilitiesHenckyAndBiot, KratosConstitutiveLawsFastSuite)
{
    // Declaration of Cauchy tensor and strain vector
    Matrix C = ZeroMatrix(3, 3);
    Vector strain_vector(6);

    // Fill Cauchy tensor
    C(0,0) = 1.1;
    C(1,1) = 1.0;
    C(2,2) = 1.0;
    C(1,0) = -1.0e-1;
    C(0,1) = -1.0e-1;

    // Compute Hencky
    AdvancedConstitutiveLawUtilities<6>::CalculateHenckyStrain(C, strain_vector);

    KRATOS_EXPECT_LE((0.0454445 - strain_vector[0])/strain_vector[0], tolerance);
    KRATOS_EXPECT_LE((-0.0023557 - strain_vector[1])/strain_vector[1], tolerance);
    KRATOS_EXPECT_LE((0.0 - strain_vector[2])/(strain_vector[2] + 1.0e-12), tolerance);
    KRATOS_EXPECT_LE((-0.0956005 - strain_vector[3])/strain_vector[3], tolerance);

    // Compute Biot
    AdvancedConstitutiveLawUtilities<6>::CalculateBiotStrain(C, strain_vector);

    KRATOS_EXPECT_LE((0.0476699 - strain_vector[0])/strain_vector[0], tolerance);
    KRATOS_EXPECT_LE((-0.00119458 - strain_vector[1])/strain_vector[1], tolerance);
    KRATOS_EXPECT_LE((0.0 - strain_vector[2])/(strain_vector[2] + 1.0e-12), tolerance);
    KRATOS_EXPECT_LE((-0.097729 - strain_vector[3])/strain_vector[3], tolerance);
}

/**
* Check the correct calculation of the CL utilities Hencky and Biot strains
*/
KRATOS_TEST_CASE_IN_SUITE(ConstitutiveLawUtilitiesAlmansiAndGreen, KratosConstitutiveLawsFastSuite)
{
    // Declaration of Cauchy tensor and strain vector
    Matrix C = ZeroMatrix(3, 3);
    Vector strain_vector(6);

    // Fill Cauchy tensor
    C(0,0) = 1.1;
    C(1,1) = 1.0;
    C(2,2) = 1.0;
    C(1,0) = -1.0e-1;
    C(0,1) = -1.0e-1;

    // Compute Green-Lagrange
    ConstitutiveLawUtilities<6>::CalculateGreenLagrangianStrain(C, strain_vector);

    KRATOS_EXPECT_LE((0.05 - strain_vector[0])/strain_vector[0], tolerance);
    KRATOS_EXPECT_LE((0.0 - strain_vector[1])/(strain_vector[1] + 1.0e-12), tolerance);
    KRATOS_EXPECT_LE((0.0 - strain_vector[2])/(strain_vector[2] + 1.0e-12), tolerance);
    KRATOS_EXPECT_LE((-0.001 - strain_vector[3])/strain_vector[3], tolerance);

    // Compute Alamansi (We use C, but it is supposed to be B, we are just using to compute and check)
    AdvancedConstitutiveLawUtilities<6>::CalculateAlmansiStrain(C, strain_vector);

    KRATOS_EXPECT_LE((0.0412844 - strain_vector[0])/strain_vector[0], tolerance);
    KRATOS_EXPECT_LE((-0.00458716 - strain_vector[1])/strain_vector[1], tolerance);
    KRATOS_EXPECT_LE((0.0 - strain_vector[2])/(strain_vector[2] + 1.0e-12), tolerance);
    KRATOS_EXPECT_LE((-0.0458716 - strain_vector[3])/strain_vector[3], tolerance);
}

/**
* Check the correct calculation of the CL utilities Hencky and Biot strains
*/
KRATOS_TEST_CASE_IN_SUITE(ConstitutiveLawUtilitiesPolarDecomposition, KratosConstitutiveLawsFastSuite)
{
    // Declaration of Cauchy tensor and strain vector
    Matrix F = ZeroMatrix(3, 3);
    Matrix R(3,3);
    Matrix U(3,3);

    // Fill Cauchy tensor
    F(0,0) = 1.1;
    F(1,1) = 1.0;
    F(2,2) = 1.0;
    F(1,0) = -1.0e-1;
    F(0,1) = -1.0e-1;

    // Compute Green-Lagrange
    ConstitutiveLawUtilities<6>::PolarDecomposition(F, R, U);

    // Reference solutions
    Matrix R_ref = ZeroMatrix(3, 3);
    Matrix U_ref = ZeroMatrix(3, 3);

    R_ref(0, 0) = 0.850651;
    R_ref(0, 1) = 0.525731;
    R_ref(1, 0) = -0.525731;
    R_ref(1, 1) = 0.850651;
    R_ref(2, 2) = 1.0;
    U_ref(0, 0) = 0.988289;
    U_ref(0, 1) = -0.610796;
    U_ref(1, 0) = 0.493239;
    U_ref(1, 1) = 0.798078;
    U_ref(2, 2) = 1.0;

    KRATOS_EXPECT_LE(norm_frobenius(R-R_ref), tolerance);
    KRATOS_EXPECT_LE(norm_frobenius(U-U_ref), tolerance);
}
/**
* Check the correct calculation of the CL utilities Hencky and Biot strains
*/
KRATOS_TEST_CASE_IN_SUITE(CalculateCharacteristicLength, KratosConstitutiveLawsFastSuite)
{
    // Model part
    Model current_model;
    ModelPart& this_model_part = current_model.CreateModelPart("Main");

    // First we create the nodes
    NodeType::Pointer p_node_1 = this_model_part.CreateNewNode(1, 0.0 , 0.0 , 0.0);
    NodeType::Pointer p_node_2 = this_model_part.CreateNewNode(2, 1.0 , 0.0 , 0.0);
    NodeType::Pointer p_node_3 = this_model_part.CreateNewNode(3, 1.0 , 1.0 , 0.0);

    // Now we create the "conditions"
    std::vector<NodeType::Pointer> element_nodes (3);
    element_nodes[0] = p_node_1;
    element_nodes[1] = p_node_2;
    element_nodes[2] = p_node_3;
    Triangle2D3 <NodeType> triangle( PointerVector<NodeType>{element_nodes} );

    array_1d<double, 3> delta_disp = ZeroVector(3);
    delta_disp[0] = 1.0e-2;
    p_node_1->Coordinates() += delta_disp;

    // Compute CalculateCharacteristicLength
    const double length = AdvancedConstitutiveLawUtilities<3>::CalculateCharacteristicLength(triangle);

    KRATOS_EXPECT_NEAR(length,  0.743871, tolerance);

    // Compute CalculateCharacteristicLengthOnReferenceConfiguration
    const double length_origin = AdvancedConstitutiveLawUtilities<3>::CalculateCharacteristicLengthOnReferenceConfiguration(triangle);

    KRATOS_EXPECT_NEAR(length_origin,  0.745356, tolerance);
}
/**
* Check the correct calculation of the shear modulus
*/
KRATOS_TEST_CASE_IN_SUITE(CalculateSearModulus, KratosConstitutiveLawsFastSuite)
{
    auto props = Properties(1);
    props.SetValue(YOUNG_MODULUS, 2.0E9);
    props.SetValue(POISSON_RATIO, 0.2);

    const double G = ConstitutiveLawUtilities<3>::CalculateShearModulus(props);

    KRATOS_EXPECT_NEAR(G,  2.0e9 / 2.4, tolerance);
}
} // namespace Testing
} // namespace Kratos
