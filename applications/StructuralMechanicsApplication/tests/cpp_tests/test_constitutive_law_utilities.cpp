// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                     license: structural_mechanics_application/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "utilities/math_utils.h"
#include "custom_utilities/constitutive_law_utilities.h"

namespace Kratos
{
namespace Testing
{

    // Tolerance
    static constexpr double tolerance = 1.0e-6;

/**
* Check the correct calculation of the CL utilities principal stress
*/
KRATOS_TEST_CASE_IN_SUITE(ConstitutiveLawUtilitiesPrincipalStresses, KratosStructuralMechanicsFastSuite)
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
    ConstitutiveLawUtilities<6>::CalculateI1Invariant(stress_vector, I1);
    ConstitutiveLawUtilities<6>::CalculateI2Invariant(stress_vector, I2);
    ConstitutiveLawUtilities<6>::CalculateI3Invariant(stress_vector, I3);

//     // DEBUG
//     std::cout << "I1: " << I1 << "\tI2: " << I2 << "\tI3: " << I3 << std::endl;

    KRATOS_CHECK_LESS_EQUAL((   -2.27e+07 - I1)/I1, tolerance);
    KRATOS_CHECK_LESS_EQUAL((-1.70812e+14 - I2)/I2, tolerance);
    KRATOS_CHECK_LESS_EQUAL(( 2.64752e+21 - I3)/I3, tolerance);

    ConstitutiveLawUtilities<6>::CalculatePrincipalStresses(principal_stress_vector, stress_vector);

//     // DEBUG
//     KRATOS_WATCH(principal_stress_vector)

    KRATOS_CHECK_LESS_EQUAL((1.16178e+07 - principal_stress_vector[0])/principal_stress_vector[0], tolerance);
    KRATOS_CHECK_LESS_EQUAL((-2.53163e+07 - principal_stress_vector[1])/principal_stress_vector[1], tolerance);
    KRATOS_CHECK_LESS_EQUAL((-9.00151e+06 - principal_stress_vector[2])/principal_stress_vector[2], tolerance);

    ConstitutiveLawUtilities<6>::CalculatePrincipalStressesWithCardano(principal_stress_vector, stress_vector);

//     // DEBUG
//     KRATOS_WATCH(principal_stress_vector)

    KRATOS_CHECK_LESS_EQUAL((1.16178e+07 - principal_stress_vector[0])/principal_stress_vector[0], tolerance);
    KRATOS_CHECK_LESS_EQUAL((-9.00151e+06 - principal_stress_vector[1])/principal_stress_vector[1], tolerance);
    KRATOS_CHECK_LESS_EQUAL((-2.53163e+07 - principal_stress_vector[2])/principal_stress_vector[2], tolerance);
}

/**
* Check the correct calculation of the CL utilities Hencky and Biot strains
*/
KRATOS_TEST_CASE_IN_SUITE(ConstitutiveLawUtilitiesHenckyAndBiot, KratosStructuralMechanicsFastSuite)
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
    ConstitutiveLawUtilities<6>::CalculateHenckyStrain(C, strain_vector);

//     // Debug
//     KRATOS_WATCH(strain_vector)

    KRATOS_CHECK_LESS_EQUAL((0.0454445 - strain_vector[0])/strain_vector[0], tolerance);
    KRATOS_CHECK_LESS_EQUAL((-0.0023557 - strain_vector[1])/strain_vector[1], tolerance);
    KRATOS_CHECK_LESS_EQUAL(strain_vector[2], tolerance);
    KRATOS_CHECK_LESS_EQUAL((-0.0956005 - strain_vector[3])/strain_vector[3], tolerance);

    // Compute Biot
    ConstitutiveLawUtilities<6>::CalculateBiotStrain(C, strain_vector);

//     // Debug
//     KRATOS_WATCH(strain_vector)

    KRATOS_CHECK_LESS_EQUAL((0.0476699 - strain_vector[0])/strain_vector[0], tolerance);
    KRATOS_CHECK_LESS_EQUAL((-0.00119458 - strain_vector[1])/strain_vector[1], tolerance);
    KRATOS_CHECK_LESS_EQUAL(strain_vector[2], tolerance);
    KRATOS_CHECK_LESS_EQUAL((-0.097729 - strain_vector[3])/strain_vector[3], tolerance);
}
} // namespace Testing
} // namespace Kratos
