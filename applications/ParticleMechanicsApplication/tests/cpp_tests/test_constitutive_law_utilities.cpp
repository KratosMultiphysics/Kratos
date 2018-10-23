//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    Bodhinanda Chandra
//


// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "utilities/math_utils.h"

namespace Kratos
{
namespace Testing
{

    // Tolerance
    static constexpr double tolerance = 1.0e-6;

/**
* Check the correct calculation of the CL utilities principal stress
*/
KRATOS_TEST_CASE_IN_SUITE(dummytestfortesting, KratosParticleMechanicsFastSuite)
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

    KRATOS_CHECK_LESS_EQUAL((   -2.27e+07 - I1)/I1, tolerance);
    KRATOS_CHECK_LESS_EQUAL((-1.70812e+14 - I2)/I2, tolerance);
    KRATOS_CHECK_LESS_EQUAL(( 2.64752e+21 - I3)/I3, tolerance);

    KRATOS_CHECK_LESS_EQUAL((1.16178e+07 - principal_stress_vector[0])/principal_stress_vector[0], tolerance);
    KRATOS_CHECK_LESS_EQUAL((-2.53163e+07 - principal_stress_vector[1])/principal_stress_vector[1], tolerance);
    KRATOS_CHECK_LESS_EQUAL((-9.00151e+06 - principal_stress_vector[2])/principal_stress_vector[2], tolerance);

    KRATOS_CHECK_LESS_EQUAL((1.16178e+07 - principal_stress_vector[0])/principal_stress_vector[0], tolerance);
    KRATOS_CHECK_LESS_EQUAL((-9.00151e+06 - principal_stress_vector[1])/principal_stress_vector[1], tolerance);
    KRATOS_CHECK_LESS_EQUAL((-2.53163e+07 - principal_stress_vector[2])/principal_stress_vector[2], tolerance);
}
} // namespace Testing
} // namespace Kratos
