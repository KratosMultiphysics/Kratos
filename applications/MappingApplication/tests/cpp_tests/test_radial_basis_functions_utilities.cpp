//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main author:    Juan Ignacio Camarotti
//

// System includes

// Project includes
#include "testing/testing.h"
#include "custom_utilities/radial_basis_functions_utilities.h"

namespace Kratos::Testing {

using namespace Kratos::RadialBasisFunctionsUtilities;

KRATOS_TEST_CASE_IN_SUITE(RadialBasisFunctionsUtilities_Gaussian, KratosMappingApplicationSerialTestSuite)
{
    Gaussian gaussian_rbf;
    gaussian_rbf.h = 2.0;

    // r = 0  -> exp(0) = 1
    const double value_at_0 = gaussian_rbf(0.0);
    KRATOS_EXPECT_NEAR(value_at_0, 1.0, 1e-12);

    // r = h -> exp(-0.5)
    const double value_at_h = gaussian_rbf(gaussian_rbf.h);
    KRATOS_EXPECT_NEAR(value_at_h, std::exp(-0.5), 1e-12);

    // symmetry in r
    const double value_at_minus_h = gaussian_rbf(-gaussian_rbf.h);
    KRATOS_EXPECT_NEAR(value_at_h, value_at_minus_h, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RadialBasisFunctionsUtilities_InverseMultiquadric, KratosMappingApplicationSerialTestSuite)
{
    InverseMultiquadric inverse_multiquadric_rbf;
    inverse_multiquadric_rbf.h = 1.5;

    const double r = 0.7;
    const double q = inverse_multiquadric_rbf.h * r;
    const double expected = 1.0 / std::sqrt(1.0 + q * q);

    const double value = inverse_multiquadric_rbf(r);
    KRATOS_EXPECT_NEAR(value, expected, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RadialBasisFunctionsUtilities_Multiquadric, KratosMappingApplicationSerialTestSuite)
{
    Multiquadric multiquadric_rbf;
    multiquadric_rbf.h = 0.5;

    const double r = 0.3;
    const double q = r / multiquadric_rbf.h;
    const double expected = std::sqrt(1.0 + q * q);

    const double value = multiquadric_rbf(r);
    KRATOS_EXPECT_NEAR(value, expected, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RadialBasisFunctionsUtilities_ThinPlateSpline, KratosMappingApplicationSerialTestSuite)
{
    ThinPlateSpline tps_rbf;

    // By definition, r = 0 must return 0
    const double value_at_zero = tps_rbf(0.0);
    KRATOS_EXPECT_NEAR(value_at_zero, 0.0, 1e-16);

    const double r = 2.0;
    const double r2 = r * r;
    const double expected = r2 * std::log(r2);

    const double value = tps_rbf(r);
    KRATOS_EXPECT_NEAR(value, expected, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RadialBasisFunctionsUtilities_WendlandC2, KratosMappingApplicationSerialTestSuite)
{
    WendlandC2 wendland_c2_rbf;
    wendland_c2_rbf.h = 1.0;

    // Outside support -> 0
    const double value_outside = wendland_c2_rbf(1.1);
    KRATOS_EXPECT_NEAR(value_outside, 0.0, 1e-16);

    // Inside support: r = 0.5, h = 1.0 => q = 0.5
    const double r = 0.5;
    const double q = r / wendland_c2_rbf.h;
    const double expected = std::pow(1.0 - q, 4) * (4.0 * q + 1.0);

    const double value_inside = wendland_c2_rbf(r);
    KRATOS_EXPECT_NEAR(value_inside, expected, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RadialBasisFunctionsUtilities_WendlandC2SupportRadius, KratosMappingApplicationSerialTestSuite)
{
    Matrix points(3, 3);
    points(0, 0) = 0.0; points(0, 1) = 0.0; points(0, 2) = 0.0;
    points(1, 0) = 1.0; points(1, 1) = 0.0; points(1, 2) = 0.0;
    points(2, 0) = 2.0; points(2, 1) = 0.0; points(2, 2) = 0.0;

    const double k = 2.5;

    const double support_radius =
        CalculateWendlandC2SupportRadius(points, k);

    // This expected value assumes the implementation uses:
    //   average_spacing = total_pairwise_distance / number_of_pairs
    //   support_radius  = k * average_spacing
    const double total_distance = 1.0 + 2.0 + 1.0; // 0-1, 0-2, 1-2
    const double count_pairs = 3.0;
    const double average_spacing = total_distance / count_pairs;
    const double expected_radius = k * average_spacing;

    KRATOS_EXPECT_NEAR(support_radius, expected_radius, 1e-12);
}

} // namespace Kratos::Testing