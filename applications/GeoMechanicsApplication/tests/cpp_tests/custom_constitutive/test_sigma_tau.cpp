// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//
//  Main authors:    Anne van de Graaf
//

#include "custom_constitutive/geo_sigma_tau.hpp"
#include "custom_utilities/ublas_utilities.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "tests/cpp_tests/test_utilities.h"

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(SigmaTau_HasZeroesAsValuesWhenDefaultConstructed, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    KRATOS_EXPECT_VECTOR_NEAR(Geo::SigmaTau().Values(), Vector(2, 0.0), Defaults::absolute_tolerance);
    KRATOS_EXPECT_VECTOR_NEAR(Geo::SigmaTau{}.Values(), Vector(2, 0.0), Defaults::absolute_tolerance);
    KRATOS_EXPECT_NEAR(Geo::SigmaTau{}.Sigma(), 0.0, Defaults::absolute_tolerance);
    KRATOS_EXPECT_NEAR(Geo::SigmaTau{}.Tau(), 0.0, Defaults::absolute_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(SigmaTau_CanBeConstructedFromAnyVectorWithSize2, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto vector_1 = UblasUtilities::CreateVector({1.0, 2.0});

    // Act
    const auto sigma_tau_1 = Geo::SigmaTau{vector_1};

    // Assert
    KRATOS_EXPECT_VECTOR_NEAR(sigma_tau_1.Values(), vector_1, Defaults::absolute_tolerance);
    KRATOS_EXPECT_NEAR(sigma_tau_1.Sigma(), vector_1[0], Defaults::absolute_tolerance);
    KRATOS_EXPECT_NEAR(sigma_tau_1.Tau(), vector_1[1], Defaults::absolute_tolerance);

    // Arrange
    auto vector_2 = BoundedVector<double, 2>{};
    vector_2[0]   = 3.0;
    vector_2[1]   = 4.0;

    // Act
    const auto sigma_tau_2 = Geo::SigmaTau{vector_2};

    // Assert
    KRATOS_EXPECT_VECTOR_NEAR(sigma_tau_2.Values(), vector_2, Defaults::absolute_tolerance);

    // Arrange
    const auto vector_3 = std::vector{5.0, 6.0};

    // Act & Assert
    KRATOS_EXPECT_VECTOR_NEAR(Geo::SigmaTau{vector_3}.Values(), vector_3, Defaults::absolute_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(SigmaTau_RaisesADebugErrorWhenAttemptingToConstructFromAVectorWithSizeUnequalTo2,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
#ifndef KRATOS_DEBUG
    GTEST_SKIP() << "This test requires a debug build";
#endif

    // Arrange
    const auto too_short = UblasUtilities::CreateVector({1.0});
    const auto too_long  = UblasUtilities::CreateVector({2.0, 3.0, 4.0});

    // Act & Assert
    EXPECT_THROW(Geo::SigmaTau{too_short}, Exception);
    EXPECT_THROW(Geo::SigmaTau{too_long}, Exception);
}

} // namespace Kratos::Testing