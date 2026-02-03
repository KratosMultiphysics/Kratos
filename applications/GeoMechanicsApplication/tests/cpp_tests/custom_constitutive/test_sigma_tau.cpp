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
#include "includes/expect.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite_without_kernel.h"
#include "tests/cpp_tests/test_utilities.h"

namespace Kratos::Testing
{

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, SigmaTau_HasZeroesAsValuesWhenDefaultConstructed)
{
    KRATOS_EXPECT_VECTOR_NEAR(Geo::SigmaTau().Values(), Vector(2, 0.0), Defaults::absolute_tolerance);
    KRATOS_EXPECT_VECTOR_NEAR(Geo::SigmaTau().Values(), Vector(2, 0.0), Defaults::absolute_tolerance);
    EXPECT_NEAR(Geo::SigmaTau().Sigma(), 0.0, Defaults::absolute_tolerance);
    EXPECT_NEAR(Geo::SigmaTau().Tau(), 0.0, Defaults::absolute_tolerance);
}

template <typename T>
class TestSigmaTauFixture : public ::testing::Test
{
};

using TestVectorTypes = ::testing::Types<Vector, BoundedVector<double, 2>, std::vector<double>>;
TYPED_TEST_SUITE(TestSigmaTauFixture, TestVectorTypes);

TYPED_TEST(TestSigmaTauFixture, SigmaTau_CanBeConstructedFromAnyVectorWithSizeOf2)
{
    // Arrange
    auto initialization_vector = TypeParam(2);
    initialization_vector[0]   = 1.0;
    initialization_vector[1]   = 2.0;

    // Act
    const auto sigma_tau = Geo::SigmaTau{initialization_vector};

    // Assert
    KRATOS_EXPECT_VECTOR_NEAR(sigma_tau.Values(), initialization_vector, Defaults::absolute_tolerance);
    EXPECT_NEAR(sigma_tau.Sigma(), initialization_vector[0], Defaults::absolute_tolerance);
    EXPECT_NEAR(sigma_tau.Tau(), initialization_vector[1], Defaults::absolute_tolerance);
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel,
       SigmaTau_RaisesADebugErrorWhenAttemptingToConstructFromAVectorWithSizeUnequalTo2)
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

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, SigmaTau_CanBeConstructedFromAStdInitializerListWithSize2)
{
    EXPECT_NEAR((Geo::SigmaTau{1.0, 2.0}.Sigma()), 1.0, Defaults::absolute_tolerance);
    EXPECT_NEAR((Geo::SigmaTau{1.0, 2.0}.Tau()), 2.0, Defaults::absolute_tolerance);
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel,
       SigmaTau_RaisesADebugErrorWhenAttemptingToConstructFromANonEmptyStdInitializerListWithSizeUnequalTo2)
{
#ifndef KRATOS_DEBUG
    GTEST_SKIP() << "This test requires a debug build";
#endif

    EXPECT_NO_THROW(Geo::SigmaTau{}); // empty list is OK
    EXPECT_THROW(Geo::SigmaTau{1.0}, Exception);
    EXPECT_THROW((Geo::SigmaTau{2.0, 3.0, 4.0}), Exception);
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, SigmaTau_ComponentsCanBeModifiedDirectly)
{
    // Arrange
    auto sigma_tau = Geo::SigmaTau{1.0, 2.0};

    // Act
    sigma_tau.Sigma() = 3.0;

    // Assert
    KRATOS_EXPECT_VECTOR_NEAR(sigma_tau.Values(), UblasUtilities::CreateVector({3.0, 2.0}),
                              Defaults::absolute_tolerance);

    // Act
    sigma_tau.Tau() = 4.0;

    // Assert
    KRATOS_EXPECT_VECTOR_NEAR(sigma_tau.Values(), UblasUtilities::CreateVector({3.0, 4.0}),
                              Defaults::absolute_tolerance);
}

TYPED_TEST(TestSigmaTauFixture, SigmaTau_CanBeCopiedToAnyVectorTypeWithSizeOf2)
{
    // Arrange
    const auto sigma_tau = Geo::SigmaTau{1.0, 2.0};

    // Act & Assert
    KRATOS_EXPECT_VECTOR_NEAR(sigma_tau.CopyTo<TypeParam>(), (std::vector{1.0, 2.0}), Defaults::absolute_tolerance);
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, SigmaTau_RaisesADebugErrorWhenAttemptingToCopyToAVectorThatIsTooSmall)
{
#ifndef KRATOS_DEBUG
    GTEST_SKIP() << "This test requires a debug build";
#endif

    // Arrange
    const auto sigma_tau = Geo::SigmaTau{1.0, 2.0};

    // Act & Assert
    EXPECT_THROW((sigma_tau.CopyTo<BoundedVector<double, 1>>()), boost::numeric::ublas::bad_size);
    EXPECT_NO_THROW((sigma_tau.CopyTo<BoundedVector<double, 3>>())); // excess element(s) are OK
}

} // namespace Kratos::Testing