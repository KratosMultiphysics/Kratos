// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Anne van de Graaf
//
#include "custom_constitutive/plane_strain.h"
#include "custom_constitutive/three_dimensional.h"
#include "custom_utilities/registration_utilities.hpp"
#include "includes/constitutive_law.h"
#include "includes/stream_serializer.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"

#include <string>

using namespace std::string_literals;

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(ConstitutiveLawDimensionPlaneStrain_CanBeSavedAndLoadedThroughInterface,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto scoped_registration =
        ScopedSerializerRegistration{std::make_pair("PlaneStrain"s, PlaneStrain{})};
    const auto p_dimension = std::unique_ptr<ConstitutiveLawDimension>{std::make_unique<PlaneStrain>()};
    auto serializer = StreamSerializer{};

    // Act
    serializer.save("test_tag"s, p_dimension);
    auto p_loaded_dimension = std::unique_ptr<ConstitutiveLawDimension>{};
    serializer.load("test_tag"s, p_loaded_dimension);

    // Assert
    ASSERT_NE(p_loaded_dimension, nullptr);
    KRATOS_EXPECT_EQ(p_loaded_dimension->GetSpatialType(), ConstitutiveLaw::PLANE_STRAIN_LAW);
}

KRATOS_TEST_CASE_IN_SUITE(ConstitutiveLawDimensionThreeDimensional_CanBeSavedAndLoadedThroughInterface,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto scoped_registration =
        ScopedSerializerRegistration{std::make_pair("ThreeDimensional"s, ThreeDimensional{})};
    const auto p_dimension =
        std::unique_ptr<ConstitutiveLawDimension>{std::make_unique<ThreeDimensional>()};
    auto serializer = StreamSerializer{};

    // Act
    serializer.save("test_tag"s, p_dimension);
    auto p_loaded_dimension = std::unique_ptr<ConstitutiveLawDimension>{};
    serializer.load("test_tag"s, p_loaded_dimension);

    // Assert
    ASSERT_NE(p_loaded_dimension, nullptr);
    KRATOS_EXPECT_EQ(p_loaded_dimension->GetSpatialType(), ConstitutiveLaw::THREE_DIMENSIONAL_LAW);
}

} // namespace Kratos::Testing