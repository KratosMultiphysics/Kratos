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

#include "custom_retention/saturated_below_phreatic_level_law.h"
#include "custom_retention/saturated_law.h"
#include "custom_retention/van_genuchten_law.h"
#include "custom_utilities/registration_utilities.hpp"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"

#include <string>

using namespace std::string_literals;

namespace Kratos::Testing
{

template <typename T>
class TestRetentionLawSerialization : public ::testing::Test
{
};

using RetentionLawTypes = ::testing::Types<SaturatedBelowPhreaticLevelLaw, SaturatedLaw, VanGenuchtenLaw>;
TYPED_TEST_SUITE(TestRetentionLawSerialization, RetentionLawTypes);

TYPED_TEST(TestRetentionLawSerialization, AnyRetentionLawCanBeSavedAndLoaded)
{
    // Arrange
    const auto scoped_registration =
        ScopedSerializerRegistration{std::make_pair(TypeParam{}.Info(), TypeParam{})};

    const auto p_retention_law = std::unique_ptr<RetentionLaw>{std::make_unique<TypeParam>()};
    auto       serializer      = StreamSerializer{};

    // Act
    serializer.save("test_tag"s, p_retention_law);
    auto p_loaded_retention_law = std::unique_ptr<RetentionLaw>{};
    serializer.load("test_tag"s, p_loaded_retention_law);

    // Assert
    KRATOS_EXPECT_NE(p_loaded_retention_law, nullptr);
    KRATOS_EXPECT_NE(dynamic_cast<TypeParam*>(p_loaded_retention_law.get()), nullptr);
}

} // namespace Kratos::Testing
