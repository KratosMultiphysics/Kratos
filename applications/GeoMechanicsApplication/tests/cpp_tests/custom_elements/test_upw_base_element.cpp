// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Richard Faasse
//

#include "custom_constitutive/incremental_linear_elastic_law.h"
#include "custom_constitutive/plane_strain.h"
#include "custom_elements/U_Pw_base_element.hpp"
#include "custom_utilities/registration_utilities.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "tests/cpp_tests/test_utilities/element_setup_utilities.h"

using namespace std::string_literals;

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(UPwBaseElement_SerializesConstitutiveLaws, KratosGeoMechanicsFastSuite)
{
    ScopedSerializerRegistration registration("PlaneStrain", PlaneStrain{});
    const auto variables = Geo::ConstVariableRefs{std::cref(DISPLACEMENT_X), std::cref(DISPLACEMENT_Y),
                                                  std::cref(DISPLACEMENT_Z), std::cref(WATER_PRESSURE)};
    auto p_element = ElementSetupUtilities::Create2D3NElement(variables);

    auto p_constitutive_law =
        std::make_shared<GeoIncrementalLinearElasticLaw>(std::make_unique<PlaneStrain>());
    p_element->GetProperties().SetValue(CONSTITUTIVE_LAW, p_constitutive_law);

    std::vector<ConstitutiveLaw::Pointer> claws;
    p_element->CalculateOnIntegrationPoints(CONSTITUTIVE_LAW, claws, ProcessInfo{});
    EXPECT_TRUE(claws.empty());

    p_element->Initialize(ProcessInfo{});
    p_element->CalculateOnIntegrationPoints(CONSTITUTIVE_LAW, claws, ProcessInfo{});
    EXPECT_EQ(claws.size(), 3);

    auto serializer = StreamSerializer{};
    serializer.save("test_tag"s, p_element);
    auto p_loaded_element = std::make_shared<UPwBaseElement>();
    serializer.load("test_tag"s, p_loaded_element);

    std::vector<ConstitutiveLaw::Pointer> loaded_claws;
    p_loaded_element->CalculateOnIntegrationPoints(CONSTITUTIVE_LAW, loaded_claws, ProcessInfo{});

    EXPECT_EQ(loaded_claws.size(), 3);
    EXPECT_EQ(loaded_claws[0]->WorkingSpaceDimension(), 2);
}

} // namespace Kratos::Testing
