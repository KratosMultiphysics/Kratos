// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Wijtze Pieter Kikstra
//                   Anne van de Graaf
//

#include "custom_constitutive/mohr_coulomb_constitutive_law.hpp"
#include "custom_constitutive/plane_strain.h"

#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "tests/cpp_tests/test_utilities.h"

#include <boost/numeric/ublas/assignment.hpp>
#include <sstream>
#include <string>

using namespace Kratos;
using namespace std::string_literals;

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(MohrCoulombConstitutiveLaw_Clone, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto original_law = MohrCoulombConstitutiveLaw(std::make_unique<PlaneStrain>());
    auto       p_cloned_law = original_law.Clone();

    KRATOS_EXPECT_NE(p_cloned_law.get(), nullptr);
    KRATOS_EXPECT_NE(p_cloned_law.get(), &original_law);
    KRATOS_EXPECT_NE(dynamic_cast<const MohrCoulombConstitutiveLaw*>(p_cloned_law.get()), nullptr);
}

KRATOS_TEST_CASE_IN_SUITE(MohrCoulombConstitutiveLaw_Check, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto                        law = MohrCoulombConstitutiveLaw(std::make_unique<PlaneStrain>());
    ConstitutiveLaw::Parameters parameters;
    Properties                  properties(3);
    parameters.SetMaterialProperties(properties);
    const auto element_geometry = Geometry<Node>{};
    const auto process_info     = ProcessInfo{};
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        law.Check(properties, element_geometry, process_info),
        "Error: GEO_COHESION is not defined or has an invalid value for property: 3")
    properties.SetValue(GEO_COHESION, 1.0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        law.Check(properties, element_geometry, process_info),
        "Error: GEO_FRICTION_ANGLE is not defined or has an invalid value for property: 3")
    properties.SetValue(GEO_FRICTION_ANGLE, 1.0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        law.Check(properties, element_geometry, process_info),
        "Error: GEO_DILATION_ANGLE is not defined or has an invalid value for property: 3")
    properties.SetValue(GEO_DILATION_ANGLE, 1.0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        law.Check(properties, element_geometry, process_info),
        "Error: GEO_TENSION_CUTOFF is not defined or has an invalid value for property: 3")
    properties.SetValue(GEO_TENSION_CUTOFF, 1.0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        law.Check(properties, element_geometry, process_info),
        "Error: YOUNG_MODULUS has Key zero, is not defined or has an invalid value for property: 3")
    properties.SetValue(YOUNG_MODULUS, 1.0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        law.Check(properties, element_geometry, process_info),
        "Error: POISSON_RATIO is not defined or has an invalid value for property: 3")
    properties.SetValue(POISSON_RATIO, 1.0);
    KRATOS_EXPECT_EQ(law.Check(properties, element_geometry, process_info), 0);
}

KRATOS_TEST_CASE_IN_SUITE(MohrCoulombConstitutiveLaw_SetValueAndGetValue, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto       law          = MohrCoulombConstitutiveLaw(std::make_unique<PlaneStrain>());
    const auto process_info = ProcessInfo{};

    constexpr auto set_value = 1.0;
    law.SetValue(STATE_VARIABLE, set_value, process_info);
    auto zero_value = 0.0;
    KRATOS_EXPECT_DOUBLE_EQ(law.GetValue(STATE_VARIABLE, zero_value), set_value);
    zero_value = 0.0;
    KRATOS_EXPECT_DOUBLE_EQ(law.GetValue(DAMAGE_VARIABLE, zero_value), set_value);
}

} // namespace Kratos::Testing
