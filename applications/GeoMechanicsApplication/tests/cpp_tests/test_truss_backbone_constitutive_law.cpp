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

#include "custom_constitutive/truss_backbone_constitutive_law.h"
#include "geo_mechanics_application_variables.h"
#include "testing/testing.h"

using namespace Kratos;

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(CheckOfBackboneLawThrowsWhenPropertiesDoesNotHaveYoungsModulus, KratosGeoMechanicsFastSuite)
{
    const auto backbone_law = TrussBackboneConstitutiveLaw{};
    const auto properties   = Properties{};
    const auto geometry     = Geometry<Node>{};
    const auto process_info = ProcessInfo{};

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        [[maybe_unused]] const auto rv = backbone_law.Check(properties, geometry, process_info),
        "Error: No YOUNGS_MODULUS found")
}

KRATOS_TEST_CASE_IN_SUITE(CheckOfBackboneLawThrowsWhenPropertiesDoesNotHaveStrains, KratosGeoMechanicsFastSuite)
{
    const auto backbone_law = TrussBackboneConstitutiveLaw{};
    auto       properties   = Properties{};
    properties.SetValue(YOUNG_MODULUS, 100.0);
    const auto geometry     = Geometry<Node>{};
    const auto process_info = ProcessInfo{};

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        [[maybe_unused]] const auto rv = backbone_law.Check(properties, geometry, process_info),
        "Error: No STRAINS_OF_PIECEWISE_LINEAR_LAW found")
}

KRATOS_TEST_CASE_IN_SUITE(CheckOfBackboneLawThrowsWhenPropertiesDoesNotHaveStresses, KratosGeoMechanicsFastSuite)
{
    const auto backbone_law = TrussBackboneConstitutiveLaw{};
    auto       properties   = Properties{};
    properties.SetValue(YOUNG_MODULUS, 100.0);
    properties.SetValue(STRAINS_OF_PIECEWISE_LINEAR_LAW, Vector{});
    const auto geometry     = Geometry<Node>{};
    const auto process_info = ProcessInfo{};

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        [[maybe_unused]] const auto rv = backbone_law.Check(properties, geometry, process_info),
        "Error: No STRESSES_OF_PIECEWISE_LINEAR_LAW found")
}

KRATOS_TEST_CASE_IN_SUITE(CheckOfBackboneLawThrowsWhenStressesAndStrainsHaveDifferentSizes, KratosGeoMechanicsFastSuite)
{
    const auto backbone_law = TrussBackboneConstitutiveLaw{};
    auto       properties   = Properties{};
    properties.SetValue(YOUNG_MODULUS, 100.0);
    properties.SetValue(STRAINS_OF_PIECEWISE_LINEAR_LAW, ScalarVector(1, 0.01));
    properties.SetValue(STRESSES_OF_PIECEWISE_LINEAR_LAW, ScalarVector(2, 2.0));
    const auto geometry     = Geometry<Node>{};
    const auto process_info = ProcessInfo{};

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        [[maybe_unused]] const auto rv = backbone_law.Check(properties, geometry, process_info),
        "Error: The number of strain components does not match the number of stress components")
}

} // namespace Kratos::Testing
