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

#include "custom_constitutive/incremental_linear_elastic_interface_law.h"
#include "custom_geometries/line_interface_geometry.h"
#include "geo_mechanics_application_variables.h"
#include "geometries/line_3d_2.h"
#include "includes/checks.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"

using namespace Kratos;

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(LinearElasticLawForInterfacesHas2DWorkingSpace, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto law = GeoIncrementalLinearElasticInterfaceLaw{};

    KRATOS_EXPECT_EQ(law.WorkingSpaceDimension(), 2);
}

KRATOS_TEST_CASE_IN_SUITE(LinearElasticLawForInterfacesHasStrainSizeOfTwo, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto law = GeoIncrementalLinearElasticInterfaceLaw{};

    KRATOS_EXPECT_EQ(law.GetStrainSize(), 2);
}

KRATOS_TEST_CASE_IN_SUITE(LinearElasticLawForInterfacesUsesInfinitesimalStrains, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto law = GeoIncrementalLinearElasticInterfaceLaw{};

    KRATOS_EXPECT_EQ(law.GetStressMeasure(), ConstitutiveLaw::StressMeasure_Cauchy);
}

KRATOS_TEST_CASE_IN_SUITE(LinearElasticLawForInterfacesIsIncremental, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto law = GeoIncrementalLinearElasticInterfaceLaw{};

    KRATOS_EXPECT_TRUE(law.IsIncremental())
}

KRATOS_TEST_CASE_IN_SUITE(LinearElasticLawForInterfacesChecksForCorrectMaterialProperties,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto law          = GeoIncrementalLinearElasticInterfaceLaw{};
    auto       properties   = Properties{};
    const auto geometry     = LineInterfaceGeometry{};
    const auto process_info = ProcessInfo{};

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(law.Check(properties, geometry, process_info),
                                      "No interface normal stiffness defined")

    properties[INTERFACE_NORMAL_STIFFNESS] = -5.0;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(law.Check(properties, geometry, process_info),
                                      "Interface normal stiffness must be positive, but got -5")

    properties[INTERFACE_NORMAL_STIFFNESS] = 0.0;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(law.Check(properties, geometry, process_info),
                                      "Interface normal stiffness must be positive, but got 0")

    properties[INTERFACE_NORMAL_STIFFNESS] = 5.0;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(law.Check(properties, geometry, process_info),
                                      "No interface shear stiffness defined")

    properties[INTERFACE_SHEAR_STIFFNESS] = -2.5;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(law.Check(properties, geometry, process_info),
                                      "Interface shear stiffness must be positive, but got -2.5")

    properties[INTERFACE_SHEAR_STIFFNESS] = 0.0;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(law.Check(properties, geometry, process_info),
                                      "Interface shear stiffness must be positive, but got 0")

    properties[INTERFACE_SHEAR_STIFFNESS] = 2.5;
    KRATOS_EXPECT_EQ(law.Check(properties, geometry, process_info), 0);
}

KRATOS_TEST_CASE_IN_SUITE(LinearElasticLawForInterfacesChecksForCorrectGeometry, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto law                         = GeoIncrementalLinearElasticInterfaceLaw{};
    auto       properties                  = Properties{};
    properties[INTERFACE_NORMAL_STIFFNESS] = 5.0;
    properties[INTERFACE_SHEAR_STIFFNESS]  = 2.5;
    const auto node1                       = make_intrusive<Node>(1, 0.0, 0.0, 0.0);
    const auto node2                       = make_intrusive<Node>(2, 5.0, 5.0, 5.0);
    const auto line_3d_geometry            = Line3D2<Node>{node1, node2};
    const auto process_info                = ProcessInfo{};

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        law.Check(properties, line_3d_geometry, process_info),
        "Expected a line interface geometry, but got 1 dimensional line with 2 nodes in 3D space")
}

} // namespace Kratos::Testing
