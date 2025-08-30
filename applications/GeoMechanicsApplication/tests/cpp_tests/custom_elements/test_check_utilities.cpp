// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Gennady Markelov
//

#include "custom_utilities/check_utilities.h"
#include "geo_mechanics_application_variables.h"
#include "includes/checks.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"

using namespace Kratos;

namespace
{

Line2D2<Node> CreatLineGeometryWithoutVariables()
{
    PointerVector<Node> nodes;
    nodes.push_back(make_intrusive<Node>(0, 0.0, 0.0, 0.0));
    nodes.push_back(make_intrusive<Node>(1, 1.0, 0.0, 0.0));
    return Line2D2<Node>(nodes);
}

Line2D2<Node> CreatLineGeometryWithVariables()
{
    Model model;
    auto& r_model_part = model.CreateModelPart("ModelPart", 1);
    r_model_part.AddNodalSolutionStepVariable(WATER_PRESSURE);
    PointerVector<Node> nodes;
    nodes.push_back(r_model_part.CreateNewNode(0, 0.0, 0.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(1, 1.0, 0.0, 0.0));
    return Line2D2<Node>(nodes);
}
} // namespace

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(CheckUtilities_CheckDomainSize, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    constexpr std::size_t id = 1;
    PointerVector<Node>   nodes;
    nodes.push_back(make_intrusive<Node>(0, 0.0, 0.0, 0.0));
    nodes.push_back(make_intrusive<Node>(1, 0.0, 0.0, 0.0));
    const auto line_with_coincident_nodes = Line2D2<Node>(nodes);
    nodes.push_back(make_intrusive<Node>(2, 0.0, 0.0, 0.0));
    const auto triangle_with_coincident_nodes = Triangle2D3<Node>(nodes);
    nodes.push_back(make_intrusive<Node>(2, 0.0, 0.0, 0.0));
    const auto tetra_with_coincident_nodes = Tetrahedra3D4<Node>(nodes);

    // Act and Assert
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        CheckUtilities::CheckDomainSize(line_with_coincident_nodes.DomainSize(), id, "Length"),
        "Length (0) is smaller than 1e-15 for element 1")

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        CheckUtilities::CheckDomainSize(triangle_with_coincident_nodes.DomainSize(), id),
        "DomainSize (0) is smaller than 1e-15 for element 1")

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        CheckUtilities::CheckDomainSize(tetra_with_coincident_nodes.DomainSize(), id),
        "DomainSize (0) is smaller than 1e-15 for element 1")
}

KRATOS_TEST_CASE_IN_SUITE(CheckUtilities_CheckNodalVariables, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto line_geometry = CreatLineGeometryWithoutVariables();

    // Act and Assert
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        CheckUtilities::CheckHasNodalSolutionStepData(
            line_geometry, {std::cref(WATER_PRESSURE), std::cref(VOLUME_ACCELERATION)}),
        "Missing variable WATER_PRESSURE on nodes 0 1")

    // Arrange 2
    line_geometry = CreatLineGeometryWithVariables();

    // Act and Assert
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        CheckUtilities::CheckHasNodalSolutionStepData(
            line_geometry, {std::cref(WATER_PRESSURE), std::cref(VOLUME_ACCELERATION)}),
        "Missing variable VOLUME_ACCELERATION on nodes 0 1")
}

KRATOS_TEST_CASE_IN_SUITE(CheckUtilities_CheckNodalDof, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto line_geometry = CreatLineGeometryWithoutVariables();

    // Act and Assert
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        CheckUtilities::CheckHasDofs(line_geometry, {std::cref(WATER_PRESSURE), std::cref(VOLUME_ACCELERATION)}),
        "Missing the DoF for the variable WATER_PRESSURE on nodes 0 1")

    // Arrange 2
    line_geometry = CreatLineGeometryWithVariables();

    for (auto& r_node : line_geometry) {
        r_node.AddDof(WATER_PRESSURE);
    }

    // Act and Assert
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        CheckUtilities::CheckHasDofs(line_geometry, {std::cref(WATER_PRESSURE), std::cref(VOLUME_ACCELERATION)}),
        "Missing the DoF for the variable VOLUME_ACCELERATION on nodes 0 1")
}

KRATOS_TEST_CASE_IN_SUITE(CheckUtilities_CheckProperties, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto properties = Properties();
    properties.SetValue(DENSITY_WATER, 1000.0);
    properties.SetValue(DENSITY, 0.0);
    CheckProperties check_properties("property", properties, CheckProperties::Bounds::AllInclusive);

    // Act and Assert
    EXPECT_NO_THROW(check_properties.Check(DENSITY_WATER));
    EXPECT_NO_THROW(check_properties.Check(DENSITY));

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(check_properties.Check(DENSITY_WATER, 500.0),
                                      "DENSITY_WATER in the property 0 has an invalid value: 1000 "
                                      "out of the range [0; 500.000000].")
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(check_properties.Check(DENSITY_WATER, 100.0, 500.0),
                                      "DENSITY_WATER in the property 0 has an invalid value: 1000 "
                                      "out of the range [100; 500.000000].")

    check_properties.SetNewBounds(CheckProperties::Bounds::AllExclusive);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(check_properties.Check(DENSITY_WATER, 1000.0),
                                      "DENSITY_WATER in the property 0 has an invalid value: 1000 "
                                      "out of the range (0; 1000.000000).")
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        check_properties.Check(DENSITY),
        "DENSITY in the property 0 has an invalid value: 0 out of the range (0; -).")
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        check_properties.SingleUseBounds(CheckProperties::Bounds::InclusiveLowerAndExclusiveUpper).Check(DENSITY_WATER, 1000.0), "DENSITY_WATER in the property 0 has an invalid value: 1000 out of the range [0; 1000.000000).")
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        check_properties.SingleUseBounds(CheckProperties::Bounds::ExclusiveLowerAndInclusiveUpper).Check(DENSITY),
        "DENSITY in the property 0 has an invalid value: 0 out of the range (0; -].")

    // Arrange 2
    auto            properties_2 = Properties();
    constexpr auto  element_Id   = 1;
    CheckProperties check_properties_2("property at element", properties_2, element_Id,
                                       CheckProperties::Bounds::AllInclusive);
    // Act and Assert
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(check_properties_2.CheckAvailabilityOnly(UDSM_NAME),
                                      " UDSM_NAME does not exist in the property at element 1.")
    properties_2.SetValue(UDSM_NAME, "");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(check_properties_2.CheckAvailabilityAndEmpty(UDSM_NAME),
                                      "UDSM_NAME is empty in the property at element 1.");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(check_properties_2.CheckAvailabilityOnly(UDSM_NUMBER),
                                      "UDSM_NUMBER does not exist in the property at element 1.")
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        check_properties_2.CheckAvailabilityOnly(IS_FORTRAN_UDSM),
        "IS_FORTRAN_UDSM does not exist in the property at element 1.")
}
} // namespace Kratos::Testing
