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
//                   Richard Faasse
//

#include "custom_utilities/check_utilities.hpp"
#include "custom_utilities/ublas_utilities.h"
#include "geo_mechanics_application_variables.h"
#include "includes/checks.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "tests/cpp_tests/stub_constitutive_law.h"

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

KRATOS_TEST_CASE_IN_SUITE(CheckUtilities_CheckPropertiesThatPrintsPropertyId, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto properties = Properties{};
    properties.SetValue(DENSITY_WATER, 1000.0);
    properties.SetValue(DENSITY, 0.0);
    const CheckProperties check_properties(properties, "property", CheckProperties::Bounds::AllInclusive);

    // Act and Assert
    EXPECT_NO_THROW(check_properties.Check(DENSITY_WATER));
    EXPECT_NO_THROW(check_properties.Check(DENSITY));

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        check_properties.Check(DENSITY_WATER, 500.0),
        "DENSITY_WATER in the property with Id 0 has an invalid value: 1000 "
        "is out of the range [0, 500].")
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        check_properties.Check(DENSITY_WATER, 100.0, 500.0),
        "DENSITY_WATER in the property with Id 0 has an invalid value: 1000 "
        "is out of the range [100, 500].")

    check_properties.SetNewRangeBounds(CheckProperties::Bounds::AllExclusive);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        check_properties.Check(DENSITY_WATER, 1000.0),
        "DENSITY_WATER in the property with Id 0 has an invalid value: 1000 "
        "is out of the range (0, 1000).")
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        check_properties.Check(DENSITY),
        "DENSITY in the property with Id 0 has an invalid value: 0 is out of the range (0, -).")
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        check_properties.SingleUseBounds(CheckProperties::Bounds::InclusiveLowerAndExclusiveUpper).Check(DENSITY_WATER, 1000.0), "DENSITY_WATER in the property with Id 0 has an invalid value: 1000 is out of the range [0, 1000).")
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        check_properties.SingleUseBounds(CheckProperties::Bounds::ExclusiveLowerAndInclusiveUpper).Check(DENSITY),
        "DENSITY in the property with Id 0 has an invalid value: 0 is out of the range (0, -).")
}

KRATOS_TEST_CASE_IN_SUITE(CheckUtilities_CheckPropertiesThatPrintsElementId, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto                  properties = Properties{};
    constexpr auto        element_Id = 1;
    const CheckProperties check_properties(properties, "property", element_Id,
                                           CheckProperties::Bounds::AllInclusive);
    // Act and Assert
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        check_properties.CheckAvailability(UDSM_NAME),
        "UDSM_NAME does not exist in the property with Id 0 at element with Id 1.")
    properties.SetValue(UDSM_NAME, "");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        check_properties.CheckAvailabilityAndNotEmpty(UDSM_NAME),
        "UDSM_NAME is empty in the property with Id 0 at element with Id 1.");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        check_properties.CheckAvailability(UDSM_NUMBER),
        "UDSM_NUMBER does not exist in the property with Id 0 at element with Id 1.")
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        check_properties.CheckAvailability(IS_FORTRAN_UDSM),
        "IS_FORTRAN_UDSM does not exist in the property with Id 0 at element with Id 1.")
}

KRATOS_TEST_CASE_IN_SUITE(CheckUtilities_CheckIntegerProperty, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto properties = Properties{};
    const CheckProperties check_properties(properties, "property", CheckProperties::Bounds::AllInclusive);
    // Act and Assert
    properties.SetValue(UDSM_NUMBER, 3);
    EXPECT_NO_THROW(check_properties.Check(UDSM_NUMBER));
    EXPECT_NO_THROW(check_properties.Check(UDSM_NUMBER, 3, 5));
    EXPECT_NO_THROW(check_properties.Check(UDSM_NUMBER, 1, 3));

    check_properties.SetNewRangeBounds(CheckProperties::Bounds::AllExclusive);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        check_properties.Check(UDSM_NUMBER, 3, 5),
        "UDSM_NUMBER in the property with Id 0 has an invalid value: 3 is out of the range (3, 5).")
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        check_properties.Check(UDSM_NUMBER, 1, 3),
        "UDSM_NUMBER in the property with Id 0 has an invalid value: 3 is out of the range (1, 3).")
}

KRATOS_TEST_CASE_IN_SUITE(CheckUtilities_CheckAvailabilityAndEquality, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto properties = Properties{};
    const CheckProperties check_properties(properties, "property", CheckProperties::Bounds::AllInclusive);
    // Act and Assert
    properties.SetValue(UDSM_NUMBER, 3);
    EXPECT_NO_THROW(check_properties.CheckAvailabilityAndEquality(UDSM_NUMBER, 3));
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        check_properties.CheckAvailabilityAndEquality(UDSM_NUMBER, 5),
        "UDSM_NUMBER has a value of 3 instead of 5 in the property with Id 0.");

    properties.SetValue(RETENTION_LAW, "SaturatedLaw");
    EXPECT_NO_THROW(check_properties.CheckAvailabilityAndEquality(RETENTION_LAW, "SaturatedLaw"));

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(check_properties.CheckAvailabilityAndEquality(RETENTION_LAW, "VanGenuchtenLaw"), "RETENTION_LAW has a value of \"SaturatedLaw\" instead of \"VanGenuchtenLaw\" in the property with Id 0.")
}

KRATOS_TEST_CASE_IN_SUITE(CheckUtilities_CheckPermeabilityPropertiesThrowsErrorsForWrongProperties,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Properties properties(2);
    const CheckProperties check_properties(properties, "properties", CheckProperties::Bounds::AllExclusive);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(check_properties.CheckPermeabilityProperties(1),
                                      "PERMEABILITY_XX does not exist in the properties with Id 2.")

    properties.SetValue(PERMEABILITY_XX, -10.0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(check_properties.CheckPermeabilityProperties(1),
                                      "PERMEABILITY_XX in the properties with Id 2 has an invalid "
                                      "value: -10 is out of the range [0, -).")

    properties.SetValue(PERMEABILITY_XX, 10.0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(check_properties.CheckPermeabilityProperties(2),
                                      "PERMEABILITY_YY does not exist in the properties with Id 2.")

    properties.SetValue(PERMEABILITY_YY, 10.0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(check_properties.CheckPermeabilityProperties(2),
                                      "PERMEABILITY_XY does not exist in the properties with Id 2.")
    properties.SetValue(PERMEABILITY_XY, 0.0);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(check_properties.CheckPermeabilityProperties(3),
                                      "PERMEABILITY_ZZ does not exist in the properties with Id 2.")
    properties.SetValue(PERMEABILITY_ZZ, 10.0);
    properties.SetValue(PERMEABILITY_YZ, 0.0);
    properties.SetValue(PERMEABILITY_ZX, -10.0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(check_properties.CheckPermeabilityProperties(3),
                                      "PERMEABILITY_ZX in the properties with Id 2 has an invalid "
                                      "value: -10 is out of the range [0, -).")
}

KRATOS_TEST_CASE_IN_SUITE(CheckUtilities_CheckPermeabilityPropertiesDoesNotThrowsErrors,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Properties properties(2);
    const CheckProperties check_properties(properties, "properties", CheckProperties::Bounds::AllExclusive);

    properties.SetValue(PERMEABILITY_XX, 10.0);
    EXPECT_NO_THROW(check_properties.CheckPermeabilityProperties(1));

    properties.SetValue(PERMEABILITY_YY, 10.0);
    properties.SetValue(PERMEABILITY_XY, 0.0);
    EXPECT_NO_THROW(check_properties.CheckPermeabilityProperties(2));

    properties.SetValue(PERMEABILITY_ZZ, 10.0);
    properties.SetValue(PERMEABILITY_YZ, 0.0);
    properties.SetValue(PERMEABILITY_ZX, 0.0);
    EXPECT_NO_THROW(check_properties.CheckPermeabilityProperties(3));
}

KRATOS_TEST_CASE_IN_SUITE(CheckUtilities_CheckForNonZeroZCoordinateIn2D, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto line_geometry = CreatLineGeometryWithVariables();

    // Act and Assert
    EXPECT_NO_THROW(CheckUtilities::CheckForNonZeroZCoordinateIn2D(line_geometry));

    // Arrange 2
    line_geometry.begin()->Z() += 1;

    // Act and Assert
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(CheckUtilities::CheckForNonZeroZCoordinateIn2D(line_geometry),
                                      "Node with Id: 0 has non-zero Z coordinate.")
}

KRATOS_TEST_CASE_IN_SUITE(CheckUtilities_CheckAvailabilityAndSpecified, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    Properties properties(1);
    const CheckProperties check_properties(properties, "properties", CheckProperties::Bounds::AllExclusive);

    // Act and Assert
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        check_properties.CheckAvailabilityAndSpecified(CONSTITUTIVE_LAW),
        " CONSTITUTIVE_LAW does not exist in the properties with Id 1.")

    properties.SetValue(CONSTITUTIVE_LAW, nullptr);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        check_properties.CheckAvailabilityAndSpecified(CONSTITUTIVE_LAW),
        "CONSTITUTIVE_LAW needs to be specified in the properties with Id 1.")

    properties.SetValue(CONSTITUTIVE_LAW, std::make_shared<StubConstitutiveLaw>());
    EXPECT_NO_THROW(check_properties.CheckAvailabilityAndSpecified(CONSTITUTIVE_LAW));
}

KRATOS_TEST_CASE_IN_SUITE(CheckValuesAreAscending_ValidNonDecreasingSequence, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto values = UblasUtilities::CreateVector({1.0, 2.5, 2.5, 5.0, 10.0, 10.0, 15.5});

    // Act & Assert - should not throw when AllowEqual=true
    EXPECT_NO_THROW(CheckUtilities::CheckValuesAreAscending(values, "test_values", true));
}

KRATOS_TEST_CASE_IN_SUITE(CheckValuesAreAscending_InvalidSequenceWithDuplicatesWhenNotAllowed,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto values = UblasUtilities::CreateVector({1.0, 2.5, 2.5, 5.0, 10.0});

    // Act & Assert - should throw when AllowEqual=false and there are duplicates
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(CheckUtilities::CheckValuesAreAscending(values, "test_values", false), "Invalid sequence in 'test_values': values must be strictly increasing. Found 2.5 at index 3, which is not greater than the previous value 2.5 at index 2.");
}

KRATOS_TEST_CASE_IN_SUITE(CheckValuesAreAscending_InvalidSequenceDecreasing, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto values = UblasUtilities::CreateVector({1.0, 5.0, 3.0, 10.0});

    // Act & Assert when duplicates are not allowed
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(CheckUtilities::CheckValuesAreAscending(values, "test_values", false), "Invalid sequence in 'test_values': values must be strictly increasing. Found 3 at index 3, which is not greater than the previous value 5 at index 2.");

    // Act & Assert when duplicates are allowed
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(CheckUtilities::CheckValuesAreAscending(values, "test_values", true), "Invalid sequence in 'test_values': values must be non-decreasing (duplicates allowed). Found 3 at index 3, which is less than the previous value 5 at index 2.");
}

KRATOS_TEST_CASE_IN_SUITE(CheckValuesAreAscending_SingleElement, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto values = UblasUtilities::CreateVector({5.0});

    // Act & Assert - single element should always pass
    EXPECT_NO_THROW(CheckUtilities::CheckValuesAreAscending(values, "test_values", false));
    EXPECT_NO_THROW(CheckUtilities::CheckValuesAreAscending(values, "test_values", true));
}

KRATOS_TEST_CASE_IN_SUITE(CheckValuesAreAscending_AscendingSequence, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange - create a large ascending sequence
    Vector values(1000);
    for (int i = 0; i < 1000; ++i) {
        values[i] = static_cast<double>(i) * 0.1 - 50.0;
    }

    // Act & Assert
    EXPECT_NO_THROW(CheckUtilities::CheckValuesAreAscending(values, "test_values", false));
    EXPECT_NO_THROW(CheckUtilities::CheckValuesAreAscending(values, "test_values", true));
}

} // namespace Kratos::Testing
