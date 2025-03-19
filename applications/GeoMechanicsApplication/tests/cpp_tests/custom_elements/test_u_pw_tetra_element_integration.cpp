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
//

#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "tests/cpp_tests/test_utilities/model_setup_utilities.h"

namespace Kratos::Testing
{
using namespace Kratos;

KRATOS_TEST_CASE_IN_SUITE(UPwSmallStrainTetrahedra3D10ElementReturnsGI_GAUS_2, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    Model model;
    const auto& r_model_part = ModelSetupUtilities::CreateModelPartWithASingle3D10NUPwDiffOrderElement(model);
    const auto& element = r_model_part.Elements().front();

    // Act & Assert
    KRATOS_EXPECT_EQ(element.GetIntegrationMethod(), GeometryData::IntegrationMethod::GI_GAUSS_2);
}

KRATOS_TEST_CASE_IN_SUITE(UPwSmallStrainTriangle2D3ElementReturnsGI_GAUS_2, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    Model       model;
    const auto& r_model_part = ModelSetupUtilities::CreateModelPartWithASingle2D3NElement(model);
    const auto& element      = r_model_part.Elements().front();

    // Act & Assert
    KRATOS_EXPECT_EQ(element.GetIntegrationMethod(), GeometryData::IntegrationMethod::GI_GAUSS_2);
}

KRATOS_TEST_CASE_IN_SUITE(UPwSmallStrainTriangle2D6ElementReturnsGI_GAUS_2, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    Model model;
    const auto& r_model_part = ModelSetupUtilities::CreateModelPartWithASingle2D6NDiffOrderElement(model);
    const auto& element = r_model_part.Elements().front();

    // Act & Assert
    KRATOS_EXPECT_EQ(element.GetIntegrationMethod(), GeometryData::IntegrationMethod::GI_GAUSS_2);
}

KRATOS_TEST_CASE_IN_SUITE(UPwSmallStrainTriangle2D10ElementReturnsGI_GAUS_4, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    Model       model;
    const auto& r_model_part = ModelSetupUtilities::CreateModelPartWithASingle2D10NElement(model);
    const auto& element      = r_model_part.Elements().front();

    // Act & Assert
    KRATOS_EXPECT_EQ(element.GetIntegrationMethod(), GeometryData::IntegrationMethod::GI_GAUSS_4);
}

KRATOS_TEST_CASE_IN_SUITE(UPwSmallStrainTriangle2D15ElementReturnsGI_GAUS_5, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    Model       model;
    const auto& r_model_part = ModelSetupUtilities::CreateModelPartWithASingle2D15NElement(model);
    const auto& element      = r_model_part.Elements().front();

    // Act & Assert
    KRATOS_EXPECT_EQ(element.GetIntegrationMethod(), GeometryData::IntegrationMethod::GI_GAUSS_5);
}

} // namespace Kratos::Testing