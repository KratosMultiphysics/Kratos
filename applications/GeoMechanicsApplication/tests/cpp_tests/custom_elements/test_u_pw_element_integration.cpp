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
//                   Richard Faasse

#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "tests/cpp_tests/test_utilities/model_setup_utilities.h"

namespace Kratos::Testing
{
using namespace Kratos;

class ParametrizedIntegrationMethodSuite
    : public ::testing::TestWithParam<std::tuple<std::string, GeometryData::IntegrationMethod>>
{
};

const Element& CreateElementWithSpecificGeometry(Model& rModel, const std::string& rGeometryDescription)
{
    if (rGeometryDescription == "Triangle2D3N") {
        return ModelSetupUtilities::CreateModelPartWithASingle2D3NElement(rModel).Elements().front();
    }
    if (rGeometryDescription == "Triangle2D6N") {
        return ModelSetupUtilities::CreateModelPartWithASingle2D6NDiffOrderElement(rModel)
            .Elements()
            .front();
    }
    if (rGeometryDescription == "Triangle2D10N") {
        return ModelSetupUtilities::CreateModelPartWithASingle2D10NElement(rModel).Elements().front();
    }
    if (rGeometryDescription == "Triangle2D15N") {
        return ModelSetupUtilities::CreateModelPartWithASingle2D15NElement(rModel).Elements().front();
    }
    if (rGeometryDescription == "Tetrahedra3D10N") {
        return ModelSetupUtilities::CreateModelPartWithASingle3D10NUPwDiffOrderElement(rModel)
            .Elements()
            .front();
    }

    KRATOS_ERROR << "Geometry description not recognized: " << rGeometryDescription << std::endl;
}

TEST_P(ParametrizedIntegrationMethodSuite, TestElementReturnsCorrectIntegrationMethod)
{
    Model model;
    const auto& [geometry_description, expected_integration_method] = GetParam();
    const auto& r_element = CreateElementWithSpecificGeometry(model, geometry_description);

    KRATOS_EXPECT_EQ(r_element.GetIntegrationMethod(), expected_integration_method)
        << "\nIntegration method is not as expected for the element with the following geometry: "
        << r_element.GetGeometry().Info();
}

INSTANTIATE_TEST_SUITE_P(
    KratosGeoMechanicsFastSuiteWithoutKernel,
    ParametrizedIntegrationMethodSuite,
    ::testing::Values(std::make_tuple("Tetrahedra3D10N", GeometryData::IntegrationMethod::GI_GAUSS_2),
                      std::make_tuple("Triangle2D3N", GeometryData::IntegrationMethod::GI_GAUSS_2),
                      std::make_tuple("Triangle2D6N", GeometryData::IntegrationMethod::GI_GAUSS_2),
                      std::make_tuple("Triangle2D10N", GeometryData::IntegrationMethod::GI_GAUSS_4),
                      std::make_tuple("Triangle2D15N", GeometryData::IntegrationMethod::GI_GAUSS_5)));

} // namespace Kratos::Testing