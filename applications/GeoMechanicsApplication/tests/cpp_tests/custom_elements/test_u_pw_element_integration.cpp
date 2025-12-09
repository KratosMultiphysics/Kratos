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

#include "containers/model.h"
#include "test_setup_utilities/element_setup_utilities.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite_without_kernel.h"

namespace Kratos::Testing
{
using namespace Kratos;

class ParametrizedIntegrationMethodSuite
    : public ::testing::TestWithParam<std::tuple<Element::Pointer, GeometryData::IntegrationMethod>>
{
};

TEST_P(ParametrizedIntegrationMethodSuite, TestElementReturnsCorrectIntegrationMethod)
{
    Model model;
    const auto& [r_element, expected_integration_method] = GetParam();

    EXPECT_EQ(r_element->GetIntegrationMethod(), expected_integration_method)
        << "\nIntegration method is not as expected for the element with the following geometry: "
        << r_element->GetGeometry().Info();
}

INSTANTIATE_TEST_SUITE_P(
    KratosGeoMechanicsFastSuiteWithoutKernel,
    ParametrizedIntegrationMethodSuite,
    ::testing::Values(
        std::make_tuple(ElementSetupUtilities::Create2D3NElement(), GeometryData::IntegrationMethod::GI_GAUSS_2),
        std::make_tuple(ElementSetupUtilities::Create2D6NElement(), GeometryData::IntegrationMethod::GI_GAUSS_2),
        std::make_tuple(ElementSetupUtilities::Create2D10NElement(), GeometryData::IntegrationMethod::GI_GAUSS_4),
        std::make_tuple(ElementSetupUtilities::Create2D15NElement(), GeometryData::IntegrationMethod::GI_GAUSS_5),
        std::make_tuple(ElementSetupUtilities::Create3D10NElement(), GeometryData::IntegrationMethod::GI_GAUSS_2)));

} // namespace Kratos::Testing