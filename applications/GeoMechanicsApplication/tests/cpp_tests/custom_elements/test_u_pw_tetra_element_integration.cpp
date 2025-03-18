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

#include "custom_elements/U_Pw_small_strain_element.hpp"
#include "custom_elements/three_dimensional_stress_state.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "tests/cpp_tests/test_utilities.h"

namespace
{

using namespace Kratos;

auto UPwSmallStrainElementWithUPwDofs(Model& rModel, const Properties::Pointer& rProperties)
{
    auto& r_model_part = rModel.CreateModelPart("Main");

    PointerVector<Node> nodes;
    nodes.push_back(r_model_part.CreateNewNode(0, 0.0, 0.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(1, 1.0, 0.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(2, 0.0, 1.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(3, 0.0, 0.0, 1.0));
    nodes.push_back(r_model_part.CreateNewNode(4, 0.5, 0.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(5, 0.5, 0.5, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(6, 0.0, 0.5, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(7, 0.0, 0.0, 0.5));
    nodes.push_back(r_model_part.CreateNewNode(8, 0.5, 0.0, 0.5));
    nodes.push_back(r_model_part.CreateNewNode(9, 0.0, 0.5, 0.5));
    const auto p_geometry = std::make_shared<Tetrahedra3D10<Node>>(nodes);
    return UPwSmallStrainElement<3, 10>(
        1, p_geometry, rProperties, std::make_unique<ThreeDimensionalStressState>());
}

} // namespace

namespace Kratos::Testing
{
using namespace Kratos;

KRATOS_TEST_CASE_IN_SUITE(UPwSmallStrainTetrahedra3D10ElementReturnsGI_GAUS_2, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto p_properties = std::make_shared<Properties>();
    auto       process_info = ProcessInfo{};
    Model      model;
    auto       element = UPwSmallStrainElementWithUPwDofs(model, p_properties);

    // Act & Assert
    KRATOS_EXPECT_EQ(element.GetIntegrationMethod(), GeometryData::IntegrationMethod::GI_GAUSS_2);
}

} // namespace Kratos::Testing