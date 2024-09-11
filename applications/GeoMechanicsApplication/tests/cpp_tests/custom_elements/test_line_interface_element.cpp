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

#include "custom_constitutive/incremental_linear_elastic_interface_law.h"
#include "custom_elements/line_interface_element.h"
#include "custom_geometries/line_interface_geometry.h"
#include "geo_mechanics_application_variables.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "tests/cpp_tests/test_utilities.h"

#include <boost/numeric/ublas/assignment.hpp>
#include <cstddef>

namespace
{

using namespace Kratos;

PointerVector<Node> CreateNodes()
{
    PointerVector<Node> result;
    result.push_back(Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0));
    result.push_back(Kratos::make_intrusive<Node>(2, 1.0, 0.0, 0.0));
    result.push_back(Kratos::make_intrusive<Node>(3, 0.0, 0.0, 0.0));
    result.push_back(Kratos::make_intrusive<Node>(4, 1.0, 0.0, 0.0));

    return result;
}

} // namespace

namespace Kratos::Testing
{

using namespace Kratos;

KRATOS_TEST_CASE_IN_SUITE(LineInterfaceElementIsAnElement, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    LineInterfaceElement element;
    auto                 casted_element = dynamic_cast<Element*>(&element);
    KRATOS_CHECK_NOT_EQUAL(casted_element, nullptr);
}

KRATOS_TEST_CASE_IN_SUITE(LineInterfaceElementCanCreateInstanceWithGeometryInput, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const LineInterfaceElement element;
    const auto geometry   = std::make_shared<LineInterfaceGeometry<Line2D2<Node>>>(CreateNodes());
    auto       properties = std::make_shared<Properties>();

    // Act
    auto created_element = element.Create(1, geometry, properties);

    // Assert
    EXPECT_NE(created_element, nullptr);
    EXPECT_EQ(created_element->Id(), 1);
    EXPECT_NE(created_element->pGetGeometry(), nullptr);
    EXPECT_NE(created_element->pGetProperties(), nullptr);
}

KRATOS_TEST_CASE_IN_SUITE(LineInterfaceElementCanCreateInstanceWithNodeInput, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto nodes      = CreateNodes();
    auto properties = std::make_shared<Properties>();

    // The source element needs to have a geometry, otherwise the version of the
    // Create method with a node input will fail.
    const auto geometry = std::make_shared<LineInterfaceGeometry<Line2D2<Node>>>(nodes);
    const LineInterfaceElement element(1, geometry, properties);

    // Act
    auto created_element = element.Create(1, nodes, properties);

    // Assert
    EXPECT_NE(created_element, nullptr);
    EXPECT_EQ(created_element->Id(), 1);
    EXPECT_NE(created_element->pGetGeometry(), nullptr);
    EXPECT_NE(created_element->pGetProperties(), nullptr);
}

KRATOS_TEST_CASE_IN_SUITE(LineInterfaceElement_ReturnsTheExpectedDoFList, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto properties = std::make_shared<Properties>();

    Model model;
    auto& model_part = model.CreateModelPart("Main");
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT);

    PointerVector<Node> result;
    result.push_back(model_part.CreateNewNode(0, 0.0, 0.0, 0.0));
    result.push_back(model_part.CreateNewNode(1, 1.0, 0.0, 0.0));
    result.push_back(model_part.CreateNewNode(2, 0.0, 0.0, 0.0));
    result.push_back(model_part.CreateNewNode(3, 1.0, 0.0, 0.0));
    auto geometry = std::make_shared<LineInterfaceGeometry<Line2D2<Node>>>(result);
    auto element  = make_intrusive<LineInterfaceElement>(1, geometry, properties);

    model_part.AddElement(element);
    for (auto& node : element->GetGeometry()) {
        node.AddDof(DISPLACEMENT_X);
        node.AddDof(DISPLACEMENT_Y);
    }

    element->GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT) = array_1d<double, 3>{1.0, 2.0, 0.0};
    element->GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT) = array_1d<double, 3>{3.0, 4.0, 0.0};
    element->GetGeometry()[2].FastGetSolutionStepValue(DISPLACEMENT) = array_1d<double, 3>{5.0, 6.0, 0.0};
    element->GetGeometry()[3].FastGetSolutionStepValue(DISPLACEMENT) = array_1d<double, 3>{7.0, 8.0, 0.0};

    // Act
    Element::DofsVectorType degrees_of_freedom;
    element->GetDofList(degrees_of_freedom, {});

    // Assert
    KRATOS_EXPECT_EQ(degrees_of_freedom.size(), 8);
    const std::vector<double> expected_dof_values = {1.0, 2.0, 3.0, 4.0,  5.0,  6.0,
                                                     7.0, 8.0, 9.0, 10.0, 11.0, 12.0};
    for (std::size_t i = 0; i < degrees_of_freedom.size(); i++) {
        KRATOS_EXPECT_DOUBLE_EQ(degrees_of_freedom[i]->GetSolutionStepValue(), expected_dof_values[i]);
    }
}

KRATOS_TEST_CASE_IN_SUITE(LineInterfaceElement_ReturnsTheExpectedEquationIdVector, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto properties = std::make_shared<Properties>();

    Model model;
    auto& model_part = model.CreateModelPart("Main");
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT);

    PointerVector<Node> result;
    result.push_back(model_part.CreateNewNode(0, 0.0, 0.0, 0.0));
    result.push_back(model_part.CreateNewNode(1, 1.0, 0.0, 0.0));
    result.push_back(model_part.CreateNewNode(2, 0.0, 0.0, 0.0));
    result.push_back(model_part.CreateNewNode(3, 1.0, 0.0, 0.0));
    auto geometry = std::make_shared<LineInterfaceGeometry<Line2D2<Node>>>(result);
    auto element  = make_intrusive<LineInterfaceElement>(1, geometry, properties);

    model_part.AddElement(element);
    int i = 0;
    for (auto& node : element->GetGeometry()) {
        node.AddDof(DISPLACEMENT_X);
        node.AddDof(DISPLACEMENT_Y);

        ++i;
        node.pGetDof(DISPLACEMENT_X)->SetEquationId(i);
        ++i;
        node.pGetDof(DISPLACEMENT_Y)->SetEquationId(i);
    }

    // Act
    Element::EquationIdVectorType equation_id_vector;
    element->EquationIdVector(equation_id_vector, {});

    // Assert
    KRATOS_EXPECT_EQ(equation_id_vector.size(), 8);
    const std::vector<int> expected_ids = {1, 2, 3, 4, 5, 6, 7, 8};
    KRATOS_EXPECT_VECTOR_EQ(equation_id_vector, expected_ids)
}

KRATOS_TEST_CASE_IN_SUITE(LineInterfaceElement_LeftHandSideContainsMaterialStiffnessContributions,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto           properties                        = std::make_shared<Properties>();
    constexpr auto normal_stiffness                  = 20.0;
    constexpr auto shear_stiffness                   = 10.0;
    properties->GetValue(INTERFACE_NORMAL_STIFFNESS) = normal_stiffness;
    properties->GetValue(INTERFACE_SHEAR_STIFFNESS)  = shear_stiffness;
    properties->GetValue(CONSTITUTIVE_LAW) = std::make_shared<GeoIncrementalLinearElasticInterfaceLaw>();

    Model model;
    auto& model_part = model.CreateModelPart("Main");
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT);

    PointerVector<Node> result;
    result.push_back(model_part.CreateNewNode(0, 0.0, 0.0, 0.0));
    result.push_back(model_part.CreateNewNode(1, 1.0, 0.0, 0.0));
    result.push_back(model_part.CreateNewNode(2, 0.0, 0.0, 0.0));
    result.push_back(model_part.CreateNewNode(3, 1.0, 0.0, 0.0));
    auto geometry = std::make_shared<LineInterfaceGeometry<Line2D2<Node>>>(result);
    auto element  = make_intrusive<LineInterfaceElement>(1, geometry, properties);

    model_part.AddElement(element);
    for (auto& node : element->GetGeometry()) {
        node.AddDof(DISPLACEMENT_X);
        node.AddDof(DISPLACEMENT_Y);
    }

    // Act
    Matrix left_hand_side;
    element->CalculateLeftHandSide(left_hand_side, {});

    // Assert
    auto expected_left_hand_side  = Matrix{IdentityMatrix{8}};
    expected_left_hand_side(0, 0) = shear_stiffness * 0.5;
    expected_left_hand_side(1, 1) = normal_stiffness * 0.5;
    expected_left_hand_side(2, 2) = shear_stiffness * 0.5;
    expected_left_hand_side(3, 3) = normal_stiffness * 0.5;
    expected_left_hand_side(4, 4) = shear_stiffness * 0.5;
    expected_left_hand_side(5, 5) = normal_stiffness * 0.5;
    expected_left_hand_side(6, 6) = shear_stiffness * 0.5;
    expected_left_hand_side(7, 7) = normal_stiffness * 0.5;
    expected_left_hand_side(0, 4) = -shear_stiffness * 0.5;
    expected_left_hand_side(1, 5) = -normal_stiffness * 0.5;
    expected_left_hand_side(2, 6) = -shear_stiffness * 0.5;
    expected_left_hand_side(3, 7) = -normal_stiffness * 0.5;
    expected_left_hand_side(4, 0) = -shear_stiffness * 0.5;
    expected_left_hand_side(5, 1) = -normal_stiffness * 0.5;
    expected_left_hand_side(6, 2) = -shear_stiffness * 0.5;
    expected_left_hand_side(7, 3) = -normal_stiffness * 0.5;
    KRATOS_EXPECT_MATRIX_RELATIVE_NEAR(left_hand_side, expected_left_hand_side, Defaults::relative_tolerance)
}

KRATOS_TEST_CASE_IN_SUITE(LineInterfaceElement_LeftHandSideContainsMaterialStiffnessContributions_Rotated,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto           properties                        = std::make_shared<Properties>();
    constexpr auto normal_stiffness                  = 20.0;
    constexpr auto shear_stiffness                   = 10.0;
    properties->GetValue(INTERFACE_NORMAL_STIFFNESS) = normal_stiffness;
    properties->GetValue(INTERFACE_SHEAR_STIFFNESS)  = shear_stiffness;
    properties->GetValue(CONSTITUTIVE_LAW) = std::make_shared<GeoIncrementalLinearElasticInterfaceLaw>();

    Model model;
    auto& model_part = model.CreateModelPart("Main");
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT);

    PointerVector<Node> result;
    result.push_back(model_part.CreateNewNode(0, 0.0, 0.0, 0.0));
    result.push_back(model_part.CreateNewNode(1, 0.5 * std::sqrt(3.0), 0.5, 0.0));
    result.push_back(model_part.CreateNewNode(2, 0.0, 0.0, 0.0));
    result.push_back(model_part.CreateNewNode(3, 0.5 * std::sqrt(3.0), 0.5, 0.0));
    auto geometry = std::make_shared<LineInterfaceGeometry<Line2D2<Node>>>(result);
    auto element  = make_intrusive<LineInterfaceElement>(1, geometry, properties);

    model_part.AddElement(element);
    for (auto& node : element->GetGeometry()) {
        node.AddDof(DISPLACEMENT_X);
        node.AddDof(DISPLACEMENT_Y);
    }

    // Act
    Matrix left_hand_side;
    element->CalculateLeftHandSide(left_hand_side, {});

    // Assert

    auto expected_left_hand_side = Matrix{8, 8};
    // clang-format off
    expected_left_hand_side <<=
         6.25,       -2.16506351,  0.0,         0.0,        -6.25,        2.16506351,  0.0,         0.0,
        -2.16506351,  8.75,        0.0,         0.0,         2.16506351, -8.75,        0.0,         0.0,
         0.0,         0.0,         6.25,       -2.16506351,  0.0,         0.0,        -6.25,        2.16506351,
         0.0,         0.0,        -2.16506351,  8.75,        0.0,         0.0,         2.16506351, -8.75,
        -6.25,        2.16506351,  0.0,         0.0,         6.25,       -2.16506351,  0.0,         0.0,
         2.16506351, -8.75,        0.0,         0.0,        -2.16506351,  8.75,        0.0,         0.0,
         0.0,         0.0,        -6.25,        2.16506351,  0.0,         0.0,         6.25,       -2.16506351,
         0.0,         0.0,         2.16506351, -8.75,        0.0,         0.0,        -2.16506351,  8.75;
    // clang-format on
    KRATOS_EXPECT_MATRIX_RELATIVE_NEAR(left_hand_side, expected_left_hand_side, Defaults::relative_tolerance)
}

KRATOS_TEST_CASE_IN_SUITE(LineInterfaceElement_SizeOfRightHandSideEqualsTheNumberOfUDofs,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    Model model;
    auto& model_part = model.CreateModelPart("Main");
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT);

    PointerVector<Node> nodes;
    nodes.push_back(model_part.CreateNewNode(0, 0.0, 0.0, 0.0));
    nodes.push_back(model_part.CreateNewNode(1, 1.0, 0.0, 0.0));
    nodes.push_back(model_part.CreateNewNode(2, 0.0, 0.0, 0.0));
    nodes.push_back(model_part.CreateNewNode(3, 1.0, 0.0, 0.0));
    auto geometry   = std::make_shared<LineInterfaceGeometry<Line2D2<Node>>>(nodes);
    auto properties = std::make_shared<Properties>();
    auto element    = make_intrusive<LineInterfaceElement>(1, geometry, properties);

    model_part.AddElement(element);
    for (auto& node : element->GetGeometry()) {
        node.AddDof(DISPLACEMENT_X);
        node.AddDof(DISPLACEMENT_Y);
    }

    // Act
    Vector actual_right_hand_side;
    element->CalculateRightHandSide(actual_right_hand_side, {});

    // Assert
    auto expected_right_hand_side = Vector{ZeroVector{8}};
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(actual_right_hand_side, expected_right_hand_side, Defaults::relative_tolerance)
}

} // namespace Kratos::Testing
