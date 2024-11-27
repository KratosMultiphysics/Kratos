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
using LineInterfaceGeometry2D2Plus2Noded = LineInterfaceGeometry<Line2D2<Node>>;
using LineInterfaceGeometry2D3Plus3Noded = LineInterfaceGeometry<Line2D3<Node>>;

PointerVector<Node> CreateNodes()
{
    PointerVector<Node> result;
    result.push_back(Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0));
    result.push_back(Kratos::make_intrusive<Node>(2, 1.0, 0.0, 0.0));
    result.push_back(Kratos::make_intrusive<Node>(3, 0.0, 0.0, 0.0));
    result.push_back(Kratos::make_intrusive<Node>(4, 1.0, 0.0, 0.0));

    return result;
}

std::shared_ptr<Properties> CreateLinearElasticMaterialProperties(double NormalStiffness, double ShearStiffness)
{
    auto result                                  = std::make_shared<Properties>();
    result->GetValue(INTERFACE_NORMAL_STIFFNESS) = NormalStiffness;
    result->GetValue(INTERFACE_SHEAR_STIFFNESS)  = ShearStiffness;
    result->GetValue(CONSTITUTIVE_LAW) = std::make_shared<GeoIncrementalLinearElasticInterfaceLaw>();

    return result;
}

ModelPart& CreateModelPartWithDisplacementVariable(Model& rModel)
{
    auto& r_result = rModel.CreateModelPart("Main");
    r_result.AddNodalSolutionStepVariable(DISPLACEMENT);

    return r_result;
}

LineInterfaceElement CreateLineInterfaceElementWithUDofs(const Properties::Pointer&     rProperties,
                                                         const Geometry<Node>::Pointer& rGeometry)
{
    auto result = LineInterfaceElement{1, rGeometry, rProperties};
    for (auto& node : result.GetGeometry()) {
        node.AddDof(DISPLACEMENT_X);
        node.AddDof(DISPLACEMENT_Y);
    }

    return result;
}

LineInterfaceElement CreateHorizontalUnitLength2Plus2NodedLineInterfaceElementWithUDofs(Model& rModel,
                                                                                        const Properties::Pointer& rProperties)
{
    auto& r_model_part = CreateModelPartWithDisplacementVariable(rModel);

    PointerVector<Node> nodes;
    nodes.push_back(r_model_part.CreateNewNode(0, 0.0, 0.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(1, 1.0, 0.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(2, 0.0, 0.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(3, 1.0, 0.0, 0.0));
    const auto p_geometry = std::make_shared<LineInterfaceGeometry2D2Plus2Noded>(nodes);
    return CreateLineInterfaceElementWithUDofs(rProperties, p_geometry);
}

LineInterfaceElement CreateHorizontalUnitLength3Plus3NodedLineInterfaceElementWithDisplacementDoF(
    Model& rModel, const Properties::Pointer& rProperties)
{
    auto& r_model_part = CreateModelPartWithDisplacementVariable(rModel);

    PointerVector<Node> nodes;
    nodes.push_back(r_model_part.CreateNewNode(0, 0.0, 0.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(1, 1.0, 0.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(2, 0.5, 0.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(3, 0.0, 0.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(4, 1.0, 0.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(5, 0.5, 0.0, 0.0));
    const auto p_geometry = std::make_shared<LineInterfaceGeometry2D3Plus3Noded>(nodes);
    return CreateLineInterfaceElementWithUDofs(rProperties, p_geometry);
}

LineInterfaceElement CreateUnitLengthLineInterfaceElementRotatedBy30DegreesWithDisplacementDoF(
    Model& rModel, const Properties::Pointer& rProperties)
{
    auto& r_model_part = CreateModelPartWithDisplacementVariable(rModel);

    PointerVector<Node> nodes;
    nodes.push_back(r_model_part.CreateNewNode(0, 0.0, 0.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(1, 0.5 * std::sqrt(3.0), 0.5, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(2, 0.0, 0.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(3, 0.5 * std::sqrt(3.0), 0.5, 0.0));
    const auto p_geometry = std::make_shared<LineInterfaceGeometry2D2Plus2Noded>(nodes);
    return CreateLineInterfaceElementWithUDofs(rProperties, p_geometry);
}

Matrix CreateExpectedStiffnessMatrixForHorizontal2Plus2NodedElement(double NormalStiffness, double ShearStiffness)
{
    auto expected_left_hand_side  = Matrix{ZeroMatrix{8, 8}};
    expected_left_hand_side(0, 0) = ShearStiffness * 0.5;
    expected_left_hand_side(1, 1) = NormalStiffness * 0.5;
    expected_left_hand_side(2, 2) = ShearStiffness * 0.5;
    expected_left_hand_side(3, 3) = NormalStiffness * 0.5;
    expected_left_hand_side(4, 4) = ShearStiffness * 0.5;
    expected_left_hand_side(5, 5) = NormalStiffness * 0.5;
    expected_left_hand_side(6, 6) = ShearStiffness * 0.5;
    expected_left_hand_side(7, 7) = NormalStiffness * 0.5;
    expected_left_hand_side(0, 4) = -ShearStiffness * 0.5;
    expected_left_hand_side(1, 5) = -NormalStiffness * 0.5;
    expected_left_hand_side(2, 6) = -ShearStiffness * 0.5;
    expected_left_hand_side(3, 7) = -NormalStiffness * 0.5;
    expected_left_hand_side(4, 0) = -ShearStiffness * 0.5;
    expected_left_hand_side(5, 1) = -NormalStiffness * 0.5;
    expected_left_hand_side(6, 2) = -ShearStiffness * 0.5;
    expected_left_hand_side(7, 3) = -NormalStiffness * 0.5;

    return expected_left_hand_side;
}

} // namespace

namespace Kratos::Testing
{

using namespace Kratos;

KRATOS_TEST_CASE_IN_SUITE(LineInterfaceElementIsAnElement, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const LineInterfaceElement element;
    auto                       p_casted_element = dynamic_cast<const Element*>(&element);
    KRATOS_CHECK_NOT_EQUAL(p_casted_element, nullptr);
}

KRATOS_TEST_CASE_IN_SUITE(LineInterfaceElementCanCreateInstanceWithGeometryInput, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const LineInterfaceElement element;
    const auto p_geometry   = std::make_shared<LineInterfaceGeometry2D2Plus2Noded>(CreateNodes());
    const auto p_properties = std::make_shared<Properties>();

    // Act
    const auto p_created_element = element.Create(1, p_geometry, p_properties);

    // Assert
    EXPECT_NE(p_created_element, nullptr);
    EXPECT_EQ(p_created_element->Id(), 1);
    EXPECT_NE(p_created_element->pGetGeometry(), nullptr);
    EXPECT_NE(p_created_element->pGetProperties(), nullptr);
}

KRATOS_TEST_CASE_IN_SUITE(LineInterfaceElementCanCreateInstanceWithNodeInput, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto nodes        = CreateNodes();
    const auto p_properties = std::make_shared<Properties>();

    // The source element needs to have a geometry, otherwise the version of the
    // Create method with a node input will fail.
    const auto p_geometry = std::make_shared<LineInterfaceGeometry2D2Plus2Noded>(nodes);
    const LineInterfaceElement element(0, p_geometry, p_properties);

    // Act
    const auto p_created_element = element.Create(1, nodes, p_properties);

    // Assert
    EXPECT_NE(p_created_element, nullptr);
    EXPECT_EQ(p_created_element->Id(), 1);
    EXPECT_NE(p_created_element->pGetGeometry(), nullptr);
    EXPECT_NE(p_created_element->pGetProperties(), nullptr);
}

KRATOS_TEST_CASE_IN_SUITE(LineInterfaceElement_ReturnsTheExpectedDoFList, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto p_properties = std::make_shared<Properties>();

    Model model;
    const auto element = CreateHorizontalUnitLength2Plus2NodedLineInterfaceElementWithUDofs(model, p_properties);

    // Act
    const auto              dummy_process_info = ProcessInfo{};
    Element::DofsVectorType degrees_of_freedom;
    element.GetDofList(degrees_of_freedom, dummy_process_info);

    // Assert
    KRATOS_EXPECT_EQ(degrees_of_freedom.size(), 8);
    for (std::size_t i = 0; i < degrees_of_freedom.size(); i += 2) {
        KRATOS_EXPECT_EQ(degrees_of_freedom[i]->GetVariable(), DISPLACEMENT_X);
        KRATOS_EXPECT_EQ(degrees_of_freedom[i + 1]->GetVariable(), DISPLACEMENT_Y);
    }
}

KRATOS_TEST_CASE_IN_SUITE(LineInterfaceElement_ReturnsTheExpectedEquationIdVector, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto p_properties = std::make_shared<Properties>();

    Model model;
    auto element = CreateHorizontalUnitLength2Plus2NodedLineInterfaceElementWithUDofs(model, p_properties);

    int i = 0;
    for (const auto& node : element.GetGeometry()) {
        ++i;
        node.pGetDof(DISPLACEMENT_X)->SetEquationId(i);

        ++i;
        node.pGetDof(DISPLACEMENT_Y)->SetEquationId(i);
    }

    // Act
    const auto                    dummy_process_info = ProcessInfo{};
    Element::EquationIdVectorType equation_id_vector;
    element.EquationIdVector(equation_id_vector, dummy_process_info);

    // Assert
    const Element::EquationIdVectorType expected_ids = {1, 2, 3, 4, 5, 6, 7, 8};
    KRATOS_EXPECT_VECTOR_EQ(equation_id_vector, expected_ids)
}

KRATOS_TEST_CASE_IN_SUITE(LineInterfaceElement_LeftHandSideContainsMaterialStiffnessContributions,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    constexpr auto normal_stiffness = 20.0;
    constexpr auto shear_stiffness  = 10.0;
    const auto p_properties = CreateLinearElasticMaterialProperties(normal_stiffness, shear_stiffness);

    Model model;
    auto element = CreateHorizontalUnitLength2Plus2NodedLineInterfaceElementWithUDofs(model, p_properties);

    const auto dummy_process_info = ProcessInfo{};
    element.Initialize(dummy_process_info);

    // Act
    Matrix actual_left_hand_side;
    element.CalculateLeftHandSide(actual_left_hand_side, dummy_process_info);

    // Assert
    const auto expected_left_hand_side =
        CreateExpectedStiffnessMatrixForHorizontal2Plus2NodedElement(normal_stiffness, shear_stiffness);
    KRATOS_EXPECT_MATRIX_RELATIVE_NEAR(actual_left_hand_side, expected_left_hand_side, Defaults::relative_tolerance)
}

KRATOS_TEST_CASE_IN_SUITE(LineInterfaceElement_LeftHandSideContainsMaterialStiffnessContributions_Rotated,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    constexpr auto normal_stiffness = 20.0;
    constexpr auto shear_stiffness  = 10.0;
    const auto p_properties = CreateLinearElasticMaterialProperties(normal_stiffness, shear_stiffness);

    Model model;
    auto element = CreateUnitLengthLineInterfaceElementRotatedBy30DegreesWithDisplacementDoF(model, p_properties);

    const auto dummy_process_info = ProcessInfo{};
    element.Initialize(dummy_process_info);

    // Act
    Matrix actual_left_hand_side;
    element.CalculateLeftHandSide(actual_left_hand_side, dummy_process_info);

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
    KRATOS_EXPECT_MATRIX_RELATIVE_NEAR(actual_left_hand_side, expected_left_hand_side, Defaults::relative_tolerance)
}

KRATOS_TEST_CASE_IN_SUITE(LineInterfaceElement_RightHandSideEqualsMinusInternalForceVector,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    constexpr auto normal_stiffness = 20.0;
    constexpr auto shear_stiffness  = 10.0;
    const auto p_properties = CreateLinearElasticMaterialProperties(normal_stiffness, shear_stiffness);

    Model model;
    auto element = CreateHorizontalUnitLength2Plus2NodedLineInterfaceElementWithUDofs(model, p_properties);

    const auto dummy_process_info = ProcessInfo{};
    element.Initialize(dummy_process_info);

    element.GetGeometry()[2].FastGetSolutionStepValue(DISPLACEMENT) = array_1d<double, 3>{0.2, 0.5, 0.0};
    element.GetGeometry()[3].FastGetSolutionStepValue(DISPLACEMENT) = array_1d<double, 3>{0.2, 0.5, 0.0};

    // Act
    Vector actual_right_hand_side;
    element.CalculateRightHandSide(actual_right_hand_side, dummy_process_info);

    // Assert
    auto expected_right_hand_side = Vector{8};
    expected_right_hand_side <<= 1.0, 5.0, 1.0, 5.0, -1.0, -5.0, -1.0, -5.0;
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(actual_right_hand_side, expected_right_hand_side, Defaults::relative_tolerance)
}

KRATOS_TEST_CASE_IN_SUITE(LineInterfaceElement_RightHandSideEqualsMinusInternalForceVector_Rotated,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    constexpr auto normal_stiffness = 20.0;
    constexpr auto shear_stiffness  = 10.0;
    const auto p_properties = CreateLinearElasticMaterialProperties(normal_stiffness, shear_stiffness);

    Model model;
    auto element = CreateUnitLengthLineInterfaceElementRotatedBy30DegreesWithDisplacementDoF(model, p_properties);

    const auto dummy_process_info = ProcessInfo{};
    element.Initialize(dummy_process_info);

    // Rotated the relative normal displacement of 0.5 and the relative shear displacement of 0.2
    element.GetGeometry()[2].FastGetSolutionStepValue(DISPLACEMENT) =
        array_1d<double, 3>{-0.07679492, 0.5330127, 0.0};
    element.GetGeometry()[3].FastGetSolutionStepValue(DISPLACEMENT) =
        array_1d<double, 3>{-0.07679492, 0.5330127, 0.0};

    // Act
    Vector actual_right_hand_side;
    element.CalculateRightHandSide(actual_right_hand_side, dummy_process_info);

    // Assert
    auto expected_right_hand_side = Vector{8};
    // clang-format off
    expected_right_hand_side <<= -1.6339746,  4.83012702, -1.6339746,  4.83012702,
                                  1.6339746, -4.83012702,  1.6339746, -4.83012702;
    // clang-format on
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(actual_right_hand_side, expected_right_hand_side, Defaults::relative_tolerance)
}

KRATOS_TEST_CASE_IN_SUITE(GetInitializedConstitutiveLawsAfterElementInitialization, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto p_properties = std::make_shared<Properties>();
    p_properties->GetValue(CONSTITUTIVE_LAW) = std::make_shared<GeoIncrementalLinearElasticInterfaceLaw>();

    Model model;
    auto element = CreateHorizontalUnitLength2Plus2NodedLineInterfaceElementWithUDofs(model, p_properties);

    const auto dummy_process_info = ProcessInfo{};
    element.Initialize(dummy_process_info);

    // Act
    auto constitutive_laws = std::vector<ConstitutiveLaw::Pointer>{};
    element.CalculateOnIntegrationPoints(CONSTITUTIVE_LAW, constitutive_laws, dummy_process_info);

    // Assert
    KRATOS_EXPECT_EQ(constitutive_laws.size(), 2);
    for (const auto& p_law : constitutive_laws) {
        KRATOS_EXPECT_NE(p_law.get(), nullptr);
        // Constitutive laws must have been cloned
        KRATOS_EXPECT_NE(p_law.get(), p_properties->GetValue(CONSTITUTIVE_LAW).get());
    }
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceElement_HasCorrectNumberOfConstitutiveLawsAfterMultipleInitializations,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto p_properties = std::make_shared<Properties>();
    p_properties->GetValue(CONSTITUTIVE_LAW) = std::make_shared<GeoIncrementalLinearElasticInterfaceLaw>();

    Model model;
    auto element = CreateHorizontalUnitLength2Plus2NodedLineInterfaceElementWithUDofs(model, p_properties);

    const auto dummy_process_info = ProcessInfo{};

    // Multiple initializations emulate a multi-stage simulation
    element.Initialize(dummy_process_info);
    element.Initialize(dummy_process_info);
    element.Initialize(dummy_process_info);

    // Assert
    auto constitutive_laws = std::vector<ConstitutiveLaw::Pointer>{};
    element.CalculateOnIntegrationPoints(CONSTITUTIVE_LAW, constitutive_laws, dummy_process_info);
    KRATOS_EXPECT_EQ(constitutive_laws.size(), 2);
}

KRATOS_TEST_CASE_IN_SUITE(LineInterfaceElement_CalculateLocalSystem_ReturnsExpectedLeftAndRightHandSide,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    constexpr auto normal_stiffness = 20.0;
    constexpr auto shear_stiffness  = 10.0;
    const auto p_properties = CreateLinearElasticMaterialProperties(normal_stiffness, shear_stiffness);

    Model model;
    auto element = CreateHorizontalUnitLength2Plus2NodedLineInterfaceElementWithUDofs(model, p_properties);

    const auto dummy_process_info = ProcessInfo{};
    element.Initialize(dummy_process_info);

    element.GetGeometry()[2].FastGetSolutionStepValue(DISPLACEMENT) = array_1d<double, 3>{0.2, 0.5, 0.0};
    element.GetGeometry()[3].FastGetSolutionStepValue(DISPLACEMENT) = array_1d<double, 3>{0.2, 0.5, 0.0};

    // Act
    Vector actual_right_hand_side;
    Matrix actual_left_hand_side;
    element.CalculateLocalSystem(actual_left_hand_side, actual_right_hand_side, dummy_process_info);

    // Assert
    const auto expected_left_hand_side =
        CreateExpectedStiffnessMatrixForHorizontal2Plus2NodedElement(normal_stiffness, shear_stiffness);
    KRATOS_EXPECT_MATRIX_RELATIVE_NEAR(actual_left_hand_side, expected_left_hand_side, Defaults::relative_tolerance)

    auto expected_right_hand_side = Vector{8};
    expected_right_hand_side <<= 1.0, 5.0, 1.0, 5.0, -1.0, -5.0, -1.0, -5.0;
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(actual_right_hand_side, expected_right_hand_side, Defaults::relative_tolerance)
}

KRATOS_TEST_CASE_IN_SUITE(LineInterfaceElement_CalculateStrain_ReturnsRelativeDisplacement,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    constexpr auto normal_stiffness = 20.0;
    constexpr auto shear_stiffness  = 10.0;
    const auto p_properties = CreateLinearElasticMaterialProperties(normal_stiffness, shear_stiffness);

    Model model;
    auto element = CreateUnitLengthLineInterfaceElementRotatedBy30DegreesWithDisplacementDoF(model, p_properties);

    const auto dummy_process_info = ProcessInfo{};
    element.Initialize(dummy_process_info);

    // Rotated the relative normal displacement of 0.5 and the relative shear displacement of 0.2
    element.GetGeometry()[2].FastGetSolutionStepValue(DISPLACEMENT) =
        array_1d<double, 3>{-0.07679492, 0.5330127, 0.0};
    element.GetGeometry()[3].FastGetSolutionStepValue(DISPLACEMENT) =
        array_1d<double, 3>{-0.07679492, 0.5330127, 0.0};

    // Act
    std::vector<Vector> relative_displacements_at_integration_points;
    element.CalculateOnIntegrationPoints(STRAIN, relative_displacements_at_integration_points, dummy_process_info);

    // Assert
    Vector expected_relative_displacement{2};
    expected_relative_displacement <<= 0.5, 0.2;
    KRATOS_EXPECT_EQ(relative_displacements_at_integration_points.size(), 2);
    for (const auto& r_relative_displacement : relative_displacements_at_integration_points) {
        KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(r_relative_displacement, expected_relative_displacement,
                                           Defaults::relative_tolerance)
    }
}

KRATOS_TEST_CASE_IN_SUITE(LineInterfaceElement_CalculateCauchyStressVector_ReturnsTraction,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    constexpr auto normal_stiffness = 20.0;
    constexpr auto shear_stiffness  = 10.0;
    const auto p_properties = CreateLinearElasticMaterialProperties(normal_stiffness, shear_stiffness);

    Model model;
    auto element = CreateUnitLengthLineInterfaceElementRotatedBy30DegreesWithDisplacementDoF(model, p_properties);

    const auto dummy_process_info = ProcessInfo{};
    element.Initialize(dummy_process_info);

    // Rotated the relative normal displacement of 0.5 and the relative shear displacement of 0.2
    element.GetGeometry()[2].FastGetSolutionStepValue(DISPLACEMENT) =
        array_1d<double, 3>{-0.07679492, 0.5330127, 0.0};
    element.GetGeometry()[3].FastGetSolutionStepValue(DISPLACEMENT) =
        array_1d<double, 3>{-0.07679492, 0.5330127, 0.0};

    // Act
    std::vector<Vector> tractions_at_integration_points;
    element.CalculateOnIntegrationPoints(CAUCHY_STRESS_VECTOR, tractions_at_integration_points, dummy_process_info);

    // Assert
    Vector expected_traction{2};
    expected_traction <<= 10.0, 2.0;
    KRATOS_EXPECT_EQ(tractions_at_integration_points.size(), 2);
    for (const auto& r_traction : tractions_at_integration_points) {
        KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(r_traction, expected_traction, Defaults::relative_tolerance)
    }
}

KRATOS_TEST_CASE_IN_SUITE(3Plus3NodedLineInterfaceElement_CalculateLocalSystem_ReturnsExpectedLeftAndRightHandSide,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    constexpr auto normal_stiffness = 20.0;
    constexpr auto shear_stiffness  = 10.0;
    const auto p_properties = CreateLinearElasticMaterialProperties(normal_stiffness, shear_stiffness);

    Model model;
    auto element = CreateHorizontalUnitLength3Plus3NodedLineInterfaceElementWithDisplacementDoF(model, p_properties);

    const auto dummy_process_info = ProcessInfo{};
    element.Initialize(dummy_process_info);

    element.GetGeometry()[3].FastGetSolutionStepValue(DISPLACEMENT) = array_1d<double, 3>{0.2, 0.5, 0.0};
    element.GetGeometry()[4].FastGetSolutionStepValue(DISPLACEMENT) = array_1d<double, 3>{0.2, 0.5, 0.0};
    element.GetGeometry()[5].FastGetSolutionStepValue(DISPLACEMENT) = array_1d<double, 3>{0.2, 0.5, 0.0};

    // Act
    Vector actual_right_hand_side;
    Matrix actual_left_hand_side;
    element.CalculateLocalSystem(actual_left_hand_side, actual_right_hand_side, dummy_process_info);

    // Assert
    auto expected_left_hand_side    = Matrix{ZeroMatrix{12, 12}};
    expected_left_hand_side(0, 0)   = shear_stiffness * (1.0 / 6.0);
    expected_left_hand_side(1, 1)   = normal_stiffness * (1.0 / 6.0);
    expected_left_hand_side(2, 2)   = shear_stiffness * (1.0 / 6.0);
    expected_left_hand_side(3, 3)   = normal_stiffness * (1.0 / 6.0);
    expected_left_hand_side(4, 4)   = shear_stiffness * (2.0 / 3.0);
    expected_left_hand_side(5, 5)   = normal_stiffness * (2.0 / 3.0);
    expected_left_hand_side(6, 6)   = shear_stiffness * (1.0 / 6.0);
    expected_left_hand_side(7, 7)   = normal_stiffness * (1.0 / 6.0);
    expected_left_hand_side(8, 8)   = shear_stiffness * (1.0 / 6.0);
    expected_left_hand_side(9, 9)   = normal_stiffness * (1.0 / 6.0);
    expected_left_hand_side(10, 10) = shear_stiffness * (2.0 / 3.0);
    expected_left_hand_side(11, 11) = normal_stiffness * (2.0 / 3.0);

    expected_left_hand_side(0, 6)  = -shear_stiffness * (1.0 / 6.0);
    expected_left_hand_side(1, 7)  = -normal_stiffness * (1.0 / 6.0);
    expected_left_hand_side(2, 8)  = -shear_stiffness * (1.0 / 6.0);
    expected_left_hand_side(3, 9)  = -normal_stiffness * (1.0 / 6.0);
    expected_left_hand_side(4, 10) = -shear_stiffness * (2.0 / 3.0);
    expected_left_hand_side(5, 11) = -normal_stiffness * (2.0 / 3.0);

    expected_left_hand_side(6, 0)  = -shear_stiffness * (1.0 / 6.0);
    expected_left_hand_side(7, 1)  = -normal_stiffness * (1.0 / 6.0);
    expected_left_hand_side(8, 2)  = -shear_stiffness * (1.0 / 6.0);
    expected_left_hand_side(9, 3)  = -normal_stiffness * (1.0 / 6.0);
    expected_left_hand_side(10, 4) = -shear_stiffness * (2.0 / 3.0);
    expected_left_hand_side(11, 5) = -normal_stiffness * (2.0 / 3.0);

    KRATOS_EXPECT_MATRIX_RELATIVE_NEAR(actual_left_hand_side, expected_left_hand_side, Defaults::relative_tolerance)

    auto expected_right_hand_side = Vector{ZeroVector{12}};
    expected_right_hand_side <<= 0.33333333, 1.66666667, 0.33333333, 1.66666667, 1.33333333,
        6.66666667, -0.33333333, -1.66666667, -0.33333333, -1.66666667, -1.33333333, -6.66666667;
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(actual_right_hand_side, expected_right_hand_side, Defaults::relative_tolerance)
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceElement_CheckThrowsWhenElementIsNotInitialized, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    constexpr auto normal_stiffness = 20.0;
    constexpr auto shear_stiffness  = 10.0;
    const auto p_properties = CreateLinearElasticMaterialProperties(normal_stiffness, shear_stiffness);

    Model model;
    auto element = CreateHorizontalUnitLength3Plus3NodedLineInterfaceElementWithDisplacementDoF(model, p_properties);

    const auto dummy_process_info = ProcessInfo{};
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        element.Check(dummy_process_info),
        "Number of integration points (3) and constitutive laws (0) do not match.\n")

    element.Initialize(dummy_process_info);

    KRATOS_EXPECT_EQ(element.Check(dummy_process_info), 0);
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceElement_CheckDoesNotThrowWhenElementIsNotActive, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    constexpr auto normal_stiffness = 20.0;
    constexpr auto shear_stiffness  = 10.0;
    const auto p_properties = CreateLinearElasticMaterialProperties(normal_stiffness, shear_stiffness);

    Model model;
    auto element = CreateHorizontalUnitLength3Plus3NodedLineInterfaceElementWithDisplacementDoF(model, p_properties);

    // In the integrated workflow, the elements are not initialized, when they are not active.
    // However, the Check method is always called on all elements, even if they are not active.
    // Therefore, the Check method should not throw an exception in this case, even though the
    // constitutive laws are not initialized.
    element.Set(ACTIVE, false);

    const auto dummy_process_info = ProcessInfo{};
    KRATOS_EXPECT_EQ(element.Check(dummy_process_info), 0);
}

} // namespace Kratos::Testing
