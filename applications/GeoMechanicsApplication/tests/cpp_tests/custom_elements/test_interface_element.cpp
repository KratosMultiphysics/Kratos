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
//                   Gennady Markelov
//

#include "custom_constitutive/incremental_linear_elastic_interface_law.h"
#include "custom_constitutive/interface_plane_strain.h"
#include "custom_constitutive/interface_three_dimensional_surface.h"
#include "custom_elements/interface_element.h"
#include "custom_elements/interface_stress_state.h"
#include "custom_geometries/interface_geometry.h"
#include "geo_mechanics_application_variables.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "tests/cpp_tests/test_utilities.h"

#include <boost/numeric/ublas/assignment.hpp>
#include <cstddef>

namespace
{

using namespace Kratos;
using LineInterfaceGeometry2D2Plus2Noded     = InterfaceGeometry<Line2D2<Node>>;
using LineInterfaceGeometry2D3Plus3Noded     = InterfaceGeometry<Line2D3<Node>>;
using TriangleInterfaceGeometry3D3Plus3Noded = InterfaceGeometry<Triangle3D3<Node>>;
using TriangleInterfaceGeometry3D6Plus6Noded = InterfaceGeometry<Triangle3D6<Node>>;
using Interface2D                            = Line2DInterfaceStressState;
using Interface3D                            = SurfaceInterfaceStressState;
using PrescribedDisplacements = std::vector<std::pair<std::size_t, array_1d<double, 3>>>;

PointerVector<Node> CreateNodesFor2Plus2LineInterfaceGeometry()
{
    PointerVector<Node> result;
    result.push_back(make_intrusive<Node>(1, 0.0, 0.0, 0.0));
    result.push_back(make_intrusive<Node>(2, 1.0, 0.0, 0.0));
    result.push_back(make_intrusive<Node>(3, 0.0, 0.0, 0.0));
    result.push_back(make_intrusive<Node>(4, 1.0, 0.0, 0.0));

    return result;
}

PointerVector<Node> CreateNodesFor3Plus3SurfaceInterfaceGeometry()
{
    PointerVector<Node> result;
    result.push_back(make_intrusive<Node>(1, 0.0, 0.0, 0.0));
    result.push_back(make_intrusive<Node>(2, 1.0, 0.0, 0.0));
    result.push_back(make_intrusive<Node>(3, 1.0, 1.0, 0.0));
    result.push_back(make_intrusive<Node>(4, 0.0, 0.0, 0.0));
    result.push_back(make_intrusive<Node>(5, 1.0, 0.0, 0.0));
    result.push_back(make_intrusive<Node>(6, 1.0, 1.0, 0.0));

    return result;
}

template <typename TConstitutiveLawDimension>
std::shared_ptr<Properties> CreateElasticMaterialProperties(double NormalStiffness, double ShearStiffness)
{
    auto result                                  = std::make_shared<Properties>();
    result->GetValue(INTERFACE_NORMAL_STIFFNESS) = NormalStiffness;
    result->GetValue(INTERFACE_SHEAR_STIFFNESS)  = ShearStiffness;
    result->GetValue(CONSTITUTIVE_LAW) = std::make_shared<GeoIncrementalLinearElasticInterfaceLaw>(
        std::make_unique<TConstitutiveLawDimension>());

    return result;
}

ModelPart& CreateModelPartWithDisplacementVariable(Model& rModel)
{
    auto& r_result = rModel.CreateModelPart("Main");
    r_result.AddNodalSolutionStepVariable(DISPLACEMENT);

    return r_result;
}

template <typename TInterfaceDimension>
InterfaceElement CreateInterfaceElementWithUDofs(const Properties::Pointer&     rProperties,
                                                 const Geometry<Node>::Pointer& rGeometry)
{
    auto result = InterfaceElement{1, rGeometry, rProperties, std::make_unique<TInterfaceDimension>()};
    for (auto& node : result.GetGeometry()) {
        node.AddDof(DISPLACEMENT_X);
        node.AddDof(DISPLACEMENT_Y);
        node.AddDof(DISPLACEMENT_Z);
    }

    return result;
}

InterfaceElement CreateHorizontalUnitLength2Plus2NodedLineInterfaceElementWithUDofs(Model& rModel,
                                                                                    const Properties::Pointer& rProperties)
{
    auto& r_model_part = CreateModelPartWithDisplacementVariable(rModel);

    PointerVector<Node> nodes;
    nodes.push_back(r_model_part.CreateNewNode(0, 0.0, 0.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(1, 1.0, 0.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(2, 0.0, 0.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(3, 1.0, 0.0, 0.0));
    const auto p_geometry = std::make_shared<LineInterfaceGeometry2D2Plus2Noded>(nodes);
    return CreateInterfaceElementWithUDofs<Interface2D>(rProperties, p_geometry);
}

InterfaceElement CreateHorizontalUnitLength3Plus3NodedLineInterfaceElementWithDisplacementDoF(
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
    return CreateInterfaceElementWithUDofs<Interface2D>(rProperties, p_geometry);
}

InterfaceElement CreateUnitLengthLineInterfaceElementRotatedBy30DegreesWithDisplacementDoF(
    Model& rModel, const Properties::Pointer& rProperties)
{
    auto& r_model_part = CreateModelPartWithDisplacementVariable(rModel);

    PointerVector<Node> nodes;
    nodes.push_back(r_model_part.CreateNewNode(0, 0.0, 0.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(1, 0.5 * std::sqrt(3.0), 0.5, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(2, 0.0, 0.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(3, 0.5 * std::sqrt(3.0), 0.5, 0.0));
    const auto p_geometry = std::make_shared<LineInterfaceGeometry2D2Plus2Noded>(nodes);
    return CreateInterfaceElementWithUDofs<Interface2D>(rProperties, p_geometry);
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

InterfaceElement CreateHorizontal3Plus3NodedTriangleInterfaceElementWithUDofs(Model& rModel,
                                                                              const Properties::Pointer& rProperties)
{
    auto& r_model_part = CreateModelPartWithDisplacementVariable(rModel);

    PointerVector<Node> nodes;
    nodes.push_back(r_model_part.CreateNewNode(0, 0.0, 0.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(1, 1.0, 0.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(2, 0.0, 1.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(3, 0.0, 0.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(4, 1.0, 0.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(5, 0.0, 1.0, 0.0));
    const auto p_geometry = std::make_shared<TriangleInterfaceGeometry3D3Plus3Noded>(nodes);
    return CreateInterfaceElementWithUDofs<Interface3D>(rProperties, p_geometry);
}

InterfaceElement CreateHorizontal3Plus3NodedTriangleInterfaceYZPlaneElementWithUDofs(Model& rModel,
                                                                                     const Properties::Pointer& rProperties)
{
    auto& r_model_part = CreateModelPartWithDisplacementVariable(rModel);

    PointerVector<Node> nodes;
    nodes.push_back(r_model_part.CreateNewNode(0, 0.0, 0.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(1, 0.0, 0.0, -1.0));
    nodes.push_back(r_model_part.CreateNewNode(2, 0.0, 1.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(3, 0.0, 0.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(4, 0.0, 0.0, -1.0));
    nodes.push_back(r_model_part.CreateNewNode(5, 0.0, 1.0, 0.0));
    const auto p_geometry = std::make_shared<TriangleInterfaceGeometry3D3Plus3Noded>(nodes);
    return CreateInterfaceElementWithUDofs<Interface3D>(rProperties, p_geometry);
}

InterfaceElement CreateHorizontal3Plus3NodedTriangleInterfaceXZPlaneElementWithUDofs(Model& rModel,
                                                                                     const Properties::Pointer& rProperties)
{
    auto& r_model_part = CreateModelPartWithDisplacementVariable(rModel);

    PointerVector<Node> nodes;
    nodes.push_back(r_model_part.CreateNewNode(0, 0.0, 0.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(1, 1.0, 0.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(2, 0.0, 0.0, 1.0));
    nodes.push_back(r_model_part.CreateNewNode(3, 0.0, 0.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(4, 1.0, 0.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(5, 0.0, 0.0, 1.0));
    const auto p_geometry = std::make_shared<TriangleInterfaceGeometry3D3Plus3Noded>(nodes);
    return CreateInterfaceElementWithUDofs<Interface3D>(rProperties, p_geometry);
}

InterfaceElement CreateHorizontal6Plus6NodedTriangleInterfaceElementWithDisplacementDoF(Model& rModel,
                                                                                        const Properties::Pointer& rProperties)
{
    auto& r_model_part = CreateModelPartWithDisplacementVariable(rModel);

    PointerVector<Node> nodes;
    nodes.push_back(r_model_part.CreateNewNode(0, 0.0, 0.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(1, 1.0, 0.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(2, 1.0, 1.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(3, 0.5, 0.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(4, 1.0, 0.5, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(5, 0.5, 0.5, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(6, 0.0, 0.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(7, 1.0, 0.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(8, 1.0, 1.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(9, 0.5, 0.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(10, 1.0, 0.5, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(11, 0.5, 0.5, 0.0));
    const auto p_geometry = std::make_shared<TriangleInterfaceGeometry3D6Plus6Noded>(nodes);
    return CreateInterfaceElementWithUDofs<Interface3D>(rProperties, p_geometry);
}

InterfaceElement CreateTriangleInterfaceElementRotatedBy30DegreesWithDisplacementDoF(Model& rModel,
                                                                                     const Properties::Pointer& rProperties)
{
    auto& r_model_part = CreateModelPartWithDisplacementVariable(rModel);

    PointerVector<Node> nodes;
    nodes.push_back(r_model_part.CreateNewNode(0, 0.0, 0.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(1, 0.5 * std::sqrt(3.0), 0.5, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(2, -0.5, 0.5 * std::sqrt(3.0), 0.0));
    nodes.push_back(r_model_part.CreateNewNode(3, 0.0, 0.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(4, 0.5 * std::sqrt(3.0), 0.5, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(5, -0.5, 0.5 * std::sqrt(3.0), 0.0));
    const auto p_geometry = std::make_shared<TriangleInterfaceGeometry3D3Plus3Noded>(nodes);
    return CreateInterfaceElementWithUDofs<Interface3D>(rProperties, p_geometry);
}

InterfaceElement CreateTriangleInterfaceElementRotatedBy30DegreesAboutYAxisWithDisplacementDoF(
    Model& rModel, const Properties::Pointer& rProperties)
{
    auto& r_model_part = CreateModelPartWithDisplacementVariable(rModel);

    PointerVector<Node> nodes;
    nodes.push_back(r_model_part.CreateNewNode(0, 0.0, 0.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(1, 0.5 * std::sqrt(3.0), 0.0, -0.5));
    nodes.push_back(r_model_part.CreateNewNode(2, 0.0, 1.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(3, 0.0, 0.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(4, 0.5 * std::sqrt(3.0), 0.0, -0.5));
    nodes.push_back(r_model_part.CreateNewNode(5, 0.0, 1.0, 0.0));
    const auto p_geometry = std::make_shared<TriangleInterfaceGeometry3D3Plus3Noded>(nodes);
    return CreateInterfaceElementWithUDofs<Interface3D>(rProperties, p_geometry);
}

template <typename TElementFactory>
InterfaceElement CreateAndInitializeElement(TElementFactory                Factory,
                                            const Properties::Pointer&     rProperties,
                                            const PrescribedDisplacements& rDisplacements = {})
{
    Model model;
    auto  element = Factory(model, rProperties);
    element.Initialize(ProcessInfo{});
    for (const auto& [idx, disp] : rDisplacements) {
        element.GetGeometry()[idx].FastGetSolutionStepValue(DISPLACEMENT) = disp;
    }
    return element;
}

Matrix ExpectedLeftHandSideForTriangleElement()
{
    // The function provides values taken from the element
    auto expected_left_hand_side = Matrix(18, 18);
    // clang-format off
    expected_left_hand_side <<=
        1.6666666666666661,0,0,0,0,0,0,0,0,-1.6666666666666661,0,0,0,0,0,0,0,0,
        0,1.6666666666666661,0,0,0,0,0,0,0,0,-1.6666666666666661,0,0,0,0,0,0,0,
        0,0,3.3333333333333321,0,0,0,0,0,0,0,0,-3.3333333333333321,0,0,0,0,0,0,
        0,0,0,1.6666666666666661,0,0,0,0,0,0,0,0,-1.6666666666666661,0,0,0,0,0,
        0,0,0,0,1.6666666666666661,0,0,0,0,0,0,0,0,-1.6666666666666661,0,0,0,0,
        0,0,0,0,0,3.3333333333333321,0,0,0,0,0,0,0,0,-3.3333333333333321,0,0,0,
        0,0,0,0,0,0,1.6666666666666661,0,0,0,0,0,0,0,0,-1.6666666666666661,0,0,
        0,0,0,0,0,0,0,1.6666666666666661,0,0,0,0,0,0,0,0,-1.6666666666666661,0,
        0,0,0,0,0,0,0,0,3.3333333333333321,0,0,0,0,0,0,0,0,-3.3333333333333321,
        -1.6666666666666661,0,0,0,0,0,0,0,0,1.6666666666666661,0,0,0,0,0,0,0,0,
        0,-1.6666666666666661,0,0,0,0,0,0,0,0,1.6666666666666661,0,0,0,0,0,0,0,
        0,0,-3.3333333333333321,0,0,0,0,0,0,0,0,3.3333333333333321,0,0,0,0,0,0,
        0,0,0,-1.6666666666666661,0,0,0,0,0,0,0,0,1.6666666666666661,0,0,0,0,0,
        0,0,0,0,-1.6666666666666661,0,0,0,0,0,0,0,0,1.6666666666666661,0,0,0,0,
        0,0,0,0,0,-3.3333333333333321,0,0,0,0,0,0,0,0,3.3333333333333321,0,0,0,
        0,0,0,0,0,0,-1.6666666666666661,0,0,0,0,0,0,0,0,1.6666666666666661,0,0,
        0,0,0,0,0,0,0,-1.6666666666666661,0,0,0,0,0,0,0,0,1.6666666666666661,0,
        0,0,0,0,0,0,0,0,-3.3333333333333321,0,0,0,0,0,0,0,0,3.3333333333333321;
    // clang-format on

    return expected_left_hand_side;
}
} // namespace

namespace Kratos::Testing
{

using namespace Kratos;

KRATOS_TEST_CASE_IN_SUITE(InterfaceElement_IsAnElement, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const InterfaceElement element(
        0, std::make_shared<LineInterfaceGeometry2D2Plus2Noded>(CreateNodesFor2Plus2LineInterfaceGeometry()),
        std::make_unique<Line2DInterfaceStressState>());
    auto p_casted_element = dynamic_cast<const Element*>(&element);
    EXPECT_NE(p_casted_element, nullptr);
}

KRATOS_TEST_CASE_IN_SUITE(LineInterfaceElement_CreatesInstanceWithGeometryInput, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const InterfaceElement element(
        0, std::make_shared<LineInterfaceGeometry2D2Plus2Noded>(CreateNodesFor2Plus2LineInterfaceGeometry()),
        std::make_unique<Line2DInterfaceStressState>());
    const auto p_geometry =
        std::make_shared<LineInterfaceGeometry2D2Plus2Noded>(CreateNodesFor2Plus2LineInterfaceGeometry());
    const auto p_properties = std::make_shared<Properties>();

    // Act
    const auto p_created_element = element.Create(1, p_geometry, p_properties);

    // Assert
    EXPECT_NE(p_created_element, nullptr);
    EXPECT_EQ(p_created_element->Id(), 1);
    EXPECT_NE(p_created_element->pGetGeometry(), nullptr);
    EXPECT_NE(p_created_element->pGetProperties(), nullptr);
}

KRATOS_TEST_CASE_IN_SUITE(LineInterfaceElement_CreatesInstanceWithNodeInput, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto nodes        = CreateNodesFor2Plus2LineInterfaceGeometry();
    const auto p_properties = std::make_shared<Properties>();

    // The source element needs to have a geometry, otherwise the version of the
    // Create method with a node input will fail.
    const auto             p_geometry = std::make_shared<LineInterfaceGeometry2D2Plus2Noded>(nodes);
    const InterfaceElement element(0, p_geometry, p_properties,
                                   std::make_unique<Line2DInterfaceStressState>());

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
    const auto     p_properties =
        CreateElasticMaterialProperties<InterfacePlaneStrain>(normal_stiffness, shear_stiffness);
    auto element = CreateAndInitializeElement(
        CreateHorizontalUnitLength2Plus2NodedLineInterfaceElementWithUDofs, p_properties);

    // Act
    Matrix actual_left_hand_side;
    element.CalculateLeftHandSide(actual_left_hand_side, ProcessInfo{});

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
    const auto     p_properties =
        CreateElasticMaterialProperties<InterfacePlaneStrain>(normal_stiffness, shear_stiffness);

    auto element = CreateAndInitializeElement(
        CreateUnitLengthLineInterfaceElementRotatedBy30DegreesWithDisplacementDoF, p_properties);

    // Act
    Matrix actual_left_hand_side;
    element.CalculateLeftHandSide(actual_left_hand_side, ProcessInfo{});

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
    KRATOS_EXPECT_MATRIX_NEAR(actual_left_hand_side, expected_left_hand_side, Defaults::relative_tolerance)
}

KRATOS_TEST_CASE_IN_SUITE(LineInterfaceElement_RightHandSideEqualsMinusInternalForceVector,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    constexpr auto normal_stiffness = 20.0;
    constexpr auto shear_stiffness  = 10.0;
    const auto     p_properties =
        CreateElasticMaterialProperties<InterfacePlaneStrain>(normal_stiffness, shear_stiffness);

    const auto prescribed_displacements = PrescribedDisplacements{
        {2, array_1d<double, 3>{0.2, 0.5, 0.0}}, {3, array_1d<double, 3>{0.2, 0.5, 0.0}}};
    auto element = CreateAndInitializeElement(CreateHorizontalUnitLength2Plus2NodedLineInterfaceElementWithUDofs,
                                              p_properties, prescribed_displacements);

    // Act
    Vector actual_right_hand_side;
    element.CalculateRightHandSide(actual_right_hand_side, ProcessInfo{});

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
    const auto     p_properties =
        CreateElasticMaterialProperties<InterfacePlaneStrain>(normal_stiffness, shear_stiffness);

    const auto prescribed_displacements =
        PrescribedDisplacements{{2, array_1d<double, 3>{-0.07679492, 0.5330127, 0.0}},
                                {3, array_1d<double, 3>{-0.07679492, 0.5330127, 0.0}}};
    auto element = CreateAndInitializeElement(CreateUnitLengthLineInterfaceElementRotatedBy30DegreesWithDisplacementDoF,
                                              p_properties, prescribed_displacements);

    // Act
    Vector actual_right_hand_side;
    element.CalculateRightHandSide(actual_right_hand_side, ProcessInfo{});

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
    p_properties->GetValue(CONSTITUTIVE_LAW) =
        std::make_shared<GeoIncrementalLinearElasticInterfaceLaw>(std::make_unique<InterfacePlaneStrain>());

    auto element = CreateAndInitializeElement(
        CreateHorizontalUnitLength2Plus2NodedLineInterfaceElementWithUDofs, p_properties);
    // Act
    auto constitutive_laws = std::vector<ConstitutiveLaw::Pointer>{};
    element.CalculateOnIntegrationPoints(CONSTITUTIVE_LAW, constitutive_laws, ProcessInfo{});

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
    p_properties->GetValue(CONSTITUTIVE_LAW) =
        std::make_shared<GeoIncrementalLinearElasticInterfaceLaw>(std::make_unique<InterfacePlaneStrain>());

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

KRATOS_TEST_CASE_IN_SUITE(LineInterfaceElement_ReturnsExpectedLeftAndRightHandSide, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    constexpr auto normal_stiffness = 20.0;
    constexpr auto shear_stiffness  = 10.0;
    const auto     p_properties =
        CreateElasticMaterialProperties<InterfacePlaneStrain>(normal_stiffness, shear_stiffness);

    const auto prescribed_displacements = PrescribedDisplacements{
        {2, array_1d<double, 3>{0.2, 0.5, 0.0}}, {3, array_1d<double, 3>{0.2, 0.5, 0.0}}};
    auto element = CreateAndInitializeElement(CreateHorizontalUnitLength2Plus2NodedLineInterfaceElementWithUDofs,
                                              p_properties, prescribed_displacements);

    // Act
    Vector actual_right_hand_side;
    Matrix actual_left_hand_side;
    element.CalculateLocalSystem(actual_left_hand_side, actual_right_hand_side, ProcessInfo{});

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
    const auto     p_properties =
        CreateElasticMaterialProperties<InterfacePlaneStrain>(normal_stiffness, shear_stiffness);

    const auto prescribed_displacements =
        PrescribedDisplacements{{2, array_1d<double, 3>{-0.07679492, 0.5330127, 0.0}},
                                {3, array_1d<double, 3>{-0.07679492, 0.5330127, 0.0}}};
    auto element = CreateAndInitializeElement(CreateUnitLengthLineInterfaceElementRotatedBy30DegreesWithDisplacementDoF,
                                              p_properties, prescribed_displacements);

    // Act
    std::vector<Vector> relative_displacements_at_integration_points;
    element.CalculateOnIntegrationPoints(STRAIN, relative_displacements_at_integration_points, ProcessInfo{});

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
    const auto     p_properties =
        CreateElasticMaterialProperties<InterfacePlaneStrain>(normal_stiffness, shear_stiffness);

    const auto prescribed_displacements =
        PrescribedDisplacements{{2, array_1d<double, 3>{-0.07679492, 0.5330127, 0.0}},
                                {3, array_1d<double, 3>{-0.07679492, 0.5330127, 0.0}}};
    auto element = CreateAndInitializeElement(CreateUnitLengthLineInterfaceElementRotatedBy30DegreesWithDisplacementDoF,
                                              p_properties, prescribed_displacements);

    // Act
    std::vector<Vector> tractions_at_integration_points;
    element.CalculateOnIntegrationPoints(CAUCHY_STRESS_VECTOR, tractions_at_integration_points, ProcessInfo{});

    // Assert
    Vector expected_traction{2};
    expected_traction <<= 10.0, 2.0;
    KRATOS_EXPECT_EQ(tractions_at_integration_points.size(), 2);
    for (const auto& r_traction : tractions_at_integration_points) {
        KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(r_traction, expected_traction, Defaults::relative_tolerance)
    }
}

KRATOS_TEST_CASE_IN_SUITE(LineInterfaceElement_3Plus3NodedElement_ReturnsExpectedLeftAndRightHandSide,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    constexpr auto normal_stiffness = 20.0;
    constexpr auto shear_stiffness  = 10.0;
    const auto     p_properties =
        CreateElasticMaterialProperties<InterfacePlaneStrain>(normal_stiffness, shear_stiffness);

    const auto prescribed_displacements =
        PrescribedDisplacements{{3, array_1d<double, 3>{0.2, 0.5, 0.0}},
                                {4, array_1d<double, 3>{0.2, 0.5, 0.0}},
                                {5, array_1d<double, 3>{0.2, 0.5, 0.0}}};
    auto element = CreateAndInitializeElement(CreateHorizontalUnitLength3Plus3NodedLineInterfaceElementWithDisplacementDoF,
                                              p_properties, prescribed_displacements);

    // Act
    Vector actual_right_hand_side;
    Matrix actual_left_hand_side;
    element.CalculateLocalSystem(actual_left_hand_side, actual_right_hand_side, ProcessInfo{});

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

KRATOS_TEST_CASE_IN_SUITE(TriangleInterfaceElement_CreatesInstanceWithGeometryInput, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const InterfaceElement element(
        0, std::make_shared<TriangleInterfaceGeometry3D3Plus3Noded>(CreateNodesFor3Plus3SurfaceInterfaceGeometry()),
        std::make_unique<SurfaceInterfaceStressState>());
    const auto p_geometry = std::make_shared<TriangleInterfaceGeometry3D3Plus3Noded>(
        CreateNodesFor3Plus3SurfaceInterfaceGeometry());
    const auto p_properties = std::make_shared<Properties>();

    // Act
    const auto p_created_element = element.Create(1, p_geometry, p_properties);

    // Assert
    EXPECT_NE(p_created_element, nullptr);
    EXPECT_EQ(p_created_element->Id(), 1);
    EXPECT_NE(p_created_element->pGetGeometry(), nullptr);
    EXPECT_NE(p_created_element->pGetProperties(), nullptr);
}

KRATOS_TEST_CASE_IN_SUITE(TriangleInterfaceElement_CreatesInstanceWithNodeInput, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto nodes        = CreateNodesFor3Plus3SurfaceInterfaceGeometry();
    const auto p_properties = std::make_shared<Properties>();

    // The source element needs to have a geometry, otherwise the version of the
    // Create method with a node input will fail.
    const auto p_geometry = std::make_shared<TriangleInterfaceGeometry3D3Plus3Noded>(nodes);
    const InterfaceElement element(0, p_geometry, p_properties,
                                   std::make_unique<SurfaceInterfaceStressState>());

    // Act
    const auto p_created_element = element.Create(1, nodes, p_properties);

    // Assert
    EXPECT_NE(p_created_element, nullptr);
    EXPECT_EQ(p_created_element->Id(), 1);
    EXPECT_NE(p_created_element->pGetGeometry(), nullptr);
    EXPECT_NE(p_created_element->pGetProperties(), nullptr);
}

KRATOS_TEST_CASE_IN_SUITE(TriangleInterfaceElement_ReturnsTheExpectedDoFList, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto p_properties = std::make_shared<Properties>();

    Model model;
    const auto element = CreateHorizontal3Plus3NodedTriangleInterfaceElementWithUDofs(model, p_properties);

    // Act
    const auto              dummy_process_info = ProcessInfo{};
    Element::DofsVectorType degrees_of_freedom;
    element.GetDofList(degrees_of_freedom, dummy_process_info);

    // Assert
    KRATOS_EXPECT_EQ(degrees_of_freedom.size(), 18);
    for (std::size_t i = 0; i < degrees_of_freedom.size(); i += 3) {
        KRATOS_EXPECT_EQ(degrees_of_freedom[i]->GetVariable(), DISPLACEMENT_X);
        KRATOS_EXPECT_EQ(degrees_of_freedom[i + 1]->GetVariable(), DISPLACEMENT_Y);
        KRATOS_EXPECT_EQ(degrees_of_freedom[i + 2]->GetVariable(), DISPLACEMENT_Z);
    }
}

KRATOS_TEST_CASE_IN_SUITE(TriangleInterfaceElement_ReturnsTheExpectedEquationIdVector, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto p_properties = std::make_shared<Properties>();

    Model model;
    auto element = CreateHorizontal3Plus3NodedTriangleInterfaceElementWithUDofs(model, p_properties);

    int i = 0;
    for (const auto& node : element.GetGeometry()) {
        ++i;
        node.pGetDof(DISPLACEMENT_X)->SetEquationId(i);

        ++i;
        node.pGetDof(DISPLACEMENT_Y)->SetEquationId(i);

        ++i;
        node.pGetDof(DISPLACEMENT_Z)->SetEquationId(i);
    }

    // Act
    const auto                    dummy_process_info = ProcessInfo{};
    Element::EquationIdVectorType equation_id_vector;
    element.EquationIdVector(equation_id_vector, dummy_process_info);

    // Assert
    const Element::EquationIdVectorType expected_ids = {1,  2,  3,  4,  5,  6,  7,  8,  9,
                                                        10, 11, 12, 13, 14, 15, 16, 17, 18};
    KRATOS_EXPECT_VECTOR_EQ(equation_id_vector, expected_ids)
}

KRATOS_TEST_CASE_IN_SUITE(TriangleInterfaceElement_LeftHandSideContainsMaterialStiffnessContributions,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    constexpr auto normal_stiffness = 20.0;
    constexpr auto shear_stiffness  = 10.0;
    const auto     p_properties = CreateElasticMaterialProperties<InterfaceThreeDimensionalSurface>(
        normal_stiffness, shear_stiffness);

    auto element = CreateAndInitializeElement(
        CreateHorizontal3Plus3NodedTriangleInterfaceElementWithUDofs, p_properties);

    // Act
    Matrix actual_left_hand_side;
    element.CalculateLeftHandSide(actual_left_hand_side, ProcessInfo{});

    // Assert
    const auto expected_left_hand_side = ExpectedLeftHandSideForTriangleElement();
    KRATOS_EXPECT_MATRIX_RELATIVE_NEAR(actual_left_hand_side, expected_left_hand_side, Defaults::relative_tolerance)
}

KRATOS_TEST_CASE_IN_SUITE(TriangleInterfaceElement_LeftHandSideContainsMaterialStiffnessContributions_Rotated,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    constexpr auto normal_stiffness = 20.0;
    constexpr auto shear_stiffness  = 10.0;
    const auto     p_properties = CreateElasticMaterialProperties<InterfaceThreeDimensionalSurface>(
        normal_stiffness, shear_stiffness);

    auto element = CreateAndInitializeElement(
        CreateTriangleInterfaceElementRotatedBy30DegreesWithDisplacementDoF, p_properties);

    // Act
    Matrix actual_left_hand_side;
    element.CalculateLeftHandSide(actual_left_hand_side, ProcessInfo{});

    // Assert
    auto expected_left_hand_side = Matrix{18, 18};
    // clang-format off

    // Since the rotation is about the z-axis (the normal of the triangle) and the two shear
    // stiffnesses are equal, the left-hand side matrix is equal to the left-hand side of a
    // non-rotated surface interface element.
    KRATOS_EXPECT_MATRIX_NEAR(actual_left_hand_side, ExpectedLeftHandSideForTriangleElement(), Defaults::relative_tolerance)
}

KRATOS_TEST_CASE_IN_SUITE(TriangleInterfaceElement_LeftHandSideContainsMaterialStiffnessContributions_RotatedAboutYAxis,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    constexpr auto normal_stiffness = 20.0;
    constexpr auto shear_stiffness  = 10.0;
    const auto     p_properties = CreateElasticMaterialProperties<InterfaceThreeDimensionalSurface>(
        normal_stiffness, shear_stiffness);

    auto element = CreateAndInitializeElement(
        CreateTriangleInterfaceElementRotatedBy30DegreesAboutYAxisWithDisplacementDoF, p_properties);

    // Act
    Matrix actual_left_hand_side;
    element.CalculateLeftHandSide(actual_left_hand_side, ProcessInfo{});

    // Assert
    auto expected_left_hand_side = Matrix{18, 18};
    expected_left_hand_side <<=
    2.083333333333333,0,0.72168783648703216,0,0,0,0,0,0,-2.083333333333333,0,-0.72168783648703216,0,0,0,0,0,0,
    0,1.6666666666666661,0,0,0,0,0,0,0,0,-1.6666666666666661,0,0,0,0,0,0,0,
    0.72168783648703216,0,2.9166666666666661,0,0,0,0,0,0,-0.72168783648703216,0,-2.9166666666666661,0,0,0,0,0,0,
    0,0,0,2.083333333333333,0,0.72168783648703216,0,0,0,0,0,0,-2.083333333333333,0,-0.72168783648703216,0,0,0,
    0,0,0,0,1.6666666666666661,0,0,0,0,0,0,0,0,-1.6666666666666661,0,0,0,0,
    0,0,0,0.72168783648703216,0,2.9166666666666661,0,0,0,0,0,0,-0.72168783648703216,0,-2.9166666666666661,0,0,0,
    0,0,0,0,0,0,2.083333333333333,0,0.72168783648703216,0,0,0,0,0,0,-2.083333333333333,0,-0.72168783648703216,
    0,0,0,0,0,0,0,1.6666666666666661,0,0,0,0,0,0,0,0,-1.6666666666666661,0,
    0,0,0,0,0,0,0.72168783648703216,0,2.9166666666666661,0,0,0,0,0,0,-0.72168783648703216,0,-2.9166666666666661,
    -2.083333333333333,0,-0.72168783648703216,0,0,0,0,0,0,2.083333333333333,0,0.72168783648703216,0,0,0,0,0,0,
    0,-1.6666666666666661,0,0,0,0,0,0,0,0,1.6666666666666661,0,0,0,0,0,0,0,
    -0.72168783648703216,0,-2.9166666666666661,0,0,0,0,0,0,0.72168783648703216,0,2.9166666666666661,0,0,0,0,0,0,
    0,0,0,-2.083333333333333,0,-0.72168783648703216,0,0,0,0,0,0,2.083333333333333,0,0.72168783648703216,0,0,0,
    0,0,0,0,-1.6666666666666661,0,0,0,0,0,0,0,0,1.6666666666666661,0,0,0,0,
    0,0,0,-0.72168783648703216,0,-2.9166666666666661,0,0,0,0,0,0,0.72168783648703216,0,2.9166666666666661,0,0,0,
    0,0,0,0,0,0,-2.083333333333333,0,-0.72168783648703216,0,0,0,0,0,0,2.083333333333333,0,0.72168783648703216,
    0,0,0,0,0,0,0,-1.6666666666666661,0,0,0,0,0,0,0,0,1.6666666666666661,0,
    0,0,0,0,0,0,-0.72168783648703216,0,-2.9166666666666661,0,0,0,0,0,0,0.72168783648703216,0,2.9166666666666661;
    // clang-format off

    KRATOS_EXPECT_MATRIX_NEAR(actual_left_hand_side, expected_left_hand_side, Defaults::relative_tolerance)
}

KRATOS_TEST_CASE_IN_SUITE(TriangleInterfaceElement_RightHandSideEqualsMinusInternalForceVector,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    constexpr auto normal_stiffness = 20.0;
    constexpr auto shear_stiffness  = 10.0;
    const auto     p_properties = CreateElasticMaterialProperties<InterfaceThreeDimensionalSurface>(
        normal_stiffness, shear_stiffness);

    const auto prescribed_displacements = PrescribedDisplacements{
        {2, array_1d<double, 3>{0.2, 0.5, 0.0}}, {3, array_1d<double, 3>{0.2, 0.5, 0.0}}};
    auto element = CreateAndInitializeElement(CreateHorizontal3Plus3NodedTriangleInterfaceElementWithUDofs,
                                              p_properties, prescribed_displacements);

    // Act
    Vector actual_right_hand_side;
    element.CalculateRightHandSide(actual_right_hand_side, ProcessInfo{});

    // Assert
    auto expected_right_hand_side = Vector{18};
    expected_right_hand_side <<= 0.333333,0.833333,0,0,0,0,-0.333333,-0.833333,0,-0.333333,-0.833333,0,0,0,0,0.333333,0.833333,0;
    constexpr auto tolerance=1e-5;
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(actual_right_hand_side, expected_right_hand_side, tolerance)
}

KRATOS_TEST_CASE_IN_SUITE(TriangleInterfaceElement_RightHandSideEqualsMinusInternalForceVector_Rotated,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    constexpr auto normal_stiffness = 20.0;
    constexpr auto shear_stiffness  = 10.0;
    const auto     p_properties = CreateElasticMaterialProperties<InterfaceThreeDimensionalSurface>(
        normal_stiffness, shear_stiffness);

    const auto prescribed_displacements =
        PrescribedDisplacements{{3, array_1d<double, 3>{1.0, 2.0, 3.0}},
                                {4, array_1d<double, 3>{1.0, 2.0, 3.0}},
                                {5, array_1d<double, 3>{1.0, 2.0, 3.0}}};
    auto element = CreateAndInitializeElement(CreateTriangleInterfaceElementRotatedBy30DegreesWithDisplacementDoF,
                                              p_properties, prescribed_displacements);

    // Act
    Vector actual_right_hand_side;
    element.CalculateRightHandSide(actual_right_hand_side, ProcessInfo{});

    // Assert
    auto expected_right_hand_side = Vector{18};
    // clang-format off
    expected_right_hand_side <<= 1.66667,3.33333,10,1.66667,3.33333,10,1.66667,3.33333,10,-1.66667,-3.33333,-10,-1.66667,-3.33333,-10,-1.66667,-3.33333,-10;
    // clang-format on
    constexpr auto tolerance = 1e-5;
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(actual_right_hand_side, expected_right_hand_side, tolerance)
}

KRATOS_TEST_CASE_IN_SUITE(TriangleIntrfaceElement_GetInitializedConstitutiveLawsAfterElementInitialization,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto p_properties = std::make_shared<Properties>();
    p_properties->GetValue(CONSTITUTIVE_LAW) =
        std::make_shared<GeoIncrementalLinearElasticInterfaceLaw>(std::make_unique<InterfacePlaneStrain>());

    auto element = CreateAndInitializeElement(
        CreateHorizontal3Plus3NodedTriangleInterfaceElementWithUDofs, p_properties);

    // Act
    auto constitutive_laws = std::vector<ConstitutiveLaw::Pointer>{};
    element.CalculateOnIntegrationPoints(CONSTITUTIVE_LAW, constitutive_laws, ProcessInfo{});

    // Assert
    KRATOS_EXPECT_EQ(constitutive_laws.size(), 3);
    for (const auto& p_law : constitutive_laws) {
        KRATOS_EXPECT_NE(p_law.get(), nullptr);
        // Constitutive laws must have been cloned
        KRATOS_EXPECT_NE(p_law.get(), p_properties->GetValue(CONSTITUTIVE_LAW).get());
    }
}

KRATOS_TEST_CASE_IN_SUITE(TriangleInterfaceElement_HasCorrectNumberOfConstitutiveLawsAfterMultipleInitializations,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto p_properties = std::make_shared<Properties>();
    p_properties->GetValue(CONSTITUTIVE_LAW) =
        std::make_shared<GeoIncrementalLinearElasticInterfaceLaw>(std::make_unique<InterfacePlaneStrain>());

    Model model;
    auto element = CreateHorizontal3Plus3NodedTriangleInterfaceElementWithUDofs(model, p_properties);

    const auto dummy_process_info = ProcessInfo{};

    // Multiple initializations emulate a multi-stage simulation
    element.Initialize(dummy_process_info);
    element.Initialize(dummy_process_info);
    element.Initialize(dummy_process_info);

    // Assert
    auto constitutive_laws = std::vector<ConstitutiveLaw::Pointer>{};
    element.CalculateOnIntegrationPoints(CONSTITUTIVE_LAW, constitutive_laws, dummy_process_info);
    KRATOS_EXPECT_EQ(constitutive_laws.size(), 3);
}

KRATOS_TEST_CASE_IN_SUITE(TriangleInterfaceElement_ReturnsExpectedLeftAndRightHandSide,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    constexpr auto normal_stiffness = 20.0;
    constexpr auto shear_stiffness  = 10.0;
    const auto     p_properties = CreateElasticMaterialProperties<InterfaceThreeDimensionalSurface>(
        normal_stiffness, shear_stiffness);

    const auto prescribed_displacements = PrescribedDisplacements{
        {2, array_1d<double, 3>{0.2, 0.5, 0.0}}, {3, array_1d<double, 3>{0.2, 0.5, 0.0}}};
    auto element = CreateAndInitializeElement(CreateHorizontal3Plus3NodedTriangleInterfaceElementWithUDofs,
                                              p_properties, prescribed_displacements);

    // Act
    Vector actual_right_hand_side;
    Matrix actual_left_hand_side;
    element.CalculateLocalSystem(actual_left_hand_side, actual_right_hand_side, ProcessInfo{});

    // Assert
    const auto expected_left_hand_side = ExpectedLeftHandSideForTriangleElement();
    KRATOS_EXPECT_MATRIX_RELATIVE_NEAR(actual_left_hand_side, expected_left_hand_side, Defaults::relative_tolerance)

    auto expected_right_hand_side = Vector{18};
    expected_right_hand_side <<= 0.333333, 0.833333, 0, 0, 0, 0, -0.333333, -0.833333, 0, -0.333333,
        -0.833333, 0, 0, 0, 0, 0.333333, 0.833333, 0;
    constexpr auto tolerance = 1e-5;
    KRATOS_EXPECT_VECTOR_NEAR(actual_right_hand_side, expected_right_hand_side, tolerance)
}

KRATOS_TEST_CASE_IN_SUITE(TriangleInterfaceElement_CalculateStrain_ReturnsRelativeDisplacement,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    constexpr auto normal_stiffness = 20.0;
    constexpr auto shear_stiffness  = 10.0;
    const auto     p_properties = CreateElasticMaterialProperties<InterfaceThreeDimensionalSurface>(
        normal_stiffness, shear_stiffness);

    const auto prescribed_displacements =
        PrescribedDisplacements{{0, array_1d<double, 3>{1.0, 2.0, 3.0}},
                                {1, array_1d<double, 3>{1.0, 2.0, 3.0}},
                                {2, array_1d<double, 3>{1.0, 2.0, 3.0}}};
    auto element = CreateAndInitializeElement(CreateTriangleInterfaceElementRotatedBy30DegreesWithDisplacementDoF,
                                              p_properties, prescribed_displacements);

    // Act
    std::vector<Vector> relative_displacements_at_integration_points;
    element.CalculateOnIntegrationPoints(STRAIN, relative_displacements_at_integration_points, ProcessInfo{});

    // Assert
    Vector expected_relative_displacement{3};
    expected_relative_displacement <<= -3.0, -0.5 * std::sqrt(3.0) - 1.0, 0.5 - std::sqrt(3);

    KRATOS_EXPECT_EQ(relative_displacements_at_integration_points.size(), 3);
    for (const auto& r_relative_displacement : relative_displacements_at_integration_points) {
        KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(r_relative_displacement, expected_relative_displacement,
                                           Defaults::relative_tolerance)
    }
}

KRATOS_TEST_CASE_IN_SUITE(TriangleInterfaceElementHorizontal_CalculateStrain_ReturnsRelativeDisplacement,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    constexpr auto normal_stiffness = 20.0;
    constexpr auto shear_stiffness  = 10.0;
    const auto     p_properties = CreateElasticMaterialProperties<InterfaceThreeDimensionalSurface>(
        normal_stiffness, shear_stiffness);

    const auto prescribed_displacements =
        PrescribedDisplacements{{0, array_1d<double, 3>{0.0, 0.0, 1.0}},
                                {1, array_1d<double, 3>{0.0, 0.0, 1.0}},
                                {2, array_1d<double, 3>{0.0, 0.0, 1.0}}};
    auto element = CreateAndInitializeElement(CreateHorizontal3Plus3NodedTriangleInterfaceElementWithUDofs,
                                              p_properties, prescribed_displacements);

    // Act
    std::vector<Vector> relative_displacements_at_integration_points;
    element.CalculateOnIntegrationPoints(STRAIN, relative_displacements_at_integration_points, ProcessInfo{});

    // Assert
    Vector expected_relative_displacement{3};
    expected_relative_displacement <<= -1.0, 0.0, 0.0;

    KRATOS_EXPECT_EQ(relative_displacements_at_integration_points.size(), 3);
    for (const auto& r_relative_displacement : relative_displacements_at_integration_points) {
        KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(r_relative_displacement, expected_relative_displacement,
                                           Defaults::relative_tolerance)
    }
}

KRATOS_TEST_CASE_IN_SUITE(TriangleInterfaceElementInYZPlane_CalculateStrain_ReturnsRelativeDisplacement,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    constexpr auto normal_stiffness = 20.0;
    constexpr auto shear_stiffness  = 10.0;
    const auto     p_properties = CreateElasticMaterialProperties<InterfaceThreeDimensionalSurface>(
        normal_stiffness, shear_stiffness);

    const auto prescribed_displacements =
        PrescribedDisplacements{{3, array_1d<double, 3>{1.0, 0.0, 0.0}},
                                {4, array_1d<double, 3>{1.0, 0.0, 0.0}},
                                {5, array_1d<double, 3>{1.0, 0.0, 0.0}}};
    auto element = CreateAndInitializeElement(CreateHorizontal3Plus3NodedTriangleInterfaceYZPlaneElementWithUDofs,
                                              p_properties, prescribed_displacements);

    // Act
    std::vector<Vector> relative_displacements_at_integration_points;
    element.CalculateOnIntegrationPoints(STRAIN, relative_displacements_at_integration_points, ProcessInfo{});

    // Assert
    Vector expected_relative_displacement{3};
    expected_relative_displacement <<= 1.0, 0.0, 0.0;

    KRATOS_EXPECT_EQ(relative_displacements_at_integration_points.size(), 3);
    for (const auto& r_relative_displacement : relative_displacements_at_integration_points) {
        KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(r_relative_displacement, expected_relative_displacement,
                                           Defaults::relative_tolerance)
    }
}

KRATOS_TEST_CASE_IN_SUITE(TriangleInterfaceElementInXZPlane_CalculateStrain_ReturnsRelativeDisplacement,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    constexpr auto normal_stiffness = 20.0;
    constexpr auto shear_stiffness  = 10.0;
    const auto     p_properties = CreateElasticMaterialProperties<InterfaceThreeDimensionalSurface>(
        normal_stiffness, shear_stiffness);

    // The prescribed displacements are in the direction of the normal of the interface
    // which is the -y direction. This results in a positive normal relative displacement,
    // since the second side of the interface is moved.
    const auto prescribed_displacements =
        PrescribedDisplacements{{3, array_1d<double, 3>{0.0, -1.0, 0.0}},
                                {4, array_1d<double, 3>{0.0, -1.0, 0.0}},
                                {5, array_1d<double, 3>{0.0, -1.0, 0.0}}};
    auto element = CreateAndInitializeElement(CreateHorizontal3Plus3NodedTriangleInterfaceXZPlaneElementWithUDofs,
                                              p_properties, prescribed_displacements);

    // Act
    std::vector<Vector> relative_displacements_at_integration_points;
    element.CalculateOnIntegrationPoints(STRAIN, relative_displacements_at_integration_points, ProcessInfo{});

    // Assert
    Vector expected_relative_displacement{3};
    expected_relative_displacement <<= 1.0, 0.0, 0.0;

    KRATOS_EXPECT_EQ(relative_displacements_at_integration_points.size(), 3);
    for (const auto& r_relative_displacement : relative_displacements_at_integration_points) {
        KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(r_relative_displacement, expected_relative_displacement,
                                           Defaults::relative_tolerance)
    }
}

KRATOS_TEST_CASE_IN_SUITE(TriangleInterfaceElement_CalculateCauchyStressVector_ReturnsTraction,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    constexpr auto normal_stiffness = 20.0;
    constexpr auto shear_stiffness  = 10.0;
    const auto     p_properties = CreateElasticMaterialProperties<InterfaceThreeDimensionalSurface>(
        normal_stiffness, shear_stiffness);

    const auto prescribed_displacements =
        PrescribedDisplacements{{0, array_1d<double, 3>{1.0, 2.0, 3.0}},
                                {1, array_1d<double, 3>{1.0, 2.0, 3.0}},
                                {2, array_1d<double, 3>{1.0, 2.0, 3.0}}};
    auto element = CreateAndInitializeElement(CreateTriangleInterfaceElementRotatedBy30DegreesWithDisplacementDoF,
                                              p_properties, prescribed_displacements);

    // Act
    std::vector<Vector> tractions_at_integration_points;
    element.CalculateOnIntegrationPoints(CAUCHY_STRESS_VECTOR, tractions_at_integration_points, ProcessInfo{});

    // Assert
    Vector expected_traction{3};
    expected_traction <<= -3.0 * normal_stiffness, (-0.5 * std::sqrt(3.0) - 1.0) * shear_stiffness,
        (0.5 - std::sqrt(3)) * shear_stiffness;
    KRATOS_EXPECT_EQ(tractions_at_integration_points.size(), 3);
    for (const auto& r_traction : tractions_at_integration_points) {
        KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(r_traction, expected_traction, Defaults::relative_tolerance)
    }
}

KRATOS_TEST_CASE_IN_SUITE(TriangleInterfaceElement_6Plus6NodedElement_ReturnsExpectedLeftAndRightHandSide,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    constexpr auto normal_stiffness = 20.0;
    constexpr auto shear_stiffness  = 10.0;
    const auto     p_properties = CreateElasticMaterialProperties<InterfaceThreeDimensionalSurface>(
        normal_stiffness, shear_stiffness);

    const auto prescribed_displacements =
        PrescribedDisplacements{{3, array_1d<double, 3>{0.2, 0.5, 0.0}},
                                {4, array_1d<double, 3>{0.2, 0.5, 0.0}},
                                {5, array_1d<double, 3>{0.2, 0.5, 0.0}}};
    auto element = CreateAndInitializeElement(CreateHorizontal6Plus6NodedTriangleInterfaceElementWithDisplacementDoF,
                                              p_properties, prescribed_displacements);

    // Act
    Vector actual_right_hand_side;
    Matrix actual_left_hand_side;
    element.CalculateLocalSystem(actual_left_hand_side, actual_right_hand_side, ProcessInfo{});

    // Assert
    auto expected_left_hand_side = Matrix(36, 36); // the values are taken from the element
    // clang-format off
    expected_left_hand_side <<=
    0.16431924882629123,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.16431924882629123,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0.16431924882629123,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.16431924882629123,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0.32863849765258246,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.32863849765258246,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0.16431924882629123,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.16431924882629123,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0.16431924882629123,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.16431924882629123,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0.32863849765258246,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.32863849765258246,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0.16431924882629123,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.16431924882629123,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0.16431924882629123,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.16431924882629123,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0.32863849765258246,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.32863849765258246,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,1.5023474178403748,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1.5023474178403748,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,1.5023474178403748,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1.5023474178403748,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,3.0046948356807497,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-3.0046948356807497,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,1.502347417840376,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1.502347417840376,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,1.502347417840376,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1.502347417840376,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.0046948356807519,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-3.0046948356807519,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1.5023474178403753,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1.5023474178403753,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1.5023474178403753,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1.5023474178403753,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.0046948356807506,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-3.0046948356807506,
    -0.16431924882629123,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.16431924882629123,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,-0.16431924882629123,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.16431924882629123,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,-0.32863849765258246,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.32863849765258246,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,-0.16431924882629123,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.16431924882629123,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,-0.16431924882629123,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.16431924882629123,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,-0.32863849765258246,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.32863849765258246,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,-0.16431924882629123,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.16431924882629123,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,-0.16431924882629123,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.16431924882629123,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,-0.32863849765258246,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.32863849765258246,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,-1.5023474178403748,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1.5023474178403748,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,-1.5023474178403748,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1.5023474178403748,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,-3.0046948356807497,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.0046948356807497,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,-1.502347417840376,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1.502347417840376,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,-1.502347417840376,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1.502347417840376,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,-3.0046948356807519,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.0046948356807519,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1.5023474178403753,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1.5023474178403753,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1.5023474178403753,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1.5023474178403753,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-3.0046948356807506,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.0046948356807506;
    // clang-format on

    KRATOS_EXPECT_MATRIX_RELATIVE_NEAR(actual_left_hand_side, expected_left_hand_side, Defaults::relative_tolerance)

    auto expected_right_hand_side = Vector(36);
    expected_right_hand_side <<= 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.300469, -0.751174, 0, -0.300469,
        -0.751174, 0, -0.300469, -0.751174, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.300469, 0.751174, 0,
        0.300469, 0.751174, 0, 0.300469, 0.751174, 0;
    constexpr auto tolerance = 1e-5;
    KRATOS_EXPECT_VECTOR_NEAR(actual_right_hand_side, expected_right_hand_side, tolerance)
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceElement_CheckThrowsWhenElementIsNotInitialized, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    constexpr auto normal_stiffness = 20.0;
    constexpr auto shear_stiffness  = 10.0;
    const auto     p_properties =
        CreateElasticMaterialProperties<InterfacePlaneStrain>(normal_stiffness, shear_stiffness);

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
    const auto     p_properties =
        CreateElasticMaterialProperties<InterfacePlaneStrain>(normal_stiffness, shear_stiffness);

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
