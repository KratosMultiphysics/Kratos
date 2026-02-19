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
#include "custom_constitutive/incremental_linear_elastic_law.h"
#include "custom_constitutive/interface_plane_strain.h"
#include "custom_constitutive/interface_three_dimensional_surface.h"
#include "custom_constitutive/plane_strain.h"
#include "custom_constitutive/three_dimensional.h"
#include "custom_elements/U_Pw_interface_element.h"
#include "custom_elements/interface_stress_state.h"
#include "custom_geometries/interface_geometry.hpp"
#include "custom_utilities/ublas_utilities.h"
#include "geo_mechanics_application_variables.h"
#include "test_setup_utilities/element_setup_utilities.hpp"
#include "test_setup_utilities/model_setup_utilities.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "tests/cpp_tests/test_utilities.h"

#include <boost/numeric/ublas/assignment.hpp>

namespace
{

using namespace Kratos;
using LineInterfaceGeometry2D2Plus2Noded          = InterfaceGeometry<Line2D2<Node>>;
using LineInterfaceGeometry2D3Plus3Noded          = InterfaceGeometry<Line2D3<Node>>;
using TriangleInterfaceGeometry3D3Plus3Noded      = InterfaceGeometry<Triangle3D3<Node>>;
using TriangleInterfaceGeometry3D6Plus6Noded      = InterfaceGeometry<Triangle3D6<Node>>;
using QuadrilateralInterfaceGeometry3D4Plus4Noded = InterfaceGeometry<Quadrilateral3D4<Node>>;
using Interface2D                                 = Line2DInterfaceStressState;
using Interface3D                                 = SurfaceInterfaceStressState;
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

ModelPart& CreateModelPartWithUPwVariables(Model& rModel)
{
    auto& r_result = rModel.CreateModelPart("Main");
    r_result.AddNodalSolutionStepVariable(DISPLACEMENT);
    r_result.AddNodalSolutionStepVariable(WATER_PRESSURE);

    return r_result;
}

template <typename TInterfaceDimension>
UPwInterfaceElement CreateInterfaceElementWithUPwDofs(const Properties::Pointer&     rpProperties,
                                                      const Geometry<Node>::Pointer& rpGeometry,
                                                      IsDiffOrderElement             IsDiffOrder,
                                                      const std::vector<CalculationContribution>& rContributions)
{
    auto result = UPwInterfaceElement{
        1,           rpGeometry,    rpProperties, std::make_unique<TInterfaceDimension>(),
        IsDiffOrder, rContributions};
    const auto solution_step_variables =
        Geo::ConstVariableDataRefs{std::cref(WATER_PRESSURE), std::cref(DISPLACEMENT)};
    const auto degrees_of_freedom =
        Geo::ConstVariableRefs{std::cref(WATER_PRESSURE), std::cref(DISPLACEMENT_X),
                               std::cref(DISPLACEMENT_Y), std::cref(DISPLACEMENT_Z)};
    auto p_interface_element = &result;
    // Note that we're a bit sloppy here, since we add the degrees of freedom to _all_ nodes (even for diff-order elements).
    // However, you'll find the same sloppiness in the fully integrated workflow. That needs to be improved later.
    Testing::ElementSetupUtilities::AddVariablesToEntity(
        p_interface_element, solution_step_variables, degrees_of_freedom);

    return result;
}

UPwInterfaceElement CreateHorizontalUnitLength2Plus2NodedLineInterfaceElementWithUPwDofs(
    Model&                                      rModel,
    const Properties::Pointer&                  rpProperties,
    IsDiffOrderElement                          IsDiffOrder,
    const std::vector<CalculationContribution>& rContributions)
{
    auto& r_model_part = CreateModelPartWithUPwVariables(rModel);

    PointerVector<Node> nodes;
    nodes.push_back(r_model_part.CreateNewNode(0, 0.0, 0.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(1, 1.0, 0.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(2, 0.0, 0.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(3, 1.0, 0.0, 0.0));
    const auto p_geometry = std::make_shared<LineInterfaceGeometry2D2Plus2Noded>(nodes);
    return CreateInterfaceElementWithUPwDofs<Interface2D>(rpProperties, p_geometry, IsDiffOrder, rContributions);
}

UPwInterfaceElement CreateHorizontalUnitLength3Plus3NodedLineInterfaceElementWithUPwDoF(
    Model&                                      rModel,
    const Properties::Pointer&                  rpProperties,
    IsDiffOrderElement                          IsDiffOrder,
    const std::vector<CalculationContribution>& rContributions)
{
    auto& r_model_part = CreateModelPartWithUPwVariables(rModel);

    PointerVector<Node> nodes;
    nodes.push_back(r_model_part.CreateNewNode(0, 0.0, 0.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(1, 1.0, 0.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(2, 0.5, 0.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(3, 0.0, 0.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(4, 1.0, 0.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(5, 0.5, 0.0, 0.0));
    const auto p_geometry = std::make_shared<LineInterfaceGeometry2D3Plus3Noded>(nodes);
    return CreateInterfaceElementWithUPwDofs<Interface2D>(rpProperties, p_geometry, IsDiffOrder, rContributions);
}

UPwInterfaceElement CreateUnitLengthLineInterfaceElementRotatedBy30DegreesWithUPwDoF(
    Model&                                      rModel,
    const Properties::Pointer&                  rpProperties,
    IsDiffOrderElement                          IsDiffOrder,
    const std::vector<CalculationContribution>& rContributions)
{
    auto& r_model_part = CreateModelPartWithUPwVariables(rModel);

    PointerVector<Node> nodes;
    nodes.push_back(r_model_part.CreateNewNode(0, 0.0, 0.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(1, 0.5 * std::sqrt(3.0), 0.5, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(2, 0.0, 0.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(3, 0.5 * std::sqrt(3.0), 0.5, 0.0));
    const auto p_geometry = std::make_shared<LineInterfaceGeometry2D2Plus2Noded>(nodes);
    return CreateInterfaceElementWithUPwDofs<Interface2D>(rpProperties, p_geometry, IsDiffOrder, rContributions);
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

UPwInterfaceElement CreateHorizontal3Plus3NodedTriangleInterfaceElementWithUPwDofs(
    Model&                                      rModel,
    const Properties::Pointer&                  rpProperties,
    IsDiffOrderElement                          IsDiffOrder,
    const std::vector<CalculationContribution>& rContributions)
{
    auto& r_model_part = CreateModelPartWithUPwVariables(rModel);

    PointerVector<Node> nodes;
    nodes.push_back(r_model_part.CreateNewNode(0, 0.0, 0.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(1, 1.0, 0.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(2, 0.0, 1.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(3, 0.0, 0.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(4, 1.0, 0.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(5, 0.0, 1.0, 0.0));
    const auto p_geometry = std::make_shared<TriangleInterfaceGeometry3D3Plus3Noded>(nodes);
    return CreateInterfaceElementWithUPwDofs<Interface3D>(rpProperties, p_geometry, IsDiffOrder, rContributions);
}

UPwInterfaceElement CreateHorizontal3Plus3NodedTriangleInterfaceYZPlaneElementWithUPwDofs(
    Model&                                      rModel,
    const Properties::Pointer&                  rpProperties,
    IsDiffOrderElement                          IsDiffOrder,
    const std::vector<CalculationContribution>& rContributions)
{
    auto& r_model_part = CreateModelPartWithUPwVariables(rModel);

    PointerVector<Node> nodes;
    nodes.push_back(r_model_part.CreateNewNode(0, 0.0, 0.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(1, 0.0, 0.0, -1.0));
    nodes.push_back(r_model_part.CreateNewNode(2, 0.0, 1.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(3, 0.0, 0.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(4, 0.0, 0.0, -1.0));
    nodes.push_back(r_model_part.CreateNewNode(5, 0.0, 1.0, 0.0));
    const auto p_geometry = std::make_shared<TriangleInterfaceGeometry3D3Plus3Noded>(nodes);
    return CreateInterfaceElementWithUPwDofs<Interface3D>(rpProperties, p_geometry, IsDiffOrder, rContributions);
}

UPwInterfaceElement CreateHorizontal3Plus3NodedTriangleInterfaceXZPlaneElementWithUPwDofs(
    Model&                                      rModel,
    const Properties::Pointer&                  rpProperties,
    IsDiffOrderElement                          IsDiffOrder,
    const std::vector<CalculationContribution>& rContributions)
{
    auto& r_model_part = CreateModelPartWithUPwVariables(rModel);

    PointerVector<Node> nodes;
    nodes.push_back(r_model_part.CreateNewNode(0, 0.0, 0.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(1, 1.0, 0.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(2, 0.0, 0.0, 1.0));
    nodes.push_back(r_model_part.CreateNewNode(3, 0.0, 0.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(4, 1.0, 0.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(5, 0.0, 0.0, 1.0));
    const auto p_geometry = std::make_shared<TriangleInterfaceGeometry3D3Plus3Noded>(nodes);
    return CreateInterfaceElementWithUPwDofs<Interface3D>(rpProperties, p_geometry, IsDiffOrder, rContributions);
}

UPwInterfaceElement CreateHorizontal6Plus6NodedTriangleInterfaceElementWithUPwDoF(
    Model&                                      rModel,
    const Properties::Pointer&                  rpProperties,
    IsDiffOrderElement                          IsDiffOrder,
    const std::vector<CalculationContribution>& rContributions)
{
    auto& r_model_part = CreateModelPartWithUPwVariables(rModel);

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
    return CreateInterfaceElementWithUPwDofs<Interface3D>(rpProperties, p_geometry, IsDiffOrder, rContributions);
}

UPwInterfaceElement CreateTriangleInterfaceElementRotatedBy30DegreesWithUPwDoF(
    Model&                                      rModel,
    const Properties::Pointer&                  rpProperties,
    IsDiffOrderElement                          IsDiffOrder,
    const std::vector<CalculationContribution>& rContributions)
{
    auto& r_model_part = CreateModelPartWithUPwVariables(rModel);

    PointerVector<Node> nodes;
    nodes.push_back(r_model_part.CreateNewNode(0, 0.0, 0.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(1, 0.5 * std::sqrt(3.0), 0.5, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(2, -0.5, 0.5 * std::sqrt(3.0), 0.0));
    nodes.push_back(r_model_part.CreateNewNode(3, 0.0, 0.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(4, 0.5 * std::sqrt(3.0), 0.5, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(5, -0.5, 0.5 * std::sqrt(3.0), 0.0));
    const auto p_geometry = std::make_shared<TriangleInterfaceGeometry3D3Plus3Noded>(nodes);
    return CreateInterfaceElementWithUPwDofs<Interface3D>(rpProperties, p_geometry, IsDiffOrder, rContributions);
}

UPwInterfaceElement CreateTriangleInterfaceElementRotatedBy30DegreesAboutYAxisWithUPwDoF(
    Model&                                      rModel,
    const Properties::Pointer&                  rpProperties,
    IsDiffOrderElement                          IsDiffOrder,
    const std::vector<CalculationContribution>& rContributions)
{
    auto& r_model_part = CreateModelPartWithUPwVariables(rModel);

    PointerVector<Node> nodes;
    nodes.push_back(r_model_part.CreateNewNode(0, 0.0, 0.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(1, 0.5 * std::sqrt(3.0), 0.0, -0.5));
    nodes.push_back(r_model_part.CreateNewNode(2, 0.0, 1.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(3, 0.0, 0.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(4, 0.5 * std::sqrt(3.0), 0.0, -0.5));
    nodes.push_back(r_model_part.CreateNewNode(5, 0.0, 1.0, 0.0));
    const auto p_geometry = std::make_shared<TriangleInterfaceGeometry3D3Plus3Noded>(nodes);
    return CreateInterfaceElementWithUPwDofs<Interface3D>(rpProperties, p_geometry, IsDiffOrder, rContributions);
}

UPwInterfaceElement CreateHorizontal4Plus4NodedQuadraliteralInterfaceElementWithUPwDoF(
    Model&,
    const Properties::Pointer&                  rpProperties,
    IsDiffOrderElement                          IsDiffOrder,
    const std::vector<CalculationContribution>& rContributions)
{
    const auto nodes = Testing::ModelSetupUtilities::CreateNodes({{1, {0.0, 0.0, 0.0}},
                                                                  {2, {1.0, 0.0, 0.0}},
                                                                  {3, {1.0, 1.0, 0.0}},
                                                                  {4, {0.0, 1.0, 0.0}},
                                                                  {5, {0.0, 0.0, 0.0}},
                                                                  {6, {1.0, 0.0, 0.0}},
                                                                  {7, {1.0, 1.0, 0.0}},
                                                                  {8, {0.0, 1.0, 0.0}}});

    auto p_geometry = std::make_shared<InterfaceGeometry<Quadrilateral3D4<Node>>>(1, nodes);
    return CreateInterfaceElementWithUPwDofs<Interface3D>(rpProperties, p_geometry, IsDiffOrder, rContributions);
}

UPwInterfaceElement CreateHorizontal8Plus8NodedQuadraliteralInterfaceElementWithUPwDoF(
    Model&,
    const Properties::Pointer&                  rpProperties,
    IsDiffOrderElement                          IsDiffOrder,
    const std::vector<CalculationContribution>& rContributions)
{
    const auto nodes = Testing::ModelSetupUtilities::CreateNodes({{1, {0.0, 0.0, 0.0}},
                                                                  {2, {1.0, 0.0, 0.0}},
                                                                  {3, {1.0, 1.0, 0.0}},
                                                                  {4, {0.0, 1.0, 0.0}},
                                                                  {5, {0.5, 0.0, 0.0}},
                                                                  {6, {1.0, 0.5, 0.0}},
                                                                  {7, {0.5, 1.0, 0.0}},
                                                                  {8, {0.0, 0.5, 0.0}},
                                                                  {9, {0.0, 0.0, 0.0}},
                                                                  {10, {1.0, 0.0, 0.0}},
                                                                  {11, {1.0, 1.0, 0.0}},
                                                                  {12, {0.0, 1.0, 0.0}},
                                                                  {13, {0.5, 0.0, 0.0}},
                                                                  {14, {1.0, 0.5, 0.0}},
                                                                  {15, {0.5, 1.0, 0.0}},
                                                                  {16, {0.0, 0.5, 0.0}}});

    auto p_geometry = std::make_shared<InterfaceGeometry<Quadrilateral3D8<Node>>>(1, nodes);
    return CreateInterfaceElementWithUPwDofs<Interface3D>(rpProperties, p_geometry, IsDiffOrder, rContributions);
}

template <typename TElementFactory>
UPwInterfaceElement CreateAndInitializeElement(TElementFactory            Factory,
                                               const Properties::Pointer& rpProperties,
                                               IsDiffOrderElement         IsDiffOrder,
                                               const std::vector<CalculationContribution>& rContributions,
                                               const PrescribedDisplacements& rDisplacements = {})
{
    Model model;
    auto  element = Factory(model, rpProperties, IsDiffOrder, rContributions);
    element.Initialize(ProcessInfo{});
    for (const auto& [idx, disp] : rDisplacements) {
        element.GetGeometry()[idx].FastGetSolutionStepValue(DISPLACEMENT) = disp;
    }
    return element;
}

Matrix ExpectedStiffnessMatrixOfLinearTriangularInterfaceElement()
{
    constexpr auto cs = 5.0 / 3.0;
    constexpr auto cn = 10.0 / 3.0;
    return UblasUtilities::CreateMatrix({
        {cs, 0, 0, 0, 0, 0, 0, 0, 0, -cs, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, cs, 0, 0, 0, 0, 0, 0, 0, 0, -cs, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, cn, 0, 0, 0, 0, 0, 0, 0, 0, -cn, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, cs, 0, 0, 0, 0, 0, 0, 0, 0, -cs, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, cs, 0, 0, 0, 0, 0, 0, 0, 0, -cs, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, cn, 0, 0, 0, 0, 0, 0, 0, 0, -cn, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, cs, 0, 0, 0, 0, 0, 0, 0, 0, -cs, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, cs, 0, 0, 0, 0, 0, 0, 0, 0, -cs, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, cn, 0, 0, 0, 0, 0, 0, 0, 0, -cn},
        {-cs, 0, 0, 0, 0, 0, 0, 0, 0, cs, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, -cs, 0, 0, 0, 0, 0, 0, 0, 0, cs, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, -cn, 0, 0, 0, 0, 0, 0, 0, 0, cn, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, -cs, 0, 0, 0, 0, 0, 0, 0, 0, cs, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, -cs, 0, 0, 0, 0, 0, 0, 0, 0, cs, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, -cn, 0, 0, 0, 0, 0, 0, 0, 0, cn, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, -cs, 0, 0, 0, 0, 0, 0, 0, 0, cs, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, -cs, 0, 0, 0, 0, 0, 0, 0, 0, cs, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, -cn, 0, 0, 0, 0, 0, 0, 0, 0, cn},
    });
}

Vector ExpectedStiffnessForceOfLinearTriangularInterfaceElement()
{
    return UblasUtilities::CreateVector({1.0 / 3.0, 5.0 / 6.0, 0.0, 0.0, 0.0, 0.0, -1.0 / 3.0,
                                         -5.0 / 6.0, 0.0, -1.0 / 3.0, -5.0 / 6.0, 0.0, 0.0, 0.0,
                                         0.0, 1.0 / 3.0, 5.0 / 6.0, 0.0});
}

class MockElementWithTotalStressVectors : public Element
{
public:
    MockElementWithTotalStressVectors(std::size_t                     ElementId,
                                      std::shared_ptr<Geometry<Node>> pGeometry,
                                      Properties::Pointer             pProperties);

    void SetValuesOnIntegrationPoints(const Variable<Vector>&    rVariable,
                                      const std::vector<Vector>& rValues,
                                      const ProcessInfo&) override;
    using Element::SetValuesOnIntegrationPoints;

    void CalculateOnIntegrationPoints(const Variable<Vector>& rVariable,
                                      std::vector<Vector>&    rOutput,
                                      const ProcessInfo&) override;
    using Element::CalculateOnIntegrationPoints;

    IntegrationMethod GetIntegrationMethod() const override;
    void              SetIntegrationMethod(IntegrationMethod CustomIntegrationMethod);

private:
    std::vector<Vector>              mTotalStressVectors;
    std::optional<IntegrationMethod> mOptionalCustomIntegrationMethod;
};

MockElementWithTotalStressVectors::MockElementWithTotalStressVectors(std::size_t ElementId,
                                                                     std::shared_ptr<Geometry<Node>> pGeometry,
                                                                     Properties::Pointer pProperties)
    : Element{ElementId, std::move(pGeometry), std::move(pProperties)}
{
}

void MockElementWithTotalStressVectors::SetValuesOnIntegrationPoints(const Variable<Vector>& rVariable,
                                                                     const std::vector<Vector>& rValues,
                                                                     const ProcessInfo&)
{
    KRATOS_DEBUG_ERROR_IF_NOT(rVariable == TOTAL_STRESS_VECTOR)
        << "This mock element can only set total stress vectors\n";

    mTotalStressVectors = rValues;
}

void MockElementWithTotalStressVectors::CalculateOnIntegrationPoints(const Variable<Vector>& rVariable,
                                                                     std::vector<Vector>& rOutput,
                                                                     const ProcessInfo&)
{
    KRATOS_DEBUG_ERROR_IF_NOT(rVariable == TOTAL_STRESS_VECTOR)
        << "This mock element can only calculate total stress vectors\n";

    rOutput = mTotalStressVectors;
}

GeometryData::IntegrationMethod MockElementWithTotalStressVectors::GetIntegrationMethod() const
{
    return mOptionalCustomIntegrationMethod.value_or(GetGeometry().GetDefaultIntegrationMethod());
}

void MockElementWithTotalStressVectors::SetIntegrationMethod(IntegrationMethod CustomIntegrationMethod)
{
    mOptionalCustomIntegrationMethod = CustomIntegrationMethod;
}

template <typename DerivedElementPtrType>
GlobalPointersVector<Element> MakeElementGlobalPtrContainerWith(const DerivedElementPtrType& rpElement)
{
    auto p_element = Element::Pointer{rpElement};
    return {GlobalPointer<Element>{p_element}};
}

template <typename TConstitutiveLawDimension, typename TElementFactory>
void GeneraizedCouplingContributionTest(TElementFactory&&                           ElementFactory,
                                        const std::vector<CalculationContribution>& rContributions,
                                        IsDiffOrderElement                          diff_order,
                                        std::size_t   number_of_u_dofs,
                                        std::size_t   number_of_pw_dofs,
                                        const Matrix& expected_up_block_matrix,
                                        const Vector& expected_up_block_vector)
{
    // Arrange
    const auto p_properties = std::make_shared<Properties>();
    p_properties->GetValue(CONSTITUTIVE_LAW) = std::make_shared<GeoIncrementalLinearElasticInterfaceLaw>(
        std::make_unique<TConstitutiveLawDimension>());

    p_properties->SetValue(RETENTION_LAW, "VanGenuchtenLaw");
    p_properties->SetValue(SATURATED_SATURATION, 1.0);
    p_properties->SetValue(RESIDUAL_SATURATION, 0.06203);
    p_properties->SetValue(VAN_GENUCHTEN_AIR_ENTRY_PRESSURE, 2.561);
    p_properties->SetValue(VAN_GENUCHTEN_GN, 1.377);
    Model model;
    auto  interface_element = ElementFactory(model, p_properties, diff_order, rContributions);

    // Set nonzero water pressure at each node to ensure coupling code is exercised
    const auto  number_of_nodes_on_side = interface_element.GetGeometry().PointsNumber() / 2;
    std::size_t counter_nodes           = 0;
    for (auto& node : interface_element.GetGeometry()) {
        node.FastGetSolutionStepValue(WATER_PRESSURE) = (counter_nodes < number_of_nodes_on_side) ? 100.0 : 1.0;
        ++counter_nodes;
    }

    interface_element.Initialize(ProcessInfo{});

    // Act
    Vector actual_right_hand_side;
    Matrix actual_left_hand_side;
    interface_element.CalculateLocalSystem(actual_left_hand_side, actual_right_hand_side, ProcessInfo{});

    // Assert
    ASSERT_EQ(actual_left_hand_side.size1(), number_of_u_dofs + number_of_pw_dofs);
    ASSERT_EQ(actual_left_hand_side.size2(), number_of_u_dofs + number_of_pw_dofs);
    ASSERT_EQ(actual_right_hand_side.size(), number_of_u_dofs + number_of_pw_dofs);

    Testing::AssertUPBlockMatrixIsNear(actual_left_hand_side, expected_up_block_matrix, number_of_u_dofs,
                                       number_of_pw_dofs, Testing::Defaults::relative_tolerance);
    KRATOS_EXPECT_VECTOR_NEAR(actual_right_hand_side, expected_up_block_vector, Testing::Defaults::relative_tolerance)
}

} // namespace

namespace Kratos::Testing
{
using namespace Kratos;

KRATOS_TEST_CASE_IN_SUITE(UPwInterfaceElement_IsAnElement, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const UPwInterfaceElement element(
        0, std::make_shared<LineInterfaceGeometry2D2Plus2Noded>(CreateNodesFor2Plus2LineInterfaceGeometry()),
        std::make_unique<Line2DInterfaceStressState>(), IsDiffOrderElement::No,
        {CalculationContribution::Stiffness});
    auto p_casted_element = dynamic_cast<const Element*>(&element);
    EXPECT_NE(p_casted_element, nullptr);
}

KRATOS_TEST_CASE_IN_SUITE(UPwLineInterfaceElement_CreatesInstanceWithGeometryInput, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const UPwInterfaceElement element(
        0, std::make_shared<LineInterfaceGeometry2D2Plus2Noded>(CreateNodesFor2Plus2LineInterfaceGeometry()),
        std::make_unique<Line2DInterfaceStressState>(), IsDiffOrderElement::No,
        {CalculationContribution::Stiffness});
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

KRATOS_TEST_CASE_IN_SUITE(UPwLineInterfaceElement_CreatesInstanceWithNodeInput, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto nodes        = CreateNodesFor2Plus2LineInterfaceGeometry();
    const auto p_properties = std::make_shared<Properties>();

    // The source element needs to have a geometry, otherwise the version of the
    // Create method with a node input will fail.
    const auto p_geometry = std::make_shared<LineInterfaceGeometry2D2Plus2Noded>(nodes);
    const UPwInterfaceElement element(0, p_geometry, p_properties,
                                      std::make_unique<Line2DInterfaceStressState>(),
                                      IsDiffOrderElement::No, {CalculationContribution::Stiffness});

    // Act
    const auto p_created_element = element.Create(1, nodes, p_properties);

    // Assert
    EXPECT_NE(p_created_element, nullptr);
    EXPECT_EQ(p_created_element->Id(), 1);
    EXPECT_NE(p_created_element->pGetGeometry(), nullptr);
    EXPECT_NE(p_created_element->pGetProperties(), nullptr);
}

KRATOS_TEST_CASE_IN_SUITE(UPwLineInterfaceElement_KeepsUDofsFirstThenPwDofs, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto p_properties = std::make_shared<Properties>();

    Model      model;
    const auto element = CreateHorizontalUnitLength2Plus2NodedLineInterfaceElementWithUPwDofs(
        model, p_properties, IsDiffOrderElement::No, {CalculationContribution::Stiffness});

    // Act
    const auto              dummy_process_info = ProcessInfo{};
    Element::DofsVectorType degrees_of_freedom;
    element.GetDofList(degrees_of_freedom, dummy_process_info);

    // Assert
    constexpr auto expected_number_of_u_dofs  = std::size_t{4 * 2};
    constexpr auto expected_number_of_pw_dofs = std::size_t{4};
    ASSERT_EQ(degrees_of_freedom.size(), expected_number_of_u_dofs + expected_number_of_pw_dofs);
    for (auto i = std::size_t{0}; i < expected_number_of_u_dofs; i += 2) {
        KRATOS_EXPECT_EQ(degrees_of_freedom[i]->GetVariable(), DISPLACEMENT_X);
        KRATOS_EXPECT_EQ(degrees_of_freedom[i + 1]->GetVariable(), DISPLACEMENT_Y);
    }
    for (auto i = std::size_t{0}; i < expected_number_of_pw_dofs; ++i) {
        KRATOS_EXPECT_EQ(degrees_of_freedom[expected_number_of_u_dofs + i]->GetVariable(), WATER_PRESSURE);
    }
}

KRATOS_TEST_CASE_IN_SUITE(UPwLineInterfaceElement_ReturnsTheExpectedEquationIdVector, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto p_properties = std::make_shared<Properties>();

    Model model;
    auto  element = CreateHorizontalUnitLength2Plus2NodedLineInterfaceElementWithUPwDofs(
        model, p_properties, IsDiffOrderElement::No, {CalculationContribution::Stiffness});

    int i = 0;
    for (const auto& node : element.GetGeometry()) {
        ++i;
        node.pGetDof(DISPLACEMENT_X)->SetEquationId(i);

        ++i;
        node.pGetDof(DISPLACEMENT_Y)->SetEquationId(i);

        ++i;
        node.pGetDof(WATER_PRESSURE)->SetEquationId(i);
    }

    // Act
    const auto                    dummy_process_info = ProcessInfo{};
    Element::EquationIdVectorType equation_id_vector;
    element.EquationIdVector(equation_id_vector, dummy_process_info);

    // Assert
    const Element::EquationIdVectorType expected_ids = {1, 2, 4, 5, 7, 8, 10, 11, 3, 6, 9, 12};
    KRATOS_EXPECT_VECTOR_EQ(equation_id_vector, expected_ids)
}

KRATOS_TEST_CASE_IN_SUITE(UPwLineInterfaceElement_LeftHandSideContainsMaterialStiffnessContributions,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    constexpr auto normal_stiffness = 20.0;
    constexpr auto shear_stiffness  = 10.0;
    const auto     p_properties =
        CreateElasticMaterialProperties<InterfacePlaneStrain>(normal_stiffness, shear_stiffness);
    auto element = CreateAndInitializeElement(
        CreateHorizontalUnitLength2Plus2NodedLineInterfaceElementWithUPwDofs, p_properties,
        IsDiffOrderElement::No, {CalculationContribution::Stiffness});

    // Act
    Matrix actual_left_hand_side;
    element.CalculateLeftHandSide(actual_left_hand_side, ProcessInfo{});

    // Assert
    constexpr auto number_of_u_dofs  = std::size_t{4 * 2};
    constexpr auto number_of_pw_dofs = std::size_t{4};
    const auto     expected_uu_block_matrix =
        CreateExpectedStiffnessMatrixForHorizontal2Plus2NodedElement(normal_stiffness, shear_stiffness);
    const auto expected_up_block_matrix = Matrix{number_of_u_dofs, number_of_pw_dofs, 0.0};
    const auto expected_pu_block_matrix = Matrix{number_of_pw_dofs, number_of_u_dofs, 0.0};
    const auto expected_pp_block_matrix = Matrix{number_of_pw_dofs, number_of_pw_dofs, 0.0};
    AssertLHSMatrixBlocksAreNear(actual_left_hand_side, expected_uu_block_matrix,
                                 expected_up_block_matrix, expected_pu_block_matrix,
                                 expected_pp_block_matrix, number_of_u_dofs, number_of_pw_dofs);
}

KRATOS_TEST_CASE_IN_SUITE(UPwLineInterfaceElement_LeftHandSideContainsMaterialStiffnessContributions_Rotated,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    constexpr auto normal_stiffness = 20.0;
    constexpr auto shear_stiffness  = 10.0;
    const auto     p_properties =
        CreateElasticMaterialProperties<InterfacePlaneStrain>(normal_stiffness, shear_stiffness);

    auto element = CreateAndInitializeElement(
        CreateUnitLengthLineInterfaceElementRotatedBy30DegreesWithUPwDoF, p_properties,
        IsDiffOrderElement::No, {CalculationContribution::Stiffness});

    // Act
    Matrix actual_left_hand_side;
    element.CalculateLeftHandSide(actual_left_hand_side, ProcessInfo{});

    // Assert
    constexpr auto number_of_u_dofs  = std::size_t{4 * 2};
    constexpr auto number_of_pw_dofs = std::size_t{4};
    ASSERT_EQ(actual_left_hand_side.size1(), number_of_u_dofs + number_of_pw_dofs);
    ASSERT_EQ(actual_left_hand_side.size2(), number_of_u_dofs + number_of_pw_dofs);

    auto expected_stiffness_matrix = Matrix{number_of_u_dofs, number_of_u_dofs};
    // clang-format off
    expected_stiffness_matrix <<=
         6.25,       -2.16506351,  0.0,         0.0,        -6.25,        2.16506351,  0.0,         0.0,
        -2.16506351,  8.75,        0.0,         0.0,         2.16506351, -8.75,        0.0,         0.0,
         0.0,         0.0,         6.25,       -2.16506351,  0.0,         0.0,        -6.25,        2.16506351,
         0.0,         0.0,        -2.16506351,  8.75,        0.0,         0.0,         2.16506351, -8.75,
        -6.25,        2.16506351,  0.0,         0.0,         6.25,       -2.16506351,  0.0,         0.0,
         2.16506351, -8.75,        0.0,         0.0,        -2.16506351,  8.75,        0.0,         0.0,
         0.0,         0.0,        -6.25,        2.16506351,  0.0,         0.0,         6.25,       -2.16506351,
         0.0,         0.0,         2.16506351, -8.75,        0.0,         0.0,        -2.16506351,  8.75;
    // clang-format on
    KRATOS_EXPECT_MATRIX_RELATIVE_NEAR(
        subrange(actual_left_hand_side, 0, 0 + number_of_u_dofs, 0, 0 + number_of_u_dofs),
        expected_stiffness_matrix, Defaults::relative_tolerance)
    const auto expected_up_block_matrix = Matrix{number_of_u_dofs, number_of_pw_dofs, 0.0};
    AssertUPBlockMatrixIsNear(actual_left_hand_side, expected_up_block_matrix, number_of_u_dofs, number_of_pw_dofs);
    const auto expected_pu_block_matrix = Matrix{number_of_pw_dofs, number_of_u_dofs, 0.0};
    AssertPUBlockMatrixIsNear(actual_left_hand_side, expected_pu_block_matrix, number_of_u_dofs, number_of_pw_dofs);
    const auto expected_pp_block_matrix = Matrix{number_of_pw_dofs, number_of_pw_dofs, 0.0};
    AssertPPBlockMatrixIsNear(actual_left_hand_side, expected_pp_block_matrix, number_of_u_dofs, number_of_pw_dofs);
}

KRATOS_TEST_CASE_IN_SUITE(UPwLineInterfaceElement_RightHandSideEqualsMinusInternalForceVector,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    constexpr auto normal_stiffness = 20.0;
    constexpr auto shear_stiffness  = 10.0;
    const auto     p_properties =
        CreateElasticMaterialProperties<InterfacePlaneStrain>(normal_stiffness, shear_stiffness);

    const auto prescribed_displacements = PrescribedDisplacements{
        {2, array_1d<double, 3>{0.2, 0.5, 0.0}}, {3, array_1d<double, 3>{0.2, 0.5, 0.0}}};
    auto element = CreateAndInitializeElement(
        CreateHorizontalUnitLength2Plus2NodedLineInterfaceElementWithUPwDofs, p_properties,
        IsDiffOrderElement::No, {CalculationContribution::Stiffness}, prescribed_displacements);

    // Act
    Vector actual_right_hand_side;
    element.CalculateRightHandSide(actual_right_hand_side, ProcessInfo{});

    // Assert
    constexpr auto number_of_u_dofs  = std::size_t{4 * 2};
    constexpr auto number_of_pw_dofs = std::size_t{4};
    const auto     expected_stiffness_force =
        UblasUtilities::CreateVector({1.0, 5.0, 1.0, 5.0, -1.0, -5.0, -1.0, -5.0});
    const auto expected_p_block_vector = Vector{number_of_pw_dofs, 0.0};
    AssertRHSVectorBlocksAreNear(actual_right_hand_side, expected_stiffness_force,
                                 expected_p_block_vector, number_of_u_dofs, number_of_pw_dofs);

    // Act
    Vector actual_external_forces_vector;
    element.Calculate(EXTERNAL_FORCES_VECTOR, actual_external_forces_vector, ProcessInfo{});

    // Assert
    auto expected_u_block_vector = Vector{number_of_u_dofs, 0.0};
    AssertRHSVectorBlocksAreNear(actual_external_forces_vector, expected_u_block_vector,
                                 expected_p_block_vector, number_of_u_dofs, number_of_pw_dofs);

    // Act
    Vector actual_internal_forces_vector;
    element.Calculate(INTERNAL_FORCES_VECTOR, actual_internal_forces_vector, ProcessInfo{});

    // Assert
    expected_u_block_vector = Vector{-1.0 * expected_stiffness_force};
    AssertRHSVectorBlocksAreNear(actual_internal_forces_vector, expected_u_block_vector,
                                 expected_p_block_vector, number_of_u_dofs, number_of_pw_dofs);
}

KRATOS_TEST_CASE_IN_SUITE(UPwLineInterfaceElement_RightHandSideEqualsMinusInternalForceVector_Rotated,
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
    auto element = CreateAndInitializeElement(
        CreateUnitLengthLineInterfaceElementRotatedBy30DegreesWithUPwDoF, p_properties,
        IsDiffOrderElement::No, {CalculationContribution::Stiffness}, prescribed_displacements);

    // Act
    Vector actual_right_hand_side;
    element.CalculateRightHandSide(actual_right_hand_side, ProcessInfo{});

    // Assert
    constexpr auto number_of_u_dofs  = std::size_t{4 * 2};
    constexpr auto number_of_pw_dofs = std::size_t{4};
    ASSERT_EQ(actual_right_hand_side.size(), number_of_u_dofs + number_of_pw_dofs);

    auto expected_stiffness_force = Vector{number_of_u_dofs};
    // clang-format off
    expected_stiffness_force <<= -1.6339746,  4.83012702, -1.6339746,  4.83012702,
                                  1.6339746, -4.83012702,  1.6339746, -4.83012702;
    // clang-format on
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(subrange(actual_right_hand_side, 0, 0 + number_of_u_dofs),
                                       expected_stiffness_force, Defaults::relative_tolerance)
    const auto expected_p_block_vector = Vector{number_of_pw_dofs, 0.0};
    AssertPBlockVectorIsNear(actual_right_hand_side, expected_p_block_vector, number_of_u_dofs, number_of_pw_dofs);
}

KRATOS_TEST_CASE_IN_SUITE(GetInitializedConstitutiveLawsAfterElementInitialization, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto p_properties = std::make_shared<Properties>();
    p_properties->GetValue(CONSTITUTIVE_LAW) =
        std::make_shared<GeoIncrementalLinearElasticInterfaceLaw>(std::make_unique<InterfacePlaneStrain>());

    auto element = CreateAndInitializeElement(
        CreateHorizontalUnitLength2Plus2NodedLineInterfaceElementWithUPwDofs, p_properties,
        IsDiffOrderElement::No, {CalculationContribution::Stiffness});
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

KRATOS_TEST_CASE_IN_SUITE(UPwInterfaceElement_HasCorrectNumberOfConstitutiveLawsAfterMultipleInitializations,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto p_properties = std::make_shared<Properties>();
    p_properties->GetValue(CONSTITUTIVE_LAW) =
        std::make_shared<GeoIncrementalLinearElasticInterfaceLaw>(std::make_unique<InterfacePlaneStrain>());

    Model model;
    auto  element = CreateHorizontalUnitLength2Plus2NodedLineInterfaceElementWithUPwDofs(
        model, p_properties, IsDiffOrderElement::No, {CalculationContribution::Stiffness});

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

KRATOS_TEST_CASE_IN_SUITE(UPwLineInterfaceElement_ReturnsExpectedLeftAndRightHandSide, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    constexpr auto normal_stiffness = 20.0;
    constexpr auto shear_stiffness  = 10.0;
    const auto     p_properties =
        CreateElasticMaterialProperties<InterfacePlaneStrain>(normal_stiffness, shear_stiffness);

    const auto prescribed_displacements = PrescribedDisplacements{
        {2, array_1d<double, 3>{0.2, 0.5, 0.0}}, {3, array_1d<double, 3>{0.2, 0.5, 0.0}}};
    auto element = CreateAndInitializeElement(
        CreateHorizontalUnitLength2Plus2NodedLineInterfaceElementWithUPwDofs, p_properties,
        IsDiffOrderElement::No, {CalculationContribution::Stiffness}, prescribed_displacements);

    // Act
    Vector actual_right_hand_side;
    Matrix actual_left_hand_side;
    element.CalculateLocalSystem(actual_left_hand_side, actual_right_hand_side, ProcessInfo{});

    // Assert
    constexpr auto number_of_u_dofs  = std::size_t{4 * 2};
    constexpr auto number_of_pw_dofs = std::size_t{4};
    const auto     expected_uu_block_matrix =
        CreateExpectedStiffnessMatrixForHorizontal2Plus2NodedElement(normal_stiffness, shear_stiffness);
    const auto expected_up_block_matrix = Matrix{number_of_u_dofs, number_of_pw_dofs, 0.0};
    const auto expected_pu_block_matrix = Matrix{number_of_pw_dofs, number_of_u_dofs, 0.0};
    const auto expected_pp_block_matrix = Matrix{number_of_pw_dofs, number_of_pw_dofs, 0.0};
    AssertLHSMatrixBlocksAreNear(actual_left_hand_side, expected_uu_block_matrix,
                                 expected_up_block_matrix, expected_pu_block_matrix,
                                 expected_pp_block_matrix, number_of_u_dofs, number_of_pw_dofs);

    const auto expected_u_block_vector =
        UblasUtilities::CreateVector({1.0, 5.0, 1.0, 5.0, -1.0, -5.0, -1.0, -5.0});
    const auto expected_p_block_vector = Vector{number_of_pw_dofs, 0.0};
    AssertRHSVectorBlocksAreNear(actual_right_hand_side, expected_u_block_vector,
                                 expected_p_block_vector, number_of_u_dofs, number_of_pw_dofs);
}

KRATOS_TEST_CASE_IN_SUITE(UPwLineInterfaceElement_CalculateRelativeDisplacementVector, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    constexpr auto normal_stiffness = 20.0;
    constexpr auto shear_stiffness  = 10.0;
    const auto     p_properties =
        CreateElasticMaterialProperties<InterfacePlaneStrain>(normal_stiffness, shear_stiffness);

    const auto prescribed_displacements =
        PrescribedDisplacements{{2, array_1d<double, 3>{-0.07679492, 0.5330127, 0.0}},
                                {3, array_1d<double, 3>{-0.07679492, 0.5330127, 0.0}}};
    auto element = CreateAndInitializeElement(
        CreateUnitLengthLineInterfaceElementRotatedBy30DegreesWithUPwDoF, p_properties,
        IsDiffOrderElement::No, {CalculationContribution::Stiffness}, prescribed_displacements);

    // Act
    std::vector<Vector> relative_displacements_at_integration_points;
    element.CalculateOnIntegrationPoints(GEO_RELATIVE_DISPLACEMENT_VECTOR,
                                         relative_displacements_at_integration_points, ProcessInfo{});

    // Assert
    Vector expected_relative_displacement{2};
    expected_relative_displacement <<= 0.5, 0.2;
    KRATOS_EXPECT_EQ(relative_displacements_at_integration_points.size(), 2);
    for (const auto& r_relative_displacement : relative_displacements_at_integration_points) {
        KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(r_relative_displacement, expected_relative_displacement,
                                           Defaults::relative_tolerance)
    }
}

KRATOS_TEST_CASE_IN_SUITE(UPwLineInterfaceElement_CalculateEffectiveTractionVector, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    constexpr auto normal_stiffness = 20.0;
    constexpr auto shear_stiffness  = 10.0;
    const auto     p_properties =
        CreateElasticMaterialProperties<InterfacePlaneStrain>(normal_stiffness, shear_stiffness);

    const auto prescribed_displacements =
        PrescribedDisplacements{{2, array_1d<double, 3>{-0.07679492, 0.5330127, 0.0}},
                                {3, array_1d<double, 3>{-0.07679492, 0.5330127, 0.0}}};
    auto element = CreateAndInitializeElement(
        CreateUnitLengthLineInterfaceElementRotatedBy30DegreesWithUPwDoF, p_properties,
        IsDiffOrderElement::No, {CalculationContribution::Stiffness}, prescribed_displacements);

    // Act
    std::vector<Vector> tractions_at_integration_points;
    element.CalculateOnIntegrationPoints(GEO_EFFECTIVE_TRACTION_VECTOR,
                                         tractions_at_integration_points, ProcessInfo{});

    // Assert
    Vector expected_traction{2};
    expected_traction <<= 10.0, 2.0;
    KRATOS_EXPECT_EQ(tractions_at_integration_points.size(), 2);
    for (const auto& r_traction : tractions_at_integration_points) {
        KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(r_traction, expected_traction, Defaults::relative_tolerance)
    }
}

KRATOS_TEST_CASE_IN_SUITE(UPwLineInterfaceElement_3Plus3NodedElement_ReturnsExpectedLeftAndRightHandSide,
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
    auto element = CreateAndInitializeElement(
        CreateHorizontalUnitLength3Plus3NodedLineInterfaceElementWithUPwDoF, p_properties,
        IsDiffOrderElement::No, {CalculationContribution::Stiffness}, prescribed_displacements);

    // Act
    Vector actual_right_hand_side;
    Matrix actual_left_hand_side;
    element.CalculateLocalSystem(actual_left_hand_side, actual_right_hand_side, ProcessInfo{});

    // Assert
    constexpr auto number_of_u_dofs  = std::size_t{6 * 2};
    constexpr auto number_of_pw_dofs = std::size_t{6};
    ASSERT_EQ(actual_left_hand_side.size1(), number_of_u_dofs + number_of_pw_dofs);
    ASSERT_EQ(actual_left_hand_side.size2(), number_of_u_dofs + number_of_pw_dofs);
    ASSERT_EQ(actual_right_hand_side.size(), number_of_u_dofs + number_of_pw_dofs);

    constexpr auto cs                       = shear_stiffness / 6.0;
    constexpr auto cn                       = normal_stiffness / 6.0;
    const auto     expected_uu_block_matrix = UblasUtilities::CreateMatrix(
        {{cs, 0.0, 0.0, 0.0, 0.0, 0.0, -cs, 0.0, 0.0, 0.0, 0.0, 0.0},
             {0.0, cn, 0.0, 0.0, 0.0, 0.0, 0.0, -cn, 0.0, 0.0, 0.0, 0.0},
             {0.0, 0.0, cs, 0.0, 0.0, 0.0, 0.0, 0.0, -cs, 0.0, 0.0, 0.0},
             {0.0, 0.0, 0.0, cn, 0.0, 0.0, 0.0, 0.0, 0.0, -cn, 0.0, 0.0},
             {0.0, 0.0, 0.0, 0.0, 4.0 * cs, 0.0, 0.0, 0.0, 0.0, 0.0, -4.0 * cs, 0.0},
             {0.0, 0.0, 0.0, 0.0, 0.0, 4.0 * cn, 0.0, 0.0, 0.0, 0.0, 0.0, -4.0 * cn},
             {-cs, 0.0, 0.0, 0.0, 0.0, 0.0, cs, 0.0, 0.0, 0.0, 0.0, 0.0},
             {0.0, -cn, 0.0, 0.0, 0.0, 0.0, 0.0, cn, 0.0, 0.0, 0.0, 0.0},
             {0.0, 0.0, -cs, 0.0, 0.0, 0.0, 0.0, 0.0, cs, 0.0, 0.0, 0.0},
             {0.0, 0.0, 0.0, -cn, 0.0, 0.0, 0.0, 0.0, 0.0, cn, 0.0, 0.0},
             {0.0, 0.0, 0.0, 0.0, -4.0 * cs, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0 * cs, 0.0},
             {0.0, 0.0, 0.0, 0.0, 0.0, -4.0 * cn, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0 * cn}});
    const auto expected_up_block_matrix = Matrix{number_of_u_dofs, number_of_pw_dofs, 0.0};
    const auto expected_pu_block_matrix = Matrix{number_of_pw_dofs, number_of_u_dofs, 0.0};
    const auto expected_pp_block_matrix = Matrix{number_of_pw_dofs, number_of_pw_dofs, 0.0};
    AssertLHSMatrixBlocksAreNear(actual_left_hand_side, expected_uu_block_matrix,
                                 expected_up_block_matrix, expected_pu_block_matrix,
                                 expected_pp_block_matrix, number_of_u_dofs, number_of_pw_dofs);

    const auto expected_u_block_vector = UblasUtilities::CreateVector(
        {1.0 / 3.0, 5.0 / 3.0, 1.0 / 3.0, 5.0 / 3.0, 4.0 / 3.0, 20.0 / 3.0, -1.0 / 3.0, -5.0 / 3.0,
         -1.0 / 3.0, -5.0 / 3.0, -4.0 / 3.0, -20.0 / 3.0});
    const auto expected_p_block_vector = Vector{number_of_pw_dofs, 0.0};
    AssertRHSVectorBlocksAreNear(actual_right_hand_side, expected_u_block_vector,
                                 expected_p_block_vector, number_of_u_dofs, number_of_pw_dofs);
}

KRATOS_TEST_CASE_IN_SUITE(UPwTriangleInterfaceElement_CreatesInstanceWithGeometryInput,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const UPwInterfaceElement element(
        0, std::make_shared<TriangleInterfaceGeometry3D3Plus3Noded>(CreateNodesFor3Plus3SurfaceInterfaceGeometry()),
        std::make_unique<SurfaceInterfaceStressState>(), IsDiffOrderElement::No,
        {CalculationContribution::Stiffness});
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

KRATOS_TEST_CASE_IN_SUITE(UPwTriangleInterfaceElement_CreatesInstanceWithNodeInput, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto nodes        = CreateNodesFor3Plus3SurfaceInterfaceGeometry();
    const auto p_properties = std::make_shared<Properties>();

    // The source element needs to have a geometry, otherwise the version of the
    // Create method with a node input will fail.
    const auto p_geometry = std::make_shared<TriangleInterfaceGeometry3D3Plus3Noded>(nodes);
    const UPwInterfaceElement element(0, p_geometry, p_properties,
                                      std::make_unique<SurfaceInterfaceStressState>(),
                                      IsDiffOrderElement::No, {CalculationContribution::Stiffness});

    // Act
    const auto p_created_element = element.Create(1, nodes, p_properties);

    // Assert
    EXPECT_NE(p_created_element, nullptr);
    EXPECT_EQ(p_created_element->Id(), 1);
    EXPECT_NE(p_created_element->pGetGeometry(), nullptr);
    EXPECT_NE(p_created_element->pGetProperties(), nullptr);
}

KRATOS_TEST_CASE_IN_SUITE(UPwTriangleInterfaceElement_ReturnsTheExpectedDoFList, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto p_properties = std::make_shared<Properties>();

    Model      model;
    const auto element = CreateHorizontal3Plus3NodedTriangleInterfaceElementWithUPwDofs(
        model, p_properties, IsDiffOrderElement::No, {CalculationContribution::Stiffness});

    // Act
    const auto              dummy_process_info = ProcessInfo{};
    Element::DofsVectorType degrees_of_freedom;
    element.GetDofList(degrees_of_freedom, dummy_process_info);

    // Assert
    constexpr auto number_of_u_dofs  = std::size_t{6 * 3};
    constexpr auto number_of_pw_dofs = std::size_t{6};
    ASSERT_EQ(degrees_of_freedom.size(), number_of_u_dofs + number_of_pw_dofs);
    for (auto i = std::size_t{0}; i < number_of_u_dofs; i += 3) {
        KRATOS_EXPECT_EQ(degrees_of_freedom[i]->GetVariable(), DISPLACEMENT_X);
        KRATOS_EXPECT_EQ(degrees_of_freedom[i + 1]->GetVariable(), DISPLACEMENT_Y);
        KRATOS_EXPECT_EQ(degrees_of_freedom[i + 2]->GetVariable(), DISPLACEMENT_Z);
    }
    for (auto i = std::size_t{0}; i < number_of_pw_dofs; ++i) {
        KRATOS_EXPECT_EQ(degrees_of_freedom[i + number_of_u_dofs]->GetVariable(), WATER_PRESSURE);
    }
}

KRATOS_TEST_CASE_IN_SUITE(UPwTriangleInterfaceElement_ReturnsTheExpectedEquationIdVector,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto p_properties = std::make_shared<Properties>();

    Model model;
    auto  element = CreateHorizontal3Plus3NodedTriangleInterfaceElementWithUPwDofs(
        model, p_properties, IsDiffOrderElement::No, {CalculationContribution::Stiffness});

    int i = 0;
    for (const auto& node : element.GetGeometry()) {
        ++i;
        node.pGetDof(DISPLACEMENT_X)->SetEquationId(i);

        ++i;
        node.pGetDof(DISPLACEMENT_Y)->SetEquationId(i);

        ++i;
        node.pGetDof(DISPLACEMENT_Z)->SetEquationId(i);

        ++i;
        node.pGetDof(WATER_PRESSURE)->SetEquationId(i);
    }

    // Act
    const auto                    dummy_process_info = ProcessInfo{};
    Element::EquationIdVectorType equation_id_vector;
    element.EquationIdVector(equation_id_vector, dummy_process_info);

    // Assert
    const Element::EquationIdVectorType expected_ids = {
        1, 2, 3, 5, 6, 7, 9, 10, 11, 13, 14, 15, 17, 18, 19, 21, 22, 23, 4, 8, 12, 16, 20, 24};
    KRATOS_EXPECT_VECTOR_EQ(equation_id_vector, expected_ids)
}

KRATOS_TEST_CASE_IN_SUITE(UPwTriangleInterfaceElement_LeftHandSideContainsMaterialStiffnessContributions,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    constexpr auto normal_stiffness = 20.0;
    constexpr auto shear_stiffness  = 10.0;
    const auto     p_properties = CreateElasticMaterialProperties<InterfaceThreeDimensionalSurface>(
        normal_stiffness, shear_stiffness);

    auto element = CreateAndInitializeElement(
        CreateHorizontal3Plus3NodedTriangleInterfaceElementWithUPwDofs, p_properties,
        IsDiffOrderElement::No, {CalculationContribution::Stiffness});

    // Act
    Matrix actual_left_hand_side;
    element.CalculateLeftHandSide(actual_left_hand_side, ProcessInfo{});

    // Assert
    constexpr auto number_of_u_dofs  = std::size_t{6 * 3};
    constexpr auto number_of_pw_dofs = std::size_t{6};
    const auto expected_uu_block_matrix = ExpectedStiffnessMatrixOfLinearTriangularInterfaceElement();
    const auto expected_up_block_matrix = Matrix{number_of_u_dofs, number_of_pw_dofs, 0.0};
    const auto expected_pu_block_matrix = Matrix{number_of_pw_dofs, number_of_u_dofs, 0.0};
    const auto expected_pp_block_matrix = Matrix{number_of_pw_dofs, number_of_pw_dofs, 0.0};
    AssertLHSMatrixBlocksAreNear(actual_left_hand_side, expected_uu_block_matrix,
                                 expected_up_block_matrix, expected_pu_block_matrix,
                                 expected_pp_block_matrix, number_of_u_dofs, number_of_pw_dofs);
}

KRATOS_TEST_CASE_IN_SUITE(UPwTriangleInterfaceElement_LeftHandSideContainsMaterialStiffnessContributions_Rotated,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    constexpr auto normal_stiffness = 20.0;
    constexpr auto shear_stiffness  = 10.0;
    const auto     p_properties = CreateElasticMaterialProperties<InterfaceThreeDimensionalSurface>(
        normal_stiffness, shear_stiffness);

    auto element = CreateAndInitializeElement(CreateTriangleInterfaceElementRotatedBy30DegreesWithUPwDoF,
                                              p_properties, IsDiffOrderElement::No,
                                              {CalculationContribution::Stiffness});

    // Act
    Matrix actual_left_hand_side;
    element.CalculateLeftHandSide(actual_left_hand_side, ProcessInfo{});

    // Assert
    constexpr auto number_of_u_dofs  = std::size_t{6 * 3};
    constexpr auto number_of_pw_dofs = std::size_t{6};
    // Since the rotation is about the z-axis (the normal of the triangle) and the two shear
    // stiffnesses are equal, the stiffness matrix is equal to the stiffness matrix of a
    // non-rotated surface interface element.
    const auto expected_uu_block_matrix = ExpectedStiffnessMatrixOfLinearTriangularInterfaceElement();
    const auto expected_up_block_matrix = Matrix{number_of_u_dofs, number_of_pw_dofs, 0.0};
    const auto expected_pu_block_matrix = Matrix{number_of_pw_dofs, number_of_u_dofs, 0.0};
    const auto expected_pp_block_matrix = Matrix{number_of_pw_dofs, number_of_pw_dofs, 0.0};
    AssertLHSMatrixBlocksAreNear(actual_left_hand_side, expected_uu_block_matrix,
                                 expected_up_block_matrix, expected_pu_block_matrix,
                                 expected_pp_block_matrix, number_of_u_dofs, number_of_pw_dofs);
}

KRATOS_TEST_CASE_IN_SUITE(UPwTriangleInterfaceElement_LeftHandSideContainsMaterialStiffnessContributions_RotatedAboutYAxis,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    constexpr auto normal_stiffness = 20.0;
    constexpr auto shear_stiffness  = 10.0;
    const auto     p_properties = CreateElasticMaterialProperties<InterfaceThreeDimensionalSurface>(
        normal_stiffness, shear_stiffness);

    auto element = CreateAndInitializeElement(
        CreateTriangleInterfaceElementRotatedBy30DegreesAboutYAxisWithUPwDoF, p_properties,
        IsDiffOrderElement::No, {CalculationContribution::Stiffness});

    // Act
    Matrix actual_left_hand_side;
    element.CalculateLeftHandSide(actual_left_hand_side, ProcessInfo{});

    // Assert
    constexpr auto number_of_u_dofs  = std::size_t{6 * 3};
    constexpr auto number_of_pw_dofs = std::size_t{6};
    ASSERT_EQ(actual_left_hand_side.size1(), number_of_u_dofs + number_of_pw_dofs);
    ASSERT_EQ(actual_left_hand_side.size2(), number_of_u_dofs + number_of_pw_dofs);

    const auto expected_uu_block_matrix = UblasUtilities::CreateMatrix(
        {{2.083333333333333, 0, 0.72168783648703216, 0, 0, 0, 0, 0, 0, -2.083333333333333, 0,
          -0.72168783648703216, 0, 0, 0, 0, 0, 0},
         {0, 1.6666666666666661, 0, 0, 0, 0, 0, 0, 0, 0, -1.6666666666666661, 0, 0, 0, 0, 0, 0, 0},
         {0.72168783648703216, 0, 2.9166666666666661, 0, 0, 0, 0, 0, 0, -0.72168783648703216, 0,
          -2.9166666666666661, 0, 0, 0, 0, 0, 0},
         {0, 0, 0, 2.083333333333333, 0, 0.72168783648703216, 0, 0, 0, 0, 0, 0, -2.083333333333333,
          0, -0.72168783648703216, 0, 0, 0},
         {0, 0, 0, 0, 1.6666666666666661, 0, 0, 0, 0, 0, 0, 0, 0, -1.6666666666666661, 0, 0, 0, 0},
         {0, 0, 0, 0.72168783648703216, 0, 2.9166666666666661, 0, 0, 0, 0, 0, 0,
          -0.72168783648703216, 0, -2.9166666666666661, 0, 0, 0},
         {0, 0, 0, 0, 0, 0, 2.083333333333333, 0, 0.72168783648703216, 0, 0, 0, 0, 0, 0,
          -2.083333333333333, 0, -0.72168783648703216},
         {0, 0, 0, 0, 0, 0, 0, 1.6666666666666661, 0, 0, 0, 0, 0, 0, 0, 0, -1.6666666666666661, 0},
         {0, 0, 0, 0, 0, 0, 0.72168783648703216, 0, 2.9166666666666661, 0, 0, 0, 0, 0, 0,
          -0.72168783648703216, 0, -2.9166666666666661},
         {-2.083333333333333, 0, -0.72168783648703216, 0, 0, 0, 0, 0, 0, 2.083333333333333, 0,
          0.72168783648703216, 0, 0, 0, 0, 0, 0},
         {0, -1.6666666666666661, 0, 0, 0, 0, 0, 0, 0, 0, 1.6666666666666661, 0, 0, 0, 0, 0, 0, 0},
         {-0.72168783648703216, 0, -2.9166666666666661, 0, 0, 0, 0, 0, 0, 0.72168783648703216, 0,
          2.9166666666666661, 0, 0, 0, 0, 0, 0},
         {0, 0, 0, -2.083333333333333, 0, -0.72168783648703216, 0, 0, 0, 0, 0, 0, 2.083333333333333,
          0, 0.72168783648703216, 0, 0, 0},
         {0, 0, 0, 0, -1.6666666666666661, 0, 0, 0, 0, 0, 0, 0, 0, 1.6666666666666661, 0, 0, 0, 0},
         {0, 0, 0, -0.72168783648703216, 0, -2.9166666666666661, 0, 0, 0, 0, 0, 0,
          0.72168783648703216, 0, 2.9166666666666661, 0, 0, 0},
         {0, 0, 0, 0, 0, 0, -2.083333333333333, 0, -0.72168783648703216, 0, 0, 0, 0, 0, 0,
          2.083333333333333, 0, 0.72168783648703216},
         {0, 0, 0, 0, 0, 0, 0, -1.6666666666666661, 0, 0, 0, 0, 0, 0, 0, 0, 1.6666666666666661, 0},
         {0, 0, 0, 0, 0, 0, -0.72168783648703216, 0, -2.9166666666666661, 0, 0, 0, 0, 0, 0,
          0.72168783648703216, 0, 2.9166666666666661}});
    KRATOS_EXPECT_MATRIX_RELATIVE_NEAR(
        subrange(actual_left_hand_side, 0, 0 + number_of_u_dofs, 0, 0 + number_of_u_dofs),
        expected_uu_block_matrix, Defaults::relative_tolerance)
    const auto expected_up_block_matrix = Matrix{number_of_u_dofs, number_of_pw_dofs, 0.0};
    AssertUPBlockMatrixIsNear(actual_left_hand_side, expected_up_block_matrix, number_of_u_dofs, number_of_pw_dofs);
    const auto expected_pu_block_matrix = Matrix{number_of_pw_dofs, number_of_u_dofs, 0.0};
    AssertPUBlockMatrixIsNear(actual_left_hand_side, expected_pu_block_matrix, number_of_u_dofs, number_of_pw_dofs);
    const auto expected_pp_block_matrix = Matrix{number_of_pw_dofs, number_of_pw_dofs, 0.0};
    AssertPPBlockMatrixIsNear(actual_left_hand_side, expected_pp_block_matrix, number_of_u_dofs, number_of_pw_dofs);
}

KRATOS_TEST_CASE_IN_SUITE(UPwTriangleInterfaceElement_RightHandSideEqualsMinusInternalForceVector,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    constexpr auto normal_stiffness = 20.0;
    constexpr auto shear_stiffness  = 10.0;
    const auto     p_properties = CreateElasticMaterialProperties<InterfaceThreeDimensionalSurface>(
        normal_stiffness, shear_stiffness);

    const auto prescribed_displacements = PrescribedDisplacements{
        {2, array_1d<double, 3>{0.2, 0.5, 0.0}}, {3, array_1d<double, 3>{0.2, 0.5, 0.0}}};
    auto element = CreateAndInitializeElement(
        CreateHorizontal3Plus3NodedTriangleInterfaceElementWithUPwDofs, p_properties,
        IsDiffOrderElement::No, {CalculationContribution::Stiffness}, prescribed_displacements);

    // Act
    Vector actual_right_hand_side;
    element.CalculateRightHandSide(actual_right_hand_side, ProcessInfo{});

    // Assert
    constexpr auto number_of_u_dofs    = std::size_t{6 * 3};
    constexpr auto number_of_pw_dofs   = std::size_t{6};
    const auto expected_u_block_vector = ExpectedStiffnessForceOfLinearTriangularInterfaceElement();
    const auto expected_p_block_vector = Vector{number_of_pw_dofs, 0.0};
    AssertRHSVectorBlocksAreNear(actual_right_hand_side, expected_u_block_vector,
                                 expected_p_block_vector, number_of_u_dofs, number_of_pw_dofs);
}

KRATOS_TEST_CASE_IN_SUITE(UPwTriangleInterfaceElement_RightHandSideEqualsMinusInternalForceVector_Rotated,
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
    auto element = CreateAndInitializeElement(
        CreateTriangleInterfaceElementRotatedBy30DegreesWithUPwDoF, p_properties,
        IsDiffOrderElement::No, {CalculationContribution::Stiffness}, prescribed_displacements);

    // Act
    Vector actual_right_hand_side;
    element.CalculateRightHandSide(actual_right_hand_side, ProcessInfo{});

    // Assert
    constexpr auto number_of_u_dofs  = std::size_t{6 * 3};
    constexpr auto number_of_pw_dofs = std::size_t{6};
    ASSERT_EQ(actual_right_hand_side.size(), number_of_u_dofs + number_of_pw_dofs);

    const auto expected_u_block_vector = UblasUtilities::CreateVector(
        {1.66667, 3.33333, 10, 1.66667, 3.33333, 10, 1.66667, 3.33333, 10, -1.66667, -3.33333, -10,
         -1.66667, -3.33333, -10, -1.66667, -3.33333, -10});
    constexpr auto tolerance = 1e-5;
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(subrange(actual_right_hand_side, 0, 0 + number_of_u_dofs),
                                       expected_u_block_vector, tolerance)
    const auto expected_p_block_vector = Vector{number_of_pw_dofs, 0.0};
    AssertPBlockVectorIsNear(actual_right_hand_side, expected_p_block_vector, number_of_u_dofs, number_of_pw_dofs);
}

KRATOS_TEST_CASE_IN_SUITE(UPwTriangleInterfaceElement_GetInitializedConstitutiveLawsAfterElementInitialization,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto p_properties = std::make_shared<Properties>();
    p_properties->GetValue(CONSTITUTIVE_LAW) =
        std::make_shared<GeoIncrementalLinearElasticInterfaceLaw>(std::make_unique<InterfacePlaneStrain>());

    auto element = CreateAndInitializeElement(
        CreateHorizontal3Plus3NodedTriangleInterfaceElementWithUPwDofs, p_properties,
        IsDiffOrderElement::No, {CalculationContribution::Stiffness});

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

KRATOS_TEST_CASE_IN_SUITE(UPwTriangleInterfaceElement_HasCorrectNumberOfConstitutiveLawsAfterMultipleInitializations,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto p_properties = std::make_shared<Properties>();
    p_properties->GetValue(CONSTITUTIVE_LAW) =
        std::make_shared<GeoIncrementalLinearElasticInterfaceLaw>(std::make_unique<InterfacePlaneStrain>());

    Model model;
    auto  element = CreateHorizontal3Plus3NodedTriangleInterfaceElementWithUPwDofs(
        model, p_properties, IsDiffOrderElement::No, {CalculationContribution::Stiffness});

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

KRATOS_TEST_CASE_IN_SUITE(UPwTriangleInterfaceElement_ReturnsExpectedLeftAndRightHandSide,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    constexpr auto normal_stiffness = 20.0;
    constexpr auto shear_stiffness  = 10.0;
    const auto     p_properties = CreateElasticMaterialProperties<InterfaceThreeDimensionalSurface>(
        normal_stiffness, shear_stiffness);

    const auto prescribed_displacements = PrescribedDisplacements{
        {2, array_1d<double, 3>{0.2, 0.5, 0.0}}, {3, array_1d<double, 3>{0.2, 0.5, 0.0}}};
    auto element = CreateAndInitializeElement(
        CreateHorizontal3Plus3NodedTriangleInterfaceElementWithUPwDofs, p_properties,
        IsDiffOrderElement::No, {CalculationContribution::Stiffness}, prescribed_displacements);

    // Act
    Vector actual_right_hand_side;
    Matrix actual_left_hand_side;
    element.CalculateLocalSystem(actual_left_hand_side, actual_right_hand_side, ProcessInfo{});

    // Assert
    constexpr auto number_of_u_dofs  = std::size_t{6 * 3};
    constexpr auto number_of_pw_dofs = std::size_t{6};
    const auto expected_uu_block_matrix = ExpectedStiffnessMatrixOfLinearTriangularInterfaceElement();
    const auto expected_up_block_matrix = Matrix{number_of_u_dofs, number_of_pw_dofs, 0.0};
    const auto expected_pu_block_matrix = Matrix{number_of_pw_dofs, number_of_u_dofs, 0.0};
    const auto expected_pp_block_matrix = Matrix{number_of_pw_dofs, number_of_pw_dofs, 0.0};
    AssertLHSMatrixBlocksAreNear(actual_left_hand_side, expected_uu_block_matrix,
                                 expected_up_block_matrix, expected_pu_block_matrix,
                                 expected_pp_block_matrix, number_of_u_dofs, number_of_pw_dofs);

    const auto expected_u_block_vector = ExpectedStiffnessForceOfLinearTriangularInterfaceElement();
    const auto expected_p_block_vector = Vector{number_of_pw_dofs, 0.0};
    AssertRHSVectorBlocksAreNear(actual_right_hand_side, expected_u_block_vector,
                                 expected_p_block_vector, number_of_u_dofs, number_of_pw_dofs);
}

KRATOS_TEST_CASE_IN_SUITE(UPwTriangleInterfaceElement_CalculateRelativeDisplacementVector,
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
    auto element = CreateAndInitializeElement(
        CreateTriangleInterfaceElementRotatedBy30DegreesWithUPwDoF, p_properties,
        IsDiffOrderElement::No, {CalculationContribution::Stiffness}, prescribed_displacements);

    // Act
    std::vector<Vector> relative_displacements_at_integration_points;
    element.CalculateOnIntegrationPoints(GEO_RELATIVE_DISPLACEMENT_VECTOR,
                                         relative_displacements_at_integration_points, ProcessInfo{});

    // Assert
    Vector expected_relative_displacement{3};
    expected_relative_displacement <<= -3.0, -0.5 * std::sqrt(3.0) - 1.0, 0.5 - std::sqrt(3);

    KRATOS_EXPECT_EQ(relative_displacements_at_integration_points.size(), 3);
    for (const auto& r_relative_displacement : relative_displacements_at_integration_points) {
        KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(r_relative_displacement, expected_relative_displacement,
                                           Defaults::relative_tolerance)
    }
}

KRATOS_TEST_CASE_IN_SUITE(UPwTriangleInterfaceElementHorizontal_CalculateRelativeDisplacementVector,
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
    auto element = CreateAndInitializeElement(
        CreateHorizontal3Plus3NodedTriangleInterfaceElementWithUPwDofs, p_properties,
        IsDiffOrderElement::No, {CalculationContribution::Stiffness}, prescribed_displacements);

    // Act
    std::vector<Vector> relative_displacements_at_integration_points;
    element.CalculateOnIntegrationPoints(GEO_RELATIVE_DISPLACEMENT_VECTOR,
                                         relative_displacements_at_integration_points, ProcessInfo{});

    // Assert
    Vector expected_relative_displacement{3};
    expected_relative_displacement <<= -1.0, 0.0, 0.0;

    KRATOS_EXPECT_EQ(relative_displacements_at_integration_points.size(), 3);
    for (const auto& r_relative_displacement : relative_displacements_at_integration_points) {
        KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(r_relative_displacement, expected_relative_displacement,
                                           Defaults::relative_tolerance)
    }
}

KRATOS_TEST_CASE_IN_SUITE(UPwTriangleInterfaceElementInYZPlane_CalculateRelativeDisplacementVector,
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
    auto element = CreateAndInitializeElement(
        CreateHorizontal3Plus3NodedTriangleInterfaceYZPlaneElementWithUPwDofs, p_properties,
        IsDiffOrderElement::No, {CalculationContribution::Stiffness}, prescribed_displacements);

    // Act
    std::vector<Vector> relative_displacements_at_integration_points;
    element.CalculateOnIntegrationPoints(GEO_RELATIVE_DISPLACEMENT_VECTOR,
                                         relative_displacements_at_integration_points, ProcessInfo{});

    // Assert
    Vector expected_relative_displacement{3};
    expected_relative_displacement <<= 1.0, 0.0, 0.0;

    KRATOS_EXPECT_EQ(relative_displacements_at_integration_points.size(), 3);
    for (const auto& r_relative_displacement : relative_displacements_at_integration_points) {
        KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(r_relative_displacement, expected_relative_displacement,
                                           Defaults::relative_tolerance)
    }
}

KRATOS_TEST_CASE_IN_SUITE(UPwTriangleInterfaceElementInXZPlane_CalculateRelativeDisplacementVector,
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
    auto element = CreateAndInitializeElement(
        CreateHorizontal3Plus3NodedTriangleInterfaceXZPlaneElementWithUPwDofs, p_properties,
        IsDiffOrderElement::No, {CalculationContribution::Stiffness}, prescribed_displacements);

    // Act
    std::vector<Vector> relative_displacements_at_integration_points;
    element.CalculateOnIntegrationPoints(GEO_RELATIVE_DISPLACEMENT_VECTOR,
                                         relative_displacements_at_integration_points, ProcessInfo{});

    // Assert
    Vector expected_relative_displacement{3};
    expected_relative_displacement <<= 1.0, 0.0, 0.0;

    KRATOS_EXPECT_EQ(relative_displacements_at_integration_points.size(), 3);
    for (const auto& r_relative_displacement : relative_displacements_at_integration_points) {
        KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(r_relative_displacement, expected_relative_displacement,
                                           Defaults::relative_tolerance)
    }
}

KRATOS_TEST_CASE_IN_SUITE(UPwTriangleInterfaceElement_CalculateEffectiveTractionVector,
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
    auto element = CreateAndInitializeElement(
        CreateTriangleInterfaceElementRotatedBy30DegreesWithUPwDoF, p_properties,
        IsDiffOrderElement::No, {CalculationContribution::Stiffness}, prescribed_displacements);

    // Act
    std::vector<Vector> tractions_at_integration_points;
    element.CalculateOnIntegrationPoints(GEO_EFFECTIVE_TRACTION_VECTOR,
                                         tractions_at_integration_points, ProcessInfo{});

    // Assert
    Vector expected_traction{3};
    expected_traction <<= -3.0 * normal_stiffness, (-0.5 * std::sqrt(3.0) - 1.0) * shear_stiffness,
        (0.5 - std::sqrt(3)) * shear_stiffness;
    KRATOS_EXPECT_EQ(tractions_at_integration_points.size(), 3);
    for (const auto& r_traction : tractions_at_integration_points) {
        KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(r_traction, expected_traction, Defaults::relative_tolerance)
    }
}

KRATOS_TEST_CASE_IN_SUITE(UPwTriangleInterfaceElement_6Plus6NodedElement_ReturnsExpectedLeftAndRightHandSide,
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
    auto element = CreateAndInitializeElement(
        CreateHorizontal6Plus6NodedTriangleInterfaceElementWithUPwDoF, p_properties,
        IsDiffOrderElement::No, {CalculationContribution::Stiffness}, prescribed_displacements);

    // Act
    Vector actual_right_hand_side;
    Matrix actual_left_hand_side;
    element.CalculateLocalSystem(actual_left_hand_side, actual_right_hand_side, ProcessInfo{});

    // Assert
    constexpr auto number_of_u_dofs  = std::size_t{12 * 3};
    constexpr auto number_of_pw_dofs = std::size_t{12};
    ASSERT_EQ(actual_left_hand_side.size1(), number_of_u_dofs + number_of_pw_dofs);
    ASSERT_EQ(actual_left_hand_side.size2(), number_of_u_dofs + number_of_pw_dofs);
    ASSERT_EQ(actual_right_hand_side.size(), number_of_u_dofs + number_of_pw_dofs);

    auto expected_uu_block_matrix =
        Matrix(number_of_u_dofs, number_of_u_dofs); // the values are taken from the element
    // clang-format off
    expected_uu_block_matrix <<=
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

    KRATOS_EXPECT_MATRIX_RELATIVE_NEAR(
        subrange(actual_left_hand_side, 0, 0 + number_of_u_dofs, 0, 0 + number_of_u_dofs),
        expected_uu_block_matrix, Defaults::relative_tolerance)
    const auto expected_up_block_matrix = Matrix{number_of_u_dofs, number_of_pw_dofs, 0.0};
    AssertUPBlockMatrixIsNear(actual_left_hand_side, expected_up_block_matrix, number_of_u_dofs, number_of_pw_dofs);
    const auto expected_pu_block_matrix = Matrix{number_of_pw_dofs, number_of_u_dofs, 0.0};
    AssertPUBlockMatrixIsNear(actual_left_hand_side, expected_pu_block_matrix, number_of_u_dofs, number_of_pw_dofs);
    const auto expected_pp_block_matrix = Matrix{number_of_pw_dofs, number_of_pw_dofs, 0.0};
    AssertPPBlockMatrixIsNear(actual_left_hand_side, expected_pp_block_matrix, number_of_u_dofs, number_of_pw_dofs);

    const auto expected_stiffness_force = UblasUtilities::CreateVector(
        {0,         0,         0, 0,         0,         0, 0,         0,         0,
         -0.300469, -0.751174, 0, -0.300469, -0.751174, 0, -0.300469, -0.751174, 0,
         0,         0,         0, 0,         0,         0, 0,         0,         0,
         0.300469,  0.751174,  0, 0.300469,  0.751174,  0, 0.300469,  0.751174,  0});
    constexpr auto tolerance = 1e-5;
    KRATOS_EXPECT_VECTOR_NEAR(subrange(actual_right_hand_side, 0, 0 + number_of_u_dofs),
                              expected_stiffness_force, tolerance)
    const auto expected_p_block_vector = Vector{number_of_pw_dofs, 0.0};
    AssertPBlockVectorIsNear(actual_right_hand_side, expected_p_block_vector, number_of_u_dofs, number_of_pw_dofs);
}

KRATOS_TEST_CASE_IN_SUITE(UPwInterfaceElement_CheckThrowsWhenElementIsNotInitialized, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    constexpr auto normal_stiffness = 20.0;
    constexpr auto shear_stiffness  = 10.0;
    const auto     p_properties =
        CreateElasticMaterialProperties<InterfacePlaneStrain>(normal_stiffness, shear_stiffness);

    Model model;
    auto  element = CreateHorizontalUnitLength3Plus3NodedLineInterfaceElementWithUPwDoF(
        model, p_properties, IsDiffOrderElement::No, {CalculationContribution::Stiffness});

    const auto dummy_process_info = ProcessInfo{};
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        element.Check(dummy_process_info),
        "Number of integration points (3) and constitutive laws (0) do not match.\n")

    element.Initialize(dummy_process_info);

    KRATOS_EXPECT_EQ(element.Check(dummy_process_info), 0);
}

KRATOS_TEST_CASE_IN_SUITE(UPwInterfaceElement_CheckDoesNotThrowWhenElementIsNotActive, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    constexpr auto normal_stiffness = 20.0;
    constexpr auto shear_stiffness  = 10.0;
    const auto     p_properties =
        CreateElasticMaterialProperties<InterfacePlaneStrain>(normal_stiffness, shear_stiffness);

    Model model;
    auto  element = CreateHorizontalUnitLength3Plus3NodedLineInterfaceElementWithUPwDoF(
        model, p_properties, IsDiffOrderElement::No, {CalculationContribution::Stiffness});

    // In the integrated workflow, the elements are not initialized, when they are not active.
    // However, the Check method is always called on all elements, even if they are not active.
    // Therefore, the Check method should not throw an exception in this case, even though the
    // constitutive laws are not initialized.
    element.Set(ACTIVE, false);

    const auto dummy_process_info = ProcessInfo{};
    KRATOS_EXPECT_EQ(element.Check(dummy_process_info), 0);
}

KRATOS_TEST_CASE_IN_SUITE(UPwLineInterfaceElement_InterpolatesNodalStresses, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    Model model;
    auto& r_model_part = model.CreateModelPart("Main");
    r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    r_model_part.AddNodalSolutionStepVariable(WATER_PRESSURE);
    PointerVector<Node> nodes;
    nodes.push_back(r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(2, 0.0, -1.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(3, 1.0, 0.0, 0.0));
    // properties for neighbour elements
    const auto p_properties = std::make_shared<Properties>();
    p_properties->SetValue(CONSTITUTIVE_LAW, std::make_shared<GeoIncrementalLinearElasticLaw>(
                                                 std::make_unique<PlaneStrain>()));
    p_properties->SetValue(YOUNG_MODULUS, 1.000000e+07);
    p_properties->SetValue(POISSON_RATIO, 0.000000e+00);
    // create a triangle neighbour element
    auto p_neighbour_element = make_intrusive<MockElementWithTotalStressVectors>(
        2, std::make_shared<Triangle2D3<Node>>(nodes), p_properties);
    p_neighbour_element->SetIntegrationMethod(GeometryData::IntegrationMethod::GI_GAUSS_2);
    ProcessInfo dummy_process_info;
    p_neighbour_element->Initialize(dummy_process_info);
    std::vector<ConstitutiveLaw::StressVectorType> total_stress_vectors;
    ConstitutiveLaw::StressVectorType              stress_vector(4);
    stress_vector <<= 3.0, 13.0 / 6.0, 4.0, 1.0;
    total_stress_vectors.emplace_back(stress_vector);
    total_stress_vectors.emplace_back(stress_vector);
    stress_vector <<= 3.0, 8.0 / 3.0, 4.0, 1.0;
    total_stress_vectors.emplace_back(stress_vector);

    p_neighbour_element->SetValuesOnIntegrationPoints(TOTAL_STRESS_VECTOR, total_stress_vectors, dummy_process_info);

    nodes.clear();
    nodes.push_back(r_model_part.pGetNode(1));
    nodes.push_back(r_model_part.pGetNode(3));
    nodes.push_back(r_model_part.CreateNewNode(4, 0.0, 1.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(5, 1.0, 1.0, 0.0));
    const auto     p_geometry       = std::make_shared<LineInterfaceGeometry2D2Plus2Noded>(nodes);
    constexpr auto normal_stiffness = 20.0;
    constexpr auto shear_stiffness  = 10.0;
    const auto     p_interface_properties =
        CreateElasticMaterialProperties<InterfacePlaneStrain>(normal_stiffness, shear_stiffness);
    auto interface_element = CreateInterfaceElementWithUPwDofs<Interface2D>(
        p_interface_properties, p_geometry, IsDiffOrderElement::No, {CalculationContribution::Stiffness});
    interface_element.SetValue(NEIGHBOUR_ELEMENTS, MakeElementGlobalPtrContainerWith(p_neighbour_element));

    // Act
    interface_element.Initialize(dummy_process_info);

    std::vector<Vector> actual_traction_vectors;
    interface_element.CalculateOnIntegrationPoints(GEO_EFFECTIVE_TRACTION_VECTOR,
                                                   actual_traction_vectors, dummy_process_info);

    // Assert
    // answers for 1st side
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(
        actual_traction_vectors[0], UblasUtilities::CreateVector({2.0, 1.0}), Defaults::relative_tolerance)
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(
        actual_traction_vectors[1], UblasUtilities::CreateVector({3.0, 1.0}), Defaults::relative_tolerance)

    // Arrange a neighbour on the other side of the interface element
    nodes.clear();
    nodes.push_back(r_model_part.pGetNode(4));
    nodes.push_back(r_model_part.pGetNode(5));
    nodes.push_back(r_model_part.CreateNewNode(6, 0.0, 2.0, 0.0));
    auto p_other_neighbour_element = make_intrusive<MockElementWithTotalStressVectors>(
        3, std::make_shared<Triangle2D3<Node>>(nodes), p_properties);
    p_other_neighbour_element->SetIntegrationMethod(GeometryData::IntegrationMethod::GI_GAUSS_2);
    p_other_neighbour_element->Initialize(dummy_process_info);
    total_stress_vectors.clear();
    stress_vector <<= 7.0, 4.0, 11.0, 5.0 / 6.0;
    total_stress_vectors.emplace_back(stress_vector);
    stress_vector <<= 7.0, 4.0, 11.0, 10.0 / 3.0;
    total_stress_vectors.emplace_back(stress_vector);
    stress_vector <<= 7.0, 4.0, 11.0, 5.0 / 6.0;
    total_stress_vectors.emplace_back(stress_vector);
    p_other_neighbour_element->SetValuesOnIntegrationPoints(
        TOTAL_STRESS_VECTOR, total_stress_vectors, dummy_process_info);
    interface_element.SetValue(NEIGHBOUR_ELEMENTS, MakeElementGlobalPtrContainerWith(p_other_neighbour_element));

    // Act
    interface_element.Initialize(dummy_process_info);
    interface_element.CalculateOnIntegrationPoints(GEO_EFFECTIVE_TRACTION_VECTOR,
                                                   actual_traction_vectors, dummy_process_info);

    // Assert
    // answers for 2nd side
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(
        actual_traction_vectors[0], UblasUtilities::CreateVector({4.0, 0.0}), Defaults::relative_tolerance)
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(
        actual_traction_vectors[1], UblasUtilities::CreateVector({4.0, 5.0}), Defaults::relative_tolerance)
}

KRATOS_TEST_CASE_IN_SUITE(UPwPlaneInterfaceElement_InterpolatesNodalStresses, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    Model model;
    auto& r_model_part = model.CreateModelPart("Main");
    r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    r_model_part.AddNodalSolutionStepVariable(WATER_PRESSURE);
    PointerVector<Node> nodes;
    nodes.push_back(r_model_part.CreateNewNode(1, 0.5, -1.5, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(2, 0.5, -0.5, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(3, -0.5, -0.5, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(4, -0.5, -1.5, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(5, 0.5, -1.5, 1.0));
    nodes.push_back(r_model_part.CreateNewNode(6, 0.5, -0.5, 1.0));
    nodes.push_back(r_model_part.CreateNewNode(7, -0.5, -0.5, 1.0));
    nodes.push_back(r_model_part.CreateNewNode(8, -0.5, -1.5, 1.0));
    // properties for neighbour elements
    const auto p_properties = std::make_shared<Properties>();
    p_properties->SetValue(CONSTITUTIVE_LAW, std::make_shared<GeoIncrementalLinearElasticLaw>(
                                                 std::make_unique<ThreeDimensional>()));
    p_properties->SetValue(YOUNG_MODULUS, 1.000000e+07);
    p_properties->SetValue(POISSON_RATIO, 0.000000e+00);
    // create a hexagonal neighbour element
    auto p_neighbour_element = make_intrusive<MockElementWithTotalStressVectors>(
        2, std::make_shared<Hexahedra3D8<Node>>(nodes), p_properties);
    p_neighbour_element->SetIntegrationMethod(GeometryData::IntegrationMethod::GI_GAUSS_2);

    ProcessInfo dummy_process_info;
    p_neighbour_element->Initialize(dummy_process_info);
    std::vector<ConstitutiveLaw::StressVectorType> total_stress_vectors;
    ConstitutiveLaw::StressVectorType              stress_vector(6);
    stress_vector <<= 3.0, (std::sqrt(3.0) - 1.0) / (2.0 * std::sqrt(3.0)), 4.0, 1.0, 2.0, 4.0;
    total_stress_vectors.emplace_back(stress_vector);
    total_stress_vectors.emplace_back(stress_vector);
    total_stress_vectors.emplace_back(stress_vector);
    total_stress_vectors.emplace_back(stress_vector);
    stress_vector <<= 3.0, (std::sqrt(3.0) + 1.0) / (2.0 * std::sqrt(3.0)), 4.0, 1.0, 2.0, 4.0;
    total_stress_vectors.emplace_back(stress_vector);
    total_stress_vectors.emplace_back(stress_vector);
    total_stress_vectors.emplace_back(stress_vector);
    total_stress_vectors.emplace_back(stress_vector);
    p_neighbour_element->SetValuesOnIntegrationPoints(TOTAL_STRESS_VECTOR, total_stress_vectors, dummy_process_info);

    nodes.clear();
    nodes.push_back(r_model_part.pGetNode(2));
    nodes.push_back(r_model_part.pGetNode(3));
    nodes.push_back(r_model_part.pGetNode(7));
    nodes.push_back(r_model_part.pGetNode(6));
    nodes.push_back(r_model_part.CreateNewNode(11, 0.5, 0.5, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(14, -0.5, 0.5, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(18, -0.5, 0.5, 1.0));
    nodes.push_back(r_model_part.CreateNewNode(15, 0.5, 0.5, 1.0));
    const auto p_geometry = std::make_shared<QuadrilateralInterfaceGeometry3D4Plus4Noded>(nodes);
    constexpr auto normal_stiffness = 20.0;
    constexpr auto shear_stiffness  = 10.0;
    const auto     p_interface_properties =
        CreateElasticMaterialProperties<InterfaceThreeDimensionalSurface>(normal_stiffness, shear_stiffness);
    auto interface_element = CreateInterfaceElementWithUPwDofs<Interface3D>(
        p_interface_properties, p_geometry, IsDiffOrderElement::No, {CalculationContribution::Stiffness});
    interface_element.SetValue(NEIGHBOUR_ELEMENTS, MakeElementGlobalPtrContainerWith(p_neighbour_element));

    // Act
    interface_element.Initialize(dummy_process_info);

    std::vector<Vector> actual_traction_vectors;
    interface_element.CalculateOnIntegrationPoints(GEO_EFFECTIVE_TRACTION_VECTOR,
                                                   actual_traction_vectors, dummy_process_info);

    // Assert
    // answers for 1st side
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(
        actual_traction_vectors[0], UblasUtilities::CreateVector({0.0, -1.0, 2.0}), Defaults::relative_tolerance)
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(
        actual_traction_vectors[1], UblasUtilities::CreateVector({0.0, -1.0, 2.0}), Defaults::relative_tolerance)
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(
        actual_traction_vectors[2], UblasUtilities::CreateVector({1.0, -1.0, 2.0}), Defaults::relative_tolerance)
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(
        actual_traction_vectors[3], UblasUtilities::CreateVector({1.0, -1.0, 2.0}), Defaults::relative_tolerance)

    // Arrange a neighbour on the other side of the interface element
    nodes.clear();
    nodes.push_back(r_model_part.pGetNode(11));
    nodes.push_back(r_model_part.CreateNewNode(12, 0.5, 1.5, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(13, -0.5, 1.5, 0.0));
    nodes.push_back(r_model_part.pGetNode(14));
    nodes.push_back(r_model_part.pGetNode(15));
    nodes.push_back(r_model_part.CreateNewNode(16, 0.5, 1.5, 1.0));
    nodes.push_back(r_model_part.CreateNewNode(17, -0.5, 1.5, 1.0));
    nodes.push_back(r_model_part.pGetNode(18));
    auto p_other_neighbour_element = make_intrusive<MockElementWithTotalStressVectors>(
        3, std::make_shared<Hexahedra3D8<Node>>(nodes), p_properties);
    p_neighbour_element->SetIntegrationMethod(GeometryData::IntegrationMethod::GI_GAUSS_2);
    p_other_neighbour_element->Initialize(dummy_process_info);
    total_stress_vectors.clear();
    stress_vector <<= 3.0, 4.0, 1.0, 2.0, (std::sqrt(3.0) - 1.0) / (2.0 * std::sqrt(3.0)), 4.0;
    total_stress_vectors.emplace_back(stress_vector);
    total_stress_vectors.emplace_back(stress_vector);
    total_stress_vectors.emplace_back(stress_vector);
    total_stress_vectors.emplace_back(stress_vector);
    stress_vector <<= 3.0, 4.0, 1.0, 2.0, (std::sqrt(3.0) + 1.0) / (2.0 * std::sqrt(3.0)), 4.0;
    total_stress_vectors.emplace_back(stress_vector);
    total_stress_vectors.emplace_back(stress_vector);
    total_stress_vectors.emplace_back(stress_vector);
    total_stress_vectors.emplace_back(stress_vector);
    p_other_neighbour_element->SetValuesOnIntegrationPoints(
        TOTAL_STRESS_VECTOR, total_stress_vectors, dummy_process_info);

    interface_element.SetValue(NEIGHBOUR_ELEMENTS, MakeElementGlobalPtrContainerWith(p_other_neighbour_element));

    // Act
    interface_element.Initialize(dummy_process_info);
    interface_element.CalculateOnIntegrationPoints(GEO_EFFECTIVE_TRACTION_VECTOR,
                                                   actual_traction_vectors, dummy_process_info);

    // Assert
    // answers for 2nd side
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(
        actual_traction_vectors[0], UblasUtilities::CreateVector({4.0, -2.0, 0.0}), Defaults::relative_tolerance)
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(
        actual_traction_vectors[1], UblasUtilities::CreateVector({4.0, -2.0, 0.0}), Defaults::relative_tolerance)
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(
        actual_traction_vectors[2], UblasUtilities::CreateVector({4.0, -2.0, 1.0}), Defaults::relative_tolerance)
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(
        actual_traction_vectors[3], UblasUtilities::CreateVector({4.0, -2.0, 1.0}), Defaults::relative_tolerance)
}

KRATOS_TEST_CASE_IN_SUITE(ThreePlusThreeUPwDiffOrderLineInterfaceElement_KeepsUDofsFirstThenPwDofs,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto nodes = ModelSetupUtilities::CreateNodes({{1, {0.0, 0.0, 0.0}},
                                                         {2, {1.0, 0.0, 0.0}},
                                                         {3, {0.5, 0.0, 0.0}},
                                                         {4, {0.0, 0.0, 0.0}},
                                                         {5, {1.0, 0.0, 0.0}},
                                                         {6, {0.5, 0.0, 0.0}}});

    auto p_line_interface_displacement_geometry = std::make_shared<InterfaceGeometry<Line2D3<Node>>>(1, nodes);

    constexpr auto normal_stiffness = 20.0;
    constexpr auto shear_stiffness  = 10.0;
    const auto     p_interface_properties =
        CreateElasticMaterialProperties<InterfacePlaneStrain>(normal_stiffness, shear_stiffness);

    auto interface_element = CreateInterfaceElementWithUPwDofs<Interface2D>(
        p_interface_properties, p_line_interface_displacement_geometry, IsDiffOrderElement::Yes,
        {CalculationContribution::Stiffness});

    // Act
    Element::DofsVectorType element_dofs;
    ProcessInfo             dummy_process_info;
    interface_element.GetDofList(element_dofs, dummy_process_info);

    // Assert
    constexpr auto expected_number_of_u_dofs  = std::size_t{6 * 2};
    constexpr auto expected_number_of_pw_dofs = std::size_t{4};
    ASSERT_EQ(element_dofs.size(), expected_number_of_u_dofs + expected_number_of_pw_dofs);
    for (auto i = std::size_t{0}; i < expected_number_of_u_dofs; i += 2) {
        KRATOS_EXPECT_EQ(element_dofs[i]->GetVariable(), DISPLACEMENT_X);
        KRATOS_EXPECT_EQ(element_dofs[i + 1]->GetVariable(), DISPLACEMENT_Y);
    }
    for (auto i = std::size_t{0}; i < expected_number_of_pw_dofs; ++i) {
        KRATOS_EXPECT_EQ(element_dofs[expected_number_of_u_dofs + i]->GetVariable(), WATER_PRESSURE);
    }
}

KRATOS_TEST_CASE_IN_SUITE(ThreePlusThreeDiffOrderLineInterface_ReturnsExpectedLeftAndRightHandSide,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto nodes = ModelSetupUtilities::CreateNodes({{1, {0.0, 0.0, 0.0}},
                                                         {2, {1.0, 0.0, 0.0}},
                                                         {3, {0.5, 0.0, 0.0}},
                                                         {4, {0.0, 0.0, 0.0}},
                                                         {5, {1.0, 0.0, 0.0}},
                                                         {6, {0.5, 0.0, 0.0}}});

    auto p_line_interface_displacement_geometry = std::make_shared<InterfaceGeometry<Line2D3<Node>>>(1, nodes);

    constexpr auto normal_stiffness = 20.0;
    constexpr auto shear_stiffness  = 10.0;
    const auto     p_interface_properties =
        CreateElasticMaterialProperties<InterfacePlaneStrain>(normal_stiffness, shear_stiffness);

    auto interface_element = CreateInterfaceElementWithUPwDofs<Interface2D>(
        p_interface_properties, p_line_interface_displacement_geometry, IsDiffOrderElement::Yes,
        {CalculationContribution::Stiffness});

    // Act
    interface_element.Initialize(ProcessInfo{});
    const auto prescribed_displacements =
        PrescribedDisplacements{{3, array_1d<double, 3>{0.2, 0.5, 0.0}},
                                {4, array_1d<double, 3>{0.2, 0.5, 0.0}},
                                {5, array_1d<double, 3>{0.2, 0.5, 0.0}}};
    for (const auto& [idx, disp] : prescribed_displacements) {
        interface_element.GetGeometry()[idx].FastGetSolutionStepValue(DISPLACEMENT) = disp;
    }

    Matrix actual_left_hand_side;
    Vector actual_right_hand_side;
    interface_element.CalculateLocalSystem(actual_left_hand_side, actual_right_hand_side, ProcessInfo{});

    // Assert
    constexpr auto number_of_u_dofs         = std::size_t{6 * 2};
    constexpr auto number_of_pw_dofs        = std::size_t{4};
    constexpr auto cs                       = shear_stiffness / 6.0;
    constexpr auto cn                       = normal_stiffness / 6.0;
    const auto     expected_uu_block_matrix = UblasUtilities::CreateMatrix(
        {{cs, 0.0, 0.0, 0.0, 0.0, 0.0, -cs, 0.0, 0.0, 0.0, 0.0, 0.0},
             {0.0, cn, 0.0, 0.0, 0.0, 0.0, 0.0, -cn, 0.0, 0.0, 0.0, 0.0},
             {0.0, 0.0, cs, 0.0, 0.0, 0.0, 0.0, 0.0, -cs, 0.0, 0.0, 0.0},
             {0.0, 0.0, 0.0, cn, 0.0, 0.0, 0.0, 0.0, 0.0, -cn, 0.0, 0.0},
             {0.0, 0.0, 0.0, 0.0, 4.0 * cs, 0.0, 0.0, 0.0, 0.0, 0.0, -4.0 * cs, 0.0},
             {0.0, 0.0, 0.0, 0.0, 0.0, 4.0 * cn, 0.0, 0.0, 0.0, 0.0, 0.0, -4.0 * cn},
             {-cs, 0.0, 0.0, 0.0, 0.0, 0.0, cs, 0.0, 0.0, 0.0, 0.0, 0.0},
             {0.0, -cn, 0.0, 0.0, 0.0, 0.0, 0.0, cn, 0.0, 0.0, 0.0, 0.0},
             {0.0, 0.0, -cs, 0.0, 0.0, 0.0, 0.0, 0.0, cs, 0.0, 0.0, 0.0},
             {0.0, 0.0, 0.0, -cn, 0.0, 0.0, 0.0, 0.0, 0.0, cn, 0.0, 0.0},
             {0.0, 0.0, 0.0, 0.0, -4.0 * cs, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0 * cs, 0.0},
             {0.0, 0.0, 0.0, 0.0, 0.0, -4.0 * cn, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0 * cn}});
    const auto expected_up_block_matrix = Matrix{number_of_u_dofs, number_of_pw_dofs, 0.0};
    const auto expected_pu_block_matrix = Matrix{number_of_pw_dofs, number_of_u_dofs, 0.0};
    const auto expected_pp_block_matrix = Matrix{number_of_pw_dofs, number_of_pw_dofs, 0.0};
    AssertLHSMatrixBlocksAreNear(actual_left_hand_side, expected_uu_block_matrix,
                                 expected_up_block_matrix, expected_pu_block_matrix,
                                 expected_pp_block_matrix, number_of_u_dofs, number_of_pw_dofs);

    const auto expected_u_block_vector = UblasUtilities::CreateVector(
        {1.0 / 3.0, 5.0 / 3.0, 1.0 / 3.0, 5.0 / 3.0, 4.0 / 3.0, 20.0 / 3.0, -1.0 / 3.0, -5.0 / 3.0,
         -1.0 / 3.0, -5.0 / 3.0, -4.0 / 3.0, -20.0 / 3.0});
    const auto expected_p_block_vector = Vector{number_of_pw_dofs, 0.0};
    AssertRHSVectorBlocksAreNear(actual_right_hand_side, expected_u_block_vector,
                                 expected_p_block_vector, number_of_u_dofs, number_of_pw_dofs);
}

KRATOS_TEST_CASE_IN_SUITE(UPwDiffOrderTriangleInterfaceElement_6Plus6NodedElement_ReturnsExpectedLeftAndRightHandSide,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto nodes = ModelSetupUtilities::CreateNodes({{1, {0.0, 0.0, 0.0}},
                                                         {2, {1.0, 0.0, 0.0}},
                                                         {3, {1.0, 1.0, 0.0}},
                                                         {4, {0.5, 0.0, 0.0}},
                                                         {5, {1.0, 0.5, 0.0}},
                                                         {6, {0.5, 0.5, 0.0}},
                                                         {7, {0.0, 0.0, 0.0}},
                                                         {8, {1.0, 0.0, 0.0}},
                                                         {9, {1.0, 1.0, 0.0}},
                                                         {10, {0.5, 0.0, 0.0}},
                                                         {11, {1.0, 0.5, 0.0}},
                                                         {12, {0.5, 0.5, 0.0}}});

    auto p_triangle_interface_displacement_geometry =
        std::make_shared<InterfaceGeometry<Triangle3D6<Node>>>(1, nodes);

    constexpr auto normal_stiffness = 20.0;
    constexpr auto shear_stiffness  = 10.0;
    const auto     p_properties = CreateElasticMaterialProperties<InterfaceThreeDimensionalSurface>(
        normal_stiffness, shear_stiffness);

    auto interface_element = CreateInterfaceElementWithUPwDofs<Interface3D>(
        p_properties, p_triangle_interface_displacement_geometry, IsDiffOrderElement::Yes,
        {CalculationContribution::Stiffness});
    interface_element.Initialize(ProcessInfo{});

    const auto prescribed_displacements = PrescribedDisplacements{
        {6, array_1d<double, 3>{0.2, 0.5, 0.0}},  {7, array_1d<double, 3>{0.2, 0.5, 0.0}},
        {8, array_1d<double, 3>{0.2, 0.5, 0.0}},  {9, array_1d<double, 3>{0.2, 0.5, 0.0}},
        {10, array_1d<double, 3>{0.2, 0.5, 0.0}}, {11, array_1d<double, 3>{0.2, 0.5, 0.0}}};
    for (const auto& [idx, disp] : prescribed_displacements) {
        interface_element.GetGeometry()[idx].FastGetSolutionStepValue(DISPLACEMENT) = disp;
    }

    // Act
    Vector actual_right_hand_side;
    Matrix actual_left_hand_side;
    interface_element.CalculateLocalSystem(actual_left_hand_side, actual_right_hand_side, ProcessInfo{});

    // Assert
    constexpr auto number_of_u_dofs  = std::size_t{12 * 3};
    constexpr auto number_of_pw_dofs = std::size_t{6};
    ASSERT_EQ(actual_left_hand_side.size1(), number_of_u_dofs + number_of_pw_dofs);
    ASSERT_EQ(actual_left_hand_side.size2(), number_of_u_dofs + number_of_pw_dofs);
    ASSERT_EQ(actual_right_hand_side.size(), number_of_u_dofs + number_of_pw_dofs);

    // clang-format off
    const auto expected_uu_block_matrix = UblasUtilities::CreateMatrix(
        {{0.16431924882629123,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.16431924882629123, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
         {0, 0.16431924882629123,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.16431924882629123, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
         {0, 0, 0.32863849765258246,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.32863849765258246, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
         {0, 0, 0, 0.16431924882629123,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.16431924882629123, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
         {0, 0, 0, 0, 0.16431924882629123,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.16431924882629123, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
         {0, 0, 0, 0, 0, 0.32863849765258246,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.32863849765258246, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
         {0, 0, 0, 0, 0, 0, 0.16431924882629123,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.16431924882629123, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
         {0, 0, 0, 0, 0, 0, 0, 0.16431924882629123,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.16431924882629123, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
         {0, 0, 0, 0, 0, 0, 0, 0, 0.32863849765258246,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.32863849765258246, 0, 0, 0, 0, 0, 0, 0, 0, 0},
         {0, 0, 0, 0, 0, 0, 0, 0, 0, 1.5023474178403748,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1.5023474178403748, 0, 0, 0, 0, 0, 0, 0, 0},
         {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.5023474178403748,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1.5023474178403748, 0, 0, 0, 0, 0, 0, 0},
         {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3.0046948356807497,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -3.0046948356807497, 0, 0, 0, 0, 0, 0},
         {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.502347417840376,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1.502347417840376, 0, 0, 0, 0, 0},
         {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.502347417840376,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1.502347417840376, 0, 0, 0, 0},
         {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3.0046948356807519,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -3.0046948356807519, 0, 0, 0},
         {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.5023474178403753,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1.5023474178403753, 0, 0},
         {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.5023474178403753,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1.5023474178403753, 0},
         {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3.0046948356807506, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -3.0046948356807506},
         {-0.16431924882629123, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.16431924882629123,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
         {0, -0.16431924882629123, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.16431924882629123,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
         {0, 0, -0.32863849765258246, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.32863849765258246,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
         {0, 0, 0, -0.16431924882629123, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.16431924882629123,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
         {0, 0, 0, 0, -0.16431924882629123, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.16431924882629123,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
         {0, 0, 0, 0, 0, -0.32863849765258246, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.32863849765258246,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
         {0, 0, 0, 0, 0, 0, -0.16431924882629123, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.16431924882629123,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
         {0, 0, 0, 0, 0, 0, 0, -0.16431924882629123, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.16431924882629123,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
         {0, 0, 0, 0, 0, 0, 0, 0, -0.32863849765258246, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.32863849765258246,  0, 0, 0, 0, 0, 0, 0, 0, 0},
         {0, 0, 0, 0, 0, 0, 0, 0, 0, -1.5023474178403748, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.5023474178403748,  0, 0, 0, 0, 0, 0, 0, 0},
         {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1.5023474178403748, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.5023474178403748,  0, 0, 0, 0, 0, 0, 0},
         {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -3.0046948356807497, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3.0046948356807497,  0, 0, 0, 0, 0, 0},
         {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1.502347417840376, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.502347417840376,  0, 0, 0, 0, 0},
         {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1.502347417840376, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.502347417840376,  0, 0, 0, 0},
         {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -3.0046948356807519, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3.0046948356807519,  0, 0, 0},
         {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1.5023474178403753, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.5023474178403753,  0, 0},
         {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1.5023474178403753, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.5023474178403753,  0},
         {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -3.0046948356807506, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3.0046948356807506}});
    // clang-format on

    KRATOS_EXPECT_MATRIX_RELATIVE_NEAR(
        subrange(actual_left_hand_side, 0, 0 + number_of_u_dofs, 0, 0 + number_of_u_dofs),
        expected_uu_block_matrix, Defaults::relative_tolerance)
    const auto expected_up_block_matrix = Matrix{number_of_u_dofs, number_of_pw_dofs, 0.0};
    AssertUPBlockMatrixIsNear(actual_left_hand_side, expected_up_block_matrix, number_of_u_dofs, number_of_pw_dofs);
    const auto expected_pu_block_matrix = Matrix{number_of_pw_dofs, number_of_u_dofs, 0.0};
    AssertPUBlockMatrixIsNear(actual_left_hand_side, expected_pu_block_matrix, number_of_u_dofs, number_of_pw_dofs);
    const auto expected_pp_block_matrix = Matrix{number_of_pw_dofs, number_of_pw_dofs, 0.0};
    AssertPPBlockMatrixIsNear(actual_left_hand_side, expected_pp_block_matrix, number_of_u_dofs, number_of_pw_dofs);

    const auto expected_u_block_vector = UblasUtilities::CreateVector(
        {0.0328638,  0.0821596,  -0, 0.0328638,  0.0821596,  -0, 0.0328638,  0.0821596,  -0,
         0.300469,   0.751174,   -0, 0.300469,   0.751174,   -0, 0.300469,   0.751174,   -0,
         -0.0328638, -0.0821596, -0, -0.0328638, -0.0821596, -0, -0.0328638, -0.0821596, -0,
         -0.300469,  -0.751174,  -0, -0.300469,  -0.751174,  -0, -0.300469,  -0.751174,  -0});
    constexpr auto tolerance = 1e-5;
    KRATOS_EXPECT_VECTOR_NEAR(subrange(actual_right_hand_side, 0, 0 + number_of_u_dofs),
                              expected_u_block_vector, tolerance)
    const auto expected_p_block_vector = Vector{number_of_pw_dofs, 0.0};
    AssertPBlockVectorIsNear(actual_right_hand_side, expected_p_block_vector, number_of_u_dofs, number_of_pw_dofs);
}

KRATOS_TEST_CASE_IN_SUITE(UPwDiffOrderQuadrilateraleInterfaceElement_8Plus8NodedElement_ReturnsExpectedLeftAndRightHandSide,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    constexpr auto normal_stiffness = 20.0;
    constexpr auto shear_stiffness  = 10.0;
    const auto     p_properties = CreateElasticMaterialProperties<InterfaceThreeDimensionalSurface>(
        normal_stiffness, shear_stiffness);
    Model model;
    auto  interface_element = CreateHorizontal8Plus8NodedQuadraliteralInterfaceElementWithUPwDoF(
        model, p_properties, IsDiffOrderElement::Yes, {CalculationContribution::Stiffness});

    interface_element.Initialize(ProcessInfo{});

    const auto prescribed_displacements = PrescribedDisplacements{
        {8, array_1d<double, 3>{0.2, 0.5, 0.0}},  {9, array_1d<double, 3>{0.2, 0.5, 0.0}},
        {10, array_1d<double, 3>{0.2, 0.5, 0.0}}, {11, array_1d<double, 3>{0.2, 0.5, 0.0}},
        {12, array_1d<double, 3>{0.2, 0.5, 0.0}}, {13, array_1d<double, 3>{0.2, 0.5, 0.0}},
        {14, array_1d<double, 3>{0.2, 0.5, 0.0}}, {15, array_1d<double, 3>{0.2, 0.5, 0.0}}};
    for (const auto& [idx, disp] : prescribed_displacements) {
        interface_element.GetGeometry()[idx].FastGetSolutionStepValue(DISPLACEMENT) = disp;
    }

    // Act
    Vector actual_right_hand_side;
    Matrix actual_left_hand_side;
    interface_element.CalculateLocalSystem(actual_left_hand_side, actual_right_hand_side, ProcessInfo{});

    // Assert
    constexpr auto number_of_u_dofs  = std::size_t{16 * 3};
    constexpr auto number_of_pw_dofs = std::size_t{8};
    ASSERT_EQ(actual_left_hand_side.size1(), number_of_u_dofs + number_of_pw_dofs);
    ASSERT_EQ(actual_left_hand_side.size2(), number_of_u_dofs + number_of_pw_dofs);
    ASSERT_EQ(actual_right_hand_side.size(), number_of_u_dofs + number_of_pw_dofs);

    // clang-format off
    const auto expected_uu_block_matrix = UblasUtilities::CreateMatrix(
           {{0.39473684210526316,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.39473684210526316,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0.39473684210526316,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.39473684210526316,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,0.78947368421052633,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.78947368421052633,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,0,0.39473684210526327,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.39473684210526327,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,0.39473684210526327,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.39473684210526327,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,0,0.78947368421052655,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.78947368421052655,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,0,0,0.39473684210526322,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.39473684210526322,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,0,0,0,0.39473684210526322,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.39473684210526322,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,0,0,0,0,0.78947368421052644,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.78947368421052644,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,0,0,0,0,0,0.39473684210526316,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.39473684210526316,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,0,0,0,0,0,0,0.39473684210526316,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.39473684210526316,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,0,0,0,0,0,0,0,0.78947368421052633,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.78947368421052633,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,0,0,0,0,0,0,0,0,2.1052631578947372,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-2.1052631578947372,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,0,0,0,0,0,0,0,0,0,2.1052631578947372,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-2.1052631578947372,0,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,4.2105263157894743,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-4.2105263157894743,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2.1052631578947367,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-2.1052631578947367,0,0,0,0,0,0,0,0},
            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2.1052631578947367,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-2.1052631578947367,0,0,0,0,0,0,0},
            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4.2105263157894735,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-4.2105263157894735,0,0,0,0,0,0},
            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2.1052631578947367,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-2.1052631578947367,0,0,0,0,0},
            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2.1052631578947367,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-2.1052631578947367,0,0,0,0},
            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4.2105263157894735,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-4.2105263157894735,0,0,0},
            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2.1052631578947372,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-2.1052631578947372,0,0},
            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2.1052631578947372,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-2.1052631578947372,0},
            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4.2105263157894743,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-4.2105263157894743},
            {-0.39473684210526316,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.39473684210526316,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,-0.39473684210526316,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.39473684210526316,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,-0.78947368421052633,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.78947368421052633,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,0,-0.39473684210526327,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.39473684210526327,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,-0.39473684210526327,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.39473684210526327,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,0,-0.78947368421052655,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.78947368421052655,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,0,0,-0.39473684210526322,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.39473684210526322,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,0,0,0,-0.39473684210526322,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.39473684210526322,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,0,0,0,0,-0.78947368421052644,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.78947368421052644,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,0,0,0,0,0,-0.39473684210526316,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.39473684210526316,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,0,0,0,0,0,0,-0.39473684210526316,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.39473684210526316,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,0,0,0,0,0,0,0,-0.78947368421052633,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.78947368421052633,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,0,0,0,0,0,0,0,0,-2.1052631578947372,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2.1052631578947372,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,0,0,0,0,0,0,0,0,0,-2.1052631578947372,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2.1052631578947372,0,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,-4.2105263157894743,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4.2105263157894743,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-2.1052631578947367,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2.1052631578947367,0,0,0,0,0,0,0,0},
            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-2.1052631578947367,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2.1052631578947367,0,0,0,0,0,0,0},
            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-4.2105263157894735,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4.2105263157894735,0,0,0,0,0,0},
            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-2.1052631578947367,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2.1052631578947367,0,0,0,0,0},
            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-2.1052631578947367,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2.1052631578947367,0,0,0,0},
            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-4.2105263157894735,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4.2105263157894735,0,0,0},
            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-2.1052631578947372,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2.1052631578947372,0,0},
            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-2.1052631578947372,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2.1052631578947372,0},
            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-4.2105263157894743,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4.2105263157894743}});
    // clang-format on

    KRATOS_EXPECT_MATRIX_RELATIVE_NEAR(
        subrange(actual_left_hand_side, 0, 0 + number_of_u_dofs, 0, 0 + number_of_u_dofs),
        expected_uu_block_matrix, Defaults::relative_tolerance)
    const auto expected_up_block_matrix = Matrix{number_of_u_dofs, number_of_pw_dofs, 0.0};
    AssertUPBlockMatrixIsNear(actual_left_hand_side, expected_up_block_matrix, number_of_u_dofs, number_of_pw_dofs);
    const auto expected_pu_block_matrix = Matrix{number_of_pw_dofs, number_of_u_dofs, 0.0};
    AssertPUBlockMatrixIsNear(actual_left_hand_side, expected_pu_block_matrix, number_of_u_dofs, number_of_pw_dofs);
    const auto expected_pp_block_matrix = Matrix{number_of_pw_dofs, number_of_pw_dofs, 0.0};
    AssertPPBlockMatrixIsNear(actual_left_hand_side, expected_pp_block_matrix, number_of_u_dofs, number_of_pw_dofs);

    const auto expected_u_block_vector = UblasUtilities::CreateVector(
        {0.0789474,  0.197368,   -0,        0.0789474,  0.197368,  -0,        0.0789474,  0.197368,
         -0,         0.0789474,  0.197368,  -0,         0.421053,  1.05263,   -0,         0.421053,
         1.05263,    -0,         0.421053,  1.05263,    -0,        0.421053,  1.05263,    -0,
         -0.0789474, -0.197368,  -0,        -0.0789474, -0.197368, -0,        -0.0789474, -0.197368,
         -0,         -0.0789474, -0.197368, -0,         -0.421053, -1.05263,  -0,         -0.421053,
         -1.05263,   -0,         -0.421053, -1.05263,   -0,        -0.421053, -1.05263,   -0});
    constexpr auto tolerance = 1e-5;
    KRATOS_EXPECT_VECTOR_NEAR(subrange(actual_right_hand_side, 0, 0 + number_of_u_dofs),
                              expected_u_block_vector, tolerance)
    const auto expected_p_block_vector = Vector{number_of_pw_dofs, 0.0};
    AssertPBlockVectorIsNear(actual_right_hand_side, expected_p_block_vector, number_of_u_dofs, number_of_pw_dofs);
}

KRATOS_TEST_CASE_IN_SUITE(UPwLineInterfaceElement_2Plus2Element_CouplingContribution, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto expected_up_block_matrix =
        UblasUtilities::CreateMatrix({{0, 0, 0, 0},
                                      {-0.080877192, 0, -0.080877192, 0},
                                      {0, 0, 0, 0},
                                      {0, -0.080877192, 0, -0.080877192},
                                      {0, 0, 0, 0},
                                      {0.080877192, 0, 0.080877192, 0},
                                      {0, 0, 0, 0},
                                      {0, 0.080877192, 0, 0.080877192}});
    const auto expected_up_block_vector = UblasUtilities::CreateVector(
        {0, -8.1685964, 0, -8.1685964, 0, 8.1685964, 0, 8.1685964, 0, 0, 0, 0});

    GeneraizedCouplingContributionTest<InterfacePlaneStrain>(
        CreateHorizontalUnitLength2Plus2NodedLineInterfaceElementWithUPwDofs,
        {CalculationContribution::UPCoupling}, IsDiffOrderElement::No, 8, 4,
        expected_up_block_matrix, expected_up_block_vector);
}

KRATOS_TEST_CASE_IN_SUITE(UPwLineInterfaceElement_3Plus3Element_CouplingContribution, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto expected_up_block_matrix =
        UblasUtilities::CreateMatrix({{0, 0, 0, 0, 0, 0},
                                      {-0.026959064, 0, 0, -0.026959064, 0, 0},
                                      {0, 0, 0, 0, 0, 0},
                                      {0, -0.026959064, 0, 0, -0.026959064, 0},
                                      {0, 0, 0, 0, 0, 0},
                                      {0, 0, -0.10783626, 0, 0, -0.10783626},
                                      {0, 0, 0, 0, 0, 0},
                                      {0.026959064, 0, 0, 0.026959064, 0, 0},
                                      {0, 0, 0, 0, 0, 0},
                                      {0, 0.026959064, 0, 0, 0.026959064, 0},
                                      {0, 0, 0, 0, 0, 0},
                                      {0, 0, 0.10783626, 0, 0, 0.10783626}});
    const auto expected_up_block_vector =
        UblasUtilities::CreateVector({0, -2.7228655, 0, -2.7228655, 0, -10.891462, 0, 2.7228655, 0,
                                      2.7228655, 0, 10.891462, 0, 0, 0, 0, 0, 0});

    GeneraizedCouplingContributionTest<InterfacePlaneStrain>(
        CreateHorizontalUnitLength3Plus3NodedLineInterfaceElementWithUPwDoF,
        {CalculationContribution::UPCoupling}, IsDiffOrderElement::No, 12, 6,
        expected_up_block_matrix, expected_up_block_vector);
}

KRATOS_TEST_CASE_IN_SUITE(UPwLineInterfaceElement_3Plus3Element_CouplingContribution_DiffOrder,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto expected_up_block_matrix =
        UblasUtilities::CreateMatrix({{0, 0, 0, 0},
                                      {-0.026959064, 0, -0.026959064, 0},
                                      {0, 0, 0, 0},
                                      {0, -0.026959064, 0, -0.026959064},
                                      {0, 0, 0, 0},
                                      {-0.053918128, -0.053918128, -0.053918128, -0.053918128},
                                      {0, 0, 0, 0},
                                      {0.026959064, 0, 0.026959064, 0},
                                      {0, 0, 0, 0},
                                      {0, 0.026959064, 0, 0.026959064},
                                      {0, 0, 0, 0},
                                      {0.053918128, 0.053918128, 0.053918128, 0.053918128}});
    const auto expected_up_block_vector = UblasUtilities::CreateVector(
        {0, -2.7228655, 0, -2.7228655, 0, -10.891462, 0, 2.7228655, 0, 2.7228655, 0, 10.891462, 0, 0, 0, 0});
    GeneraizedCouplingContributionTest<InterfacePlaneStrain>(
        CreateHorizontalUnitLength3Plus3NodedLineInterfaceElementWithUPwDoF,
        {CalculationContribution::UPCoupling}, IsDiffOrderElement::Yes, 12, 4,
        expected_up_block_matrix, expected_up_block_vector);
}

KRATOS_TEST_CASE_IN_SUITE(UPwTriangleInterfaceElement_3Plus3Element_CouplingContribution,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto expected_up_block_matrix =
        UblasUtilities::CreateMatrix({{0, 0, 0, 0, 0, 0},
                                      {0, 0, 0, 0, 0, 0},
                                      {-0.026959064, 0, 0, -0.026959064, 0, 0},
                                      {0, 0, 0, 0, 0, 0},
                                      {0, 0, 0, 0, 0, 0},
                                      {0, -0.026959064, 0, 0, -0.026959064, 0},
                                      {0, 0, 0, 0, 0, 0},
                                      {0, 0, 0, 0, 0, 0},
                                      {0, 0, -0.026959064, 0, 0, -0.026959064},
                                      {0, 0, 0, 0, 0, 0},
                                      {0, 0, 0, 0, 0, 0},
                                      {0.026959064, 0, 0, 0.026959064, 0, 0},
                                      {0, 0, 0, 0, 0, 0},
                                      {0, 0, 0, 0, 0, 0},
                                      {0, 0.026959064, 0, 0, 0.026959064, 0},
                                      {0, 0, 0, 0, 0, 0},
                                      {0, 0, 0, 0, 0, 0},
                                      {0, 0, 0.026959064, 0, 0, 0.026959064}});
    const auto expected_up_block_vector = UblasUtilities::CreateVector(
        {0, 0, -2.7228655, 0, 0, -2.7228655, 0, 0, -2.7228655, 0, 0, 2.7228655,
         0, 0, 2.7228655,  0, 0, 2.7228655,  0, 0, 0,          0, 0, 0});
    GeneraizedCouplingContributionTest<InterfaceThreeDimensionalSurface>(
        CreateHorizontal3Plus3NodedTriangleInterfaceElementWithUPwDofs, {CalculationContribution::UPCoupling},
        IsDiffOrderElement::No, 18, 6, expected_up_block_matrix, expected_up_block_vector);
}

KRATOS_TEST_CASE_IN_SUITE(UPwTriangleInterfaceElement_6Plus6Element_CouplingContribution,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto expected_up_block_matrix =
        UblasUtilities::CreateMatrix({{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                      {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                      {-0.0026579359, 0, 0, 0, 0, 0, -0.0026579359, 0, 0, 0, 0, 0},
                                      {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                      {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                      {0, -0.0026579359, 0, 0, 0, 0, 0, -0.0026579359, 0, 0, 0, 0},
                                      {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                      {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                      {0, 0, -0.0026579359, 0, 0, 0, 0, 0, -0.0026579359, 0, 0, 0},
                                      {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                      {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                      {0, 0, 0, -0.024301128, 0, 0, 0, 0, 0, -0.024301128, 0, 0},
                                      {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                      {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                      {0, 0, 0, 0, -0.024301128, 0, 0, 0, 0, 0, -0.024301128, 0},
                                      {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                      {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                      {0, 0, 0, 0, 0, -0.024301128, 0, 0, 0, 0, 0, -0.024301128},
                                      {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                      {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                      {0.0026579359, 0, 0, 0, 0, 0, 0.0026579359, 0, 0, 0, 0, 0},
                                      {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                      {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                      {0, 0.0026579359, 0, 0, 0, 0, 0, 0.0026579359, 0, 0, 0, 0},
                                      {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                      {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                      {0, 0, 0.0026579359, 0, 0, 0, 0, 0, 0.0026579359, 0, 0, 0},
                                      {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                      {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                      {0, 0, 0, 0.024301128, 0, 0, 0, 0, 0, 0.024301128, 0, 0},
                                      {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                      {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                      {0, 0, 0, 0, 0.024301128, 0, 0, 0, 0, 0, 0.024301128, 0},
                                      {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                      {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                      {0, 0, 0, 0, 0, 0.024301128, 0, 0, 0, 0, 0, 0.024301128}});
    const auto expected_up_block_vector = UblasUtilities::CreateVector(
        {0, 0, -0.26845153, 0, 0, -0.26845153, 0, 0, -0.26845153, 0, 0, -2.454414,
         0, 0, -2.454414,   0, 0, -2.454414,   0, 0, 0.26845153,  0, 0, 0.26845153,
         0, 0, 0.26845153,  0, 0, 2.454414,    0, 0, 2.454414,    0, 0, 2.454414,
         0, 0, 0,           0, 0, 0,           0, 0, 0,           0, 0, 0});
    GeneraizedCouplingContributionTest<InterfaceThreeDimensionalSurface>(
        CreateHorizontal6Plus6NodedTriangleInterfaceElementWithUPwDoF, {CalculationContribution::UPCoupling},
        IsDiffOrderElement::No, 36, 12, expected_up_block_matrix, expected_up_block_vector);
}

KRATOS_TEST_CASE_IN_SUITE(UPwTriangleInterfaceElement_6Plus6Element_CouplingContribution_DiffOrder,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto expected_up_block_matrix =
        UblasUtilities::CreateMatrix({{0, 0, 0, 0, 0, 0},
                                      {0, 0, 0, 0, 0, 0},
                                      {-0.0026579359, 0, 0, -0.0026579359, 0, 0},
                                      {0, 0, 0, 0, 0, 0},
                                      {0, 0, 0, 0, 0, 0},
                                      {0, -0.0026579359, 0, 0, -0.0026579359, 0},
                                      {0, 0, 0, 0, 0, 0},
                                      {0, 0, 0, 0, 0, 0},
                                      {0, 0, -0.0026579359, 0, 0, -0.0026579359},
                                      {0, 0, 0, 0, 0, 0},
                                      {0, 0, 0, 0, 0, 0},
                                      {-0.012150564, -0.012150564, 0, -0.012150564, -0.012150564, 0},
                                      {0, 0, 0, 0, 0, 0},
                                      {0, 0, 0, 0, 0, 0},
                                      {0, -0.012150564, -0.012150564, 0, -0.012150564, -0.012150564},
                                      {0, 0, 0, 0, 0, 0},
                                      {0, 0, 0, 0, 0, 0},
                                      {-0.012150564, 0, -0.012150564, -0.012150564, 0, -0.012150564},
                                      {0, 0, 0, 0, 0, 0},
                                      {0, 0, 0, 0, 0, 0},
                                      {0.0026579359, 0, 0, 0.0026579359, 0, 0},
                                      {0, 0, 0, 0, 0, 0},
                                      {0, 0, 0, 0, 0, 0},
                                      {0, 0.0026579359, 0, 0, 0.0026579359, 0},
                                      {0, 0, 0, 0, 0, 0},
                                      {0, 0, 0, 0, 0, 0},
                                      {0, 0, 0.0026579359, 0, 0, 0.0026579359},
                                      {0, 0, 0, 0, 0, 0},
                                      {0, 0, 0, 0, 0, 0},
                                      {0.012150564, 0.012150564, 0, 0.012150564, 0.012150564, 0},
                                      {0, 0, 0, 0, 0, 0},
                                      {0, 0, 0, 0, 0, 0},
                                      {0, 0.012150564, 0.012150564, 0, 0.012150564, 0.012150564},
                                      {0, 0, 0, 0, 0, 0},
                                      {0, 0, 0, 0, 0, 0},
                                      {0.012150564, 0, 0.012150564, 0.012150564, 0, 0.012150564}});
    const auto expected_up_block_vector = UblasUtilities::CreateVector(
        {0, 0, -0.26845153, 0, 0, -0.26845153, 0, 0, -0.26845153, 0, 0, -2.454414,  0, 0, -2.454414,
         0, 0, -2.454414,   0, 0, 0.26845153,  0, 0, 0.26845153,  0, 0, 0.26845153, 0, 0, 2.454414,
         0, 0, 2.454414,    0, 0, 2.454414,    0, 0, 0,           0, 0, 0});
    GeneraizedCouplingContributionTest<InterfaceThreeDimensionalSurface>(
        CreateHorizontal6Plus6NodedTriangleInterfaceElementWithUPwDoF, {CalculationContribution::UPCoupling},
        IsDiffOrderElement::Yes, 36, 6, expected_up_block_matrix, expected_up_block_vector);
}

KRATOS_TEST_CASE_IN_SUITE(UPwQuadraliteralInterfaceElement_4Plus4Element_CouplingContribution,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto expected_up_block_matrix =
        UblasUtilities::CreateMatrix({{0, 0, 0, 0, 0, 0, 0, 0},
                                      {0, 0, 0, 0, 0, 0, 0, 0},
                                      {-0.040438596, 0, 0, 0, -0.040438596, 0, 0, 0},
                                      {0, 0, 0, 0, 0, 0, 0, 0},
                                      {0, 0, 0, 0, 0, 0, 0, 0},
                                      {0, -0.040438596, 0, 0, 0, -0.040438596, 0, 0},
                                      {0, 0, 0, 0, 0, 0, 0, 0},
                                      {0, 0, 0, 0, 0, 0, 0, 0},
                                      {0, 0, -0.040438596, 0, 0, 0, -0.040438596, 0},
                                      {0, 0, 0, 0, 0, 0, 0, 0},
                                      {0, 0, 0, 0, 0, 0, 0, 0},
                                      {0, 0, 0, -0.040438596, 0, 0, 0, -0.040438596},
                                      {0, 0, 0, 0, 0, 0, 0, 0},
                                      {0, 0, 0, 0, 0, 0, 0, 0},
                                      {0.040438596, 0, 0, 0, 0.040438596, 0, 0, 0},
                                      {0, 0, 0, 0, 0, 0, 0, 0},
                                      {0, 0, 0, 0, 0, 0, 0, 0},
                                      {0, 0.040438596, 0, 0, 0, 0.040438596, 0, 0},
                                      {0, 0, 0, 0, 0, 0, 0, 0},
                                      {0, 0, 0, 0, 0, 0, 0, 0},
                                      {0, 0, 0.040438596, 0, 0, 0, 0.040438596, 0},
                                      {0, 0, 0, 0, 0, 0, 0, 0},
                                      {0, 0, 0, 0, 0, 0, 0, 0},
                                      {0, 0, 0, 0.040438596, 0, 0, 0, 0.040438596}});
    const auto expected_up_block_vector = UblasUtilities::CreateVector(
        {0, 0, -4.0842982, 0, 0, -4.0842982, 0, 0, -4.0842982, 0, 0, -4.0842982,
         0, 0, 4.0842982,  0, 0, 4.0842982,  0, 0, 4.0842982,  0, 0, 4.0842982,
         0, 0, 0,          0, 0, 0,          0, 0});
    GeneraizedCouplingContributionTest<InterfaceThreeDimensionalSurface>(
        CreateHorizontal4Plus4NodedQuadraliteralInterfaceElementWithUPwDoF,
        {CalculationContribution::UPCoupling}, IsDiffOrderElement::No, 24, 8,
        expected_up_block_matrix, expected_up_block_vector);
}

KRATOS_TEST_CASE_IN_SUITE(UPwQuadraliteralInterfaceElement_8Plus8Element_CouplingContribution,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto expected_up_block_matrix = UblasUtilities::CreateMatrix(
        {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
         {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
         {-0.0063850415, 0, 0, 0, 0, 0, 0, 0, -0.0063850415, 0, 0, 0, 0, 0, 0, 0},
         {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
         {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
         {0, -0.0063850415, 0, 0, 0, 0, 0, 0, 0, -0.0063850415, 0, 0, 0, 0, 0, 0},
         {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
         {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
         {0, 0, -0.0063850415, 0, 0, 0, 0, 0, 0, 0, -0.0063850415, 0, 0, 0, 0, 0},
         {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
         {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
         {0, 0, 0, -0.0063850415, 0, 0, 0, 0, 0, 0, 0, -0.0063850415, 0, 0, 0, 0},
         {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
         {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
         {0, 0, 0, 0, -0.034053555, 0, 0, 0, 0, 0, 0, 0, -0.034053555, 0, 0, 0},
         {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
         {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
         {0, 0, 0, 0, 0, -0.034053555, 0, 0, 0, 0, 0, 0, 0, -0.034053555, 0, 0},
         {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
         {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
         {0, 0, 0, 0, 0, 0, -0.034053555, 0, 0, 0, 0, 0, 0, 0, -0.034053555, 0},
         {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
         {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
         {0, 0, 0, 0, 0, 0, 0, -0.034053555, 0, 0, 0, 0, 0, 0, 0, -0.034053555},
         {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
         {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
         {0.0063850415, 0, 0, 0, 0, 0, 0, 0, 0.0063850415, 0, 0, 0, 0, 0, 0, 0},
         {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
         {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
         {0, 0.0063850415, 0, 0, 0, 0, 0, 0, 0, 0.0063850415, 0, 0, 0, 0, 0, 0},
         {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
         {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
         {0, 0, 0.0063850415, 0, 0, 0, 0, 0, 0, 0, 0.0063850415, 0, 0, 0, 0, 0},
         {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
         {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
         {0, 0, 0, 0.0063850415, 0, 0, 0, 0, 0, 0, 0, 0.0063850415, 0, 0, 0, 0},
         {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
         {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
         {0, 0, 0, 0, 0.034053555, 0, 0, 0, 0, 0, 0, 0, 0.034053555, 0, 0, 0},
         {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
         {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
         {0, 0, 0, 0, 0, 0.034053555, 0, 0, 0, 0, 0, 0, 0, 0.034053555, 0, 0},
         {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
         {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
         {0, 0, 0, 0, 0, 0, 0.034053555, 0, 0, 0, 0, 0, 0, 0, 0.034053555, 0},
         {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
         {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
         {0, 0, 0, 0, 0, 0, 0, 0.034053555, 0, 0, 0, 0, 0, 0, 0, 0.034053555}});
    const auto expected_up_block_vector = UblasUtilities::CreateVector(
        {0, 0, -0.64488919, 0, 0, -0.64488919, 0, 0, -0.64488919, 0, 0, -0.64488919,
         0, 0, -3.439409,   0, 0, -3.439409,   0, 0, -3.439409,   0, 0, -3.439409,
         0, 0, 0.64488919,  0, 0, 0.64488919,  0, 0, 0.64488919,  0, 0, 0.64488919,
         0, 0, 3.439409,    0, 0, 3.439409,    0, 0, 3.439409,    0, 0, 3.439409,
         0, 0, 0,           0, 0, 0,           0, 0, 0,           0, 0, 0,
         0, 0, 0,           0});
    GeneraizedCouplingContributionTest<InterfaceThreeDimensionalSurface>(
        CreateHorizontal8Plus8NodedQuadraliteralInterfaceElementWithUPwDoF,
        {CalculationContribution::UPCoupling}, IsDiffOrderElement::No, 48, 16,
        expected_up_block_matrix, expected_up_block_vector);
}

KRATOS_TEST_CASE_IN_SUITE(UPwQuadraliteralInterfaceElement_8Plus8Element_CouplingContribution_DiffOrder,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto expected_up_block_matrix = UblasUtilities::CreateMatrix(
        {{0, 0, 0, 0, 0, 0, 0, 0},
         {0, 0, 0, 0, 0, 0, 0, 0},
         {-0.0063850415, 0, 0, 0, -0.0063850415, 0, 0, 0},
         {0, 0, 0, 0, 0, 0, 0, 0},
         {0, 0, 0, 0, 0, 0, 0, 0},
         {0, -0.0063850415, 0, 0, 0, -0.0063850415, 0, 0},
         {0, 0, 0, 0, 0, 0, 0, 0},
         {0, 0, 0, 0, 0, 0, 0, 0},
         {0, 0, -0.0063850415, 0, 0, 0, -0.0063850415, 0},
         {0, 0, 0, 0, 0, 0, 0, 0},
         {0, 0, 0, 0, 0, 0, 0, 0},
         {0, 0, 0, -0.0063850415, 0, 0, 0, -0.0063850415},
         {0, 0, 0, 0, 0, 0, 0, 0},
         {0, 0, 0, 0, 0, 0, 0, 0},
         {-0.017026777, -0.017026777, 0, 0, -0.017026777, -0.017026777, 0, 0},
         {0, 0, 0, 0, 0, 0, 0, 0},
         {0, 0, 0, 0, 0, 0, 0, 0},
         {0, -0.017026777, -0.017026777, 0, 0, -0.017026777, -0.017026777, 0},
         {0, 0, 0, 0, 0, 0, 0, 0},
         {0, 0, 0, 0, 0, 0, 0, 0},
         {0, 0, -0.017026777, -0.017026777, 0, 0, -0.017026777, -0.017026777},
         {0, 0, 0, 0, 0, 0, 0, 0},
         {0, 0, 0, 0, 0, 0, 0, 0},
         {-0.017026777, 0, 0, -0.017026777, -0.017026777, 0, 0, -0.017026777},
         {0, 0, 0, 0, 0, 0, 0, 0},
         {0, 0, 0, 0, 0, 0, 0, 0},
         {0.0063850415, 0, 0, 0, 0.0063850415, 0, 0, 0},
         {0, 0, 0, 0, 0, 0, 0, 0},
         {0, 0, 0, 0, 0, 0, 0, 0},
         {0, 0.0063850415, 0, 0, 0, 0.0063850415, 0, 0},
         {0, 0, 0, 0, 0, 0, 0, 0},
         {0, 0, 0, 0, 0, 0, 0, 0},
         {0, 0, 0.0063850415, 0, 0, 0, 0.0063850415, 0},
         {0, 0, 0, 0, 0, 0, 0, 0},
         {0, 0, 0, 0, 0, 0, 0, 0},
         {0, 0, 0, 0.0063850415, 0, 0, 0, 0.0063850415},
         {0, 0, 0, 0, 0, 0, 0, 0},
         {0, 0, 0, 0, 0, 0, 0, 0},
         {0.017026777, 0.017026777, 0, 0, 0.017026777, 0.017026777, 0, 0},
         {0, 0, 0, 0, 0, 0, 0, 0},
         {0, 0, 0, 0, 0, 0, 0, 0},
         {0, 0.017026777, 0.017026777, 0, 0, 0.017026777, 0.017026777, 0},
         {0, 0, 0, 0, 0, 0, 0, 0},
         {0, 0, 0, 0, 0, 0, 0, 0},
         {0, 0, 0.017026777, 0.017026777, 0, 0, 0.017026777, 0.017026777},
         {0, 0, 0, 0, 0, 0, 0, 0},
         {0, 0, 0, 0, 0, 0, 0, 0},
         {0.017026777, 0, 0, 0.017026777, 0.017026777, 0, 0, 0.017026777}});
    const auto expected_up_block_vector = UblasUtilities::CreateVector(
        {0, 0, -0.64488919, 0, 0, -0.64488919, 0, 0, -0.64488919, 0, 0, -0.64488919,
         0, 0, -3.439409,   0, 0, -3.439409,   0, 0, -3.439409,   0, 0, -3.439409,
         0, 0, 0.64488919,  0, 0, 0.64488919,  0, 0, 0.64488919,  0, 0, 0.64488919,
         0, 0, 3.439409,    0, 0, 3.439409,    0, 0, 3.439409,    0, 0, 3.439409,
         0, 0, 0,           0, 0, 0,           0, 0});

    GeneraizedCouplingContributionTest<InterfaceThreeDimensionalSurface>(
        CreateHorizontal8Plus8NodedQuadraliteralInterfaceElementWithUPwDoF,
        {CalculationContribution::UPCoupling}, IsDiffOrderElement::Yes, 48, 8,
        expected_up_block_matrix, expected_up_block_vector);
}

} // namespace Kratos::Testing
