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

#include "element_setup_utilities.h"
#include "custom_elements/U_Pw_small_strain_element.hpp"
#include "custom_elements/plane_strain_stress_state.h"
#include "custom_elements/small_strain_U_Pw_diff_order_element.hpp"
#include "custom_elements/three_dimensional_stress_state.h"
#include "geometries/tetrahedra_3d_10.h"
#include "geometries/triangle_2d_10.h"
#include "geometries/triangle_2d_15.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_2d_6.h"

namespace
{

using namespace Kratos;

PointerVector<Node> GenerateNodes(const std::vector<Point>& rPoints)
{
    PointerVector<Node> nodes;
    nodes.reserve(rPoints.size());
    for (const auto& r_point : rPoints) {
        nodes.push_back(make_intrusive<Node>(nodes.size() + 1, r_point.X(), r_point.Y(), r_point.Z()));
    }
    return nodes;
}

} // namespace

namespace Kratos::Testing
{

std::vector<Kratos::Point> ElementSetupUtilities::CreatePointsFor2D3NElement()
{
    return {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}};
}

std::vector<Kratos::Point> ElementSetupUtilities::CreatePointsFor2D6NElement()
{
    return {{0.0, 0.0, 0.0},    {0.0, -0.05, 0.0},    {0.05, 0.0, 0.0},
            {0.0, -0.025, 0.0}, {0.025, -0.025, 0.0}, {0.025, 0.0, 0.0}};
}

std::vector<Kratos::Point> ElementSetupUtilities::CreatePointsFor2D10NElement()
{
    return {{0.0, 0.0, 0.0},
            {1.0, 0.0, 0.0},
            {0.0, 1.0, 0.0},
            {1.0 / 3.0, 0.0, 0.0},
            {2.0 / 3.0, 0.0, 0.0},
            {2.0 / 3.0, 1.0 / 3.0, 0.0},
            {1.0 / 3.0, 2.0 / 3.0, 0.0},
            {0.0, 2.0 / 3.0, 0.0},
            {0.0, 1.0 / 3.0, 0.0},
            {1.0 / 3.0, 1.0 / 3.0, 0.0}};
}

std::vector<Kratos::Point> ElementSetupUtilities::CreatePointsFor3D10NElement()
{
    return {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}, {0.5, 0.0, 0.0},
            {0.5, 0.5, 0.0}, {0.0, 0.5, 0.0}, {0.0, 0.0, 0.5}, {0.5, 0.0, 0.5}, {0.0, 0.5, 0.5}};
}

std::vector<Kratos::Point> ElementSetupUtilities::CreatePointsFor2D15NElement()
{
    return {{0.00, 0.00, 0.0}, {1.00, 0.00, 0.0}, {0.00, 1.00, 0.0}, {0.25, 0.00, 0.0},
            {0.50, 0.00, 0.0}, {0.75, 0.00, 0.0}, {0.75, 0.25, 0.0}, {0.50, 0.50, 0.0},
            {0.25, 0.75, 0.0}, {0.00, 0.75, 0.0}, {0.00, 0.50, 0.0}, {0.00, 0.25, 0.0},
            {0.25, 0.25, 0.0}, {0.50, 0.25, 0.0}, {0.25, 0.50, 0.0}};
}

Element::Pointer ElementSetupUtilities::Create2D3NElement(const PointerVector<Node>& rNodes,
                                                          const Properties::Pointer& rProperties)
{
    return make_intrusive<UPwSmallStrainElement<2, 3>>(1, Kratos::make_shared<Triangle2D3<Node>>(rNodes), rProperties,
                                                       std::make_unique<PlaneStrainStressState>());
}

Element::Pointer ElementSetupUtilities::Create2D3NElement()
{
    return Create2D3NElement(GenerateNodes(CreatePointsFor2D3NElement()), std::make_shared<Properties>(0));
}

Element::Pointer ElementSetupUtilities::Create2D3NElement(const Geo::ConstVariableRefs& rVariableRefs)
{
    auto result          = Create2D3NElement();
    auto p_variable_list = make_intrusive<VariablesList>();
    for (const auto& r_variable_ref : rVariableRefs) {
        p_variable_list->Add(r_variable_ref);
    }

    for (auto& r_node : result->GetGeometry()) {
        r_node.SetSolutionStepVariablesList(p_variable_list);
    }

    return result;
}

Element::Pointer ElementSetupUtilities::Create2D6NElement(const PointerVector<Node>& rNodes,
                                                          const Properties::Pointer& rProperties)
{
    return make_intrusive<UPwSmallStrainElement<2, 6>>(1, Kratos::make_shared<Triangle2D6<Node>>(rNodes), rProperties,
                                                       std::make_unique<PlaneStrainStressState>());
}

Element::Pointer ElementSetupUtilities::Create2D6NElement()
{
    return Create2D6NElement(GenerateNodes(CreatePointsFor2D6NElement()), std::make_shared<Properties>(0));
}

Element::Pointer ElementSetupUtilities::Create2D6NDiffOrderElement(const PointerVector<Node>& rNodes,
                                                                   const Properties::Pointer& rProperties)
{
    return make_intrusive<SmallStrainUPwDiffOrderElement>(1, std::make_shared<Triangle2D6<Node>>(rNodes), rProperties,
                                                          std::make_unique<PlaneStrainStressState>());
}

Element::Pointer ElementSetupUtilities::Create2D6NDiffOrderElement()
{
    return Create2D6NDiffOrderElement(GenerateNodes(CreatePointsFor2D6NElement()),
                                      std::make_shared<Properties>(0));
}

Element::Pointer ElementSetupUtilities::Create2D10NElement(const PointerVector<Node>& rNodes,
                                                           const Properties::Pointer& rProperties)
{
    return make_intrusive<UPwSmallStrainElement<2, 10>>(
        1, Kratos::make_shared<Triangle2D10<Node>>(rNodes), rProperties,
        std::make_unique<PlaneStrainStressState>());
}

Element::Pointer ElementSetupUtilities::Create2D10NElement()
{
    return Create2D10NElement(GenerateNodes(CreatePointsFor2D10NElement()), std::make_shared<Properties>(0));
}

Element::Pointer ElementSetupUtilities::Create2D15NElement(const PointerVector<Node>& rNodes,
                                                           const Properties::Pointer& rProperties)
{
    return make_intrusive<UPwSmallStrainElement<2, 15>>(
        1, Kratos::make_shared<Triangle2D15<Node>>(rNodes), rProperties,
        std::make_unique<PlaneStrainStressState>());
}

Element::Pointer ElementSetupUtilities::Create2D15NElement()
{
    return Create2D15NElement(GenerateNodes(CreatePointsFor2D15NElement()), std::make_shared<Properties>(0));
}

Element::Pointer ElementSetupUtilities::Create3D10NElement(const PointerVector<Node>& rNodes,
                                                           const Properties::Pointer& rProperties)
{
    return make_intrusive<UPwSmallStrainElement<3, 10>>(
        1, Kratos::make_shared<Tetrahedra3D10<Node>>(rNodes), rProperties,
        std::make_unique<ThreeDimensionalStressState>());
}

Element::Pointer ElementSetupUtilities::Create3D10NElement()
{
    return Create3D10NElement(GenerateNodes(CreatePointsFor3D10NElement()), std::make_shared<Properties>(0));
}

} // namespace Kratos::Testing