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

#include "element_setup_utilities.hpp"
#include "custom_conditions/Pw_point_flux_condition.h"
#include "custom_conditions/U_Pw_normal_face_load_condition.h"
#include "custom_conditions/line_load_2D_diff_order_condition.h"
#include "custom_elements/Pw_element.hpp"
#include "custom_elements/U_Pw_small_strain_element.h"
#include "custom_elements/calculation_contribution.h"
#include "custom_elements/interface_element.h"
#include "custom_elements/interface_stress_state.h"
#include "custom_elements/plane_strain_stress_state.h"
#include "custom_elements/small_strain_U_Pw_diff_order_element.h"
#include "custom_elements/three_dimensional_stress_state.h"
#include "custom_geometries/interface_geometry.hpp"
#include "geometries/hexahedra_3d_20.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/line_2d_2.h"
#include "geometries/point_3d.h"
#include "geometries/quadrilateral_2d_4.h"
#include "geometries/quadrilateral_2d_8.h"
#include "geometries/quadrilateral_3d_4.h"
#include "geometries/quadrilateral_3d_8.h"
#include "geometries/tetrahedra_3d_10.h"
#include "geometries/triangle_2d_10.h"
#include "geometries/triangle_2d_15.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_2d_6.h"

namespace Kratos::Testing
{

Condition::Pointer ElementSetupUtilities::CreateCondition(const std::string&         rType,
                                                          const PointerVector<Node>& rNodes)
{
    if (rType == "3D3NCondition") return Create3D3NCondition(rNodes);
    if (rType == "3D4NCondition") return Create3D4NCondition(rNodes);
    if (rType == "3D6NCondition") return Create3D6NCondition(rNodes);
    if (rType == "3D8NCondition") return Create3D8NCondition(rNodes);
    if (rType == "3D1NCondition") return Create3D1NCondition(rNodes);
    if (rType == "2D2NCondition") return Create2D2NCondition(rNodes);
    if (rType == "3D3NLineCondition") return Create3D3NLineCondition(rNodes);

    KRATOS_ERROR << "Condition type " << rType << " not recognized.";
}

std::vector<Point> ElementSetupUtilities::CreatePointsFor2D2NElement()
{
    return {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}};
}

std::vector<Kratos::Point> ElementSetupUtilities::CreatePointsFor2D3NElement()
{
    return {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}};
}

std::vector<Kratos::Point> ElementSetupUtilities::CreatePointsFor2D3NLineEntity()
{
    return {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.5, 0.0, 0.0}};
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

std::vector<Point> ElementSetupUtilities::CreatePointsFor3D20NElement()
{
    return {{0.0, 1.0, 1.0}, {0.0, 1.0, 0.5}, {0.5, 1.0, 1.0}, {0.0, 0.5, 1.0}, {0.0, 0.0, 1.0},
            {1.0, 1.0, 1.0}, {0.0, 1.0, 0.0}, {1.0, 0.5, 1.0}, {1.0, 1.0, 0.5}, {0.0, 0.0, 0.5},
            {0.5, 1.0, 0.0}, {0.0, 0.5, 0.0}, {0.5, 0.0, 1.0}, {1.0, 1.0, 0.0}, {1.0, 0.0, 1.0},
            {0.0, 0.0, 0.0}, {0.5, 0.0, 0.0}, {1.0, 0.5, 0.0}, {1.0, 0.0, 0.5}, {1.0, 0.0, 0.0}};
}

std::vector<Kratos::Point> ElementSetupUtilities::CreatePointsFor3D6NInterfaceElement()
{
    return {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0},
            {0.0, 0.0, 0.1}, {1.0, 0.0, 0.1}, {0.0, 1.0, 0.1}};
}

PointerVector<Node> ElementSetupUtilities::GenerateNodes(const std::vector<Point>& rPoints)
{
    PointerVector<Node> result;
    result.reserve(rPoints.size());
    for (const auto& r_point : rPoints) {
        result.push_back(make_intrusive<Node>(result.size() + 1, r_point.X(), r_point.Y(), r_point.Z()));
    }
    return result;
}

std::vector<Kratos::Point> ElementSetupUtilities::CreatePointsFor2D15NElement()
{
    return {{0.00, 0.00, 0.0}, {1.00, 0.00, 0.0}, {0.00, 1.00, 0.0}, {0.25, 0.00, 0.0},
            {0.50, 0.00, 0.0}, {0.75, 0.00, 0.0}, {0.75, 0.25, 0.0}, {0.50, 0.50, 0.0},
            {0.25, 0.75, 0.0}, {0.00, 0.75, 0.0}, {0.00, 0.50, 0.0}, {0.00, 0.25, 0.0},
            {0.25, 0.25, 0.0}, {0.50, 0.25, 0.0}, {0.25, 0.50, 0.0}};
}

std::vector<Point> ElementSetupUtilities::CreatePointsFor3D4NElement()
{
    return {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
}

std::vector<Point> ElementSetupUtilities::CreatePointsFor3D8NElement()
{
    return {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}, {0.0, 1.0, 0.0},
            {0.0, 0.0, 1.0}, {1.0, 0.0, 1.0}, {1.0, 1.0, 1.0}, {0.0, 1.0, 1.0}};
}

Element::Pointer ElementSetupUtilities::Create2D3NElement(const PointerVector<Node>& rNodes,
                                                          const Properties::Pointer& rProperties)
{
    return Create2D3NElement(rNodes, rProperties, 1);
}

Element::Pointer ElementSetupUtilities::Create2D3NElement(const PointerVector<Node>& rNodes,
                                                          const Properties::Pointer& rProperties,
                                                          std::size_t                Id)
{
    return make_intrusive<UPwSmallStrainElement<2, 3>>(
        Id, std::make_shared<Triangle2D3<Node>>(rNodes), rProperties,
        std::make_unique<PlaneStrainStressState>(), nullptr);
}

Element::Pointer ElementSetupUtilities::Create2D3NElement()
{
    return Create2D3NElement(GenerateNodes(CreatePointsFor2D3NElement()), std::make_shared<Properties>(0));
}

Element::Pointer ElementSetupUtilities::Create2D4NElement(const PointerVector<Node>& rNodes,
                                                          const Properties::Pointer& rProperties)
{
    return make_intrusive<UPwSmallStrainElement<2, 4>>(
        1, std::make_shared<Quadrilateral2D4<Node>>(rNodes), rProperties,
        std::make_unique<PlaneStrainStressState>(), nullptr);
}

Element::Pointer ElementSetupUtilities::Create2D8NElement(const PointerVector<Node>& rNodes,
                                                          const Properties::Pointer& rProperties)
{
    return make_intrusive<UPwSmallStrainElement<2, 8>>(
        1, std::make_shared<Quadrilateral2D8<Node>>(rNodes), rProperties,
        std::make_unique<PlaneStrainStressState>(), nullptr);
}

Element::Pointer ElementSetupUtilities::Create2D3NLineElement(const PointerVector<Node>& rNodes,
                                                              const Properties::Pointer& rProperties)
{
    return make_intrusive<UPwSmallStrainElement<2, 3>>(1, std::make_shared<Line2D3<Node>>(rNodes), rProperties,
                                                       std::make_unique<PlaneStrainStressState>(), nullptr);
}

Element::Pointer ElementSetupUtilities::Create2D2NElement(const PointerVector<Node>& rNodes,
                                                          const Properties::Pointer& rProperties)
{
    using enum CalculationContribution;
    const auto contributions = {Permeability, Compressibility, FluidBodyFlow};

    return make_intrusive<PwElement<2, 2>>(1, std::make_shared<Line2D2<Node>>(rNodes), rProperties,
                                           contributions, nullptr);
}

Condition::Pointer ElementSetupUtilities::Create3D3NCondition(const PointerVector<Node>& rNodes)
{
    return make_intrusive<UPwNormalFaceLoadCondition<3, 3>>(1, std::make_shared<Triangle3D3<Node>>(rNodes));
}

Condition::Pointer ElementSetupUtilities::Create3D4NCondition(const PointerVector<Node>& rNodes)
{
    return make_intrusive<UPwNormalFaceLoadCondition<3, 4>>(1, std::make_shared<Quadrilateral3D4<Node>>(rNodes));
}

Condition::Pointer ElementSetupUtilities::Create3D6NCondition(const PointerVector<Node>& rNodes)
{
    return make_intrusive<UPwNormalFaceLoadCondition<3, 6>>(1, std::make_shared<Triangle3D6<Node>>(rNodes));
}

Condition::Pointer ElementSetupUtilities::Create3D8NCondition(const PointerVector<Node>& rNodes)
{
    return make_intrusive<UPwNormalFaceLoadCondition<3, 8>>(1, std::make_shared<Quadrilateral3D8<Node>>(rNodes));
}

Condition::Pointer ElementSetupUtilities::Create3D1NCondition(const PointerVector<Node>& rNodes)
{
    return make_intrusive<PwPointFluxCondition<3, 1>>(1, std::make_shared<Point3D<Node>>(rNodes));
}

Condition::Pointer ElementSetupUtilities::Create2D2NCondition(const PointerVector<Node>& rNodes)
{
    return make_intrusive<UPwNormalFaceLoadCondition<2, 2>>(1, std::make_shared<Line2D2<Node>>(rNodes));
}

Condition::Pointer ElementSetupUtilities::Create3D3NLineCondition(const PointerVector<Node>& rNodes)
{
    return make_intrusive<UPwNormalFaceLoadCondition<3, 3>>(1, std::make_shared<Line3D3<Node>>(rNodes));
}

Element::Pointer ElementSetupUtilities::Create2D6NElement(const PointerVector<Node>& rNodes,
                                                          const Properties::Pointer& rProperties)
{
    return make_intrusive<UPwSmallStrainElement<2, 6>>(1, std::make_shared<Triangle2D6<Node>>(rNodes), rProperties,
                                                       std::make_unique<PlaneStrainStressState>(), nullptr);
}

Element::Pointer ElementSetupUtilities::Create2D6NElement()
{
    return Create2D6NElement(GenerateNodes(CreatePointsFor2D6NElement()), std::make_shared<Properties>(0));
}

Element::Pointer ElementSetupUtilities::Create2D6NDiffOrderElement(const PointerVector<Node>& rNodes,
                                                                   const Properties::Pointer& rProperties)
{
    return make_intrusive<SmallStrainUPwDiffOrderElement>(
        1, std::make_shared<Triangle2D6<Node>>(rNodes), rProperties,
        std::make_unique<PlaneStrainStressState>(), nullptr);
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
        1, std::make_shared<Triangle2D10<Node>>(rNodes), rProperties,
        std::make_unique<PlaneStrainStressState>(), nullptr);
}

Element::Pointer ElementSetupUtilities::Create2D10NElement()
{
    return Create2D10NElement(GenerateNodes(CreatePointsFor2D10NElement()), std::make_shared<Properties>(0));
}

Element::Pointer ElementSetupUtilities::Create2D15NElement(const PointerVector<Node>& rNodes,
                                                           const Properties::Pointer& rProperties)
{
    return make_intrusive<UPwSmallStrainElement<2, 15>>(
        1, std::make_shared<Triangle2D15<Node>>(rNodes), rProperties,
        std::make_unique<PlaneStrainStressState>(), nullptr);
}

Element::Pointer ElementSetupUtilities::Create2D15NElement()
{
    return Create2D15NElement(GenerateNodes(CreatePointsFor2D15NElement()), std::make_shared<Properties>(0));
}

Element::Pointer ElementSetupUtilities::Create2D4NInterfaceElement(const PointerVector<Node>& rNodes,
                                                                   const Properties::Pointer& rProperties)
{
    return make_intrusive<InterfaceElement>(1, std::make_shared<InterfaceGeometry<Line2D2<Node>>>(rNodes),
                                            rProperties, std::make_unique<Line2DInterfaceStressState>());
}

Element::Pointer ElementSetupUtilities::Create2D6NInterfaceElement(const PointerVector<Node>& rNodes,
                                                                   const Properties::Pointer& rProperties)
{
    return make_intrusive<InterfaceElement>(1, std::make_shared<InterfaceGeometry<Line2D3<Node>>>(rNodes),
                                            rProperties, std::make_unique<Line2DInterfaceStressState>());
}

Element::Pointer ElementSetupUtilities::Create3D6NInterfaceElement(const PointerVector<Node>& rNodes,
                                                                   const Properties::Pointer& rProperties)
{
    return make_intrusive<InterfaceElement>(1, std::make_shared<InterfaceGeometry<Triangle3D3<Node>>>(rNodes),
                                            rProperties, std::make_unique<SurfaceInterfaceStressState>());
}

Element::Pointer ElementSetupUtilities::Create3D12NInterfaceElement(const PointerVector<Node>& rNodes,
                                                                    const Properties::Pointer& rProperties)
{
    return make_intrusive<InterfaceElement>(1, std::make_shared<InterfaceGeometry<Triangle3D6<Node>>>(rNodes),
                                            rProperties, std::make_unique<SurfaceInterfaceStressState>());
}

Element::Pointer ElementSetupUtilities::Create3D8NInterfaceElement(const PointerVector<Node>& rNodes,
                                                                   const Properties::Pointer& rProperties,
                                                                   std::size_t Id)
{
    return make_intrusive<InterfaceElement>(
        Id, std::make_shared<InterfaceGeometry<Quadrilateral3D4<Node>>>(rNodes), rProperties,
        std::make_unique<SurfaceInterfaceStressState>());
}

Element::Pointer ElementSetupUtilities::Create3D4NElement(const PointerVector<Node>& rNodes,
                                                          const Properties::Pointer& rProperties)
{
    return make_intrusive<UPwSmallStrainElement<3, 4>>(1, std::make_shared<Tetrahedra3D4<Node>>(rNodes), rProperties,
                                                       std::make_unique<ThreeDimensionalStressState>());
}

Element::Pointer ElementSetupUtilities::Create3D10NElement(const PointerVector<Node>& rNodes,
                                                           const Properties::Pointer& rProperties)
{
    return make_intrusive<UPwSmallStrainElement<3, 10>>(
        1, std::make_shared<Tetrahedra3D10<Node>>(rNodes), rProperties,
        std::make_unique<ThreeDimensionalStressState>(), nullptr);
}

Element::Pointer ElementSetupUtilities::Create3D10NElement()
{
    return Create3D10NElement(GenerateNodes(CreatePointsFor3D10NElement()), std::make_shared<Properties>(0));
}

Element::Pointer ElementSetupUtilities::Create3D8NElement(const PointerVector<Node>& rNodes,
                                                          const Properties::Pointer& rProperties)
{
    return Create3D8NElement(rNodes, rProperties, 1);
}

Element::Pointer ElementSetupUtilities::Create3D8NElement(const PointerVector<Node>& rNodes,
                                                          const Properties::Pointer& rProperties,
                                                          std::size_t                Id)
{
    return make_intrusive<UPwSmallStrainElement<3, 8>>(Id, std::make_shared<Hexahedra3D8<Node>>(rNodes), rProperties,
                                                       std::make_unique<ThreeDimensionalStressState>());
}

Element::Pointer ElementSetupUtilities::Create3D20NElement(const PointerVector<Node>& rNodes,
                                                           const Properties::Pointer& rProperties)
{
    return make_intrusive<UPwSmallStrainElement<3, 20>>(
        1, std::make_shared<Hexahedra3D20<Node>>(rNodes), rProperties,
        std::make_unique<ThreeDimensionalStressState>());
}

Condition::Pointer ElementSetupUtilities::Create2D3NLineCondition(const PointerVector<Node>& rNodes,
                                                                  const Properties::Pointer& rProperties)
{
    return make_intrusive<LineLoad2DDiffOrderCondition>(1, std::make_shared<Line2D3<Node>>(rNodes), rProperties);
}

Condition::Pointer ElementSetupUtilities::Create2D3NLineCondition()
{
    return Create2D3NLineCondition(GenerateNodes(CreatePointsFor2D3NLineEntity()),
                                   std::make_shared<Properties>(0));
}

} // namespace Kratos::Testing