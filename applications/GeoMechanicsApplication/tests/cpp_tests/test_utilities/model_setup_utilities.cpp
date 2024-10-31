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
#include "model_setup_utilities.h"
#include "containers/model.h"
#include "includes/model_part.h"
#include "includes/properties.h"

#include <custom_elements/U_Pw_small_strain_element.hpp>
#include <custom_elements/plane_strain_stress_state.h>
#include <custom_elements/small_strain_U_Pw_diff_order_element.hpp>
#include <custom_elements/three_dimensional_stress_state.h>
#include <geometries/tetrahedra_3d_10.h>
#include <geometries/tetrahedra_3d_4.h>
#include <geometries/triangle_2d_3.h>
#include <geometries/triangle_2d_6.h>

namespace
{

using namespace Kratos;

void AddNodalVariablesToModelPart(ModelPart& rModelPart, const Geo::ConstVariableRefs& rNodalVariables)
{
    for (const auto& r_variable : rNodalVariables) {
        rModelPart.AddNodalSolutionStepVariable(r_variable.get());
    }
}

PointerVector<Node> CreateNewNodes(ModelPart& rModelPart, const std::vector<Point>& rPoints)
{
    PointerVector<Node> nodes;
    auto                NodeIndex = rModelPart.NumberOfNodes();
    for (const auto& r_point : rPoints) {
        ++NodeIndex;
        nodes.push_back(rModelPart.CreateNewNode(NodeIndex, r_point.X(), r_point.Y(), r_point.Z()));
    }

    return nodes;
}

template <typename InputIt>
void AddDofsToNodes(InputIt NodeRangeBegin, InputIt NodeRangeEnd, const Geo::ConstVariableRefs& rNodalVariables)
{
    for (const auto& r_variable : rNodalVariables) {
        for (auto it = NodeRangeBegin; it != NodeRangeEnd; ++it) {
            it->AddDof(r_variable.get());
        }
    }
}

template <typename NodeRange>
void AddDofsToNodes(const NodeRange& rNodeRange, const Geo::ConstVariableRefs& rNodalVariables)
{
    AddDofsToNodes(std::begin(rNodeRange), std::end(rNodeRange), rNodalVariables);
}

} // namespace

namespace Kratos::Testing::ModelSetupUtilities
{

ModelPart& CreateModelPartWithASingle2D3NElement(Model& rModel, const Geo::ConstVariableRefs& rNodalVariables)
{
    ModelPart& result = rModel.CreateModelPart("Main");
    AddNodalVariablesToModelPart(result, rNodalVariables);

    auto nodes = CreateNewNodes(result, {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}});
    AddDofsToNodes(result.Nodes(), rNodalVariables);

    auto element = make_intrusive<UPwSmallStrainElement<2, 3>>(
        1, Kratos::make_shared<Triangle2D3<Node>>(nodes), result.CreateNewProperties(0),
        std::make_unique<PlaneStrainStressState>());

    result.AddElement(element);

    return result;
}

ModelPart& CreateModelPartWithASingle3D4NElement(Model& rModel, const Geo::ConstVariableRefs& rNodalVariables)
{
    ModelPart& result = rModel.CreateModelPart("Main");
    AddNodalVariablesToModelPart(result, rNodalVariables);

    auto nodes =
        CreateNewNodes(result, {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}});
    AddDofsToNodes(result.Nodes(), rNodalVariables);

    auto element = make_intrusive<UPwSmallStrainElement<3, 4>>(
        1, Kratos::make_shared<Tetrahedra3D4<Node>>(nodes), result.CreateNewProperties(0),
        std::make_unique<ThreeDimensionalStressState>());

    result.AddElement(element);

    return result;
}

ModelPart& CreateModelPartWithASingle2D6NDiffOrderElement(Model& rModel)
{
    // similar to a 1D-consolidation test
    auto& result = rModel.CreateModelPart("Main");
    const auto variables = Geo::ConstVariableRefs{std::cref(DISPLACEMENT_X), std::cref(DISPLACEMENT_Y),
                                                  std::cref(DISPLACEMENT_Z), std::cref(WATER_PRESSURE)};
    AddNodalVariablesToModelPart(result, variables);
    auto nodes = CreateNewNodes(
        result,
        {{0.0, 0.0, 0.0}, {0.0, -0.05, 0.0}, {0.05, 0.0, 0.0}, {0.0, -0.025, 0.0}, {0.025, -0.025, 0.0}, {0.025, 0.0, 0.0}});

    AddDofsToNodes(result.Nodes(), variables);

    auto element = make_intrusive<SmallStrainUPwDiffOrderElement>(
        1, Kratos::make_shared<Triangle2D6<Node>>(nodes), result.CreateNewProperties(0),
        std::make_unique<PlaneStrainStressState>());

    result.AddElement(element);
    return result;
}

ModelPart& CreateModelPartWithASingle2D6NUPwDiffOrderElement(Model& rModel)
{
    auto&      r_result = rModel.CreateModelPart("Main");
    const auto second_order_variables =
        Geo::ConstVariableRefs{std::cref(DISPLACEMENT_X), std::cref(DISPLACEMENT_Y)};
    AddNodalVariablesToModelPart(r_result, second_order_variables);
    const auto first_order_variables = Geo::ConstVariableRefs{std::cref(WATER_PRESSURE)};
    AddNodalVariablesToModelPart(r_result, first_order_variables);

    auto node_pointers = CreateNewNodes(
        r_result,
        {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.5, 0.0, 0.0}, {0.5, 0.5, 0.0}, {0.0, 0.5, 0.0}});

    const auto nodes = r_result.Nodes();
    AddDofsToNodes(nodes, second_order_variables);
    AddDofsToNodes(nodes.begin(), nodes.begin() + 3, first_order_variables);

    const std::vector<ModelPart::IndexType> node_ids{1, 2, 3, 4, 5, 6};

    auto element = make_intrusive<SmallStrainUPwDiffOrderElement>(
        1, Kratos::make_shared<Triangle2D6<Node>>(node_pointers), r_result.CreateNewProperties(0),
        std::make_unique<PlaneStrainStressState>());

    r_result.AddElement(element);

    return r_result;
}

ModelPart& CreateModelPartWithASingle3D10NUPwDiffOrderElement(Model& rModel)
{
    auto&      r_result               = rModel.CreateModelPart("Main");
    const auto second_order_variables = Geo::ConstVariableRefs{
        std::cref(DISPLACEMENT_X), std::cref(DISPLACEMENT_Y), std::cref(DISPLACEMENT_Z)};
    AddNodalVariablesToModelPart(r_result, second_order_variables);
    const auto first_order_variables = Geo::ConstVariableRefs{std::cref(WATER_PRESSURE)};
    AddNodalVariablesToModelPart(r_result, first_order_variables);

    auto node_pointers = CreateNewNodes(r_result, {{0.0, 0.0, 0.0},
                                                   {1.0, 0.0, 0.0},
                                                   {0.0, 1.0, 0.0},
                                                   {0.0, 0.0, 1.0},
                                                   {0.5, 0.0, 0.0},
                                                   {0.5, 0.5, 0.0},
                                                   {0.0, 0.5, 0.0},
                                                   {0.0, 0.0, 0.5},
                                                   {0.5, 0.0, 0.5},
                                                   {0.0, 0.5, 0.5}});

    const auto nodes = r_result.Nodes();
    AddDofsToNodes(nodes, second_order_variables);
    AddDofsToNodes(nodes.begin(), nodes.begin() + 4, first_order_variables);

    auto element = make_intrusive<SmallStrainUPwDiffOrderElement>(
        1, Kratos::make_shared<Tetrahedra3D10<Node>>(node_pointers),
        r_result.CreateNewProperties(0), std::make_unique<ThreeDimensionalStressState>());

    r_result.AddElement(element);

    return r_result;
}

Triangle2D3<Node> Create2D3NTriangleGeometry()
{
    const auto node_1 = make_intrusive<Node>(1, 0.0, 0.0, 0.0);
    const auto node_2 = make_intrusive<Node>(2, 1.0, 0.0, 0.0);
    const auto node_3 = make_intrusive<Node>(3, 1.0, 1.0, 0.0);
    return {node_1, node_2, node_3};
}

Tetrahedra3D4<Node> Create3D4NTetrahedraGeometry()
{
    const auto node_1 = make_intrusive<Node>(1, 0.0, 0.0, 0.0);
    const auto node_2 = make_intrusive<Node>(2, 1.0, 0.0, 0.0);
    const auto node_3 = make_intrusive<Node>(3, 0.0, 1.0, 0.0);
    const auto node_4 = make_intrusive<Node>(3, 0.0, 0.0, 1.0);
    return {node_1, node_2, node_3, node_4};
}

} // namespace Kratos::Testing::ModelSetupUtilities
