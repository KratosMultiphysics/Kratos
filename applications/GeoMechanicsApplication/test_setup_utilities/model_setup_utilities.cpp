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

#include "custom_elements/plane_strain_stress_state.h"
#include "custom_elements/small_strain_U_Pw_diff_order_element.hpp"
#include "custom_elements/three_dimensional_stress_state.h"
#include "element_setup_utilities.h"
#include "geometries/tetrahedra_3d_10.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_2d_6.h"

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

namespace Kratos::Testing
{

PointerVector<Node> ModelSetupUtilities::CreateNodes(ModelPart& rModelPart, const NodeDefinitionVector& rNodeDefinitions)
{
    auto result = PointerVector<Node>{};
    result.reserve(rNodeDefinitions.size());
    for (const auto& r_node_definition : rNodeDefinitions) {
        result.push_back(rModelPart.CreateNewNode(r_node_definition.id, r_node_definition.position.X(),
                                                  r_node_definition.position.Y(),
                                                  r_node_definition.position.Z()));
    }
    return result;
}

ModelPart& ModelSetupUtilities::CreateModelPartWithASingle2D3NElement(Model& rModel,
                                                                      const Geo::ConstVariableRefs& rNodalVariables)
{
    auto& r_result = rModel.CreateModelPart("Main");
    AddNodalVariablesToModelPart(r_result, rNodalVariables);

    auto nodes = CreateNewNodes(r_result, ElementSetupUtilities::CreatePointsFor2D3NElement());
    AddDofsToNodes(r_result.Nodes(), rNodalVariables);

    auto element = ElementSetupUtilities::Create2D3NElement(nodes, r_result.CreateNewProperties(0));

    r_result.AddElement(element);

    return r_result;
}

ModelPart& ModelSetupUtilities::CreateModelPartWithASingle3D4NElement(Model& rModel,
                                                                      const Geo::ConstVariableRefs& rNodalVariables)
{
    auto& r_result = rModel.CreateModelPart("Main");
    AddNodalVariablesToModelPart(r_result, rNodalVariables);

    auto nodes = CreateNewNodes(r_result, {ElementSetupUtilities::CreatePointsFor3D4NElement()});
    AddDofsToNodes(r_result.Nodes(), rNodalVariables);

    auto element = ElementSetupUtilities::Create3D4NElement(nodes, r_result.CreateNewProperties(0));
    r_result.AddElement(element);

    return r_result;
}

ModelPart& ModelSetupUtilities::CreateModelPartWithASingle3D8NElement(Model& rModel,
                                                                      const Geo::ConstVariableRefs& rNodalVariables)
{
    auto& r_result = rModel.CreateModelPart("Main");
    AddNodalVariablesToModelPart(r_result, rNodalVariables);

    auto nodes = CreateNewNodes(r_result, ElementSetupUtilities::CreatePointsFor3D8NElement());
    AddDofsToNodes(r_result.Nodes(), rNodalVariables);

    auto element = ElementSetupUtilities::Create3D8NElement(nodes, r_result.CreateNewProperties(0));
    r_result.AddElement(element);

    return r_result;
}

ModelPart& ModelSetupUtilities::CreateModelPartWithASingle3D20NElement(Model& rModel,
                                                                       const Geo::ConstVariableRefs& rNodalVariables)
{
    auto& r_result = rModel.CreateModelPart("Main");
    AddNodalVariablesToModelPart(r_result, rNodalVariables);

    auto nodes = CreateNewNodes(r_result, ElementSetupUtilities::CreatePointsFor3D20NElement());
    AddDofsToNodes(r_result.Nodes(), rNodalVariables);

    auto element = ElementSetupUtilities::Create3D20NElement(nodes, r_result.CreateNewProperties(0));
    r_result.AddElement(element);

    return r_result;
}

ModelPart& ModelSetupUtilities::CreateModelPartWithASingle2D2NElement(Model& rModel,
                                                                      const Geo::ConstVariableRefs& rNodalVariables)
{
    auto& r_result = rModel.CreateModelPart("Main");
    AddNodalVariablesToModelPart(r_result, rNodalVariables);

    auto nodes = CreateNewNodes(r_result, ElementSetupUtilities::CreatePointsFor2D2NElement());
    AddDofsToNodes(r_result.Nodes(), rNodalVariables);

    auto element = ElementSetupUtilities::Create2D2NElement(nodes, r_result.CreateNewProperties(0));
    r_result.AddElement(element);

    return r_result;
}

ModelPart& ModelSetupUtilities::CreateModelPartWithASingle2D10NElement(Model& rModel,
                                                                       const Geo::ConstVariableRefs& rNodalVariables)
{
    auto& r_result = rModel.CreateModelPart("Main");
    AddNodalVariablesToModelPart(r_result, rNodalVariables);

    auto nodes = CreateNewNodes(r_result, ElementSetupUtilities::CreatePointsFor2D10NElement());
    AddDofsToNodes(r_result.Nodes(), rNodalVariables);

    auto element = ElementSetupUtilities::Create2D10NElement(nodes, r_result.CreateNewProperties(0));

    r_result.AddElement(element);

    return r_result;
}

ModelPart& ModelSetupUtilities::CreateModelPartWithASingle2D15NElement(Model& rModel,
                                                                       const Geo::ConstVariableRefs& rNodalVariables)
{
    auto& r_result = rModel.CreateModelPart("Main");
    AddNodalVariablesToModelPart(r_result, rNodalVariables);

    auto nodes = CreateNewNodes(r_result, ElementSetupUtilities::CreatePointsFor2D15NElement());
    AddDofsToNodes(r_result.Nodes(), rNodalVariables);

    auto element = ElementSetupUtilities::Create2D15NElement(nodes, r_result.CreateNewProperties(0));

    r_result.AddElement(element);

    return r_result;
}

ModelPart& ModelSetupUtilities::CreateModelPartWithASingle3D6NInterfaceElement(Model& rModel,
                                                                               const Geo::ConstVariableRefs& rNodalVariables)
{
    auto& r_result = rModel.CreateModelPart("Main");
    AddNodalVariablesToModelPart(r_result, rNodalVariables);

    auto nodes = CreateNewNodes(r_result, ElementSetupUtilities::CreatePointsFor3D6NInterfaceElement());
    AddDofsToNodes(r_result.Nodes(), rNodalVariables);

    auto element =
        ElementSetupUtilities::Create3D6NInterfaceElement(nodes, r_result.CreateNewProperties(0));

    r_result.AddElement(element);

    return r_result;
}

ModelPart& ModelSetupUtilities::CreateModelPartWithASingle2D6NDiffOrderElement(Model& rModel)
{
    // similar to a 1D-consolidation test
    auto& r_result = rModel.CreateModelPart("Main");
    const auto variables = Geo::ConstVariableRefs{std::cref(DISPLACEMENT_X), std::cref(DISPLACEMENT_Y),
                                                  std::cref(DISPLACEMENT_Z), std::cref(WATER_PRESSURE)};
    AddNodalVariablesToModelPart(r_result, variables);
    auto nodes = CreateNewNodes(r_result, ElementSetupUtilities::CreatePointsFor2D6NElement());

    AddDofsToNodes(r_result.Nodes(), variables);

    auto element = make_intrusive<SmallStrainUPwDiffOrderElement>(
        1, Kratos::make_shared<Triangle2D6<Node>>(nodes), r_result.CreateNewProperties(0),
        std::make_unique<PlaneStrainStressState>());

    r_result.AddElement(element);
    return r_result;
}

ModelPart& ModelSetupUtilities::CreateModelPartWithASingle2D6NUPwDiffOrderElement(Model& rModel)
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
        std::make_unique<PlaneStrainStressState>(), nullptr);

    r_result.AddElement(element);

    return r_result;
}

ModelPart& ModelSetupUtilities::CreateModelPartWithASingle3D10NUPwDiffOrderElement(Model& rModel)
{
    auto&      r_result               = rModel.CreateModelPart("Main");
    const auto second_order_variables = Geo::ConstVariableRefs{
        std::cref(DISPLACEMENT_X), std::cref(DISPLACEMENT_Y), std::cref(DISPLACEMENT_Z)};
    AddNodalVariablesToModelPart(r_result, second_order_variables);
    const auto first_order_variables = Geo::ConstVariableRefs{std::cref(WATER_PRESSURE)};
    AddNodalVariablesToModelPart(r_result, first_order_variables);

    auto node_pointers = CreateNewNodes(r_result, {ElementSetupUtilities::CreatePointsFor3D10NElement()});

    const auto nodes = r_result.Nodes();
    AddDofsToNodes(nodes, second_order_variables);
    AddDofsToNodes(nodes.begin(), nodes.begin() + 4, first_order_variables);

    auto element = make_intrusive<SmallStrainUPwDiffOrderElement>(
        1, Kratos::make_shared<Tetrahedra3D10<Node>>(node_pointers),
        r_result.CreateNewProperties(0), std::make_unique<ThreeDimensionalStressState>(), nullptr);

    r_result.AddElement(element);

    return r_result;
}

Triangle2D3<Node> ModelSetupUtilities::Create2D3NTriangleGeometry()
{
    const auto node_1 = make_intrusive<Node>(1, 0.0, 0.0, 0.0);
    const auto node_2 = make_intrusive<Node>(2, 1.0, 0.0, 0.0);
    const auto node_3 = make_intrusive<Node>(3, 1.0, 1.0, 0.0);
    return {node_1, node_2, node_3};
}

Tetrahedra3D4<Node> ModelSetupUtilities::Create3D4NTetrahedraGeometry()
{
    const auto node_1 = make_intrusive<Node>(1, 0.0, 0.0, 0.0);
    const auto node_2 = make_intrusive<Node>(2, 1.0, 0.0, 0.0);
    const auto node_3 = make_intrusive<Node>(3, 0.0, 1.0, 0.0);
    const auto node_4 = make_intrusive<Node>(3, 0.0, 0.0, 1.0);
    return {node_1, node_2, node_3, node_4};
}

void ModelSetupUtilities::CreateNumberOfNewNodes(ModelPart& rModelPart, std::size_t NumberOfNodes)
{
    for (std::size_t i = 0; i < NumberOfNodes; ++i) {
        rModelPart.CreateNewNode(i + 1, 0.0, 0.0, 0.0);
    }
}

PointerVector<Node> ModelSetupUtilities::GetNodesFromIds(ModelPart&                      rModelPart,
                                                         const std::vector<std::size_t>& rNodeIds)
{
    PointerVector<Node> result(rNodeIds.size());
    std::ranges::transform(rNodeIds, result.ptr_begin(),
                           [&rModelPart](auto Id) { return rModelPart.pGetNode(Id); });
    return result;
}

} // namespace Kratos::Testing
