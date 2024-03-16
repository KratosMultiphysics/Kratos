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

namespace
{

using namespace Kratos;

void AddNodalVariablesToModelPart(ModelPart& rModelPart, const Geo::ConstVariableRefs& rNodalVariables)
{
    for (const auto& r_variable : rNodalVariables) {
        rModelPart.AddNodalSolutionStepVariable(r_variable.get());
    }
}

void CreateNewNodes(ModelPart& rModelPart, const std::vector<Point>& rPoints)
{
    auto NodeIndex = rModelPart.NumberOfNodes();
    for (const auto& r_point : rPoints) {
        ++NodeIndex;
        rModelPart.CreateNewNode(NodeIndex, r_point.X(), r_point.Y(), r_point.Z());
    }
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

    CreateNewNodes(result, {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}});
    AddDofsToNodes(result.Nodes(), rNodalVariables);

    const std::vector<ModelPart::IndexType> node_ids{1, 2, 3};
    result.CreateNewElement("UPwSmallStrainElement2D3N", 1, node_ids, result.CreateNewProperties(0));

    return result;
}

ModelPart& CreateModelPartWithASingle3D4NElement(Model& rModel, const Geo::ConstVariableRefs& rNodalVariables)
{
    ModelPart& result = rModel.CreateModelPart("Main");
    AddNodalVariablesToModelPart(result, rNodalVariables);

    CreateNewNodes(result, {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}});
    AddDofsToNodes(result.Nodes(), rNodalVariables);

    const std::vector<ModelPart::IndexType> node_ids{1, 2, 3, 4};
    result.CreateNewElement("UPwSmallStrainElement3D4N", 1, node_ids, result.CreateNewProperties(0));

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

    CreateNewNodes(
        r_result,
        {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.5, 0.0, 0.0}, {0.5, 0.5, 0.0}, {0.0, 0.5, 0.0}});

    const auto nodes = r_result.Nodes();
    AddDofsToNodes(nodes, second_order_variables);
    AddDofsToNodes(nodes.begin(), nodes.begin() + 3, first_order_variables);

    const std::vector<ModelPart::IndexType> node_ids{1, 2, 3, 4, 5, 6};
    r_result.CreateNewElement("SmallStrainUPwDiffOrderElement2D6N", 1, node_ids,
                              r_result.CreateNewProperties(0));

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

    CreateNewNodes(r_result, {{0.0, 0.0, 0.0},
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

    const std::vector<ModelPart::IndexType> node_ids{1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    r_result.CreateNewElement("SmallStrainUPwDiffOrderElement3D10N", 1, node_ids,
                              r_result.CreateNewProperties(0));

    return r_result;
}

} // namespace Kratos::Testing::ModelSetupUtilities
