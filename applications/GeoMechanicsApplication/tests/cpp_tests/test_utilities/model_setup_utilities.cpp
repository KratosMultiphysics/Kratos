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

ModelPart& CreateEmptyModelPart(Model& rModel, const Geo::ConstVariableRefs& rNodalVariables)
{
    auto& r_result = rModel.CreateModelPart("Main");
    for (const auto& r_variable : rNodalVariables) {
        r_result.AddNodalSolutionStepVariable(r_variable);
    }
    return r_result;
}

} // namespace

namespace Kratos::Testing::ModelSetupUtilities
{

ModelPart& CreateModelPartWithASingle2D3NElement(Model& rModel, const Geo::ConstVariableRefs& rNodalVariables)
{
    ModelPart& result = CreateEmptyModelPart(rModel, rNodalVariables);

    result.CreateNewNode(1, 0.0, 0.0, 0.0);
    result.CreateNewNode(2, 1.0, 0.0, 0.0);
    result.CreateNewNode(3, 1.0, 1.0, 0.0);

    for (const auto& r_variable : rNodalVariables) {
        for (auto& r_node : result.Nodes()) {
            r_node.AddDof(r_variable.get());
        }
    }

    const std::vector<ModelPart::IndexType> node_ids{1, 2, 3};
    result.CreateNewElement("UPwSmallStrainElement2D3N", 1, node_ids, result.CreateNewProperties(0));

    return result;
}

ModelPart& CreateModelPartWithASingle3D4NElement(Model& rModel, const Geo::ConstVariableRefs& rNodalVariables)
{
    ModelPart& result = CreateEmptyModelPart(rModel, rNodalVariables);

    result.CreateNewNode(1, 0.0, 0.0, 0.0);
    result.CreateNewNode(2, 1.0, 0.0, 0.0);
    result.CreateNewNode(3, 0.0, 1.0, 0.0);
    result.CreateNewNode(4, 0.0, 0.0, 1.0);

    for (const auto& r_variable : rNodalVariables) {
        for (auto& r_node : result.Nodes()) {
            r_node.AddDof(r_variable.get());
        }
    }

    const std::vector<ModelPart::IndexType> node_ids{1, 2, 3, 4};
    result.CreateNewElement("UPwSmallStrainElement3D4N", 1, node_ids, result.CreateNewProperties(0));

    return result;
}

ModelPart& CreateModelPartWithASingle2D6NUPwDiffOrderElement(Model& rModel)
{
    const auto rSecondOrderVariables =
        Geo::ConstVariableRefs{std::cref(DISPLACEMENT_X), std::cref(DISPLACEMENT_Y)};
    const auto rFirstOrderVariables = Geo::ConstVariableRefs{std::cref(WATER_PRESSURE)};

    auto& r_result = CreateEmptyModelPart(rModel, rSecondOrderVariables);
    for (const auto& r_variable : rFirstOrderVariables) {
        r_result.AddNodalSolutionStepVariable(r_variable.get());
    }

    auto p_node1 = r_result.CreateNewNode(1, 0.0, 0.0, 0.0);
    auto p_node2 = r_result.CreateNewNode(2, 1.0, 0.0, 0.0);
    auto p_node3 = r_result.CreateNewNode(3, 0.0, 1.0, 0.0);
    auto p_node4 = r_result.CreateNewNode(4, 0.5, 0.0, 0.0);
    auto p_node5 = r_result.CreateNewNode(5, 0.5, 0.5, 0.0);
    auto p_node6 = r_result.CreateNewNode(6, 0.0, 0.5, 0.0);
    auto nodes   = std::vector<Node*>{p_node1.get(), p_node2.get(), p_node3.get(),
                                      p_node4.get(), p_node5.get(), p_node6.get()};

    for (const auto& r_variable : rSecondOrderVariables) {
        for (auto p_node : nodes) {
            p_node->AddDof(r_variable.get());
        }
    }

    for (const auto& r_variable : rFirstOrderVariables) {
        for (auto it = nodes.begin(); it != nodes.begin() + 3; ++it) {
            (*it)->AddDof(r_variable.get());
        }
    }

    const std::vector<ModelPart::IndexType> node_ids{1, 2, 3, 4, 5, 6};
    r_result.CreateNewElement("SmallStrainUPwDiffOrderElement2D6N", 1, node_ids,
                              r_result.CreateNewProperties(0));

    return r_result;
}

} // namespace Kratos::Testing::ModelSetupUtilities
