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

namespace Kratos::Testing::ModelSetupUtilities
{

ModelPart& CreateModelPartWithASingle2D3NElement(Model& rModel, const Geo::ConstVariableRefs& rNodalVariables)
{
    ModelPart& result = rModel.CreateModelPart("Main");
    for (const auto& r_variable : rNodalVariables) {
        result.AddNodalSolutionStepVariable(r_variable);
    }

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

ModelPart& CreateModelPartWithASingle3D4NElement(Model& rModel)
{
    ModelPart& result = rModel.CreateModelPart("Main");

    result.CreateNewNode(1, 0.0, 0.0, 0.0);
    result.CreateNewNode(2, 1.0, 0.0, 0.0);
    result.CreateNewNode(3, 0.0, 1.0, 0.0);
    result.CreateNewNode(4, 0.0, 0.0, 1.0);

    const std::vector<ModelPart::IndexType> node_ids{1, 2, 3, 4};
    result.CreateNewElement("UPwSmallStrainElement3D4N", 1, node_ids, result.CreateNewProperties(0));

    return result;
}

} // namespace Kratos::Testing::ModelSetupUtilities
