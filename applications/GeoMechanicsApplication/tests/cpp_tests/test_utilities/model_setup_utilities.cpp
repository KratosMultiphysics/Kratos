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
#include "includes/model_part.h"
#include "containers/model.h"

namespace Kratos::Testing::ModelSetupUtilities
{

ModelPart& CreateModelPartWithASingle2D3NElement(Model& rModel)
{
    ModelPart& result = rModel.CreateModelPart("Main");

    result.CreateNewNode(1, 0.0, 0.0, 0.0);
    result.CreateNewNode(2, 1.0, 0.0, 0.0);
    result.CreateNewNode(3, 1.0, 1.0, 0.0);

    std::vector<ModelPart::IndexType> node_ids{1, 2, 3};
    result.CreateNewElement("UPwSmallStrainElement2D3N", 1, node_ids, result.CreateNewProperties(0));

    return result;
}

}
