//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "utilities/model_utilities.h"

namespace Kratos
{
namespace ModelUtilities
{
ModelPart& GetModelPartFromModelAndSettings(
    Model& rModel,
    const Parameters& rParameters
    )
{
    const std::string& r_model_part_name = rParameters["model_part_name"].GetString();
    const std::string sub_model_part_name = rParameters.Has("sub_model_part_name") ? rParameters["sub_model_part_name"].GetString() : "";

    if (sub_model_part_name == "") {
        return rModel.GetModelPart(r_model_part_name);
    } else {
        return rModel.GetModelPart(r_model_part_name).GetSubModelPart(sub_model_part_name);
    }
}

/***********************************************************************************/
/***********************************************************************************/

} // namespace ModelUtilities
} // namespace Kratos
