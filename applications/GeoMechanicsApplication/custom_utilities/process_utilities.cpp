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
//                   Markelov Gennady
//

// Project includes
#include "process_utilities.h"
#include "containers/model.h"
#include "includes/kratos_parameters.h"

namespace Kratos
{
std::vector<std::reference_wrapper<ModelPart>> ProcessUtilities::GetModelPartsFromSettings(
    Model& rModel, const Parameters& rProcessSettings, const std::string& rProcessInfo)
{
    KRATOS_ERROR_IF_NOT(rProcessSettings.Has("model_part_name") ||
                        rProcessSettings.Has("model_part_name_list"))
        << "Please specify 'model_part_name' or 'model_part_name_list' for " << rProcessInfo;

    KRATOS_ERROR_IF(rProcessSettings.Has("model_part_name") &&
                    rProcessSettings.Has("model_part_name_list"))
        << "The parameters 'model_part_name' and 'model_part_name_list' are mutually exclusive for "
        << rProcessInfo;

    const auto model_part_names = rProcessSettings.Has("model_part_name_list")
                                      ? rProcessSettings["model_part_name_list"].GetStringArray()
                                      : std::vector{rProcessSettings["model_part_name"].GetString()};
    KRATOS_ERROR_IF(model_part_names.empty()) << "The parameters 'model_part_name_list' needs "
                                                 "to contain at least one model part name for "
                                              << rProcessInfo;

    std::vector<std::reference_wrapper<ModelPart>> result;
    result.reserve(model_part_names.size());
    std::ranges::transform(
        model_part_names, std::back_inserter(result),
        [&rModel](const auto& rName) -> ModelPart& { return rModel.GetModelPart(rName); });
    return result;
}

}; // namespace Kratos
