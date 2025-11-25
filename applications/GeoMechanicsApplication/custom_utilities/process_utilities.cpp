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

namespace
{
std::set<std::string> ExtractModelPartNames(const auto& rProcessList, std::string_view RootName, std::string_view Prefix)
{
    std::set<std::string> result;
    for (const auto& r_process : rProcessList) {
        std::set<std::string> model_part_names;
        if (r_process.Has("Parameters")) {
            if (r_process["Parameters"].Has("model_part_name")) {
                model_part_names.insert(r_process["Parameters"]["model_part_name"].GetString());
            } else if (r_process["Parameters"].Has("model_part_name_list")) {
                for (const auto& name : r_process["Parameters"]["model_part_name_list"].GetStringArray()) {
                    model_part_names.insert(name);
                }
            }
        }
        for (auto model_part_name : model_part_names) {
            if (model_part_name == RootName) continue;
            if (model_part_name.starts_with(Prefix)) model_part_name.erase(0, Prefix.size());
            result.insert(model_part_name);
        }
    }
    return result;
};
} // namespace

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

void ProcessUtilities::AddProcessesSubModelPartListToSolverSettings(const Parameters& rProjectParameters,
                                                                    Parameters& rSolverSettings)
{
    std::set<std::string> domain_condition_names;
    const auto            root_name = rSolverSettings["model_part_name"].GetString();
    const auto            prefix    = root_name + ".";

    if (rProjectParameters.Has("processes")) {
        for (const auto& r_process : rProjectParameters["processes"]) {
            const auto modelpart_names = ExtractModelPartNames(r_process, root_name, prefix);
            domain_condition_names.insert(modelpart_names.begin(), modelpart_names.end());
        }
    }
    if (rSolverSettings.Has("processes_sub_model_part_list")) {
        KRATOS_WARNING("ProcessUtilities")
            << "'processes_sub_model_part_list' is deprecated. This list is built automatically "
               "from the model parts used in all processes."
            << std::endl;
        rSolverSettings.RemoveValue("processes_sub_model_part_list");
    }
    rSolverSettings.AddEmptyArray("processes_sub_model_part_list");

    for (const auto& r_name : domain_condition_names) {
        rSolverSettings["processes_sub_model_part_list"].Append(r_name);
    }
}

}; // namespace Kratos
