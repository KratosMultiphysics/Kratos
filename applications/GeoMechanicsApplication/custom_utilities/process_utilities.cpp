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
#include "custom_utilities/string_utilities.h"
#include "includes/kratos_parameters.h"

using namespace std::string_literals;

namespace
{

std::vector<std::string> GetProcessModelPartNames(const Kratos::Parameters& rProcessSettings,
                                                  const std::string&        rProcessInfo,
                                                  const std::vector<std::string>& rModelPartNameKeys)
{
    auto has_name_key = [&rProcessSettings](const std::string& rKey) {
        return rProcessSettings.Has(rKey);
    };
    KRATOS_ERROR_IF_NOT(std::ranges::any_of(rModelPartNameKeys, has_name_key))
        << "Please specify any of " << Kratos::GeoStringUtilities::Join(rModelPartNameKeys, ", ")
        << " for " << rProcessInfo;

    KRATOS_ERROR_IF(rModelPartNameKeys.size() > 1 && std::ranges::all_of(rModelPartNameKeys, has_name_key))
        << "The parameters " << Kratos::GeoStringUtilities::Join(rModelPartNameKeys, ", ")
        << " are mutually exclusive for " << rProcessInfo;

    const auto& r_name_key = *std::find_if(rModelPartNameKeys.begin(), rModelPartNameKeys.end(), has_name_key);
    const auto name_or_names = rProcessSettings[r_name_key];
    return name_or_names.IsStringArray() ? name_or_names.GetStringArray()
                                         : std::vector{name_or_names.GetString()};
}

std::set<std::string, std::less<>> ExtractModelPartNames(const auto&      rProcessList,
                                                         std::string_view RootName,
                                                         std::string_view Prefix)
{
    const std::set master_slave_process_names = {"AssignAverageMasterSlaveConstraintsProcess"s,
                                                 "ApplyPeriodicConditionProcess"s, "SkinDetectionProcess"s};

    std::set<std::string, std::less<>> result;
    for (const auto& r_process : rProcessList) {
        if (!r_process.Has("Parameters")) continue;

        const auto model_part_name_keys =
            (r_process.Has("process_name") &&
             master_slave_process_names.contains(r_process["process_name"].GetString()))
                ? std::vector{"computing_model_part_name"s}
                : std::vector{"model_part_name"s, "model_part_name_list"s};
        const auto model_part_names =
            GetProcessModelPartNames(r_process["Parameters"], {}, model_part_name_keys);
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
    Model&                          rModel,
    const Parameters&               rProcessSettings,
    const std::string&              rProcessInfo,
    const std::vector<std::string>& rModelPartNameKeys)
{
    const auto model_part_names = GetProcessModelPartNames(rProcessSettings, rProcessInfo, rModelPartNameKeys);
    KRATOS_ERROR_IF(model_part_names.empty()) << "The parameters 'model_part_name_list' needs "
                                                 "to contain at least one model part name for "
                                              << rProcessInfo << ".";

    std::vector<std::reference_wrapper<ModelPart>> result;
    result.reserve(model_part_names.size());
    std::ranges::transform(
        model_part_names, std::back_inserter(result),
        [&rModel](const auto& rName) -> ModelPart& { return rModel.GetModelPart(rName); });

    const std::set<std::string, std::less<>> unique_names(model_part_names.begin(), model_part_names.end());
    KRATOS_ERROR_IF_NOT(unique_names.size() == model_part_names.size())
        << "model_part_name_list has duplicated names for " << rProcessInfo << "." << std::endl;

    return result;
}

void ProcessUtilities::AddProcessesSubModelPartListToSolverSettings(const Parameters& rProjectParameters,
                                                                    Parameters& rSolverSettings)
{
    std::set<std::string, std::less<>> domain_condition_names;
    const auto                         root_name = rSolverSettings["model_part_name"].GetString();
    const auto                         prefix    = root_name + ".";

    if (rProjectParameters.Has("processes")) {
        for (const auto& r_process_list : rProjectParameters["processes"]) {
            const auto modelpart_names = ExtractModelPartNames(r_process_list, root_name, prefix);
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
