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
struct StringHash {
    using is_transparent = void; // Enables heterogeneous operations.

    std::size_t operator()(std::string_view sv) const
    {
        std::hash<std::string_view> hasher;
        return hasher(sv);
    }
};

void ExtractModelPartNames(const auto& rProcessList,
                           std::unordered_set<std::string, StringHash, std::equal_to<>>& rDomainConditionNames,
                           std::string_view RootName,
                           std::string_view Prefix)
{
    for (const auto& r_process : rProcessList) {
        if (r_process.Has("Parameters") && r_process["Parameters"].Has("model_part_name")) {
            auto model_part_name = r_process["Parameters"]["model_part_name"].GetString();
            if (model_part_name == RootName) continue;
            if (model_part_name.rfind(Prefix, 0) == 0) model_part_name.erase(0, Prefix.size());
            rDomainConditionNames.insert(model_part_name);
        }
    }
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

void ProcessUtilities::AddProcessesSubModelPartList(const Parameters& rProjectParameters, Parameters& rSolverSettings)
{
    std::unordered_set<std::string, StringHash, std::equal_to<>> domain_condition_names;
    const auto root_name = rSolverSettings["model_part_name"].GetString();
    const auto prefix    = root_name + ".";

    const std::vector<std::string> process_lists_to_be_checked = {
        "constraints_process_list", "loads_process_list", "auxiliary_process_list"};
    if (rProjectParameters.Has("processes")) {
        const auto processes = rProjectParameters["processes"];
        for (const auto& r_process_list_name : process_lists_to_be_checked) {
            if (processes.Has(r_process_list_name))
                ExtractModelPartNames(processes[r_process_list_name], domain_condition_names, root_name, prefix);
        }
    }
    if (rSolverSettings.Has("processes_sub_model_part_list")) {
        KRATOS_INFO("ProcessUtilities") << "processes_sub_model_part_list is deprecated" << std::endl;
        rSolverSettings.RemoveValue("processes_sub_model_part_list");
    }
    rSolverSettings.AddEmptyArray("processes_sub_model_part_list");

    for (const auto& r_name : domain_condition_names) {
        rSolverSettings["processes_sub_model_part_list"].Append(Kratos::Parameters("\"" + r_name + "\""));
    }
}

}; // namespace Kratos
