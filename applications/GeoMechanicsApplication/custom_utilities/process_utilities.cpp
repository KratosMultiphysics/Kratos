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

void ExtractModelPartNames(const auto& process_list,
                           std::unordered_set<std::string, StringHash, std::equal_to<>>& domain_condition_names,
                           std::string_view root_name,
                           std::string_view prefix)
{
    for (const auto& process : process_list) {
        if (process.Has("Parameters") && process["Parameters"].Has("model_part_name")) {
            auto model_part_name = process["Parameters"]["model_part_name"].GetString();
            if (model_part_name == root_name) continue;
            if (model_part_name.rfind(prefix, 0) == 0) model_part_name.erase(0, prefix.size());
            domain_condition_names.insert(model_part_name);
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
        for (const auto& process_list_name : process_lists_to_be_checked) {
            if (processes.Has(process_list_name))
                ExtractModelPartNames(processes[process_list_name], domain_condition_names, root_name, prefix);
        }
    }
    if (rSolverSettings.Has("processes_sub_model_part_list")) {
        KRATOS_INFO("KratosGeoSettlement") << "processes_sub_model_part_list is deprecated" << std::endl;
        rSolverSettings.RemoveValue("processes_sub_model_part_list");
    }
    rSolverSettings.AddEmptyArray("processes_sub_model_part_list");

    for (const auto& name : domain_condition_names) {
        rSolverSettings["processes_sub_model_part_list"].Append(Kratos::Parameters("\"" + name + "\""));
    }
}

}; // namespace Kratos
