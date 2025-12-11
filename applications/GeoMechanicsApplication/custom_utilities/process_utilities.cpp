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

std::vector<std::string> GetProcessModelPartNames(const Kratos::Parameters& rProcessSettings,
                                                  const std::string&        rProcessInfo)
{
    KRATOS_ERROR_IF_NOT(rProcessSettings.Has("model_part_name") ||
                        rProcessSettings.Has("model_part_name_list"))
        << "Please specify 'model_part_name' or 'model_part_name_list' for " << rProcessInfo;

    KRATOS_ERROR_IF(rProcessSettings.Has("model_part_name") &&
                    rProcessSettings.Has("model_part_name_list"))
        << "The parameters 'model_part_name' and 'model_part_name_list' are mutually exclusive for "
        << rProcessInfo;

    return rProcessSettings.Has("model_part_name_list")
               ? rProcessSettings["model_part_name_list"].GetStringArray()
               : std::vector{rProcessSettings["model_part_name"].GetString()};
}

std::set<std::string, std::less<>> ExtractModelPartNames(const auto&      rProcessList,
                                                         std::string_view RootName,
                                                         std::string_view Prefix)
{
    std::set<std::string, std::less<>> result;
    for (const auto& r_process : rProcessList) {
        if (!r_process.Has("Parameters")) continue;
        const auto model_part_names = GetProcessModelPartNames(r_process["Parameters"], {});
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
    const auto model_part_names = GetProcessModelPartNames(rProcessSettings, rProcessInfo);
    KRATOS_ERROR_IF(model_part_names.empty()) << "The parameters 'model_part_name_list' needs "
                                                 "to contain at least one model part name for "
                                              << rProcessInfo << ".";

    std::vector<std::reference_wrapper<ModelPart>> result;
    result.reserve(model_part_names.size());
    std::ranges::transform(
        model_part_names, std::back_inserter(result),
        [&rModel](const auto& rName) -> ModelPart& { return rModel.GetModelPart(rName); });

    // check for duplicated names
    std::set<std::string> unique_names;
    for (auto& mp_ref : result) {
        unique_names.insert(mp_ref.get().Name());
    }
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
