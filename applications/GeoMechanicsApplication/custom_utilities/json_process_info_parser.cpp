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
#include "json_process_info_parser.h"

namespace Kratos
{

std::vector<ProcessParameters> JsonProcessInfoParser::GetProcessList(const Parameters& rProcessParameters) const
{
    std::vector<ProcessParameters> result;
    const std::vector<std::string> process_list_names = {
        "constraints_process_list", "loads_process_list", "auxiliary_process_list"};

    for (const auto& process_list_name : process_list_names) {
        const auto processes_for_list = AddProcessesForList(process_list_name, rProcessParameters);
        result.insert(result.end(), processes_for_list.begin(), processes_for_list.end());
    }

    return result;
}

std::vector<ProcessParameters> JsonProcessInfoParser::AddProcessesForList(const std::string& rProcessListName,
                                                                          const Parameters& rProcessParameters)
{
    if (!rProcessParameters.Has(rProcessListName)) {
        return {};
    }

    std::vector<ProcessParameters> result;
    for (Parameters process : rProcessParameters[rProcessListName]) {
        const std::string process_name_entry = "process_name";
        if (process.Has(process_name_entry)) {
            result.emplace_back(process[process_name_entry].GetString(), process["Parameters"]);
        }
    }

    return result;
}

} // namespace Kratos
