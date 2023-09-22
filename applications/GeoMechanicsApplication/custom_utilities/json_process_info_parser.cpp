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

std::vector<ProcessParameters> JsonProcessInfoParser::GetProcessList(const Parameters& rProcessParameters)
{
    mProcessParameters = rProcessParameters;
    const std::vector<std::string> process_list_names = {"constraints_process_list",
                                                         "loads_process_list",
                                                         "auxiliar_process_list"};

    for (const auto& process_list_name : process_list_names)
    {
        JsonProcessInfoParser::AddProcessesForList(process_list_name);
    }

    return mProcessNames;
}

void JsonProcessInfoParser::AddProcessesForList(const std::string& rProcessListName) {
    if (!mProcessParameters.Has(rProcessListName))
    {
        return;
    }

    const auto processes_in_list = mProcessParameters[rProcessListName];

    for (Parameters process : processes_in_list)
    {
        const std::string process_name_entry = "process_name";
        if (process.Has(process_name_entry))
        {
            mProcessNames.emplace_back(process["Parameters"], process[process_name_entry].GetString());
        }
    }
}

}

