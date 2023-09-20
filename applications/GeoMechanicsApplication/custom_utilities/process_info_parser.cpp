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
#include "process_info_parser.h"

namespace Kratos
{

ProcessInfoParser::ProcessInfoParser(const Parameters& rProcessParameters) :
    mProcessParameters{rProcessParameters}
{
}

std::vector<ProcessInfo> ProcessInfoParser::GetProcessList()
{
    const std::vector<std::string> process_list_names = {"constraints_process_list",
                                                         "loads_process_list",
                                                         "auxiliar_process_list"};

    for (const auto& process_list_name : process_list_names)
    {
        AddProcessesForList(process_list_name);
    }

    return mProcessNames;
}

void ProcessInfoParser::AddProcessesForList(const std::string& rProcessListName) {
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
            mProcessNames.push_back({process["Parameters"], process[process_name_entry].GetString()});
        }
    }
}

}

