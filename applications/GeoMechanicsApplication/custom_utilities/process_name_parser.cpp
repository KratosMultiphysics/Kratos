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
#include "process_name_parser.h"

namespace Kratos
{

std::vector<std::string> ProcessNameParser::GetProcessNames(const Kratos::Parameters& processParameters)
{
    AddConstraintsProcesses(processParameters);

    return mProcessNames;
}

void ProcessNameParser::AddConstraintsProcesses(const Parameters& processParameters) {
    const auto constraints_processes = processParameters["constraints_process_list"];

    for (Parameters process : constraints_processes)
    {
        const std::string process_name_entry = "process_name";
        if (process.Has(process_name_entry))
        {
            mProcessNames.push_back(process[process_name_entry].GetString());
        }
    }
}

}

