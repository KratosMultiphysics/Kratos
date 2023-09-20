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

#pragma once

#include "process_info.h"
#include <vector>

namespace Kratos {

class ProcessInfoParser {
public:
    explicit ProcessInfoParser(const Parameters& rProcessParameters);
    std::vector<ProcessInfo> GetProcessList();

private:
    std::vector<ProcessInfo> mProcessNames;
    const Parameters mProcessParameters;

    void AddProcessesForList(const std::string& rProcessListName);
};

}