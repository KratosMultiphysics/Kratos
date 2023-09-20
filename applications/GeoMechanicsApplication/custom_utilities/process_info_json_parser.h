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

#include "process_parameters.h"
#include "process_info_parser.h"
#include <vector>

namespace Kratos {

class ProcessInfoJsonParser : public ProcessInfoParser {
public:
    std::vector<ProcessParameters> GetProcessList(const Parameters& rProcessParameters) override;

private:
    std::vector<ProcessParameters> mProcessNames;
    Parameters mProcessParameters;

    void AddProcessesForList(const std::string& rProcessListName);
};

}