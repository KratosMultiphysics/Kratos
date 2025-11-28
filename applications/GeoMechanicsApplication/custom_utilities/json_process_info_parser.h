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

#include "process_info_parser.h"
#include "process_parameters.h"

#include <vector>

namespace Kratos
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) JsonProcessInfoParser : public ProcessInfoParser
{
public:
    std::vector<ProcessParameters> GetProcessList(const Parameters& rProcessParameters) const override;

private:
    static std::vector<ProcessParameters> AddProcessesForList(const std::string& rProcessListName,
                                                              const Parameters& rProcessParameters);
};

} // namespace Kratos