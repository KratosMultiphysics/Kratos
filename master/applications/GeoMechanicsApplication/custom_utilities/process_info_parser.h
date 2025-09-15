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

#include <vector>

namespace Kratos
{

class ProcessInfoParser
{
public:
    virtual ~ProcessInfoParser() = default;
    virtual std::vector<ProcessParameters> GetProcessList(const Parameters& rProcessParameters) const = 0;
};

} // namespace Kratos