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

namespace Kratos {

class ProcessInfoStubParser : public ProcessInfoParser {

public:
    std::vector<ProcessParameters> GetProcessList(const Parameters& rProcessParameters) override;
};

} // Kratos