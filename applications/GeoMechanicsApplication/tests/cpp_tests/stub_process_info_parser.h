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

#include "custom_utilities/process_info_parser.h"

namespace Kratos
{

class StubProcessInfoParser : public ProcessInfoParser
{
public:
    std::vector<ProcessParameters> GetProcessList(const Parameters& rProcessParameters) const override;

    static std::size_t NumberOfProcesses();
};

} // namespace Kratos