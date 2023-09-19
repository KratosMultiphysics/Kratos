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

#include "includes/kratos_parameters.h"

#include <string>
#include <vector>

namespace Kratos {

class ProcessNameParser {
public:
    std::vector<std::string> GetProcessNames(const Parameters& processParameters);


private:
    std::vector<std::string> mProcessNames;


    void AddConstraintsProcesses(const Parameters& processParameters);
};

}