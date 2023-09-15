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

#include <string>
#include <includes/kratos_parameters.h>

namespace Kratos {

class InterfaceInputUtility {
public:
    virtual Parameters ProjectParametersFrom(const std::string &rProjectFilePath) const = 0;
};

}
