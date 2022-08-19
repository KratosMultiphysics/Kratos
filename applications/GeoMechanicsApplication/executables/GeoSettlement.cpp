// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Jonathan Nuttall
//

#include "GeoSettlement.h"
#include <iostream>

int main(int argc, char **argv)
{
    std::string workingDirectory = argv[1];
    std::string projectName = argv[2];

    auto execute = Kratos::KratosExecute();
    execute.geosettlement(workingDirectory, projectName);
}
