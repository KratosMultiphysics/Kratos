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

#pragma once

#include "dgeosettlement.h"


namespace Kratos
{
    KratosGeoSettlement::KratosGeoSettlement()
        : KratosGeoApplication()
    {
    }

    int execute_application(std::string workingDirectory, std::string parameterName,
        std::function<void(char*)> logCallback,
        std::function<void(double)> reportProgress,
        std::function<void(char*)> reportTextualProgress,
        std::function<bool()> shouldCancel)
    {
    	execute_kratos_calculation(model_part, processes, p_solving_strategy);
        return 0;
    }

}
