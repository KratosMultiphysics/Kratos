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

// System includes

/* External includes */

#include <geo_mechanics_application.h>
#include "dgeoapplication.h"

namespace Kratos
{
    class KRATOS_API(GEO_MECHANICS_APPLICATION) KratosGeoSettlement : public KratosGeoApplication
    {
    public:
        KratosGeoSettlement();
        ~KratosGeoSettlement(){}

        int execute_application(std::string workingDirectory, std::string parameterName,
            std::function<void(char*)> logCallback,
            std::function<void(double)> reportProgress,
            std::function<void(char*)> reportTextualProgress,
            std::function<bool()> shouldCancel);
    };
}