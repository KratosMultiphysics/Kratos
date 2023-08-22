// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Anne van de Graaf
//

#pragma once

#include <string>
#include <functional>

#include "includes/kernel.h"
#include "includes/kratos_export_api.h"

#include "geo_mechanics_application.h"


namespace Kratos
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) KratosGeoSettlement
{
public:
    KratosGeoSettlement();

    int RunStage(const std::string&                      rWorkingDirectory,
                 const std::string&                      rProjectParametersFileName,
                 const std::function<void(const char*)>& rLogCallback,
                 const std::function<void(double)>&      rReportProgress,
                 const std::function<void(const char*)>& rReportTextualProgress,
                 const std::function<bool()>&            rShouldCancel);

private:
    static void AddNodalSolutionStepVariablesTo(ModelPart& rModelPart);
    static void AddDegreesOfFreedomTo(ModelPart& rModelPart);

    Kernel mKernel;
    Model mModel;
    KratosGeoMechanicsApplication::Pointer mpGeoApp;
};

}
