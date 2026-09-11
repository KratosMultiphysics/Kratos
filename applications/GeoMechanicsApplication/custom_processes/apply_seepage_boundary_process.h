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

#pragma once

#include "includes/kratos_export_api.h"
#include "processes/process.h"

namespace Kratos
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) ApplySeepageBoundaryProcess : public Process
{
public:
    ApplySeepageBoundaryProcess();
    ~ApplySeepageBoundaryProcess() override;
};

} // namespace Kratos