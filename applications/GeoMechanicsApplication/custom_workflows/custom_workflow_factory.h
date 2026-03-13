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

#include "includes/kratos_export_api.h"

namespace Kratos
{

class KratosGeoSettlement;
class KratosExecute;

class KRATOS_API(GEO_MECHANICS_APPLICATION) CustomWorkflowFactory
{
public:
    static KratosGeoSettlement* CreateKratosGeoSettlement();
    static KratosExecute*       CreateKratosGeoFlow();
};

} // namespace Kratos
