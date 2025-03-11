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

#include "includes/kratos_export_api.h"

namespace Kratos
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) RegistrationUtilities
{
public:
    static void RegisterStressStatePolicies();
    static void DeregisterStressStatePolicies();
};

} // namespace Kratos