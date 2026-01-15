// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Mohamed Nabi
//                   Anne van de Graaf
//

#pragma once

#include "geo_aliases.h"
#include "includes/kratos_parameters.h"

namespace Kratos
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) FunctionObjectUtilities
{
public:
    static Geo::KappaDependentFunction MakeConstantFunction(double Value);
    static Geo::KappaDependentFunction MakeLinearFunction(double Value, double Coefficient);
};

} // namespace Kratos
