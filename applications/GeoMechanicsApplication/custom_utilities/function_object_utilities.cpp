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

#include "function_object_utilities.h"

namespace Kratos
{

Geo::KappaDependentFunction FunctionObjectUtilities::MakeConstantFunction(double Value)
{
    return [Value](double /* unused kappa */) { return Value; };
}

Geo::KappaDependentFunction FunctionObjectUtilities::MakeLinearFunction(double Value, double Coefficient)
{
    return [Value, Coefficient](double Kappa) { return Value + Coefficient * Kappa; };
}

} // namespace Kratos