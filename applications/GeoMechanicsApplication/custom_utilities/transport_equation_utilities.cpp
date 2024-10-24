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

#include "transport_equation_utilities.hpp"

namespace Kratos
{

double GeoTransportEquationUtilities::CalculateParticleDiameter(const Properties& rProperties)
{
    double diameter;

    if (rProperties[PIPE_MODIFIED_D])
        diameter = 2.08e-4 * pow((rProperties[PIPE_D_70] / 2.08e-4), 0.4);
    else diameter = rProperties[PIPE_D_70];
    return diameter;
}

} // namespace Kratos
