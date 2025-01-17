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

Matrix GeoTransportEquationUtilities::CalculatePermeabilityMatrix(const Matrix& rGradNpT,
                                                                  double DynamicViscosityInverse,
                                                                  const Matrix& rMaterialPermeabilityMatrix,
                                                                  double RelativePermeability,
                                                                  double IntegrationCoefficient)
{
    return -PORE_PRESSURE_SIGN_FACTOR * DynamicViscosityInverse *
           prod(rGradNpT, Matrix(prod(rMaterialPermeabilityMatrix, trans(rGradNpT)))) *
           RelativePermeability * IntegrationCoefficient;
}

double GeoTransportEquationUtilities::CalculateParticleDiameter(const Properties& rProperties)
{
    return rProperties.Has(PIPE_MODIFIED_D) && rProperties[PIPE_MODIFIED_D]
               ? 2.08e-4 * std::pow((rProperties[PIPE_D_70] / 2.08e-4), 0.4)
               : rProperties[PIPE_D_70];
}

} // namespace Kratos
