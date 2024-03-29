// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Gennady Markelov
//

#pragma once

// Project includes
#include <optional>

// Application includes
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

class GeoTransportEquationUtilities
{
public:
    template <unsigned int TDim, unsigned int TNumNodes>
    static inline BoundedMatrix<double, TNumNodes, TNumNodes> CalculatePermeabilityMatrix(
        const Matrix&                            rGradNpT,
        double                                   DynamicViscosityInverse,
        const BoundedMatrix<double, TDim, TDim>& rMaterialPermeabilityMatrix,
        double                                   RelativePermeability,
        double                                   PermeabilityUpdateFactor,
        double                                   IntegrationCoefficient)
    {
        return CalculatePermeabilityMatrix(rGradNpT, DynamicViscosityInverse,
                                           rMaterialPermeabilityMatrix, RelativePermeability,
                                           PermeabilityUpdateFactor, IntegrationCoefficient);
    }

    static inline Matrix CalculatePermeabilityMatrix(const Matrix& rGradNpT,
                                                     double        DynamicViscosityInverse,
                                                     const Matrix& rMaterialPermeabilityMatrix,
                                                     double        RelativePermeability,
                                                     double        PermeabilityUpdateFactor,
                                                     double        IntegrationCoefficient)
    {
        return -PORE_PRESSURE_SIGN_FACTOR * DynamicViscosityInverse *
               prod(rGradNpT, Matrix(prod(rMaterialPermeabilityMatrix, trans(rGradNpT)))) *
               RelativePermeability * PermeabilityUpdateFactor * IntegrationCoefficient;
    }

    template <unsigned int TNumNodes>
    static inline BoundedMatrix<double, TNumNodes, TNumNodes> CalculateCompressibilityMatrix(
        const Vector&         rNp,
        double                BiotModulusInverse,
        double                IntegrationCoefficient,
        std::optional<double> DtPressureCoefficient = std::nullopt)
    {
        return CalculateCompressibilityMatrix(rNp, BiotModulusInverse, IntegrationCoefficient, DtPressureCoefficient);
    }

    static inline Matrix CalculateCompressibilityMatrix(const Vector& rNp,
                                                        double        BiotModulusInverse,
                                                        double        IntegrationCoefficient,
                                                        std::optional<double> DtPressureCoefficient = std::nullopt)
    {
        auto result = -PORE_PRESSURE_SIGN_FACTOR * BiotModulusInverse * outer_prod(rNp, rNp) * IntegrationCoefficient;
        if (DtPressureCoefficient.has_value()) return result * DtPressureCoefficient.value();
        else return result;
    }
}; /* Class GeoTransportEquationUtilities*/
} /* namespace Kratos.*/
