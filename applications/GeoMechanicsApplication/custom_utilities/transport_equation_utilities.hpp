// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Ignasi de Pouplana,
//                   Vahid Galavi
//

#pragma once

// Project includes
#include "includes/element.h"
#include "utilities/math_utils.h"

// Application includes
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

class GeoTransportEquationUtilities
{
public:
    using VectorType = Vector;
    using MatrixType = Matrix;

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    template <unsigned int TDim, unsigned int TNumNodes>
    static inline void CalculatePermeabilityMatrixH(BoundedMatrix<double, TNumNodes, TDim>& rPDimMatrix,
                                                    BoundedMatrix<double, TNumNodes, TNumNodes>& rPMatrix,
                                                    const Matrix& GradNpT,
                                                    double        DynamicViscosityInverse,
                                                    const BoundedMatrix<double, TDim, TDim>& PermeabilityMatrix,
                                                    double IntegrationCoefficient)
    {
        noalias(rPDimMatrix) = prod(GradNpT, PermeabilityMatrix);
        noalias(rPMatrix) = DynamicViscosityInverse * prod(rPDimMatrix, trans(GradNpT)) * IntegrationCoefficient;
    }

    template <unsigned int TNumNodes>
    static inline void PreparePermeabilityMatrixHForIntegration(BoundedMatrix<double, TNumNodes, TNumNodes>& rPMatrix,
                                                                double RelativePermeability,
                                                                double PermeabilityUpdateFactor)
    {
        rPMatrix *= (-PORE_PRESSURE_SIGN_FACTOR * RelativePermeability * PermeabilityUpdateFactor);
    }

    static inline Matrix CalculatePermeabilityMatrixH(const Matrix& GradNpT,
                                                      double        DynamicViscosityInverse,
                                                      const Matrix& PermeabilityMatrix,
                                                      double        IntegrationCoefficient)
    {
        return DynamicViscosityInverse *
               prod(GradNpT, Matrix(prod(PermeabilityMatrix, trans(GradNpT)))) * IntegrationCoefficient;
    }

    static inline void PreparePermeabilityMatrixHForIntegration(Matrix& rPMatrix,
                                                                double  RelativePermeability,
                                                                double  PermeabilityUpdateFactor)
    {
        rPMatrix *= (-PORE_PRESSURE_SIGN_FACTOR * RelativePermeability * PermeabilityUpdateFactor);
    }

}; /* Class GeoTransportEquationUtilities*/
} /* namespace Kratos.*/
