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
//                   Gennady Markelov
//

#pragma once

// Project includes

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
    static inline void CalculatePermeabilityMatrixH(BoundedMatrix<double, TNumNodes, TNumNodes>& rPMatrix,
                                                    const Matrix& GradNpT,
                                                    double        DynamicViscosityInverse,
                                                    const BoundedMatrix<double, TDim, TDim>& PermeabilityMatrix,
                                                    double IntegrationCoefficient)
    {
        noalias(rPMatrix) = DynamicViscosityInverse *
                            prod(Matrix(prod(GradNpT, PermeabilityMatrix)), trans(GradNpT)) * IntegrationCoefficient;
    }

    static inline Matrix CalculatePermeabilityMatrixH(const Matrix& GradNpT,
                                                      double        DynamicViscosityInverse,
                                                      const Matrix& PermeabilityMatrix,
                                                      double        IntegrationCoefficient)
    {
        return DynamicViscosityInverse *
               prod(GradNpT, Matrix(prod(PermeabilityMatrix, trans(GradNpT)))) * IntegrationCoefficient;
    }

    template <unsigned int TNumNodes>
    static inline void PreparePermeabilityMatrixHForIntegration(BoundedMatrix<double, TNumNodes, TNumNodes>& rPMatrix,
                                                                double RelativePermeability,
                                                                double PermeabilityUpdateFactor)
    {
        rPMatrix *= (-PORE_PRESSURE_SIGN_FACTOR * RelativePermeability * PermeabilityUpdateFactor);
    }

    static inline void PreparePermeabilityMatrixHForIntegration(Matrix& rPMatrix,
                                                                double  RelativePermeability,
                                                                double  PermeabilityUpdateFactor)
    {
        rPMatrix *= (-PORE_PRESSURE_SIGN_FACTOR * RelativePermeability * PermeabilityUpdateFactor);
    }

}; /* Class GeoTransportEquationUtilities*/
} /* namespace Kratos.*/
