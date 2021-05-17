//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Ruben Zorrilla
//                   Zhiming Guo
//

#if !defined(KRATOS_MLS_SHAPE_FUNCTIONS_UTILITY_H_INCLUDED)
#define  KRATOS_MLS_SHAPE_FUNCTIONS_UTILITY_H_INCLUDED


// System includes

// External includes

// Project includes
#include "containers/array_1d.h"
#include "includes/define.h"

namespace Kratos
{

/**
 * @brief Moving Least-Squares utility to calculate shape function values
 * This class uses a linear MLS kernel to calculate shape function values for a given point
 */
class MLSShapeFunctionsUtility
{

public:

    static double CalculateKernel(
        const array_1d<double,3>& rRadVect,
        const double h);

    template<std::size_t TDim>
    static void CalculateKernelDerivative(
        const array_1d<double,3>& rRadVect,
        const double h,
        array_1d<double,TDim>& rKernelDerivative);

    template<std::size_t TDim>
    static void EvaluateLinearPolynomialBasis(
        const array_1d<double,3>& rX,
        array_1d<double,TDim+1>& rBasis);

    template<std::size_t TDim>
    static void CalculateShapeFunctions(
        const Matrix& rPoints,
        const array_1d<double,3>& rX,
        const double h,
        Vector& rN);

    template<std::size_t TDim>
    static void CalculateShapeFunctionsAndGradients(
        const Matrix& rPoints,
        const array_1d<double,3>& rX,
        const double h,
        Vector& rN,
        Matrix& rDNDX);
};

}  // namespace Kratos.

#endif // KRATOS_MLS_SHAPE_FUNCTIONS_UTILITY_H_INCLUDED  defined
