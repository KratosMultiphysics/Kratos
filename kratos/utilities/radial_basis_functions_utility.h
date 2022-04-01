//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Sebastian Ares de Parga
//

#if !defined(KRATOS_RADIAL_BASIS_FUNCTIONS_UTILITY_H_INCLUDED)
#define  KRATOS_RADIAL_BASIS_FUNCTIONS_UTILITY_H_INCLUDED


// System includes

// External includes

// Project includes
#include "containers/array_1d.h"
#include "includes/define.h"

namespace Kratos
{

/**
 * @brief Radial Basis Functions utility to calculate shape function values
 * This class uses Gaussian raidal basis functions to calculate
 * the shape function values for a given value (i.e. norm of a point) with partition of unity.
 */
class KRATOS_API(KRATOS_CORE) RadialBasisFunctionsUtility
{

public:

    /**
     * @brief Calculate the RBF value
     * This function evaluates the Gaussian RBF for a norm
     * @param x Norm of RBF argument (i.e. norm of radial vector)
     * @param h Gaussian radius
     * @return double The RBF value
     */
    static double EvaluateRBF(
        const double x,
        const double h);

    /**
     * @brief Calculates the RBF shape function values
     * This method calculates the RBF shape function values in one point using as
     * support the given cloud of points.
     * @param rPoints Matrix containing the coordinates of the support cloud of points
     * @param rX Coordinates where the shape functions are to be computed
     * @param h RBF shape parameter
     * @param rN Shape functions container
     * @param Y Function to interpolate (if values of RHS are known apriori)
     */
    static void CalculateShapeFunctions(
        const Matrix& rPoints,
        const array_1d<double,3>& rX,
        const double h,
        Vector& rN);

    static void CalculateShapeFunctions(
        const Matrix& rPoints,
        const array_1d<double,3>& rX,
        Vector& rN);
    
    static double CalculateShapeFunctionsAndInterpolation(
        const Matrix& rPoints,
        const array_1d<double,3>& rX,
        const double h,
        Vector& rN,
        Vector& Y);
};

}  // namespace Kratos.

#endif // KRATOS_RADIAL_BASIS_FUNCTIONS_UTILITY_H_INCLUDED  defined