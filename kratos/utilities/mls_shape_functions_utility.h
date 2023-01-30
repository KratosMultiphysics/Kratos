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

#pragma once

// System includes

// External includes

// Project includes
#include "containers/array_1d.h"
#include "includes/define.h"
#include "includes/ublas_interface.h"

namespace Kratos
{

/**
 * @brief Moving Least-Squares utility to calculate shape function values
 * This class uses a linear polynomial basis and an exponential kernel to calculate
 * the shape function values and gradients in a given point using a Moving Least-Squares minimization
 */
class MLSShapeFunctionsUtility
{

public:

    /**
     * @brief Calculate the kernel value
     * This function calculates the exponential kernel value in one point
     * @param rRadVect Vector from the current point to the point of interest
     * @param h Kernel radius
     * @return double The kernel value
     */
    static double CalculateKernel(
        const array_1d<double,3>& rRadVect,
        const double h);

    /**
     * @brief Calculate the kernel derivative values
     * This function calculates the exponential kernel derivatives in one point
     * @tparam TDim Dimension (2 in 2D and 3 in 3D)
     * @param rRadVect Vector from the current point to the point of interest
     * @param h Kernel radius
     * @param rKernelDerivative The kernel derivative value
     */
    template<std::size_t TDim>
    static void CalculateKernelDerivative(
        const array_1d<double,3>& rRadVect,
        const double h,
        array_1d<double,TDim>& rKernelDerivative);

    /**
     * @brief Evaluates the linear polynomial basis
     * This method evaluates the 2D linear polynommial basis in one point
     * @param rX Coordinates where the basis is evaluated
     * @param rBasis 2D linear polynomial basis values
     */
    static void EvaluatePolynomialBasis(
        const array_1d<double,3>& rX,
        array_1d<double, 3>& rBasis);

    /**
     * @brief Evaluates the linear polynomial basis
     * This method evaluates the 3D linear polynommial basis in one point
     * @param rX Coordinates where the basis is evaluated
     * @param rBasis 3D linear polynomial basis values
     */
    static void EvaluatePolynomialBasis(
        const array_1d<double,3>& rX,
        array_1d<double, 4>& rBasis);

    /**
     * @brief Evaluates the quadratic polynomial basis
     * This method evaluates the 2D quadratic polynommial basis in one point
     * @param rX Coordinates where the basis is evaluated
     * @param rBasis 2D linear polynomial basis values
     */
    static void EvaluatePolynomialBasis(
        const array_1d<double,3>& rX,
        array_1d<double, 6>& rBasis);

    /**
     * @brief Evaluates the quadratic polynomial basis
     * This method evaluates the 3D quadratic polynommial basis in one point
     * @param rX Coordinates where the basis is evaluated
     * @param rBasis 3D quadratic polynomial basis values
     */
    static void EvaluatePolynomialBasis(
        const array_1d<double,3>& rX,
        array_1d<double, 10>& rBasis);

    /**
     * @brief Evaluates the linear polynomial basis derivatives
     * This method evaluates the 2D linear polynommial basis derivatives in one point
     * @param rX Coordinates where the basis is evaluated
     * @param rBasis 2D linear polynomial basis derivatives values
     */
    static void EvaluatePolynomialBasisDerivatives(
        const array_1d<double,3>& rX,
        BoundedMatrix<double, 2, 3>& rBasisDerivatives);

    /**
     * @brief Evaluates the linear polynomial basis derivatives
     * This method evaluates the 3D linear polynommial basis derivatives in one point
     * @param rX Coordinates where the basis is evaluated
     * @param rBasis 3D linear polynomial basis derivatives values
     */
    static void EvaluatePolynomialBasisDerivatives(
        const array_1d<double,3>& rX,
        BoundedMatrix<double, 3, 4>& rBasisDerivatives);

    /**
     * @brief Evaluates the quadratic polynomial basis derivatives
     * This method evaluates the 2D quadratic polynommial basis derivatives in one point
     * @param rX Coordinates where the basis is evaluated
     * @param rBasis 2D quadratic polynomial basis derivatives values
     */
    static void EvaluatePolynomialBasisDerivatives(
        const array_1d<double,3>& rX,
        BoundedMatrix<double, 2, 6>& rBasisDerivatives);

    /**
     * @brief Evaluates the quadratic polynomial basis derivatives
     * This method evaluates the 3D quadratic polynommial basis derivatives in one point
     * @param rX Coordinates where the basis is evaluated
     * @param rBasis 3D quadratic polynomial basis derivatives values
     */
    static void EvaluatePolynomialBasisDerivatives(
        const array_1d<double,3>& rX,
        BoundedMatrix<double, 3, 10>& rBasisDerivatives);

    /**
     * @brief Calculates the MLS shape function values
     * This method calculates the shape function values in one point using as
     * support the given cloud of points.
     * @tparam TDim Dimension (2 in 2D and 3 in 3D)
     * @param rPoints Matrix containing the coordinates of the support cloud of points
     * @param rX Coordinates where the shape functions are to be computed
     * @param h Kernel radius
     * @param rN Shape functions container
     */
    template<std::size_t TDim, std::size_t TOrder>
    static void CalculateShapeFunctions(
        const Matrix& rPoints,
        const array_1d<double,3>& rX,
        const double h,
        Vector& rN);

    /**
     * @brief Calculates the MLS shape function values
     * This method calculates the shape function values and their gradients in one point using as
     * support the given cloud of points.
     * @tparam TDim Dimension (2 in 2D and 3 in 3D)
     * @param rPoints Matrix containing the coordinates of the support cloud of points
     * @param rX Coordinates where the shape functions are to be computed
     * @param h Kernel radius
     * @param rN Shape functions container
     * @param rDNDX Shape functions gradients container
     */
    template<std::size_t TDim, std::size_t TOrder>
    static void CalculateShapeFunctionsAndGradients(
        const Matrix& rPoints,
        const array_1d<double,3>& rX,
        const double h,
        Vector& rN,
        Matrix& rDNDX);
};

}  // namespace Kratos.
