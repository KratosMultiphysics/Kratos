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
#include "includes/node.h"
#include "includes/variables.h"
#include "geometries/point_3d.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"

#define PI 3.14159265

namespace Kratos
{

/**
 * @brief Moving Least-Squares utility to calculate shape function values
 * This class uses a linear MLS kernel to calculate shape function values for a given point
 */
class MLSShapeFunctionsUtility
{

public:

    static void CalculateKernel(
        const array_1d<double,3>& rRadVect,
        const double h,
        double& rKernel);

    template<std::size_t TDim>
    static void CalculateKernelDerivative(
        const array_1d<double,TDim>& rRadVect,
        const double h,
        array_1d<double,TDim>& rKernelDerivative);

    ///A is the area we associate to the point
    ///N is a vector such that N(j) is the shape function associated to node j
    ///DN_DX(i,k) is the derivative of the shape function of node i, component k
    ///coordinates, is a input matrix with the coordinates of all of the points in the cloud, Coordinates(i,k) is the component k of the i-th node in the cloud
    ///nn, is number of gauss points in the cloud
    // only 2d now

    static void ComputeMLSKernel(
        Vector& Ng,
        Matrix& DN_DX,
        const Matrix& Coordinates,
        const array_1d<double,3>& x_size3,
        const double& h);
};

}  // namespace Kratos.

#endif // KRATOS_MLS_SHAPE_FUNCTIONS_UTILITY_H_INCLUDED  defined
