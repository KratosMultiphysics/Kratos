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

// System includes

// External includes

// Project includes
#include "includes/global_variables.h"
#include "utilities/math_utils.h"
#include "rbf_shape_functions_utility.h"

namespace Kratos
{

    double RBFShapeFunctionsUtility::EvaluateRBF(
        const double x,
        const double h)
    {
        // Evaluate Inverse multiquadric
        // const double q = x*h;
        // return 1/std::sqrt(1+std::pow(q,2));
        
        // // Evaluate linear
        return x;

        // Evaluate cubic
        // return std::pow(x, 3);

         // Evaluate multiquadric
        // return std::sqrt(std::pow(x,2) + std::pow(h,2));

        // // Evaluate Gaussian function
        // const double q = x/2.0;
        // return std::exp(-0.5*std::pow(q,2));
    }

    void RBFShapeFunctionsUtility::CalculateShapeFunctions(
        const Matrix& rPoints,
        const array_1d<double,3>& rX,
        const double h,
        Vector& rN,
        DenseQRPointerType pDenseQR)
    {
        KRATOS_TRY;

        KRATOS_ERROR_IF(h < 1.0e-12) << "Reference distance close to zero." << std::endl;

        // Set RBF shape functions containers
        const std::size_t n_points = rPoints.size1();
        if (rN.size() != n_points) {
            rN.resize(n_points, false);
        }

        // Initialize the RBF interpolation matrix and RBF interpolated vector
        Matrix A = ZeroMatrix(n_points,n_points);
        Vector Phi = ZeroVector(n_points);

        // Build the RBF interpolation matrix and RBF interpolated vector
        for (std::size_t i_pt = 0; i_pt < n_points; ++i_pt) {
            for (std::size_t j_pt = 0; j_pt < n_points; ++j_pt) {
                const double norm_xij = norm_2(row(rPoints,i_pt)-row(rPoints,j_pt));
                A(i_pt,j_pt) = EvaluateRBF(norm_xij,h);
            }
            const double norm_X =  norm_2(rX-row(rPoints,i_pt));
            Phi[i_pt] = EvaluateRBF(norm_X,h);
        }

        // Obtain the RBF shape functions (N)
        // Note that we do a QR solve if a QR decomposition pointer is provided
        // Otherwise we solve the system with standard LU factorization
        if (pDenseQR) {
            pDenseQR->Compute(A);
            pDenseQR->Solve(Phi, rN);
        } else {
            MathUtils<double>::Solve(A, rN, Phi);
        }

        // Partition of unity
        noalias(rN) = rN/sum(rN);

        KRATOS_CATCH("");
    }

    void RBFShapeFunctionsUtility::CalculateShapeFunctions(
        const Matrix& rPoints,
        const array_1d<double,3>& rX,
        Vector& rN,
        DenseQRPointerType pDenseQR)
    {
        KRATOS_TRY;

        const double h = CalculateInverseMultiquadricShapeParameter(rPoints);
        CalculateShapeFunctions(rPoints, rX, h, rN, pDenseQR);

        KRATOS_CATCH("");
    }

    double RBFShapeFunctionsUtility::CalculateShapeFunctionsAndInterpolation(
        const Matrix& rPoints,
        const array_1d<double,3>& rX,
        const double h,
        Vector& rN,
        Vector& rY)
    {
        double interpolation = 0;
        KRATOS_TRY;

        KRATOS_ERROR_IF(h < 1.0e-12) << "Reference distance close to zero." << std::endl;

        // Set RBF shape functions containers
        const std::size_t n_points = rPoints.size1();
        if (rN.size() != n_points) {
            rN.resize(n_points, false);
        }

        // Initialize the RBF interpolation matrix and RBF interpolated vector
        Matrix A = ZeroMatrix(n_points,n_points);
        double norm_xij;

        // Build the RBF interpolation matrix and RBF interpolated vector
        for (std::size_t i_pt = 0; i_pt < n_points; ++i_pt) {
            for (std::size_t j_pt = 0; j_pt < n_points; ++j_pt) {
                norm_xij = norm_2(row(rPoints,i_pt)-row(rPoints,j_pt));
                A(i_pt,j_pt) = EvaluateRBF(norm_xij,h);
            }
        }

        // Obtain the RBF shape functions (N)
        DenseHouseholderQRDecomposition<DenseSpace> qr_decomposition;
        qr_decomposition.Compute(A);
        qr_decomposition.Solve(rY, rN);

        // Interpolate solution
        for (std::size_t i_pt = 0; i_pt < n_points; ++i_pt) {
            const double norm_xi = norm_2(rX-row(rPoints,i_pt));
            const double Phi = EvaluateRBF(norm_xi,h);
            interpolation += Phi*rN[i_pt];
        }

        KRATOS_CATCH("");

        return interpolation;
    }

    double RBFShapeFunctionsUtility::CalculateInverseMultiquadricShapeParameter(const Matrix& rPoints)
    {
        // Find shape parameter http://www.math.iit.edu/~fass/Dolomites.pdf [Hardy]
        double d = 0; // Total distance of nearest x_i neighbors
        const std::size_t n_points = rPoints.size1();
        for (std::size_t i_pt = 0; i_pt < n_points; ++i_pt) {
            double d_nearest_neighbor = std::numeric_limits<double>::max(); // Initialize distance to nearest x_i neighbor
            for (std::size_t j_pt = 0; j_pt < n_points; ++j_pt) {
                if (i_pt!=j_pt){ // Avoid measuring distance between same points
                    const double norm_xij = norm_2(row(rPoints,i_pt)-row(rPoints,j_pt));
                    if (norm_xij < d_nearest_neighbor){
                        d_nearest_neighbor = norm_xij;
                    }
                }

            }
            d += d_nearest_neighbor;
        }
        d /= n_points;
        KRATOS_ERROR_IF(d < 1.0e-12) << "Nearest neighbours distance is close to zero. Check that the cloud points are not overlapping." << std::endl;
        return 1/(0.815*d);// Shape parameter (Only for inverted multiquadratic)
    }

}  // namespace Kratos.