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
#include "radial_basis_functions_utility.h"

namespace Kratos
{

    double RadialBasisFunctionsUtility::EvaluateRBF(
        double x,
        const double h)
    {   
        // Evaluate Inverse multiquadric
        const double q = x*h;
        return 1/std::sqrt(1+std::pow(q,2));

        // Evaluate Gaussian function
        // const double q = x/h;
        // return std::exp(-0.5*std::pow(q,2));
    }

    void RadialBasisFunctionsUtility::CalculateShapeFunctions(
        const Matrix& rPoints,
        const array_1d<double,3>& rX,
        const double h,
        Vector& rN)
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
        double norm_xij;
        double norm_X;

        // Build the RBF interpolation matrix and RBF interpolated vector
        for (std::size_t i_pt = 0; i_pt < n_points; ++i_pt) {
            for (std::size_t j_pt = 0; j_pt < n_points; ++j_pt) {
                norm_xij = norm_2(row(rPoints,i_pt)-row(rPoints,j_pt));
                A(i_pt,j_pt) = EvaluateRBF(norm_xij,h);
            }
            norm_X =  norm_2(rX-row(rPoints,i_pt));
            Phi[i_pt] = EvaluateRBF(norm_X,h);
        }

        // Invert the RBF interpolation matrix
        double A_det;
        Matrix A_inv(n_points,n_points);
        MathUtils<double>::InvertMatrix(A, A_inv, A_det);

        // Obtain the RBF shape functions 
        noalias(rN) = prod(A_inv,Phi);
        

        // Partition of unity
        noalias(rN) =  rN/sum(rN);

        KRATOS_CATCH("");
    }

    void RadialBasisFunctionsUtility::CalculateShapeFunctions(
        const Matrix& rPoints,
        const array_1d<double,3>& rX,
        Vector& rN)
    {
        KRATOS_TRY;

        // Set RBF shape functions containers
        const std::size_t n_points = rPoints.size1();
        if (rN.size() != n_points) {
            rN.resize(n_points, false);
        }

        // Initialize the RBF interpolation matrix and RBF interpolated vector
        Matrix A = ZeroMatrix(n_points,n_points);
        Vector Phi = ZeroVector(n_points);
        double norm_xij;
        double norm_X;
        

        // Find shape parameter http://www.math.iit.edu/~fass/Dolomites.pdf [Hardy]
        // TODO: Build a matrix with distance norm values to recycle for RBFs evaluation.
        // Make a function for this 
        double d = 0; // Total distance of nearest x_i neighbors 
        double h; // Shape parameter
        for (std::size_t i_pt = 0; i_pt < n_points; ++i_pt) {
            double d_nearest_neighbor = 1e30; // Initialize distance to nearest x_i neighbor 
            for (std::size_t j_pt = 0; j_pt < n_points; ++j_pt) {
                if (i_pt!=j_pt){ // Avoid measuring distance between same points
                    norm_xij = norm_2(row(rPoints,i_pt)-row(rPoints,j_pt)); 
                    if (norm_xij < d_nearest_neighbor){
                        d_nearest_neighbor = norm_xij;
                    }
                }
                
            }
            d += d_nearest_neighbor;
        }
        d /= n_points;
        h = 1/(0.815*d);// Only for inverted multiquadratic
        

        // Build the RBF interpolation matrix and RBF interpolated vector
        // TODO: Make a function for this, get the norm values from the matrix created for finding the shape parameter
        for (std::size_t i_pt = 0; i_pt < n_points; ++i_pt) {
            for (std::size_t j_pt = 0; j_pt < n_points; ++j_pt) {
                norm_xij = norm_2(row(rPoints,i_pt)-row(rPoints,j_pt));
                A(i_pt,j_pt) = EvaluateRBF(norm_xij,h);
            }
            norm_X =  norm_2(rX-row(rPoints,i_pt));
            Phi[i_pt] = EvaluateRBF(norm_X,h);
        }

        // Invert the RBF interpolation matrix
        double A_det;
        Matrix A_inv(n_points,n_points);
        MathUtils<double>::InvertMatrix(A, A_inv, A_det);

        // Obtain the RBF shape functions 
        noalias(rN) = prod(Phi,A_inv);

        // Partition of unity
        noalias(rN) =  rN/sum(rN);

        KRATOS_CATCH("");
    }

    double RadialBasisFunctionsUtility::CalculateShapeFunctionsAndInterpolation(
        const Matrix& rPoints,
        const array_1d<double,3>& rX,
        const double h,
        Vector& rN,
        Vector Y)
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

        // Invert the RBF interpolation matrix
        double A_det;
        Matrix A_inv(n_points,n_points);
        MathUtils<double>::InvertMatrix(A, A_inv, A_det);

        // Obtain the RBF shape functions 
        noalias(rN) = prod(A_inv,Y);

        // Interpolate solution
        double Phi;
        double norm_xi;
        for (std::size_t i_pt = 0; i_pt < n_points; ++i_pt) {
            norm_xi = norm_2(rX-row(rPoints,i_pt));
            Phi = EvaluateRBF(norm_xi,h);
            interpolation += Phi*rN[i_pt];
        }


        // Partition of unity
        // noalias(rN) =  rN/sum(rN);

        KRATOS_CATCH("");

        return interpolation;
    }

}  // namespace Kratos.