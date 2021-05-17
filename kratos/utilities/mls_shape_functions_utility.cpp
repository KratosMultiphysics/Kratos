//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//                   Riccardo Rossi
//                   Zhiming Guo
//

// System includes

// External includes

// Project includes
#include "includes/global_variables.h"
#include "utilities/math_utils.h"
#include "mls_shape_functions_utility.h"

namespace Kratos
{

    double MLSShapeFunctionsUtility::CalculateKernel(
        const array_1d<double,3>& rRadVect,
        const double h)
    {
        const double q = norm_2(rRadVect) / h;
        return std::exp(-std::pow(q,2)) /(Globals::Pi*std::pow(h,2));
    }

    template<>
    void MLSShapeFunctionsUtility::CalculateKernelDerivative<2>(
        const array_1d<double,3>& rRadVect,
        const double h,
        array_1d<double,2>& rKernelDerivative)
    {
        const double rad = norm_2(rRadVect);

        if (rad < 1.0e-12) {
            noalias(rKernelDerivative) = ZeroVector(2);
        } else {
            const double q = rad / h;
            const double kernel_value = (std::exp(-std::pow(q,2)) * (-2.0*q/h)) / (Globals::Pi*std::pow(h,2));
            rKernelDerivative[0] = (kernel_value / rad) * rRadVect[0];
            rKernelDerivative[1] = (kernel_value / rad) * rRadVect[1];
        }
    }

    // template<std::size_t TDim>
    // void MLSShapeFunctionsUtility::CalculateKernelDerivative(
    //     const array_1d<double,3>& rRadVect,
    //     const double h,
    //     array_1d<double,TDim>& rKernelDerivative)
    // {
    //     const double rad = norm_2(rRadVect);

    //     if (rad < 1.0e-12) {
    //         noalias(rKernelDerivative) = ZeroVector(TDim);
    //     } else {
    //         const double q = rad / h;
    //         const double kernel_value = (1.0/(Globals::Pi*h*h)) * std::exp(-1.0*q*q) * (-2.0*q*(1.0/h));
    //         noalias(rKernelDerivative) = (kernel_value / rad) * rRadVect;
    //     }
    // }

    void MLSShapeFunctionsUtility::CalculateShapeFunctionsAndGradients(
        const Matrix& rPoints,
        const array_1d<double,3>& rX,
        const double h,
        Vector& rN,
        Matrix& rDNDX)
    {
        KRATOS_TRY;

        //TODO: Write a more meaningful message.
        KRATOS_ERROR_IF(h < 1.0e-12) << "Reference distance close to zero." << std::endl;

        // Set MLS shape functions containers
        const std::size_t n_points = rPoints.size1();
        if (rN.size() != n_points) {
            rN.resize(n_points, false);
        }
        if (rDNDX.size1() != n_points || rDNDX.size2() != 2) {
            rDNDX.resize(n_points, 2, false);
        }

        // Set the auxiliary arrays for the L2-norm problem minimization
        std::vector<array_1d<double,3>> B_vect(n_points);
        std::vector<array_1d<double,3>> DB_Dx_vect(n_points);
        std::vector<array_1d<double,3>> DB_Dy_vect(n_points);

        BoundedMatrix<double,3,3> M = ZeroMatrix(3,3);
        BoundedMatrix<double,3,3> DM_Dx = ZeroMatrix(3,3);
        BoundedMatrix<double,3,3> DM_Dy = ZeroMatrix(3,3);

        // Evaluate the L2-norm minimization problem
        array_1d<double,2> w;
        array_1d<double,3> p;
        array_1d<double,3> rad_vect;
        BoundedMatrix<double,3,3> p_outer_mat;
        for (std::size_t i_pt = 0; i_pt < n_points; ++i_pt) {
            // Set current point data
            const array_1d<double,3>& r_i_pt_coords = row(rPoints, i_pt);
            noalias(rad_vect) = rX - r_i_pt_coords;

            // Calculate kernel values
            const double kernel = MLSShapeFunctionsUtility::CalculateKernel(rad_vect, h);
            MLSShapeFunctionsUtility::CalculateKernelDerivative<2>(rad_vect, h, w);

            // Evaluate the current point basis
            p[0] = 1.0;
            p[1] = r_i_pt_coords[0];
            p[2] = r_i_pt_coords[1];
            noalias(p_outer_mat) = outer_prod(p,p);

            // Add shape functions data
            noalias(M) += kernel*p_outer_mat;
            B_vect[i_pt] = kernel * p;

            // Add shape functions gradients data
            noalias(DM_Dx) += w[0]*p_outer_mat;
            noalias(DM_Dy) += w[1]*p_outer_mat;
            DB_Dx_vect[i_pt] = w[0]*p;
            DB_Dy_vect[i_pt] = w[1]*p;
        }

        // Least-Squares problem resolution
        double M_det;
        BoundedMatrix<double,3,3> M_inv;
        MathUtils<double>::InvertMatrix(M, M_inv, M_det);

        // Set the polynomial basis values at the point of interest
        array_1d<double,3> p0;
        p0[0] = 1;
        p0[1] = rX[0];
        p0[2] = rX[1];

        // Set the polynomial basis x-derivative values at the point of interest
        array_1d<double,3> Dp0_Dx;
        Dp0_Dx[0] = 0;
        Dp0_Dx[1] = 1;
        Dp0_Dx[2] = 0;

        // Set the polynomial basis y-derivative values at the point of interest
        array_1d<double,3> Dp0_Dy;
        Dp0_Dy[0] = 0;
        Dp0_Dy[1] = 0;
        Dp0_Dy[2] = 1;

        // MLS shape function
        array_1d<double,3> aux_prod;
        for (std::size_t i_pt = 0; i_pt < n_points; ++i_pt) {
            aux_prod = prod(M_inv, B_vect[i_pt]);
            rN[i_pt] = inner_prod(p0, aux_prod);
        }

        // MLS shape function gradients
        const array_1d<double,3> alpha = prod(M_inv, p0);
        const array_1d<double,3> Dalpha_Dx = prod(Dp0_Dx - Vector(prod(alpha,DM_Dx)), M_inv);
        const array_1d<double,3> Dalpha_Dy = prod(Dp0_Dy - Vector(prod(alpha,DM_Dy)), M_inv);
        for(std::size_t i_pt = 0; i_pt < n_points; ++i_pt) {
            rDNDX(i_pt,0) = inner_prod(DB_Dx_vect[i_pt], alpha)  + inner_prod(B_vect[i_pt], Dalpha_Dx);
            rDNDX(i_pt,1) = inner_prod(DB_Dy_vect[i_pt], alpha)  + inner_prod(B_vect[i_pt], Dalpha_Dy);
        }

        KRATOS_CATCH("");
    }

}  // namespace Kratos.
