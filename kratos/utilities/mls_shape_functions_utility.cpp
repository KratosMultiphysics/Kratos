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
        Matrix Bg(3, n_points);
        Matrix DBg_Dx(3, n_points);
        Matrix DBg_Dy(3, n_points);
        std::vector<array_1d<double,3>> B_vect(n_points);

        BoundedMatrix<double,3,3> M = ZeroMatrix(3,3);
        BoundedMatrix<double,3,3> DM_Dx = ZeroMatrix(3,3);
        BoundedMatrix<double,3,3> DM_Dy = ZeroMatrix(3,3);

        // Set the polynomial basis values at the point of interest
        array_1d<double,3> p0;
        p0[0] = 1;
        p0[1] = rX[0];
        p0[2] = rX[1];

        array_1d<double,2> w;
        array_1d<double,3> p;
        array_1d<double,3> rad_vect;
        BoundedMatrix<double,3,3> p_outer_mat;
        for (unsigned int i_pt = 0; i_pt < n_points; ++i_pt) {
            // Set current point data
            const array_1d<double,3>& r_i_pt_coords = row(rPoints, i_pt);
            noalias(rad_vect) = rX - r_i_pt_coords;

            // Calculate kernel values
            const double kernel = MLSShapeFunctionsUtility::CalculateKernel(rad_vect, h);
            MLSShapeFunctionsUtility::CalculateKernelDerivative<2>(rad_vect, h, w);

            // Set current point basis
            p[0] = 1.0;
            p[1] = r_i_pt_coords[0];
            p[2] = r_i_pt_coords[1];
            noalias(p_outer_mat) = outer_prod(p,p);

            // Add shape functions data
            noalias(M) += kernel*p_outer_mat;
            B_vect[i_pt] = kernel * p;
            //TODO: REMOVE THIS --> ARRANGE THE GRADIENTS FIRST
            Bg(0,i_pt) = kernel*p[0];
            Bg(1,i_pt) = kernel*p[1];
            Bg(2,i_pt) = kernel*p[2];

            // Add shape functions gradients data
            noalias(DM_Dx) += w[0]*p_outer_mat;
            noalias(DM_Dy) += w[1]*p_outer_mat;

            DBg_Dx(0,i_pt) = w[0]*p[0];
            DBg_Dx(1,i_pt) = w[0]*p[1];
            DBg_Dx(2,i_pt) = w[0]*p[2];

            DBg_Dy(0,i_pt) = w[1]*p[0];
            DBg_Dy(1,i_pt) = w[1]*p[1];
            DBg_Dy(2,i_pt) = w[1]*p[2];
        }

        // Least-Squares problem resolution
        double det_M;
        BoundedMatrix<double,3,3> Inverted_M;
        MathUtils<double>::InvertMatrix(M, Inverted_M, det_M);

        //MLS shape function
        // rN = prod( p0, prod(Inverted_M , Bg));
        array_1d<double,3> aux_prod;
        for (std::size_t i_pt = 0; i_pt < n_points; ++i_pt) {
            aux_prod = prod(Inverted_M, B_vect[i_pt]);
            rN[i_pt] = inner_prod(p0, aux_prod);
        }

        //// first derivative of MLS shape function
        array_1d<double,3> Dp0_Dx;
        array_1d<double,3> Dp0_Dy;

        /// only for 2d
        Dp0_Dx[0] = 0;
        Dp0_Dx[1] = 1;
        Dp0_Dx[2] = 0;

        Dp0_Dy[0] = 0;
        Dp0_Dy[1] = 0;
        Dp0_Dy[2] = 1;

        const array_1d<double,3> Alpha = prod( Inverted_M, p0 );
        const array_1d<double,3> DAlpha_Dx = prod(Dp0_Dx, Inverted_M) - prod(Alpha, prod(DM_Dx, Inverted_M));
        const array_1d<double,3> DAlpha_Dy = prod(Dp0_Dy, Inverted_M) - prod(Alpha, prod(DM_Dy, Inverted_M));
        const Vector DN_Dxg = prod(trans(DBg_Dx), Alpha)  + prod(trans(Bg), DAlpha_Dx);
        const Vector DN_Dyg = prod(trans(DBg_Dy), Alpha)  + prod(trans(Bg), DAlpha_Dy);

        for(std::size_t i_pt = 0; i_pt < n_points; ++i_pt) {
            rDNDX(i_pt,0) = DN_Dxg(i_pt);
            rDNDX(i_pt,1) = DN_Dyg(i_pt);
        }

        KRATOS_CATCH("");
    }

}  // namespace Kratos.
