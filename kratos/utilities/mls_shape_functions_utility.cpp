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

// System includes

// External includes

// Project includes
#include "mls_shape_functions_utility.h"

namespace Kratos
{

    void MLSShapeFunctionsUtility::CalculateKernel(
        const array_1d<double,3>& rRadVect,
        const double h,
        double& rKernel)
    {
        const double q = norm_2(rRadVect) / h;
        rKernel = (1.0/(Globals::Pi*h*h)) * std::exp(-1.0*q*q);
    }

    template<std::size_t TDim>
    void MLSShapeFunctionsUtility::CalculateKernelDerivative(
        const array_1d<double,TDim>& rRadVect,
        const double h,
        array_1d<double,TDim>& rKernelDerivative)
    {
        const double rad = norm_2(rRadVect);

        if (rad < 1.0e-12) {
            noalias(rKernelDerivative) = ZeroVector(TDim);
        } else {
            const double q = rad / h;
            const double kernel_value = (1.0/(Globals::Pi*h*h)) * std::exp(-1.0*q*q) * (-2.0*q*(1.0/h));
            noalias(rKernelDerivative) = (kernel_value / rad) * rRadVect;
        }
    }

    void MLSShapeFunctionsUtility::ComputeMLSKernel(
        Vector& Ng,
        Matrix& DN_DX,
        const Matrix& Coordinates,
        const array_1d<double,3>& x_size3,
        const double& h)
    {
        KRATOS_TRY;

        //TODO: Write a more meaningful message.
        KRATOS_ERROR_IF(h < 1.0e-12) << "Reference distance close to zero." << std::endl;

        // Set MLS shape functions containers
        const std::size_t n_points = Coordinates.size1();
        if (Ng.size() != n_points) {
            Ng.resize(n_points, false);
        }
        if (DN_DX.size1() != n_points || DN_DX.size2() != 2) {
            DN_DX.resize(n_points, 2, false);
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
        p0[1] = x_size3[0]; // Coordinates(0,0);
        p0[2] = x_size3[1]; // this->GetGeometry()(0)->Y();

        double kernel;
        array_1d<double,2> w;
        array_1d<double,3> p;
        array_1d<double,3> rad_vect;
        BoundedMatrix<double,3,3> p_outer_mat;
        for (unsigned int i_point = 0; i_point < n_points; ++i_point) {
            // Set current point data
            const array_1d<double,3>& r_i_pt_coords = row(Coordinates, i_point);
            noalias(rad_vect) = x_size3 - r_i_pt_coords;

            // Calculate kernel values
            MLSShapeFunctionsUtility::CalculateKernel(rad_vect, h, kernel);
            MLSShapeFunctionsUtility::CalculateKernelDerivative<2>(rad_vect, h, w);

            // Set current point basis
            p[0] = 1.0;
            p[1] = r_i_pt_coords[0];
            p[2] = r_i_pt_coords[1];
            noalias(p_outer_mat) = outer_prod(p,p);

            // Add shape functions data
            noalias(M) += kernel*p_outer_mat;
            B_vect[i_point] = kernel * p;
            //TODO: REMOVE THIS
            Bg(0,i_point) = kernel*p[0];
            Bg(1,i_point) = kernel*p[1];
            Bg(2,i_point) = kernel*p[2];

            // Add shape functions gradients data
            noalias(DM_Dx) += w[0]*p_outer_mat;
            noalias(DM_Dy) += w[1]*p_outer_mat;

            DBg_Dx(0,i_point) = w[0]*p[0];
            DBg_Dx(1,i_point) = w[0]*p[1];
            DBg_Dx(2,i_point) = w[0]*p[2];

            DBg_Dy(0,i_point) = w[1]*p[0];
            DBg_Dy(1,i_point) = w[1]*p[1];
            DBg_Dy(2,i_point) = w[1]*p[2];
        }

        // Least-Squares problem resolution
        const double determinant_M = MathUtils<double>::Det(M);
        if (std::abs(determinant_M) < std::numeric_limits<double>::epsilon()) {

            double M1 = 0.0;
            double DM1_Dx = 0.0;
            double DM1_Dy = 0.0;
            Vector Bg1(n_points);
            Vector DBg1_Dx(n_points);
            Vector DBg1_Dy(n_points);

            for (std::size_t i_point = 0; i_point < n_points; ++i_point) {
                // Calculate kernel values
                double kernel;
                noalias(rad_vect) = x_size3 - row(Coordinates, i_point);
                MLSShapeFunctionsUtility::CalculateKernel(rad_vect, h, kernel);
                MLSShapeFunctionsUtility::CalculateKernelDerivative<2>(rad_vect, h, w);

                // Add shape functions data
                M1 += kernel;
                Bg1(i_point) = kernel;

                // Add shape functions gradient data
                DM1_Dx += w[0];
                DM1_Dy += w[1];
                DBg1_Dx(i_point) = w[0];
                DBg1_Dy(i_point) = w[1];
            }

            Ng = Bg1/M1;

            /// only for 2d
            const double M1_pow = std::pow(M1, 2);
            const double DAlpha_Dx = -DM1_Dx/M1_pow;
            const double DAlpha_Dy = -DM1_Dy/M1_pow;
            const vector<double> DN_Dxg = (DBg1_Dx/M1)  + (Bg1*DAlpha_Dx);
            const vector<double> DN_Dyg = (DBg1_Dy/M1)  + (Bg1*DAlpha_Dy);

            for(unsigned int ii = 0; ii < Coordinates.size1(); ii++) {
                DN_DX(ii,0) = DN_Dxg(ii);
                DN_DX(ii,1) = DN_Dyg(ii);
            }

        } else {
            double det_M;
            BoundedMatrix<double,3,3> Inverted_M;
            MathUtils<double>::InvertMatrix(M, Inverted_M, det_M);

            //MLS shape function
            // Ng = prod( p0, prod(Inverted_M , Bg));
            array_1d<double,3> aux_prod;
            for (std::size_t i_point = 0; i_point < n_points; ++i_point) {
                aux_prod = prod(Inverted_M, B_vect[i_point]);
                Ng[i_point] = inner_prod(p0, aux_prod);
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

            for(std::size_t i_point = 0; i_point < n_points; ++i_point) {
                DN_DX(i_point,0) = DN_Dxg(i_point);
                DN_DX(i_point,1) = DN_Dyg(i_point);
            }
        }

        KRATOS_CATCH("");
    }

}  // namespace Kratos.
