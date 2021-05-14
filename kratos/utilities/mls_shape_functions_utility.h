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
        double& rKernel)
    {
        const double q = norm_2(rRadVect) / h;
        rKernel = (1.0/(Globals::Pi*h*h)) * std::exp(-1.0*q*q);
    }

    template<std::size_t TDim>
    static void CalculateKernelDerivative(
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
        const double& h)
    {
        KRATOS_TRY;

        //TODO: Write a more meaningful message.
        KRATOS_ERROR_IF(h < 1.0e-12) << "Reference distance close to zero." << std::endl;



        /// for MlS shape function

        Matrix Bg(3,Coordinates.size1());
        Matrix DBg_Dx(3,Coordinates.size1());
        Matrix DBg_Dy(3,Coordinates.size1());

        BoundedMatrix<double,3,3> M = ZeroMatrix(3,3);
        BoundedMatrix<double,3,3> DM_Dx = ZeroMatrix(3,3);
        BoundedMatrix<double,3,3> DM_Dy = ZeroMatrix(3,3);

        array_1d<double,3> p0;
        p0[0] = 1;
        p0[1] = x_size3[0]; // Coordinates(0,0);
        p0[2] = x_size3[1]; // this->GetGeometry()(0)->Y();

        array_1d<double,2> w;
        array_1d<double,3> p;
        array_1d<double,3> rad_vect;
        BoundedMatrix<double,3,3> p_outer_mat;
        const std::size_t n_points = Coordinates.size();
        for (unsigned int i_point = 0; i_point < n_points; ++i_point) {
            // Calculate kernel values
            double kernel;
            const array_1d<double,3>& r_i_pt_coords = row(Coordinates, i_point);
            noalias(rad_vect) = x_size3 - r_i_pt_coords;
            MLSShapeFunctionsUtility::CalculateKernel(rad_vect, h, kernel);
            MLSShapeFunctionsUtility::CalculateKernelDerivative<2>(rad_vect, h, w);

            p[0] = 1.0;
            p[1] = r_i_pt_coords[0];
            p[2] = r_i_pt_coords[1];
            noalias(p_outer_mat) = outer_prod(p,p);

            noalias(M) += kernel*p_outer_mat;


            DM_Dx += w[0]*p_outer_mat;
            DM_Dy += w[1]*p_outer_mat;

            Bg(0,ii) = kernel*p[0];
            Bg(1,ii) = kernel*p[1];
            Bg(2,ii) = kernel*p[2];


            /// the derivatives of Bg
            DBg_Dx(0,ii) = w[0]*p[0];
            DBg_Dx(1,ii) = w[0]*p[1];
            DBg_Dx(2,ii) = w[0]*p[2];

            DBg_Dy(0,ii) = w[1]*p[0];
            DBg_Dy(1,ii) = w[1]*p[1];
            DBg_Dy(2,ii) = w[1]*p[2];


        }
        //// moment matrix inverse
        const double tol = 1e-16;

        double determinant_M = MathUtils<double>::Det3(M);

        ///return false

        Matrix Inverted_M = ZeroMatrix(3,3);

        if (std::abs(determinant_M) < tol)
        {


            double M1 = 0.0;
            double DM1_Dx = 0.0;
            double DM1_Dy = 0.0;
            Vector Bg1(Coordinates.size1());
            Vector DBg1_Dx(Coordinates.size1());
            Vector DBg1_Dy(Coordinates.size1());

            for (std::size_t i_point = 0; i_point < n_points; ++i_point) {
                // Calculate kernel values
                double kernel;
                noalias(rad_vect) = x_size3 - row(Coordinates, i_point);
                MLSShapeFunctionsUtility::CalculateKernel(rad_vect, h, kernel);
                MLSShapeFunctionsUtility::CalculateKernelDerivative(rad_vec, h, w);

                M1 += kernel;

                DM1_Dx += w[0];
                DM1_Dy += w[1];

                Bg1(ii) = kernel;

                /// the derivatives of Bg
                DBg1_Dx(ii) = w[0];
                DBg1_Dy(ii) = w[1];
            }

            Ng = Bg1/M1;

            double DAlpha_Dx;
            double DAlpha_Dy;
            vector<double> DN_Dxg;
            vector<double> DN_Dyg;
            DN_DX = ZeroMatrix(Coordinates.size1(),Coordinates.size2());


            /// only for 2d

            DAlpha_Dx = - DM1_Dx/(M1*M1);
            DAlpha_Dy =  - DM1_Dy/(M1*M1);

            DN_Dxg = (DBg1_Dx/M1)  + (Bg1*DAlpha_Dx  );
            DN_Dyg = (DBg1_Dy/M1)  + (Bg1*DAlpha_Dy  );

            for(unsigned int ii = 0; ii < Coordinates.size1(); ii++)
            {
                DN_DX(ii,0) = DN_Dxg(ii);
                DN_DX(ii,1) = DN_Dyg(ii);

            }



        }
        else
        {

            MathUtils<double>::InvertMatrix3(M,Inverted_M,determinant_M);



            //MLS shape function
           // Ng = prod( prod( p0, Inverted_M ), Bg);
            Ng = prod( p0, prod(Inverted_M , Bg));


            //// first derivative of MLS shape function
            array_1d<double,3> Dp0_Dx;
            array_1d<double,3> Dp0_Dy;
            array_1d<double,3> Alpha;
            array_1d<double,3> DAlpha_Dx;
            array_1d<double,3> DAlpha_Dy;
            vector<double> DN_Dxg;
            vector<double> DN_Dyg;
            DN_DX = ZeroMatrix(Coordinates.size1(),Coordinates.size2());


            /// only for 2d
            Dp0_Dx[0] = 0;
            Dp0_Dx[1] = 1;
            Dp0_Dx[2] = 0;

            Dp0_Dy[0] = 0;
            Dp0_Dy[1] = 0;
            Dp0_Dy[2] = 1;

            Alpha = prod( Inverted_M, p0 );
            DAlpha_Dx = prod( Dp0_Dx, Inverted_M ) - prod(  Alpha, prod(DM_Dx, Inverted_M) );
            DAlpha_Dy = prod( Dp0_Dy, Inverted_M ) - prod(  Alpha, prod(DM_Dy, Inverted_M) );



            DN_Dxg = prod( trans(DBg_Dx),Alpha  )  + prod( trans(Bg),DAlpha_Dx  );
            DN_Dyg = prod( trans(DBg_Dy),Alpha  )  + prod( trans(Bg),DAlpha_Dy  );

            for(unsigned int ii = 0; ii < Coordinates.size1(); ii++)
            {
                DN_DX(ii,0) = DN_Dxg(ii);
                DN_DX(ii,1) = DN_Dyg(ii);

            }


        }

        KRATOS_CATCH("");

    }
};

}  // namespace Kratos.

#endif // KRATOS_MLS_SHAPE_FUNCTIONS_UTILITY_H_INCLUDED  defined
