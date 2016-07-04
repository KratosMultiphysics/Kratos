//
//   Project Name:        Kratos
//   Last modified by:    $Author: zhiming guo, Riccardo Rossi
//   Date:                $Date: 2014-05-20  $
//   Revision:            $Revision: 1.0 $
//https://mail.google.com/mail/u/0/#inbox
//

#if !defined(KRATOS_MLS_SHAPE_FUNCTIONS_UTILITY_H_INCLUDED)
#define  KRATOS_MLS_SHAPE_FUNCTIONS_UTILITY_H_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <fstream>
//using namespace std;


// Project includes
#include "lagrangian_mpm_application_variables.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"

#include "includes/node.h"
#include "geometries/point_3d.h"
#include "includes/variables.h"

// External includes
#include "boost/smart_ptr.hpp"


// Project includes

#define PI 3.14159265

namespace Kratos
{



///////////////////////////////////////////////////////////////////////////////////////////


class LinearMLSKernel
{
public:

    static void ComputeKernel(double r, double h, double &kernel)
    {

        double q;
        q=r/h;
        //double Kernel;
        kernel = ( 1.0/ (PI*h*h) ) * exp(-1.0*q*q);
        //return Kernel;
    }

    static void ComputeKernelDerivative(array_1d<double,2> rvec, double h, array_1d<double,2>& w)
    {
        double r = norm_2(rvec);
        double q;
        q=r/h;

        double KernelValue;

        KernelValue=( 1.0/ (PI*h*h) ) * exp(-1.0*q*q) * (-2.0 * q * (1.0/h) );


        if (r==0.0)
        {
            w.clear();    // To avoid having nan
        }
        else noalias(w) = (KernelValue/r)*rvec;

        //return w;
    }





    ///A is the area we associate to the point
    ///N is a vector such that N(j) is the shape function associated to node j
    ///DN_DX(i,k) is the derivative of the shape function of node i, component k
    ///coordinates, is a input matrix with the coordinates of all of the points in the cloud, Coordinates(i,k) is the component k of the i-th node in the cloud
    ///nn, is number of gauss points in the cloud
    // only 2d now

    static void ComputeMLSKernel( Vector& Ng, Matrix& DN_DX, const Matrix& Coordinates, const array_1d<double,3>& x_size3, const double& h)
    {
        KRATOS_TRY;
        
        array_1d<double,2> x;
        x[0] = x_size3[0];
        x[1] = x_size3[1];

        /// for kenel and derivative of kernel computation
        //const double h=this->GetGeometry()(0)->FastGetSolutionStepValue(EFFECTIVE_RADIUS);
        array_1d<double,2> rvec;
        array_1d<double,2> w;

        /// for MlS shape function
        array_1d<double,3> p0;

        //array_1d<double,3> Ng;
       // boost::numeric::ublas::bounded_matrix<double,3,3> M = ZeroMatrix(3,3);
        Matrix M = ZeroMatrix(3,3);

        Matrix Bg(3,Coordinates.size1());
        Matrix DBg_Dx(3,Coordinates.size1());
        Matrix DBg_Dy(3,Coordinates.size1());
        array_1d<double,3> p;

        boost::numeric::ublas::bounded_matrix<double,3,3> DM_Dx = ZeroMatrix(3,3);
        boost::numeric::ublas::bounded_matrix<double,3,3> DM_Dy = ZeroMatrix(3,3);


        p0[0] = 1;
        p0[1] = x[0]; // Coordinates(0,0);
        p0[2] = x[1]; // this->GetGeometry()(0)->Y();

        for(unsigned int ii = 0; ii < Coordinates.size1(); ii++)
        {

            noalias(rvec) = x; //this->GetGeometry()(0)->Coordinates();
            noalias(rvec) -= row(Coordinates,ii); //Coordinates[ii];
            //const double r = norm_2(rvec);

            p[0] = 1.0;
            p[1] = Coordinates(ii,0);
            p[2] = Coordinates(ii,1);

            double kernel;
            double r = norm_2(rvec);
            LinearMLSKernel::ComputeKernel(r,h,kernel);
            LinearMLSKernel::ComputeKernelDerivative(rvec,h,w);
            //const double kernel = ComputeKernel(rvec,h);
            //w = ComputeKernelDerivative(rvec,h);

            noalias(M) += kernel*outer_prod( p,p );



            DM_Dx += w[0]*outer_prod( p,p );
            DM_Dy += w[1]*outer_prod( p,p );



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

        if ( fabs(determinant_M) < tol)
        {


            double M1 = 0.0;
            double DM1_Dx = 0.0;
            double DM1_Dy = 0.0;
            Vector Bg1(Coordinates.size1());
            Vector DBg1_Dx(Coordinates.size1());
            Vector DBg1_Dy(Coordinates.size1());

            for(unsigned int ii = 0; ii < Coordinates.size1(); ii++)
            {

                noalias(rvec) = x; //this->GetGeometry()(0)->Coordinates();
                noalias(rvec) -= row(Coordinates,ii); //Coordinates[ii];

                double kernel;
                double r = norm_2(rvec);
                LinearMLSKernel::ComputeKernel(r,h,kernel);
                LinearMLSKernel::ComputeKernelDerivative(rvec,h,w);
                //const double kernel = ComputeKernel(rvec,h);
                //w = ComputeKernelDerivative(rvec,h);



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
