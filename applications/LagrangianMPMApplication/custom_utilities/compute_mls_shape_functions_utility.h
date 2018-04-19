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

    static void ComputeKernel(double r, double h, double& kernel)
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

    //KernelC2
    /*static void ComputeKernel(double r, double h, double& kernel)
    {
        kernel=0.0;
        double q;
        q= r/h;
        if (q<0.0) {std::cout<<"WARNING: q is negative";}
        else if (q>=0.0 && q<1.0) {kernel=(10.0/(7.0*PI*h*h)) * (1.0-1.5*pow(q,2)+0.75*pow(q,3)); }
        else if (q>=1.0 && q<2.0) {kernel=(10.0/(7.0*PI*h*h)) * (0.25 * pow(2.0-q,3));}
        else if (q>=2.0) {kernel=0.0;}


        //return Kernel;
    }

    static void ComputeKernelDerivative(array_1d<double,2> rvec, double h, array_1d<double,2>& w)
    {
        double r = norm_2(rvec);

        double q;
        q=r/h;
        double KernelValue=0.0;

        if (q<0.0) {std::cout<<"WARNING: q is negative (GRADIENT) ";}
        else if (q>=0.0 && q<1.0) {KernelValue=(10.0/(7.0*PI*h*h)) * (-3.0 * q * (1.0/h) + 2.25 * q * q * (1.0/h)); }
        else if (q>=1.0 && q<2.0) {KernelValue=(10.0/(7.0*PI*h*h)) * (-0.75 * (1.0/h) * pow(2.0-q,2) );}
        else if (q>=2.0) {KernelValue=0;}

        //array_1d<double,3> w;
        if (r==0.0) {noalias(w)=ZeroVector(3); } // To avoid having nan
        else noalias(w) = (KernelValue/r)*rvec;

        //return w;
    }*/

    //KernelQuintic
    /*static void ComputeKernel(double r, double h, double& Kernel)
    {

        double q;
        q=r/h;

        //Calculate the Value
        Kernel=0.0;



        if (q<0.0) {std::cout<<"WARNING: q is negative";}
        else if (q>=0.0 && q<1.0) {Kernel=(7.0/(478.0*PI*h*h)) * (pow(3.0-q,5) - 6.0*pow(2.0-q,5) + 15.0*pow(1.0-q,5)); }
        else if (q>=1.0 && q<2.0) {Kernel=(7.0/(478.0*PI*h*h)) * (pow(3.0-q,5) - 6.0*pow(2.0-q,5)); }
        else if (q>=2.0 && q<3.0) {Kernel=(7.0/(478.0*PI*h*h)) * (pow(3.0-q,5)); }
        else if (q>=3.0) {Kernel=0.0;}

        //return Kernel;
    }

    static void ComputeKernelDerivative(array_1d<double,2> rvec, double h, array_1d<double,2>& w)
    {
        double r = norm_2(rvec);

        double KernelValue=0.0;
        double q;
        q=r/h;

        if (q<0.0) {std::cout<<"WARNING: q is negative (GRADIENT) ";}
        else if (q>=0.0 && q<1.0) {KernelValue=(7.0/(478.0*PI*h*h)) * ( ( 5.0 * pow(3.0-q,4) * (-1.0/h) ) - ( 30.0 * pow(2.0-q,4) * (-1.0/h) ) + ( 75.0 * pow(1.0-q,4) * (-1.0/h) ) ); }
        else if (q>=1.0 && q<2.0) {KernelValue=(7.0/(478.0*PI*h*h)) * ( ( 5.0 * pow(3.0-q,4) * (-1.0/h) ) - ( 30.0 * pow(2.0-q,4) * (-1.0/h) ) ) ; }
        else if (q>=2.0 && q<3.0) {KernelValue=(7.0/(478.0*PI*h*h)) * ( ( 5.0 * pow(3.0-q,4) * (-1.0/h) ) ) ;}
        else if (q>=3.0) {KernelValue=0;}

        //array_1d<double,3> w;
        if (r==0.0) {noalias(w)=ZeroVector(3); } // To avoid having nan
        else noalias(w) = (KernelValue/r)*rvec;

        //return w;
    }*/


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

/*class LinearULFLME
{
public:

    static void ComputeGamma(Vector& p_a, const Matrix& Coordinates, Matrix& hgam,Vector& dgam,const array_1d<double,2>& x,const double& beta, Vector& lam)
    {

        KRATOS_TRY;
        const unsigned int dim = 2;

        Vector sum1,sum2,temp;
        //array_1d<double,dim> dgam;

        double Z,gam;
        Matrix xx(Coordinates.size1(),dim),xx1(Coordinates.size1(),dim),xx2(Coordinates.size1(),dim),hgam0(dim,dim);
        Vector dgam0;



        Z = 0;
        gam = 0;
        sum1 = zero_vector<double> (Coordinates.size1());
        sum2 = zero_vector<double> (Coordinates.size1());
        temp = zero_vector<double> (Coordinates.size1());
        dgam = zero_vector<double> (dim);
        dgam0 = zero_vector<double> (dim);
        xx = ZeroMatrix(Coordinates.size1(),dim);
        xx1 = ZeroMatrix(Coordinates.size1(),dim);
        xx2 = ZeroMatrix(Coordinates.size1(),dim);
        hgam = ZeroMatrix(dim,dim);
        hgam0 = ZeroMatrix(dim,dim);

        //lam = zero_vector<double> (dim);

        for(unsigned int ii = 0; ii < Coordinates.size1(); ii++)
        {
           for(unsigned int d = 0; d < dim; d++)
           {
                //sum1(ii) += pow((x(d)-Coordinates(ii,d)),2);
                sum1(ii) += ((x(d)-Coordinates(ii,d))*(x(d)-Coordinates(ii,d)));
                sum2(ii) += lam(d)*(x(d)-Coordinates(ii,d));
            }
        }
        //KRATOS_WATCH(sum1);
        //KRATOS_WATCH(sum2);


        for(unsigned int ii = 0; ii < Coordinates.size1(); ii++)
        {
            temp(ii) = exp(-beta*sum1(ii) + sum2(ii));
            Z += temp(ii);
        }



        for(unsigned int ii = 0; ii < Coordinates.size1(); ii++)
        {
            p_a(ii) = temp(ii)/Z;
        }

        gam = log(Z);

        //KRATOS_WATCH(p_a);

        for(unsigned int d = 0; d < dim; d++)
        {
           for(unsigned int ii = 0; ii < Coordinates.size1(); ii++)
           {

                dgam0(d) += (x(d)-Coordinates(ii,d));

            }

        }

        dgam0 /= Coordinates.size1();

        for(unsigned int d = 0; d < dim; d++)
        {
           for(unsigned int ii = 0; ii < Coordinates.size1(); ii++)
           {
                //dgam(d) += (x(d)-Coordinates(ii,d))*p_a(ii);
                dgam(d) += (x(d)-Coordinates(ii,d)-dgam0(d))*p_a(ii);

            }


        }

        dgam += dgam0;

        //dgam =  prod(trans(xx),p_a);
        //KRATOS_WATCH(dgam);

        for(unsigned int d1 = 0; d1 < dim; d1++)
        {
           for(unsigned int d2 = 0; d2 < dim; d2++)
           {

               for(unsigned int ii = 0; ii < Coordinates.size1(); ii++)
               {
                   //xx1(ii,d1) = x(d1)-Coordinates(ii,d1);
                   //xx2(ii,d2) = x(d2)-Coordinates(ii,d2);

                   //hgam0(ii,d1) = p_a(ii)*xx1(ii,d1);



                   //hgam0(d1,d2) += p_a(ii)*(x(d1)-Coordinates(ii,d1))*(x(d2)-Coordinates(ii,d2));
                   hgam0(d1,d2) += p_a(ii)*(x(d1)-Coordinates(ii,d1)-dgam0(d1))*(x(d2)-Coordinates(ii,d2)-dgam0(d2));

               }


               //hgam(d1,d2) =  hgam0(d1,d2) - dgam(d1)*dgam(d2);
               hgam(d1,d2) =  hgam0(d1,d2) - (dgam(d1)-dgam0(d1))*(dgam(d2)-dgam0(d2));

            }
        }




        //hgam = hgam0 - dgam0;
        //hgam = prod(trans(hgam0), xx2) - dgam0;

        //KRATOS_WATCH(hgam);


        KRATOS_CATCH("");


    }







    ///A is the area we associate to the point
    ///N is a vector such that N(j) is the shape function associated to node j
    ///DN_DX(i,k) is the derivative of the shape function of node i, component k
    ///coordinates, is a input matrix with the coordinates of all of the points in the cloud, Coordinates(i,k) is the component k of the i-th node in the cloud
    ///nn, is number of gauss points in the cloud
    // only 2d now

    static void ComputeLMEShapef( Vector& Ng, Matrix& DN_DX, const Matrix& Coordinates, const array_1d<double,2>& x, const double& h)
    {
        KRATOS_TRY;


        const unsigned int dim = 2;
        double gamma = 1.8;
        //This sets the numerical threshold for the support of the shape functions
        double target_zero = 1.e-5;
        //This is the tolerance for the Newton iteration to compute the shape functions
        double TolLag = 1.e-6;

        double blk_spacing = h;

        double beta = gamma/(blk_spacing*blk_spacing);
        //beta=beta_(i);

        Vector lam,R,dlam;

        lam = zero_vector<double> (dim);
        R = 10*unit_vector<double> ( dim );
        dlam = 10*unit_vector<double> ( dim );

        unsigned int niter,niter_mean;
        niter = 0;
        niter_mean = 0;

        //Newton iteration
        bool iflag = 0;


        Matrix M = ZeroMatrix(2,2);
        Matrix Inverted_M = ZeroMatrix(2,2);

        while (norm_2(R)>TolLag)
        {
            LinearULFLME::ComputeGamma(Ng,Coordinates,M,R,x,beta,lam);

            //KRATOS_WATCH(norm_2(R));


            //KRATOS_WATCH(Ng);
            //KRATOS_WATCH(J);
            //KRATOS_WATCH(R);

            double determinant_M = MathUtils<double>::Det2(M);
            MathUtils<double>::InvertMatrix(M,Inverted_M,determinant_M);
            if(fabs(determinant_M) < 1e-8)
            {
                iflag = 1;
                //std::cout << "Newton Failed, near to singular matrix" << std::endl;
            }




            dlam = -prod(Inverted_M,R);

            for(unsigned int d = 0; d < dim; d++)
            {
                lam(d) += dlam(d);
            }
            //lam += dlam;
            niter++;



            //KRATOS_WATCH(determinant_M);

            if (niter>100)
            {
                iflag = 1;
                std::cout << "Newton Failed 2, no convergence in 100 iterations" << std::endl;
            }


        }



        //Spacial Gradients
        DN_DX = ZeroMatrix(Coordinates.size1(),dim);
        Matrix xx3;
        Vector pa_0;

        xx3 = ZeroMatrix(Coordinates.size1(),dim);
        pa_0 = zero_vector<double> (dim);

        for (unsigned int ia=0;ia<Coordinates.size1(); ia++)
        {
            for(unsigned int d = 0; d < dim; d++)
            {
                 //xx3(ia,d) = x(d)-Coordinates(ia,d);
                 pa_0(d) = Ng(ia)*(x(d)-Coordinates(ia,d));

                 //DN_DX(ia,d) = Ng(ia)*(x(d)-Coordinates(ia,d));
            }

            row(DN_DX,ia) = -prod(Inverted_M,pa_0);


        }


        //DN_DX = -prod(pa_0,Inverted_M);


        //niter_mean = niter_mean + niter;

        //KRATOS_WATCH(DN_DX);
        //KRATOS_WATCH(pa_0);

        KRATOS_CATCH("");

    }
};*/

}  // namespace Kratos.

#endif // KRATOS_MLS_SHAPE_FUNCTIONS_UTILITY_H_INCLUDED  defined
