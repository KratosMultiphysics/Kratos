#if !defined(KRATOS_LME_SHAPE_FUNCTIONS_UTILITY_H_INCLUDED)
#define  KRATOS_LME_SHAPE_FUNCTIONS_UTILITY_H_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <fstream>
//using namespace std;




/* Project includes */
#include "includes/define.h"

#include "includes/ublas_interface.h"
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
class LinearULFLME //LinearLMEKernel
{
public:

	    template<class T>
        static bool InvertMatrix(T& input, T& inverse)
    {
        typedef permutation_matrix<std::size_t> pmatrix;

        // create a working copy of the input
        T A(input);

        // create a permutation matrix for the LU-factorization
        pmatrix pm(A.size1());

        // perform LU-factorization
        int res = lu_factorize(A, pm);
        if (res != 0)
            return false;

        // create identity matrix of "inverse"
        inverse.assign(identity_matrix<double> (A.size1()));

        // backsubstitute to get the inverse
        lu_substitute(A, pm, inverse);

        return true;
    }


    static void ComputeGamma(Vector& p_a, const Matrix& Coordinates, Matrix& hgam,Vector& dgam,const array_1d<double,2>& x,Vector&  beta, Vector& lam, double& gam,int& el_id)
    //static void ComputeGamma(Vector& p_a, const Matrix& Coordinates, Matrix& hgam,Vector& dgam,const array_1d<double,2>& x,const double& beta, Vector& lam)
    {

        KRATOS_TRY;
        const unsigned int dim = 2;

        Vector sum1,sum2,temp;
        //array_1d<double,dim> dgam;

        double Z;
        Matrix hgam0(dim,dim);
 



        Z = 0;
        //gam = 0;
        sum1 = zero_vector<double> (Coordinates.size1());
        sum2 = zero_vector<double> (Coordinates.size1());
        temp = zero_vector<double> (Coordinates.size1());
        dgam = zero_vector<double> (dim);
  
        hgam = ZeroMatrix(dim,dim);
        hgam0 = ZeroMatrix(dim,dim);



        for(unsigned int ii = 0; ii < Coordinates.size1(); ii++)
        {
           for(unsigned int d = 0; d < dim; d++)
           {
         
                sum1(ii) += ((x(d)-Coordinates(ii,d))*(x(d)-Coordinates(ii,d)));
                sum2(ii) += lam(d)*(x(d)-Coordinates(ii,d));
            }
        }
   


        for(unsigned int ii = 0; ii < Coordinates.size1(); ii++)
        {
            temp(ii) = exp(-beta(ii)*sum1(ii) + sum2(ii));
            Z += temp(ii);
        }



        for(unsigned int ii = 0; ii < Coordinates.size1(); ii++)
        {
            p_a(ii) = temp(ii)/Z;
        }

        gam = log(Z);
		

        for(unsigned int d = 0; d < dim; d++)
        {
           for(unsigned int ii = 0; ii < Coordinates.size1(); ii++)
           {
      
                dgam(d) += (x(d)-Coordinates(ii,d))*p_a(ii);

            }
        }
   

        for(unsigned int d1 = 0; d1 < dim; d1++)
        {
           for(unsigned int d2 = 0; d2 < dim; d2++)
           {

               for(unsigned int ii = 0; ii < Coordinates.size1(); ii++)
               {
  
                   hgam0(d1,d2) += p_a(ii)*(x(d1)-Coordinates(ii,d1))*(x(d2)-Coordinates(ii,d2));

               }


               hgam(d1,d2) =  hgam0(d1,d2) - dgam(d1)*dgam(d2);

            }
        }


 


        KRATOS_CATCH("");


    }

    
    static void ComputeGamma(Vector& p_a, const Matrix& Coordinates, Matrix& hgam,Vector& dgam,const array_1d<double,2>& x,const double& beta, Vector& lam)
    {

        KRATOS_TRY;
        const unsigned int dim = 2;

        Vector sum1,sum2,temp;
        //array_1d<double,dim> dgam;

        double Z;//,gam;
        Matrix hgam0(dim,dim);
 



        Z = 0;
        //gam = 0;
        sum1 = zero_vector<double> (Coordinates.size1());
        sum2 = zero_vector<double> (Coordinates.size1());
        temp = zero_vector<double> (Coordinates.size1());
        dgam = zero_vector<double> (dim);
  
        hgam = ZeroMatrix(dim,dim);
        hgam0 = ZeroMatrix(dim,dim);



        for(unsigned int ii = 0; ii < Coordinates.size1(); ii++)
        {
           for(unsigned int d = 0; d < dim; d++)
           {
         
                sum1(ii) += ((x(d)-Coordinates(ii,d))*(x(d)-Coordinates(ii,d)));
                sum2(ii) += lam(d)*(x(d)-Coordinates(ii,d));
            }
        }
   


        for(unsigned int ii = 0; ii < Coordinates.size1(); ii++)
        {
            temp(ii) = exp(-beta*sum1(ii) + sum2(ii));
            Z += temp(ii);
        }



        for(unsigned int ii = 0; ii < Coordinates.size1(); ii++)
        {
            p_a(ii) = temp(ii)/Z;
        }

        //double gam = log(Z);

  

        for(unsigned int d = 0; d < dim; d++)
        {
           for(unsigned int ii = 0; ii < Coordinates.size1(); ii++)
           {
      
                dgam(d) += (x(d)-Coordinates(ii,d))*p_a(ii);

            }
        }
   

        for(unsigned int d1 = 0; d1 < dim; d1++)
        {
           for(unsigned int d2 = 0; d2 < dim; d2++)
           {

               for(unsigned int ii = 0; ii < Coordinates.size1(); ii++)
               {
  
                   hgam0(d1,d2) += p_a(ii)*(x(d1)-Coordinates(ii,d1))*(x(d2)-Coordinates(ii,d2));

               }


               hgam(d1,d2) =  hgam0(d1,d2) - dgam(d1)*dgam(d2);

            }
        }


 


        KRATOS_CATCH("");


    }




    static void ComputeLMEShapef( Vector& Ng, Matrix& DN_DX, const Matrix& Coordinates, const array_1d<double,2>& x, Vector& h, int& el_id, bool& iflag)
    //static void ComputeLMEShapef( Vector& Ng, Matrix& DN_DX, const Matrix& Coordinates, const array_1d<double,2>& x, const double& h)
    {
        KRATOS_TRY;


        const unsigned int dim = 2;
        double gamma = 1.0;
        //This sets the numerical threshold for the support of the shape functions
        //double target_zero = 1.e-5;
        //This is the tolerance for the Newton iteration to compute the shape functions
        double TolLag = 1e-10;
       
        double log_Z = 1.0;
        
        //beta=beta_(i);

        Vector lam,R,dlam, lam_aux;

        lam = zero_vector<double> (dim);
        R = 10*unit_vector<double> ( dim );
        dlam = 10*unit_vector<double> ( dim );


		Vector beta;
		beta = zero_vector<double> (Coordinates.size1());

		for(unsigned int i = 0; i<Coordinates.size1();i++)
		{
			
			beta(i) = gamma /( h(i)*h(i));
			
		}


        unsigned int niter;//,niter_mean;
        niter = 0;
        //niter_mean = 0;

        //Newton iteration
        //bool iflag = 0;


        Matrix M = ZeroMatrix(dim, dim);
        Matrix Inverted_M = ZeroMatrix(dim,dim);

//************************************************************************************************			
//************************************************************************************************			        
		////LINE SEARCH
		
		Vector R0, R1_2, R1;//, R_old, R_new, DeltaR;
		Vector lam0, lam1_2, lam1;
		R0 = 10*unit_vector<double> ( dim );
		R1_2 = 10*unit_vector<double> ( dim );
		R1 = 10*unit_vector<double> ( dim );
		
		//R_old = 10*unit_vector<double> ( dim );
		//R_new = 10*unit_vector<double> ( dim );
		//DeltaR = 10*unit_vector<double> ( dim );
			
		lam0 = zero_vector<double> (dim);
		lam1_2 = zero_vector<double> (dim);
		lam1 = zero_vector<double> (dim);
			
		Matrix M0 = ZeroMatrix(dim,dim);	
		Matrix M1_2 = ZeroMatrix(dim,dim);
		Matrix M1 = ZeroMatrix(dim,dim);
			
		Vector Ng0, Ng1_2, Ng1;
		Ng0 = zero_vector<double> (Coordinates.size1());	
		Ng1_2 = zero_vector<double> (Coordinates.size1());	
		Ng1 = zero_vector<double> (Coordinates.size1());	
			
		//double norm_R0, norm_R1_2, norm_R1;	
		
		while (norm_2(R)>TolLag)
        {
            niter++;
            
            //R_old = R;
            
            LinearULFLME::ComputeGamma(Ng,Coordinates,M,R,x,beta,lam, log_Z,el_id);

			//R_new = R;
			//DeltaR = R_new - R_old;
            
            //double determinant_M = MathUtils<double>::Det2(M);
            //MathUtils<double>::InvertMatrix(M,Inverted_M,determinant_M);
            
            bool inversion_successful = LinearULFLME::InvertMatrix(M, Inverted_M);

            if (inversion_successful == false)
            {
                iflag = 1;
                //KRATOS_ERROR << "Matrix inversion error within the Jacobian approximation computation!!";
            }
            
            
            double InvCondNumber = 0.0;
            double normM = 0.0;
            double normInvM = 0.0;
            for(unsigned int i = 0; i<dim; i++)
            {
				for(unsigned int j = 0;j<dim;j++)
				{
					normM += M(i,j) * M(i,j);
					normInvM += Inverted_M(i,j) * Inverted_M(i,j);
				}
			}
			
            normM = sqrt(normM);
            normInvM = sqrt(normInvM);
            InvCondNumber = 1/(normM * normInvM);
            
            if(InvCondNumber < 1e-9 || norm_2(Ng) == 0.0)//(fabs(determinant_M) < 1e-9 || norm_2(Ng) == 0.0)
            {
                iflag = 1;
               
                
				break;
                
            }  
			
            dlam = -prod(R,Inverted_M);
            
         
			
			double gTdlam = R(0) * dlam(0) + R(1) * dlam(1);
            if(gTdlam > 0)
            {
				for(unsigned int d = 0; d < dim; d++)
				{
				dlam(d) *= (-1);
				}
			}
			
            
            //std::cout << "dlam after" << dlam<< std::endl;
            
	
			//double alfa0 = 0.0;
			double log_Z0 = 0.0;
			lam0 = lam;
			LinearULFLME::ComputeGamma(Ng0,Coordinates,M0,R0,x,beta,lam0,log_Z0,el_id);
			//norm_R0 = norm_2(R0);
			
////**********************************************************************************************
//BACKTRACKING-ARMIJO LINESEARCH

//initial value of alfa
		    //double alfa_v = 0.0;
		    double alfa_it = 1.0;
		    //double alfa_v0 = 0.0;
		    double B1 = 0.001;
		    double tau = 0.5;
		    
		    double log_Zaux = log_Z0 + alfa_it * B1 * gTdlam;
		    
		    
		    lam_aux = zero_vector<double> (dim);
		    
		    
		    
		    lam_aux = lam + dlam * alfa_it;
		    LinearULFLME::ComputeGamma(Ng,Coordinates,M,R,x,beta,lam_aux, log_Z,el_id);
		    int lcount = 0;
		    
		    //double determinant_M = MathUtils<double>::Det2(M);
		    
		    
		   
			while(log_Z > log_Zaux || norm_2(Ng) == 0.0)
			{
				
				alfa_it = tau * alfa_it;
				
				
				lam_aux = lam + dlam * alfa_it;
				LinearULFLME::ComputeGamma(Ng,Coordinates,M,R,x,beta,lam_aux, log_Z,el_id);
				
				lcount++;
				
				log_Zaux = log_Z0 + alfa_it * B1 * gTdlam;
				
				
				if(alfa_it< 1e-4)//norm_2(lam_aux) == 0.0)
				    {
					iflag = 1;
					break; //I don't know if I have to use a break here..
					}  
			}            //I saw there is the case that I can get a vector of nan  
			
			//alfa_v = alfa_it;
			lam = lam_aux;
			
			
		
////**********************************************************************************************		
//************************************************************************************************		

			
            
            
           
            if (niter>20)
            {
                iflag = 1;
                //std::cout << "norm_2(R) " << norm_2(R)<<std::endl;
                //if(norm_2(R)>1e-9)
                //{
                //std::cout << "Newton Failed 2, no convergence in 20 iterations" << std::endl;
                //std::cout << "norm_2(R) " << norm_2(R)<<std::endl;
                //std::cout<<"ID ELEMENT"<<el_id<<std::endl;
				//}
                break;
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
 
                 pa_0(d) = Ng(ia)*(x(d)-Coordinates(ia,d));

 
            }

            row(DN_DX,ia) = -prod(pa_0,Inverted_M);


        }
//************************************************************************************************		
//************************************************************************************************		


        KRATOS_CATCH("");

    }


    
    static void ComputeLMEShapef( Vector& Ng, Matrix& DN_DX, const Matrix& Coordinates, const array_1d<double,2>& x, const double& h)
    {
        KRATOS_TRY;


        const unsigned int dim = 2;
        double gamma = 0.8;
        //This sets the numerical threshold for the support of the shape functions
        //double target_zero = 1.e-5;
        //This is the tolerance for the Newton iteration to compute the shape functions
        double TolLag = 1.e-8;

//**********************************************************************
        double blk_spacing = h;

        double beta = gamma/(blk_spacing*blk_spacing);
//**********************************************************************
        
        
        
        //beta=beta_(i);

        Vector lam,R,dlam;

        lam = zero_vector<double> (dim);
        R = 10*unit_vector<double> ( dim );
        dlam = 10*unit_vector<double> ( dim );


        unsigned int niter;//,niter_mean;
        niter = 0;
        //niter_mean = 0;

        //Newton iteration
        bool iflag = 0;


        Matrix M = ZeroMatrix(2,2);
        Matrix Inverted_M = ZeroMatrix(2,2);

        //NEWTON-RAPHSON
        while (norm_2(R)>TolLag)
        {
            niter++;
            LinearULFLME::ComputeGamma(Ng,Coordinates,M,R,x,beta,lam);



            double determinant_M = MathUtils<double>::Det2(M);
            MathUtils<double>::InvertMatrix(M,Inverted_M,determinant_M);
            if(fabs(determinant_M) < 1e-8)
            {
                iflag = 1;
            }  
            



            dlam = -prod(R,Inverted_M);

			for(unsigned int d = 0; d < dim; d++)
            {
                lam(d) += dlam(d);
            }
            //lam += dlam;
            

			if(iflag == 1)
            { 
                //std::cout << "niter " << niter<<std::endl;
                //std::cout << "Ng " << Ng<<std::endl;
                //std::cout << "Coordinates " << Coordinates<<std::endl;
                //std::cout << "x " << x<<std::endl;
                //std::cout << "beta " << beta<<std::endl;
                //std::cout << "lam " << lam<<std::endl;
                //std::cout << "norm_2(R) " << norm_2(R)<<std::endl;
                //std::cout << "M " << M<<std::endl;
                //std::cout << "determinant_M " << determinant_M<<std::endl;
                
                std::cout << "Newton Failed, near to singular matrix" << std::endl;
            }


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
 
                 pa_0(d) = Ng(ia)*(x(d)-Coordinates(ia,d));

 
            }

            row(DN_DX,ia) = -prod(pa_0,Inverted_M);


        }



        KRATOS_CATCH("");

    }
};
//class LinearULFLME
//{
//public:

    //static void ComputeGamma(Vector& p_a, const Matrix& Coordinates, Matrix& hgam,Vector& dgam,const array_1d<double,2>& x,const double& beta, Vector& lam)
    //{

        //KRATOS_TRY;
        //const unsigned int dim = 2;

        //Vector sum1,sum2,temp;
        ////array_1d<double,dim> dgam;

        //double Z,gam;
        //Matrix xx(Coordinates.size1(),dim),xx1(Coordinates.size1(),dim),xx2(Coordinates.size1(),dim),hgam0(dim,dim);
        //Vector dgam0;



        //Z = 0.0;
        //gam = 0.0;
        //sum1 = zero_vector<double> (Coordinates.size1());
        //sum2 = zero_vector<double> (Coordinates.size1());
        //temp = zero_vector<double> (Coordinates.size1());
        //dgam = zero_vector<double> (dim);
        //dgam0 = zero_vector<double> (dim);
        //xx = ZeroMatrix(Coordinates.size1(),dim);
        //xx1 = ZeroMatrix(Coordinates.size1(),dim);
        //xx2 = ZeroMatrix(Coordinates.size1(),dim);
        //hgam = ZeroMatrix(dim,dim);
        //hgam0 = ZeroMatrix(dim,dim);

        ////lam = zero_vector<double> (dim);

        //for(unsigned int ii = 0; ii < Coordinates.size1(); ii++)
        //{
           //for(unsigned int d = 0; d < dim; d++)
           //{

                //sum1(ii) += ((x(d)-Coordinates(ii,d))*(x(d)-Coordinates(ii,d)));
                //sum2(ii) += lam(d)*(x(d)-Coordinates(ii,d));
            //}
        //}



        //for(unsigned int ii = 0; ii < Coordinates.size1(); ii++)
        //{
            //temp(ii) = exp(-beta*sum1(ii) + sum2(ii));
            //Z += temp(ii);
        //}

		////Z += 1e-30;

        //for(unsigned int ii = 0; ii < Coordinates.size1(); ii++)
        //{
            //p_a(ii) = temp(ii)/(Z);
        //}

        //gam = log(Z);

        ////KRATOS_WATCH(p_a);

        //for(unsigned int d = 0; d < dim; d++)
        //{
           //for(unsigned int ii = 0; ii < Coordinates.size1(); ii++)
           //{

                //dgam0(d) += (x(d)-Coordinates(ii,d));

            //}

        //}

        //dgam0 /= Coordinates.size1();   //this is like an average distance

        //for(unsigned int d = 0; d < dim; d++)
        //{
           //for(unsigned int ii = 0; ii < Coordinates.size1(); ii++)
           //{

                //dgam(d) += (x(d)-Coordinates(ii,d)-dgam0(d))*p_a(ii);

            //}


        //}

        //dgam += dgam0;

        //for(unsigned int d1 = 0; d1 < dim; d1++)
        //{
           //for(unsigned int d2 = 0; d2 < dim; d2++)
           //{

               //for(unsigned int ii = 0; ii < Coordinates.size1(); ii++)
               //{

                   //hgam0(d1,d2) += p_a(ii)*(x(d1)-Coordinates(ii,d1)-dgam0(d1))*(x(d2)-Coordinates(ii,d2)-dgam0(d2));

               //}



               //hgam(d1,d2) =  hgam0(d1,d2) - (dgam(d1)-dgam0(d1))*(dgam(d2)-dgam0(d2));

            //}
        //}



        //KRATOS_CATCH("");


    //}


    //static void ComputeLMEShapef( Vector& Ng, Matrix& DN_DX, const Matrix& Coordinates, const array_1d<double,2>& x, const double& h)
    //{
        //KRATOS_TRY;


        //const unsigned int dim = 2;
        //double gamma = 0.8;
        ////This sets the numerical threshold for the support of the shape functions
        //double target_zero = 1.e-5;
        ////This is the tolerance for the Newton iteration to compute the shape functions
        //double TolLag = 1.e-8;

        //double blk_spacing = h;

        //double beta = gamma/(blk_spacing*blk_spacing);
        ////beta=beta_(i);

        //Vector lam,R,dlam;

		////lam: Lagrange multipliers vector
        //lam = zero_vector<double> (dim);
        ////R: residual vector
        //R = 10*unit_vector<double> ( dim );
        //dlam = 10*unit_vector<double> ( dim );

        //unsigned int niter,niter_mean;
        //niter = 0;
        //niter_mean = 0;

        ////Newton iteration
        //bool iflag = 0;

		////M: Hessain matrix
        //Matrix M = ZeroMatrix(2,2);
        //Matrix Inverted_M = ZeroMatrix(2,2);

        //while (norm_2(R)>TolLag)
        //{
            //LinearULFLME::ComputeGamma(Ng,Coordinates,M,R,x,beta,lam);

            ////KRATOS_WATCH(norm_2(R));


            ////KRATOS_WATCH(Ng);
            ////KRATOS_WATCH(J);
            ////KRATOS_WATCH(R);

            //double determinant_M = MathUtils<double>::Det2(M);
            //MathUtils<double>::InvertMatrix(M,Inverted_M,determinant_M);
            //if(fabs(determinant_M) < 1e-8)
            //{
                //iflag = 1;
                //std::cout << "niter " << niter<<std::endl;
                //std::cout << "Ng " << Ng<<std::endl;
                //std::cout << "Coordinates " << Coordinates<<std::endl;
                //std::cout << "x " << x<<std::endl;
                //std::cout << "beta " << beta<<std::endl;
                //std::cout << "lam " << lam<<std::endl;
                //std::cout << "M " << M<<std::endl;
                //std::cout << "determinant_M " << determinant_M<<std::endl;
                //std::cout << "Newton Failed, near to singular matrix" << std::endl;
                
            //}


			////Update Lagrange multipliers

            //dlam = -prod(Inverted_M,R);

            //for(unsigned int d = 0; d < dim; d++)
            //{
                //lam(d) += dlam(d);
            //}
            ////lam += dlam;
            //niter++;
			


            ////KRATOS_WATCH(determinant_M);

            //if (niter>100)
            //{
                //iflag = 1;
                //std::cout << "Newton Failed 2, no convergence in 100 iterations" << std::endl;
            //}

		//for(unsigned int i = 0; i< Coordinates.size1(); i++)
		//{
			//if(Ng(i) < 1e-16)
			//{
				//Ng(i) = 0.0;
			//}
			
		//}
        //}

		

        ////Spacial Gradients
        //DN_DX = ZeroMatrix(Coordinates.size1(),dim);
        //Matrix xx3;
        //Vector pa_0;

        //xx3 = ZeroMatrix(Coordinates.size1(),dim);
        //pa_0 = zero_vector<double> (dim);

        //for (unsigned int ia=0;ia<Coordinates.size1(); ia++)
        //{
            //for(unsigned int d = 0; d < dim; d++)
            //{
                 ////xx3(ia,d) = x(d)-Coordinates(ia,d);
                 //pa_0(d) = Ng(ia)*(x(d)-Coordinates(ia,d));

                 ////DN_DX(ia,d) = Ng(ia)*(x(d)-Coordinates(ia,d));
            //}

            //row(DN_DX,ia) = -prod(Inverted_M,pa_0);


        //}




        //KRATOS_CATCH("");

    //}
//};
}
#endif
