  /*
==============================================================================
KratosStructuralApplication 
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel 
pooyan@cimne.upc.edu 
rrossi@cimne.upc.edu
janosch.stascheit@rub.de
nagel@sd.rub.de
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-University Bochum, Institute for Structural Mechanics, Germany


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/


//#include "utilities/math_utils.h"
//#include "custom_utilities/sd_math_utils.h"
#include "custom_utilities/tensor_utils.h"
#include "fluency_criteria/energy_yield_function.h"
#include <cmath>



namespace Kratos
  {

  
	    typedef boost::numeric::ublas::vector<Vector> Second_Order_Tensor; // dos opciones: un tensor de segundo orden y/o un vector que almacena un vector		  
	    typedef boost::numeric::ublas::vector<Second_Order_Tensor> Third_Order_Tensor;
			  
            typedef boost::numeric::ublas::vector<boost::numeric::ublas::vector<Matrix> > Fourth_Order_Tensor;
			  
	    typedef matrix<Second_Order_Tensor> Matrix_Second_Tensor; // Acumulo un tensor de 2 orden en una matri
           
            Energy_Yield_Function::Energy_Yield_Function(myState State )
	    :FluencyCriteria()
	    {
               mState = State; 
	    }


             Energy_Yield_Function::~Energy_Yield_Function() {}

    
//***********************************************************************
//***********************************************************************
		  // Energy_Criteria
		  // Diferent limits in traccion and compresion
      
		    void  Energy_Yield_Function::InitializeMaterial(const Properties& props) {
		      mprops = &props;
                      mSigma_y = (*mprops)[FC]; 
		      mSigma_o =  mSigma_y; 
		      }
		     

                     /// WARNING = Revisar por un posible bag
		    void  Energy_Yield_Function::CalculateEquivalentUniaxialStress(
		    const Vector& StressVector, const Vector& StrainVector, double& Result)

			    {
			      
			    int    iter      = 50;
			    double zero      = 1.0E-9;
			    double teta_a    = 0.00;
			    double teta_b    = 0.00;
			    double teta      = 0.00;
                            double n         = (*mprops)[FC]/(*mprops)[FT];
			    Matrix StressTensor    = ZeroMatrix(3,3);
                            Matrix EigenVectors    = ZeroMatrix(3,3);
			    Vector PrincipalStress = ZeroVector(3);
		            this->State_Tensor(StressVector,StressTensor);
		            this->Comprobate_State_Tensor(StressTensor, StressVector); 
                            SD_MathUtils<double>::EigenVectors(StressTensor, EigenVectors,PrincipalStress, zero, iter);			       
			    teta_a =  Tensor_Utils<double>::Mc_aully(PrincipalStress);
			    teta_b =  norm_1(PrincipalStress);  
			    
			    // Condicion para no tener divison de 0/0 
			    if (teta_a==0.00 && teta_b==0.00)
			    {teta = 1.00;}
			    else
			    {teta = teta_a/teta_b;}
			      
                            Result              = inner_prod(StressVector,StrainVector); 
			    Result              = sqrt(Result);
			    Result              = (teta + (1.00-teta)/n)*Result;
			    mSigma_e = Result;
                            Result  -= mSigma_y; 
			    mElasticDomain = Result;
			    }


		    void  Energy_Yield_Function::CalculateEquivalentUniaxialStressViaInvariants( 
		    const Vector& StressVector,double& Result)
			{
			  
		          double iter  = 50;
                          double zero  = 1E-10;  
		          double u_a   = 0.00;
                          double u_b   = 0.00;
                          double u_c   = 0.00;      
                          double r     = 0.00;

			  Matrix StressTensor     =  ZeroMatrix(3,3);
			  Vector PrincipalStress  =  ZeroVector(3);
                          mP_Stress               =  ZeroVector(3);
                          Matrix EigenVectors     =  ZeroMatrix(3,3);
		          double n                = (*mprops)[FC]/(*mprops)[FT];

			  
			  this->State_Tensor(StressVector,StressTensor);
                          this->Comprobate_State_Tensor(StressTensor, StressVector);
		          
			  {SD_MathUtils<double>::EigenVectors(StressTensor, EigenVectors,PrincipalStress, zero, iter);}
                          if(norm_2(PrincipalStress)<1E-5){PrincipalStress = ZeroVector(3);}
                          mu_a = u_a = norm_2(PrincipalStress);
                          mu_b = u_b = sum(PrincipalStress);
			  mu_c = u_c = norm_1(PrincipalStress); 

                          noalias(mP_Stress) = PrincipalStress;     

                          mr = r = 0.5 + 0.5*(u_b/u_c);
			  if (u_b== 0.00 && u_c==0.00)
			     {mr = r = 1.00; }

			  //Result = (r +(1.00-r)/n)*u_a/sqrt((*mprops)[YOUNG_MODULUS]);
			  Result = (1.00 + r*(n-1.00))*u_a; //sqrt((*mprops)[YOUNG_MODULUS]);
                          mSigma_e = Result;
                          Result -= mSigma_y;
                          mElasticDomain = Result;
			}


                    void Energy_Yield_Function::GetValue(const Variable<double>& rVariable, double& Result)
                    {
		      if(rVariable==YIELD_STRESS)
			Result = mSigma_y;
		      else if(rVariable==YIELD_SURFACE)
			Result =  mSigma_e; 
		      return;
		    }

                    void Energy_Yield_Function::FinalizeSolutionStep()
                    {
                         mSigma_y = mSigma_e;
		    }
		    
		    void Energy_Yield_Function::CalculateDerivateFluencyCriteria(const Vector& StressVector,Vector& DerivateFluencyCriteria)
			{
			  
                          Second_Order_Tensor a;  
                          Second_Order_Tensor C;
                          Vector Dev_sigma_a = ZeroVector(3);
                          Vector Unit        = ZeroVector(3);
                          DerivateFluencyCriteria = ZeroVector(6);
                          C.resize(3,false);
                          C[0].resize(3,false);
                          C[1].resize(3,false); 
                          C[2].resize(3,false);


			  double tetha_Lode = 0.00;  
                          double n          = (*mprops)[FC]/(*mprops)[FT];
			  Vector I          = ZeroVector(3);
			  Vector J          = ZeroVector(3);
			  Vector J_des      = ZeroVector(3);		      

			  Matrix StressTensor     = ZeroMatrix(3,3);
			  Vector PrincipalStress  = ZeroVector(3);
				  
                          Unit(0) = 1.00;
                          Unit(1) = 1.00;
                          Unit(2) = 1.00;
			  
			  this->State_Tensor(StressVector,StressTensor);

			  this->Comprobate_State_Tensor(StressTensor, StressVector); // funcion definida en clase base;
			  Tensor_Utils<double>::TensorialInvariants(StressTensor, I, J, J_des);
		          
			  //KRATOS_WATCH(mP_Stress)
			  //KRATOS_WATCH(J_des)


			  if (J_des(1)==0.00 && J_des(2)==0.00) 
			    {
				tetha_Lode = PI/6.00;                          
			    }
			  else
			    {  
			    tetha_Lode = -(3.00*sqrt(3.00)*J_des(2))/(2.00*pow(J_des(1), 1.50));
                            //KRATOS_WATCH(J_des)
			    tetha_Lode = asin(tetha_Lode)/3.00;
			    if(fabs(tetha_Lode) > PI/6) {std::cout<<"Warning: Angle must be less than PI/6"<<std::endl;}   
			    }

 
			    double par_a = (1.00 - 1.00/n)*mu_a/sqrt((*mprops)[YOUNG_MODULUS]);
                            double par_b = (mr - (1.00-mr)/n)/sqrt((*mprops)[YOUNG_MODULUS]);

			    if(mu_a==0.00 && mu_b==0.00 && mu_c==0.00)
                             {
                               Dev_sigma_a = ZeroVector(3);
                             }
                            else
                             {
                             noalias(Dev_sigma_a) = (par_b/mu_a)*(mP_Stress);
                             noalias(Dev_sigma_a) += (par_a*((1/(2.00*mu_c))*Unit -(mu_b/(2.00*mu_c*mu_c))*Tensor_Utils<double>::Sign(mP_Stress)));
                            }

                        
                           /*KRATOS_WATCH(mP_Stress)
			   KRATOS_WATCH(mu_a)
                           KRATOS_WATCH(mu_b)
                           KRATOS_WATCH(mu_c)*/
                           //KRATOS_WATCH(Dev_sigma_a)                          
                            

			    Vector Vector_Sin_Teta = ZeroVector(3);
                            Vector_Sin_Teta(0) =  sin(tetha_Lode + 2.00*PI/3.00);
                            Vector_Sin_Teta(1) =  sin(tetha_Lode);
                            Vector_Sin_Teta(2) =  sin(tetha_Lode - 2.00*PI/3.00);  

                            Vector Vector_Cos_Teta = ZeroVector(3);
                            Vector_Cos_Teta(0) =  cos(tetha_Lode + 2.00*PI/3.00);
                            Vector_Cos_Teta(1) =  cos(tetha_Lode);
                            Vector_Cos_Teta(2) =  cos(tetha_Lode - 2.00*PI/3.00);

                           C[0] = Unit/3.00;  

                           // evitando puntos singulares
                           if (fabs(PI/6 -fabs(tetha_Lode)) < 1.00E-2)
                            {
                                  Vector Singlular(3);
				  Singlular(0) = 0.500; 
				  Singlular(1) = 0.500;
				  Singlular(2) = -1.00;

                                  C[1] = (2.00/(sqrt(3.00)))*Singlular;
                                  C[2] = ZeroVector(3);         
                            }
                            else
                            {               
				C[1] = (2.00/(sqrt(3.00)))*(Vector_Sin_Teta - tan(3.00*tetha_Lode)*Vector_Cos_Teta);
				C[2] = -(1.00/(J_des(1)*cos(3.00*tetha_Lode)))*Vector_Cos_Teta; 
			    }
                            this->CalculateVectorFlowDerivate(StressVector, a);

                            //KRATOS_WATCH(a)
                            //KRATOS_WATCH(C)
                            Matrix aux_a = ZeroMatrix(3,6);
                            Matrix aux_b = ZeroMatrix(3,3);
                            
                            for(unsigned int i=0; i<3; i++){
                                for(unsigned int j=0; j<3; j++){
                                   aux_b(i,j)  =  C[i](j);
                              }
                            }

                            for(unsigned int i=0; i<3; i++){
                                for(unsigned int j=0; j<6; j++){
                                   aux_a(i,j)  =  a[i](j);
                              }
                            }


                            DerivateFluencyCriteria =  prod(Dev_sigma_a, Matrix(prod(aux_b,aux_a)));
                           // KRATOS_WATCH(DerivateFluencyCriteria)

			}

  
  
}



