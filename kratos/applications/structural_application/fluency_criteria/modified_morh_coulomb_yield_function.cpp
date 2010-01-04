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
//#include "custom_utilities/tensor_utils.h"
#include "fluency_criteria/modified_morh_coulomb_yield_function.h"
//#include <cmath>



namespace Kratos
  {

  
	    typedef boost::numeric::ublas::vector<Vector> Second_Order_Tensor; // dos opciones: un tensor de segundo orden y/o un vector que almacena un vector		  
	    typedef boost::numeric::ublas::vector<Second_Order_Tensor> Third_Order_Tensor;
			  
            typedef boost::numeric::ublas::vector<boost::numeric::ublas::vector<Matrix> > Fourth_Order_Tensor;
			  
	    typedef matrix<Second_Order_Tensor> Matrix_Second_Tensor; // Acumulo un tensor de 2 orden en una matri
           
            Modified_Morh_Coulomb_Yield_Function::Modified_Morh_Coulomb_Yield_Function(myState State, myPotencialPlastic PotencialPlastic )
	    :FluencyCriteria()
	    {
              mState            = State;
              mPotencialPlastic = PotencialPlastic;
	    }

             Modified_Morh_Coulomb_Yield_Function::~Modified_Morh_Coulomb_Yield_Function() {}

    
//***********************************************************************
//***********************************************************************

		    void Modified_Morh_Coulomb_Yield_Function::InitializeMaterial(const Properties& props) 
                   {mprops = &props;
                   double maxfriction_angle  = (*mprops)[MAX_FRICTION_INTERNAL_ANGLE]*PI/180.00;
		   mEta = 2.00*tan(maxfriction_angle/2.00 + PI/4.00); 
                   }
		     

		    void Modified_Morh_Coulomb_Yield_Function:: CalculateEquivalentUniaxialStress(
		    const Vector& StressVector,double& Result)  {}


		    void Modified_Morh_Coulomb_Yield_Function::CalculateEquivalentUniaxialStressViaPrincipalStress(
		    const Vector& StressVector,double& Result)
                       {
                       }



		    void Modified_Morh_Coulomb_Yield_Function::CalculateEquivalentUniaxialStressViaInvariants(
		    const Vector& StressVector,double& Result)
			   {
		
		      // Nota: la resistencia de comparacion de esta superficie de fluencia es FC.
                      double tetha_Lode          = 0.00;
		      double frictional_internal = 0.00;
                      double cohesion            = 0.00;
                      double R_morh              = 0.00;
                      double Alfa                = 0.00;
		      double K_one               = 0.00;
                      double K_two               = 0.00;
                      double K_three             = 0.00; 
		      double sin_fric_inter      = 0.00;
		      double cos_fric_inter      = 0.00;
		      double sin_tetha_lode      = 0.00;
		      double cos_tetha_lode      = 0.00;
		      
                      Vector I        = ZeroVector(3);
		      Vector J        = ZeroVector(3);
                      Vector J_des    = ZeroVector(3);		      

		      Matrix StressTensor     = ZeroMatrix(3,3);
		      Vector PrincipalStress  = ZeroVector(3);
			      
		      
		      this->State_Tensor(StressVector,StressTensor);
		      this->Comprobate_State_Tensor(StressTensor, StressVector); // funcion definida en clase base;
		      Tensor_Utils<double>::TensorialInvariants(StressTensor, I, J, J_des);
                      //KRATOS_WATCH(I)
                      //KRATOS_WATCH(J_des)
		
		      if (J_des[1]==0.00) 
                        {
			    tetha_Lode = 0.00;                              
			}
                      else
		        {  
			tetha_Lode = -(3.00*sqrt(3.00)*J_des(2))/(2.00*pow(J_des[1], 1.50));
                        if(fabs(tetha_Lode) > 1.00){tetha_Lode = 1.00; }
			tetha_Lode = asin(tetha_Lode)/3.00; 
		        }

                        
			frictional_internal = mfriction_angle; // (*mprops)[FRICTION_INTERNAL_ANGLE]*PI/180.00;
                        cohesion = mcohesion;
			
                       
                        R_morh = tan(frictional_internal/2.00 + PI/4.00);
                        R_morh = R_morh*R_morh;
			
			Alfa = ((*mprops)[FC]/(*mprops)[FT])/R_morh;
			sin_fric_inter = sin(frictional_internal);
                        cos_fric_inter = cos(frictional_internal);
                        sin_tetha_lode = sin(tetha_Lode);
                        cos_tetha_lode = cos(tetha_Lode);

			K_one   = 0.50*((1.00 + Alfa)-(1.00-Alfa)*sin_fric_inter);
                        K_two   = 1E15;
                        if(sin_fric_inter!=0.00)  
			      {K_two   = 0.50*((1.00 + Alfa)-(1.00-Alfa)/sin_fric_inter);}
                        K_three = 0.50*((1.00 + Alfa)*sin_fric_inter-(1.00-Alfa));
	  
                        //KRATOS_WATCH(K_one)
                        //KRATOS_WATCH(K_two)
                        //KRATOS_WATCH(K_three)

			Result   = I(0)*K_three/3.00 + sqrt(J_des(1))*(K_one*cos_tetha_lode - K_two*sin_tetha_lode*sin_fric_inter/sqrt(3.00));
                        Result   = fabs(Result);
		        mSigma_e = Result;
                        
			mSigma_y = fabs(cohesion*cos_fric_inter); //(*mprops)[COHESION]*cos(frictional_internal);
                        //KRATOS_WATCH(cos_fric_inter)
                        Result   -=mSigma_y;  
                           }


		    void Modified_Morh_Coulomb_Yield_Function::CalculateEquivalentUniaxialStressViaCilindricalCoordinate(
		    const Vector& StressVector,double& Result){}

                    // solo respecto de las tensiones
                    void Modified_Morh_Coulomb_Yield_Function::CalculateDerivateFluencyCriteria(const Vector& StressVector, Vector& DerivateFluencyCriteria)
		    {

			  KRATOS_TRY  
		      	  
			  double tetha_Lode      = 0.00;
			  double R_morh          = 0.00;
			  double Alfa            = 0.00;
			  double K_one           = 0.00;
			  double K_two           = 0.00;
			  double K_three         = 0.00; 
			  double sin_fric_inter  = 0.00;
                          double cos_fric_inter  = 0.00;
                          double sin_tetha_lode  = 0.00;
                          double cos_tetha_lode  = 0.00;
                          double frictional_internal = 0.00;

                          Second_Order_Tensor a;  
			  Matrix StressTensor     = ZeroMatrix(3,3);
			  Vector PrincipalStress  = ZeroVector(3);

                          array_1d<double,3> C; C = ZeroVector(3);  
                          Vector A(6); A.resize(6,false); A = ZeroVector(6);    
                          DerivateFluencyCriteria.resize(6, false); 
			  DerivateFluencyCriteria = ZeroVector(6);
			  a.resize(3, false);
			  a[0].resize(6, false); a[0] = ZeroVector(6);
                          a[1].resize(6, false); a[1] = ZeroVector(6);
                          a[2].resize(6, false); a[2] = ZeroVector(6);

                         
                          
                          frictional_internal = mfriction_angle;                                     
			    
			  Vector I          = ZeroVector(3);
			  Vector J          = ZeroVector(3);
			  Vector J_des      = ZeroVector(3);      				  
			  
		      this->State_Tensor(StressVector,StressTensor);
		      this->Comprobate_State_Tensor(StressTensor, StressVector); // funcion definida en clase base;
		      Tensor_Utils<double>::TensorialInvariants(StressTensor, I, J, J_des);
                      //KRATOS_WATCH(StressVector)
                      //KRATOS_WATCH(J_des)
		
		      if (J_des(1)==0.00) 
                        {
			    tetha_Lode = 0.00;                              
			}
                      else
		        {  
			tetha_Lode = -(3.00*sqrt(3.00)*J_des(2))/(2.00*pow(J_des(1), 1.50));
                        if(fabs(tetha_Lode) > 1.00){tetha_Lode = 1.00; }
			tetha_Lode = asin(tetha_Lode)/3.00; 
		        }
                        
			//frictional_internal = mfriction_angle; // (*mprops)[FRICTION_INTERNAL_ANGLE]*PI/180.00;
			
                       
                        R_morh = tan(frictional_internal/2.00 + PI/4.00);
                        R_morh = R_morh*R_morh;
			
			Alfa = ((*mprops)[FC]/(*mprops)[FT])/R_morh;
			sin_fric_inter = sin(frictional_internal);
                        cos_fric_inter = cos(frictional_internal);
                        sin_tetha_lode = sin(tetha_Lode);
                        cos_tetha_lode = cos(tetha_Lode);

			K_one   = 0.50*((1.00 + Alfa)-(1.00-Alfa)*sin_fric_inter);
                        K_two   = 1E15;
                        if(sin_fric_inter!=0.00)  
			      {K_two   = 0.50*((1.00 + Alfa)-(1.00-Alfa)/sin_fric_inter);}
                        K_three = 0.50*((1.00 + Alfa)*sin_fric_inter-(1.00-Alfa));
                             
                         if(fabs((fabs(tetha_Lode) - 0.523599)) < 0.0001)  // Angulo de  +-30.00
	                    {
                              C[0] =  K_three/3.00;  
                              C[1] =  0.50*(sqrt(3.00)*K_one - K_two*sin_fric_inter/sqrt(3.00));
                              if(tetha_Lode<0.00) {C(1) =  0.50*(sqrt(3.00)*K_one + K_two*sin_fric_inter/sqrt(3.00)); } 
                              C[2]  =  0.00;
                            }
                            else
                           {
                            C[0] =  K_three/3.00; 
                            C[1] = (K_one*cos_tetha_lode - K_two*sin_tetha_lode*sin_fric_inter/sqrt(3.00));  
                            C[1] = C[1] + tan(3.00*tetha_Lode)*(K_one*sin_tetha_lode + K_two*cos_tetha_lode*sin_fric_inter/sqrt(3.00));  
                            C[2] = (sqrt(3.00)*K_one*sin_tetha_lode + K_two*cos_tetha_lode*sin_fric_inter)/(2.00*cos(3.00*tetha_Lode)*J_des[1]); 
                           } 

			    this->CalculateVectorFlowDerivate(StressVector, a);
			    for(unsigned  int i=0; i<3; i++)
                               {
                                 noalias(A) += a[i]*C(i); 
                               }  
                           noalias(DerivateFluencyCriteria) = A;
                           //KRATOS_WATCH(DerivateFluencyCriteria)
                           KRATOS_CATCH("")
                         }

                   // solo respecto de las tensiones
                   void Modified_Morh_Coulomb_Yield_Function::CalculateDerivatePotencialFlowCriteria(const Vector& StressVector, Vector& DerivateFluencyCriteria)
                    {
			  KRATOS_TRY  
		      	  
			  double tetha_Lode      = 0.00;
			  double R_morh          = 0.00;
			  double Alfa            = 0.00;
			  double K_one           = 0.00;
			  double K_two           = 0.00;
			  double K_three         = 0.00; 
			  double sin_fric_inter  = 0.00;
                          double cos_fric_inter  = 0.00;
                          double sin_tetha_lode  = 0.00;
                          double cos_tetha_lode  = 0.00;
                          double frictional_internal = 0.00;

                          Second_Order_Tensor a;  
			  Matrix StressTensor     = ZeroMatrix(3,3);
			  Vector PrincipalStress  = ZeroVector(3);

                          array_1d<double,3> C; C = ZeroVector(3);  
                          Vector A(6); A.resize(6,false); A = ZeroVector(6);    
                          DerivateFluencyCriteria.resize(6, false); 
			  DerivateFluencyCriteria = ZeroVector(6);
			  a.resize(3, false);
			  a[0].resize(6, false); a[0] = ZeroVector(6);
                          a[1].resize(6, false); a[1] = ZeroVector(6);
                          a[2].resize(6, false); a[2] = ZeroVector(6);

                         
                          frictional_internal = mdilatancy_angle;
                                    
			    
			  Vector I          = ZeroVector(3);
			  Vector J          = ZeroVector(3);
			  Vector J_des      = ZeroVector(3);      				  
			  
		      this->State_Tensor(StressVector,StressTensor);
		      this->Comprobate_State_Tensor(StressTensor, StressVector); // funcion definida en clase base;
		      Tensor_Utils<double>::TensorialInvariants(StressTensor, I, J, J_des);
                      //KRATOS_WATCH(I)
                      //KRATOS_WATCH(J_des)
		
		      if (J_des[1]==0.00) 
                        {
			    tetha_Lode = 0.00;                              
			}
                      else
		        {  
			tetha_Lode = -(3.00*sqrt(3.00)*J_des(2))/(2.00*pow(J_des(1), 1.50));
                        if(fabs(tetha_Lode) > 1.00){tetha_Lode = 1.00; }
			tetha_Lode = asin(tetha_Lode)/3.00; 
		        }
                      
                       
                        R_morh = tan(frictional_internal/2.00 + PI/4.00);
                        R_morh = R_morh*R_morh;
			
			Alfa = ((*mprops)[FC]/(*mprops)[FT])/R_morh;
			sin_fric_inter = sin(frictional_internal);
                        cos_fric_inter = cos(frictional_internal);
                        sin_tetha_lode = sin(tetha_Lode);
                        cos_tetha_lode = cos(tetha_Lode);

			K_one   = 0.50*((1.00 + Alfa)-(1.00-Alfa)*sin_fric_inter);
                        K_two   = 1E15;
                        if(sin_fric_inter!=0.00)  
			      {K_two   = 0.50*((1.00 + Alfa)-(1.00-Alfa)/sin_fric_inter);}
                        K_three = 0.50*((1.00 + Alfa)*sin_fric_inter-(1.00-Alfa));
                             
                         if(fabs((fabs(tetha_Lode) - 0.523599)) < 0.0001)  // Angulo de  +-30.00
	                    {
                              C[0] =  K_three/3.00;  
                              C[1] =  0.50*(sqrt(3.00)*K_one - K_two*sin_fric_inter/sqrt(3.00));
                              if(tetha_Lode<0.00) {C(1) =  0.50*(sqrt(3.00)*K_one + K_two*sin_fric_inter/sqrt(3.00)); } 
                              C[2]  =  0.00;
                            }
                            else
                           {
                            C[0] =  K_three/3.00; 
                            C[1] = (K_one*cos_tetha_lode - K_two*sin_tetha_lode*sin_fric_inter/sqrt(3.00));  
                            C[1] = C[1] + tan(3.00*tetha_Lode)*(K_one*sin_tetha_lode + K_two*cos_tetha_lode*sin_fric_inter/sqrt(3.00));  
                            C[2] = (sqrt(3.00)*K_one*sin_tetha_lode + K_two*cos_tetha_lode*sin_fric_inter)/(2.00*cos(3.00*tetha_Lode)*J_des[1]); 
                           } 

			    this->CalculateVectorFlowDerivate(StressVector, a);
			    for(unsigned  int i=0; i<3; i++)
                               {
                                 noalias(A) += a[i]*C(i); 
                               }  
                           noalias(DerivateFluencyCriteria) = A;
                           KRATOS_CATCH("")

                    }


                  void Modified_Morh_Coulomb_Yield_Function::UpdateVariables(const Vector& Variables)
                     {
                          mcohesion        = Variables(0);
                          mfriction_angle  = Variables(1);
                          mdilatancy_angle = Variables(2); 

                     }

    }


