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
#include "fluency_criteria/isotropic_rankine_yield_function.h"
#include <cmath>



namespace Kratos
  {

  
	    typedef boost::numeric::ublas::vector<Vector> Second_Order_Tensor; // dos opciones: un tensor de segundo orden y/o un vector que almacena un vector		  
	    typedef boost::numeric::ublas::vector<Second_Order_Tensor> Third_Order_Tensor;
			  
            typedef boost::numeric::ublas::vector<boost::numeric::ublas::vector<Matrix> > Fourth_Order_Tensor;
			  
	    typedef matrix<Second_Order_Tensor> Matrix_Second_Tensor; // Acumulo un tensor de 2 orden en una matri
           
            Isotropic_Rankine_Yield_Function::Isotropic_Rankine_Yield_Function(myState State)
	    :FluencyCriteria()
	    {
              mState = State;
	    }

             Isotropic_Rankine_Yield_Function::~Isotropic_Rankine_Yield_Function() {}

    
//***********************************************************************
//***********************************************************************

	  void Isotropic_Rankine_Yield_Function::InitializeMaterial(const Properties& props) 
	  {   
           
           mprops      = &props;
	   mFt         =  (*mprops)[FT];
	   mcurrent_Ft = mFt;

           minitialize = false; 
           maccumulated_plastic_strain_current = 0.00;  
           maccumulated_plastic_strain_old     = 0.00;  
	  }

          
          void Isotropic_Rankine_Yield_Function::CalculateEquivalentUniaxialStressViaPrincipalStress(
	  const Vector& StressVector,double& Result){}


	  void Isotropic_Rankine_Yield_Function:: CalculateEquivalentUniaxialStress(
	  const Vector& StressVector,double& Result)
	  {

 	  array_1d<double,3> Trial_Stress_Vector = ZeroVector(3);

          CalculatePrincipalStressVector(StressVector, Trial_Stress_Vector);  

          mMultisurface_Platicity_Sigma       = ZeroVector(3);
          mMultisurface_Platicity_Yield       = ZeroVector(3); 
	  ///* Multisurface Representation 
	  mMultisurface_Platicity_Sigma[0]    =   Trial_Stress_Vector[0]; 
	  mMultisurface_Platicity_Yield[0]    =   mMultisurface_Platicity_Sigma[0] - mcurrent_Ft;

	  mMultisurface_Platicity_Sigma[1]    =   Trial_Stress_Vector[1]; 
	  mMultisurface_Platicity_Yield[1]    =   mMultisurface_Platicity_Sigma[1] - mcurrent_Ft;

	  mMultisurface_Platicity_Sigma[2]    =   Trial_Stress_Vector[2]; 
	  mMultisurface_Platicity_Yield[2]    =   mMultisurface_Platicity_Sigma[2] - mcurrent_Ft;
 
          Result = (*max_element(mMultisurface_Platicity_Yield.begin(), mMultisurface_Platicity_Yield.end()));  
           
	  }

	  void Isotropic_Rankine_Yield_Function::UpdateVariables(const Vector& Variables)
	  {
	   mFt = Variables[0];
           mcurrent_Ft = mFt;  
           maccumulated_plastic_strain_current = maccumulated_plastic_strain_old; 
 
           if(minitialize == false)
               { 
                 const double& ft     =  (*mprops)[FT];
                 const double& gt     =  (*mprops)[FRACTURE_ENERGY]; 
                 const double& length = Variables[3];      
	         mH            = length * ft * ft / ( 2.00 * gt); 
                 minitialize = true;
               }
	  }



	  void Isotropic_Rankine_Yield_Function::CalculateEquivalentUniaxialStressViaInvariants(
	  const Vector& StressVector,double& Result)
	  {

	  }


	  void Isotropic_Rankine_Yield_Function::CalculateEquivalentUniaxialStressViaCilindricalCoordinate(
	  const Vector& StressVector,double& Result){}



	  void Isotropic_Rankine_Yield_Function::CalculateDerivateFluencyCriteria(const Vector& StressVector, Vector& DerivateFluencyCriteria)
	  {
	  }

          void Isotropic_Rankine_Yield_Function::ReturnMapping(const Vector& StressVector, 
            const Vector& StrainVector,   
            Vector& delta_lamda,
            array_1d<double,3>& Result)
            {               

               
                array_1d<double,3> Trial_Stress_Vector     = ZeroVector(3); 
                array_1d<double,3> Aux_Trial_Stress_Vector = ZeroVector(3); 
                CalculatePrincipalStressVector(StressVector, Trial_Stress_Vector);   
                noalias(Aux_Trial_Stress_Vector) = Trial_Stress_Vector;

                ///* return to main plane
                One_Vector_Return_Mapping_To_Main_Plane(StressVector, delta_lamda,  Trial_Stress_Vector);
                  
                ///*check validty
                bool check = false;   
                check = CheckValidity(Trial_Stress_Vector);
                if( check==false)
                 {
                     ///*return to corner
                     noalias(Trial_Stress_Vector) =  Aux_Trial_Stress_Vector;
                     Two_Vector_Return_Mapping_To_Corner(StressVector,delta_lamda,  Trial_Stress_Vector);  
                     check = CheckValidity(Trial_Stress_Vector);
                     if (check==false)
                          {
                             ///return to apex
                              noalias(Trial_Stress_Vector) =  Aux_Trial_Stress_Vector;
                              Three_Vector_Return_Mapping_To_Apex (StressVector, delta_lamda, Trial_Stress_Vector);   
                          }    

                      
                  }

                   Result[0] =  Trial_Stress_Vector[0];
                   Result[1] =  Trial_Stress_Vector[1];
                   Result[2] =  Trial_Stress_Vector[2];
            }

         void Isotropic_Rankine_Yield_Function::GetValue(double Result)
                 {
                    Result = mcurrent_Ft; 
                 } 

        void Isotropic_Rankine_Yield_Function::GetValue(Vector& Result)
                 {
                    Result.resize(3, false);   
                    Result[0] = mcurrent_Ft;
                    Result[1] = mcurrent_Ft;
                    Result[2] = mcurrent_Ft;  
                 }  

         void Isotropic_Rankine_Yield_Function::Finalize()
              {
                   mFt = mcurrent_Ft;  
                   maccumulated_plastic_strain_old     =  maccumulated_plastic_strain_current;     
                   maccumulated_plastic_strain_current = 0.00;    
              }

void Isotropic_Rankine_Yield_Function::One_Vector_Return_Mapping_To_Main_Plane(const Vector& StressVector, Vector& delta_lamda,  array_1d<double,3>& Result)
              {

              unsigned int iter    = 0;    
              double norma         = 1.00;   
	      double delta_lamda_a = 0.00; 
	      double E             = (*mprops)[YOUNG_MODULUS];
	      double NU            = (*mprops)[POISSON_RATIO];             
              double G             = 0.5*E / (1.00 + NU);
              double K             =  E / (3.00 * (1.00-2.00*NU));
              double H             =  mH;
	      double d             = 0.00;  
	      double residual      = 0.00;
              const double toler   = 1E-6; 
     

              mcurrent_Ft = mFt;
              delta_lamda = ZeroVector(1); 
              CalculateEquivalentUniaxialStress(StressVector, norma);   
 
              while(iter++<=100 && norma>= toler) 
	          {  
                    if(iter>=100){KRATOS_WATCH("WARNING = DO NOT CONVERGENCE FOR ONE ACTIVE SURFACE RANKINE" )}
                    d = 4.00 * G /3.00 + K - H;
                    delta_lamda[0]  += (mMultisurface_Platicity_Yield[0]) / d;; 
                           
                    ///* Updatinf mFt
                    maccumulated_plastic_strain_current +=  delta_lamda[0];
                    mcurrent_Ft = mFt - H * delta_lamda[0]; 
                    delta_lamda_a    = delta_lamda[0];   

                    ///* comprobando si mft se cumplio   
		    if(mcurrent_Ft <= toler) 
                       {
		           mcurrent_Ft = 0.00;
                           break;         
                       } 
                    else
                    { 
                    ///* update teh current value                       
                    CalculateEquivalentUniaxialStress(StressVector, norma); 
                    residual    =  norma - delta_lamda_a * (4.00  * G / 3.00 + K ) ;
                    norma       =  fabs(residual);   
                    
                    }
                  }   
                     ///* Updating Stress  
                    //KRATOS_WATCH(Trial_Stress_Vector) 
		    if(mcurrent_Ft <=toler) 
                    {
		    Result[0] = 0.00;
		    Result[1] = 0.00; 
		    Result[2] = 0.00;  
		    } 
                    else
                    { 
	            Result[0] = Result[0] - delta_lamda_a*(4.00  * G / 3.00 + K );    
		    Result[1] = Result[1] - delta_lamda_a*(-2.00 * G / 3.00 + K ); 
		    Result[2] = Result[2] - delta_lamda_a*(-2.00 * G / 3.00 + K );    
                    }      
                                   
 
              }


    void Isotropic_Rankine_Yield_Function::Two_Vector_Return_Mapping_To_Corner (const Vector& StressVector, Vector& delta_lamda ,array_1d<double,3>& Result)
       {
              unsigned int iter    = 0;    
              int singular         = 0;  
              double norma         = 1.00;   
	      double delta_lamda_a = 0.00; 
              double delta_lamda_b = 0.00; 
              const  double toler  = 1E-6;
	      double E             = (*mprops)[YOUNG_MODULUS];
	      double NU            = (*mprops)[POISSON_RATIO];             
              double G             = 0.5*E / (1.00 + NU);
              double K             =  E / (3.00 * (1.00-2.00*NU));
              double H             =  mH;
	      Matrix d             = ZeroMatrix(2,2);  
              Matrix d_inv         = ZeroMatrix(2,2); 
	      Vector residual      = ZeroVector(2);
              delta_lamda          = ZeroVector(2); 
            
              mcurrent_Ft = mFt;
              CalculateEquivalentUniaxialStress(StressVector, norma);  
              
                
	      residual[0] = mMultisurface_Platicity_Yield[0];
              residual[1] = mMultisurface_Platicity_Yield[1];
              //KRATOS_WATCH(residual)  
 
              d.resize(2,2);
	      d_inv.resize(2,2);
	      d(0,0) = -4.00 * G / 3.00 - K + H;  d(0,1)   =  2.00 * G / 3.00 - K  + H;  
	      d(1,0) = 2.00 * G / 3.00 - K  + H;  d(1,1)   = -4.00 * G / 3.00 - K  + H; 
              singular =  SD_MathUtils<double>::InvertMatrix(d, d_inv);    
               

              while(iter++<=100 && norma>= toler) 
		  {
                      if(iter>=100){KRATOS_WATCH("WARNING = DO NOT CONVERGENCE FOR TWO ACTIVE SURFACE RANKINE" )}
		      noalias(delta_lamda) =  delta_lamda - Vector(prod(d_inv, residual)); 
		
                      //if(delta_lamda[0] < 0.00) {delta_lamda[0] = 0.00; } //std::cout<<"WARNING = GAMMA NEGATIVE FOR TWO ACTIVE SURFACE "<<std::endl; }  
		      //if(delta_lamda[1] < 0.00) {delta_lamda[1] = 0.00; } // std::cout<<"WARNING = GAMMA NEGATIVE FOR TWO ACTIVE SURFACE "<<std::endl; }     
                      //KRATOS_WATCH(residual)
                      //KRATOS_WATCH(delta_lamda) Catholic
                        
		      delta_lamda_a = delta_lamda[0];
		      delta_lamda_b = delta_lamda[1];
                      
                       /// Updatinf mFt
                      mcurrent_Ft =  mFt - H * (delta_lamda_a + delta_lamda_b);  
                      
                       /// ft se anulan
                      if(mcurrent_Ft<=toler)
                        {
                            mcurrent_Ft = 0.00;
                            delta_lamda_a  = delta_lamda[0];  
                            delta_lamda_b  = delta_lamda[1];    
                            break;
                        }
                       else 
                       {
                                            
                          CalculateEquivalentUniaxialStress(StressVector, norma);        
                  
                          residual[0] = mMultisurface_Platicity_Yield[0];
                          residual[1] = mMultisurface_Platicity_Yield[1];
                          //KRATOS_WATCH(mMultisurface_Platicity_Yield) 
                       
                          residual[0]=  residual[0] - delta_lamda_a*( 4.00 * G / 3.00 + K ) - delta_lamda_b*( -2.00  * G / 3.00 + K );   
                          residual[1]=  residual[1] - delta_lamda_a*(-2.00 * G / 3.00 + K ) - delta_lamda_b*(  4.00  * G / 3.00 + K );   
                          norma      =  norm_2(residual); 
                          //KRATOS_WATCH(residual)  
                           //KRATOS_WATCH(norma)
                        }
                      }

		     ///* Updating Stress  
                     //if( delta_lamda[0]< 0.00 ) { delta_lamda[0]=0.00; }
                     //if( delta_lamda[1]< 0.00 ) { delta_lamda[1]=0.00; }  

                     if(mcurrent_Ft <= 0.00) 
                        {
                           Result[0] = 0.00;
                           Result[1] = 0.00;
                           Result[2] = 0.00;
                        }
                     else 
                      {
                             Result[0] = Result[0] - delta_lamda_a*( 4.00  * G / 3.00 + K ) - delta_lamda_b*( -2.00  * G / 3.00 + K );
                             Result[1] = Result[1] - delta_lamda_a*(-2.00 * G  / 3.00 + K ) - delta_lamda_b*( 4.00  * G / 3.00 + K ); 
                             Result[2] = Result[2] - (delta_lamda_a + delta_lamda_b) * (-2.00 * G / 3.00 + K );  
                       }

              }


      

     void Isotropic_Rankine_Yield_Function::Three_Vector_Return_Mapping_To_Apex (const Vector& StressVector, Vector& delta_lamda ,array_1d<double,3>& Result)
      {
          
              unsigned int iter    = 0;    
              int singular         = 0;  
              double norma         = 1.00;   
	      double delta_lamda_a = 0.00; 
              double delta_lamda_b = 0.00; 
              double delta_lamda_c = 0.00; 
              const  double toler  = 1E-6;
	      double E             = (*mprops)[YOUNG_MODULUS];
	      double NU            = (*mprops)[POISSON_RATIO];             
              double G             = 0.5*E / (1.00 + NU);
              double K             =  E / (3.00 * (1.00-2.00*NU));
              double H             =  mH;
	      Matrix d             = ZeroMatrix(3,3);  
              Matrix d_inv         = ZeroMatrix(3,3); 
	      Vector residual      = ZeroVector(3);
              delta_lamda          = ZeroVector(3); 
            
              mcurrent_Ft = mFt;
              CalculateEquivalentUniaxialStress(StressVector, norma); 
              
                
	      residual[0] = mMultisurface_Platicity_Yield[0];
              residual[1] = mMultisurface_Platicity_Yield[1];
              residual[2] = mMultisurface_Platicity_Yield[2]; 
              //KRATOS_WATCH(residual) 
 
              d.resize(3,3);
	      d_inv.resize(3,3);
	      d(0,0) = -4.00 * G / 3.00 - K + H;  d(0,1)  =  2.00 * G / 3.00 - K + H;  d(0,2)   =   2.00 * G / 3.00 - K  + H; 
	      d(1,0) =  2.00 * G / 3.00 - K + H;  d(1,1)  = -4.00 * G / 3.00 - K + H;  d(1,2)   =   2.00 * G / 3.00 - K  + H; 
              d(2,0) =  2.00 * G / 3.00 - K + H;  d(2,1)  =  2.00 * G / 3.00 - K + H;  d(2,2)   =  -4.00 * G / 3.00 - K  + H ;  
              //KRATOS_WATCH(d)     
              singular =  SD_MathUtils<double>::InvertMatrix(d, d_inv);    
              //KRATOS_WATCH(d_inv)   

              while(iter++<=100 && norma>= toler) 
		  {
                      if(iter>=100){KRATOS_WATCH("WARNING = DO NOT CONVERGENCE FOR TWO ACTIVE SURFACE RANKINE" )}
		      noalias(delta_lamda) =  delta_lamda - Vector(prod(d_inv, residual)); 
		
                      //if(delta_lamda[0] < 0.00) {delta_lamda[0] = 0.00; } //std::cout<<"WARNING = GAMMA NEGATIVE FOR TWO ACTIVE SURFACE "<<std::endl; }  
		      //if(delta_lamda[1] < 0.00) {delta_lamda[1] = 0.00; } // std::cout<<"WARNING = GAMMA NEGATIVE FOR TWO ACTIVE SURFACE "<<std::endl; }     
                      //KRATOS_WATCH(residual)
                      //KRATOS_WATCH(delta_lamda) Catholic
                        
		      delta_lamda_a = delta_lamda[0];
		      delta_lamda_b = delta_lamda[1];
                      delta_lamda_c = delta_lamda[2];
                      
                       /// Updatinf mFt
                      mcurrent_Ft =  mFt - H * (delta_lamda_a + delta_lamda_b + delta_lamda_c);  
                      
                       /// ft se anulan
                      if(mcurrent_Ft<=toler)
                        {
                            mcurrent_Ft = 0.00;
                            delta_lamda_a  = delta_lamda[0];  
                            delta_lamda_b  = delta_lamda[1];    
                            delta_lamda_c  = delta_lamda[2];                             
                            break;
                        }
                       else 
                       {
                                            
                          CalculateEquivalentUniaxialStress(StressVector, norma);        
                  
                          residual[0] = mMultisurface_Platicity_Yield[0];
                          residual[1] = mMultisurface_Platicity_Yield[1];
                          residual[2] = mMultisurface_Platicity_Yield[2];
                          //KRATOS_WATCH(mMultisurface_Platicity_Yield) 
                       
         residual[0] =  residual[0] - delta_lamda_a*( 4.00 * G / 3.00 + K ) - delta_lamda_b*( -2.00 * G / 3.00 + K ) - delta_lamda_c*( -2.00 * G / 3.00 + K );   
         residual[1] =  residual[1] - delta_lamda_a*(-2.00 * G / 3.00 + K ) - delta_lamda_b*( 4.00  * G / 3.00 + K ) - delta_lamda_c*( -2.00 * G / 3.00 + K );   
         residual[2] =  residual[2] - delta_lamda_a*(-2.00 * G / 3.00 + K ) - delta_lamda_b*( -2.00 * G / 3.00 + K ) - delta_lamda_c*( 4.00  * G / 3.00 + K );         
         norma      =  norm_2(residual); 
         //KRATOS_WATCH(residual)  
                           //KRATOS_WATCH(norma)
                        }
                      }

		     ///* Updating Stress  
                     //if( delta_lamda[0]< 0.00 ) { delta_lamda[0]=0.00; }
                     //if( delta_lamda[1]< 0.00 ) { delta_lamda[1]=0.00; }  

                     if(mcurrent_Ft <= 0.00) 
                        {
                           Result[0] = 0.00;
                           Result[1] = 0.00;
                           Result[2] = 0.00;
                        }
                     else 
                      {
          Result[0] = Result[0] - delta_lamda_a*( 4.00  * G / 3.00 + K ) - delta_lamda_b*( -2.00  * G / 3.00 + K ) - delta_lamda_c*( -2.00  * G / 3.00 + K );
          Result[1] = Result[1] - delta_lamda_a*(-2.00 * G  / 3.00 + K ) - delta_lamda_b*( 4.00  * G / 3.00 + K )  - delta_lamda_c*( -2.00  * G / 3.00 + K );
          Result[2] = Result[2] - delta_lamda_a*(-2.00 * G  / 3.00 + K ) - delta_lamda_b*( -2.00  * G / 3.00 + K )  - delta_lamda_c*( 4.00  * G / 3.00 + K ); 
                       }


      }


     bool Isotropic_Rankine_Yield_Function::CheckValidity( array_1d<double,3>& Sigma)
           {
                 bool check   = false;
                 array_1d<double,3> Aux_Sigma;
                 Aux_Sigma[0] = fabs(Sigma[0]);
                 Aux_Sigma[1] = fabs(Sigma[1]);
                 Aux_Sigma[2] = fabs(Sigma[2]);
                 double delta = (*max_element(Aux_Sigma.begin(), Aux_Sigma.end())) * 1.00E-6; 
                 if( (Sigma[0] + delta) >= Sigma[1] && (Sigma[1] + delta ) >= Sigma[2]){ check = true;}
                 return check;
           }  




}


