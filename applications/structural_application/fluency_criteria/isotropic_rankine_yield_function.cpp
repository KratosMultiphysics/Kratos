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

            Isotropic_Rankine_Yield_Function::Isotropic_Rankine_Yield_Function() {}       
    
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
           
           mprops      =  &props;
	   mFt         =  (*mprops)[FT];
	   mcurrent_Ft =  mFt;

           minitialize                                 = false; 
           mrankine_accumulated_plastic_strain_current = 0.00;  
           mrankine_accumulated_plastic_strain_old     = 0.00;
	   mplastic_strain                             = ZeroVector(4);
	   mplastic_strain_old                         = ZeroVector(4);
	   mPrincipalPlasticStrain_current             = ZeroVector(3);
	   mPrincipalPlasticStrain_old                 = ZeroVector(3);

	  }

	  void Isotropic_Rankine_Yield_Function::UpdateMaterial()
	  {
	   
           mcurrent_Ft                              = mFt;  
           maccumulated_plastic_strain_current      = maccumulated_plastic_strain_old; 
	   noalias(mPrincipalPlasticStrain_current) = mPrincipalPlasticStrain_old;
	   noalias(mplastic_strain)                 = mplastic_strain_old; 
	   
           if(minitialize == false)
               { 
                 //const double& ft     =  (*mprops)[FT];
                 //const double& gt     =  (*mprops)[FRACTURE_ENERGY]; 
                 //const double& length =  0.00;
	         //mH                   =  0.00; //length * ft * ft / ( 2.00 * gt); 
                 minitialize = true;
               }
	  }

   
          bool Isotropic_Rankine_Yield_Function::CheckPlasticAdmisibility(const Vector& Stress)
          { 
	     array_1d<double, 3> check;
	     double toler = 1E-8;
	     check[0] = Stress[0] - mcurrent_Ft;
	     check[1] = Stress[1] - mcurrent_Ft;
	     check[2] = Stress[2] - mcurrent_Ft;
	     
	     if(check[0]>= toler || check[1]>= toler ||  check[2]>= toler)
	        return true;
	     else 
	        return false;
	  }

          void Isotropic_Rankine_Yield_Function::ReturnMapping(const Vector& StrainVector, Vector& StressVector)
            {               

		array_1d<double,3> PrincipalStress;
		array_1d<array_1d < double,3 > ,3> EigenVectors;
		array_1d<unsigned int,3> Order;
		array_1d<double,3> Sigma;
		
		SpectralDecomposition(StressVector, PrincipalStress, EigenVectors);
		IdentifyMaximunAndMinumumPrincipalStres_CalculateOrder(PrincipalStress, Order);
		

                if(CheckPlasticAdmisibility(PrincipalStress))
		{	  
                // return to main plane
		 Vector delta_lamda; 
                 bool check = One_Vector_Return_Mapping_To_Main_Plane(PrincipalStress, delta_lamda,  Sigma);
                   
                 if(check==false)
                  {
                     //return to corner
                     check =  Two_Vector_Return_Mapping_To_Corner(PrincipalStress, delta_lamda, Sigma);  
                     if (check==false)
                          {
                            //return to apex
                              Three_Vector_Return_Mapping_To_Apex (PrincipalStress, delta_lamda, Sigma);   
                          }    
                       
                      
                  }

                 AssembleUpdateStressAndStrainTensor(Sigma,  EigenVectors,  Order, StrainVector, StressVector);
		} 
		
	       // Updating the plastic strain tensor
                CalculatePlasticStrain(StrainVector,StressVector);
            
	    }
	    

void Isotropic_Rankine_Yield_Function::FinalizeSolutionStep()
    {
	  mFt                                     =  mcurrent_Ft;  
	  mrankine_accumulated_plastic_strain_old =  mrankine_accumulated_plastic_strain_current;     
	  noalias(mPrincipalPlasticStrain_old)    =  mPrincipalPlasticStrain_current; 
	  noalias(mplastic_strain_old)            =  mplastic_strain; 
    }


bool Isotropic_Rankine_Yield_Function::One_Vector_Return_Mapping_To_Main_Plane(const array_1d<double,3>& PrincipalStress, Vector& delta_lamda,  array_1d<double,3>& Sigma)
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
	      norma =  PrincipalStress[0] - mcurrent_Ft;
              norma =  norma/mcurrent_Ft;
	      
              while(iter++<=100 && norma>= toler) 
	          {  
                    d = 4.00 * G /3.00 + K - H;
                    delta_lamda[0]  +=  (PrincipalStress[0] - mcurrent_Ft) / d;; 
                           
                    ///normal acuulated
                    ///maccumulated_plastic_strain_current = maccumulated_plastic_strain_old + delta_lamda[0];
		    
		    /// von mises
		    maccumulated_plastic_strain_current = maccumulated_plastic_strain_old + 0.81649658092773 * delta_lamda[0];    
		    
		    ///* Updatinf mFt
                    mcurrent_Ft   =  mFt - H * (maccumulated_plastic_strain_current - maccumulated_plastic_strain_old);
                    delta_lamda_a  = delta_lamda[0];   
		    
                    ///* comprobando si mft se cumplio   
		    if(mcurrent_Ft <= toler) 
                       {
		           mcurrent_Ft = 0.00;
                           break;         
                       } 
                    else
                    { 
                    ///* update teh current value                       
                    norma       =  PrincipalStress[0] - mcurrent_Ft; 
                    residual    =  norma - delta_lamda_a * (4.00  * G / 3.00 + K ) ;
                    residual    =  residual / norma;
                    norma       =  fabs(residual);   
                    
                    }
                  }   
                     ///* Updating Stress  
		    if(mcurrent_Ft <=toler) 
                    {
		    Sigma[0] = 0.00;
		    Sigma[1] = 0.00; 
		    Sigma[2] = 0.00;  
		    } 
                    else
                    { 
	            Sigma[0] = PrincipalStress[0] - delta_lamda_a*(4.00  * G / 3.00 + K );    
		    Sigma[1] = PrincipalStress[1] - delta_lamda_a*(-2.00 * G / 3.00 + K ); 
		    Sigma[2] = PrincipalStress[2] - delta_lamda_a*(-2.00 * G / 3.00 + K );    
                    }  
                    
                    
                    bool check = CheckValidity(Sigma);
		    if(check==true) 
		     {
		    //updating the correct principal pastic strain
		    mPrincipalPlasticStrain_current[0] = mPrincipalPlasticStrain_old[0] +  delta_lamda_a;
		    mPrincipalPlasticStrain_current[1] = mPrincipalPlasticStrain_old[1]; 
		    mPrincipalPlasticStrain_current[2] = mPrincipalPlasticStrain_old[2];
		    return true;
		   }
		else
		  return false;
		
		}


    
    bool Isotropic_Rankine_Yield_Function::Two_Vector_Return_Mapping_To_Corner (const array_1d<double,3>& PrincipalStress, Vector& delta_lamda ,array_1d<double,3>& Sigma)
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
                
	      residual[0] = PrincipalStress[0] - mcurrent_Ft;
              residual[1] = PrincipalStress[1] - mcurrent_Ft;

 
              d.resize(2,2);
	      d_inv.resize(2,2);
	      d(0,0) = -4.00 * G / 3.00 - K + H;  d(0,1)   =  2.00 * G / 3.00 - K  + H;  
	      d(1,0) = 2.00 * G / 3.00 - K  + H;  d(1,1)   = -4.00 * G / 3.00 - K  + H; 
              singular =  SD_MathUtils<double>::InvertMatrix(d, d_inv);    
               
              while(iter++<=100 && norma>= toler) 
		  {
		      delta_lamda =  delta_lamda - Vector(prod(d_inv, residual)); 
		      delta_lamda_a = delta_lamda[0];
		      delta_lamda_b = delta_lamda[1];
                      
		      // normal
		      //maccumulated_plastic_strain_current = maccumulated_plastic_strain_old + delta_lamda_a + delta_lamda_b;
		      
		      // von mises
		      maccumulated_plastic_strain_current = maccumulated_plastic_strain_old + 0.81649658092773 * norm_2(delta_lamda);    
		      
                       // Updatinf mFt
                      mcurrent_Ft =  mFt - H * (maccumulated_plastic_strain_current - maccumulated_plastic_strain_old); 
                      
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
                                                 
                  
                 	  residual[0] = PrincipalStress[0] - mcurrent_Ft;
                          residual[1] = PrincipalStress[1] - mcurrent_Ft;
                       
                          residual[0]=  residual[0] - delta_lamda_a*( 4.00 * G / 3.00 + K ) - delta_lamda_b*( -2.00  * G / 3.00 + K );   
                          residual[1]=  residual[1] - delta_lamda_a*(-2.00 * G / 3.00 + K ) - delta_lamda_b*(  4.00  * G / 3.00 + K );   
                          norma      =  norm_2(residual); 

                        }
                      }



                     if(mcurrent_Ft <= 0.00) 
                        {
                           Sigma[0] = 0.00;
                           Sigma[1] = 0.00;
                           Sigma[2] = 0.00;
                        }
                     else 
                      {
                             Sigma[0] = PrincipalStress[0] - delta_lamda_a*( 4.00  * G / 3.00 + K ) - delta_lamda_b*( -2.00  * G / 3.00 + K );
                             Sigma[1] = PrincipalStress[1] - delta_lamda_a*(-2.00 * G  / 3.00 + K ) - delta_lamda_b*( 4.00  * G / 3.00 + K ); 
                             Sigma[2] = PrincipalStress[2] - (delta_lamda_a + delta_lamda_b) * (-2.00 * G / 3.00 + K );  
                       }
                       
                      bool check = CheckValidity(Sigma);
		      if(check==true) 
		      {
			//updating the correct principal pastic strain
			mPrincipalPlasticStrain_current[0] = mPrincipalPlasticStrain_old[0] + delta_lamda_a;
			mPrincipalPlasticStrain_current[1] = mPrincipalPlasticStrain_old[1] + delta_lamda_b; 
			mPrincipalPlasticStrain_current[2] = mPrincipalPlasticStrain_old[2];
			return true;
		      }
		      else
		        return false;
		              

              }

     
     void Isotropic_Rankine_Yield_Function::Three_Vector_Return_Mapping_To_Apex(const array_1d<double,3>& PrincipalStress, Vector& delta_lamda ,array_1d<double,3>& Sigma)
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
              
                
	      residual[0] = PrincipalStress[0] - mcurrent_Ft;
              residual[1] = PrincipalStress[1] - mcurrent_Ft;
	      residual[2] = PrincipalStress[2] - mcurrent_Ft;
	      
	      
              d.resize(3,3);
	      d_inv.resize(3,3);
	      d(0,0) = -4.00 * G / 3.00 - K + H;  d(0,1)  =  2.00 * G / 3.00 - K + H;  d(0,2)   =   2.00 * G / 3.00 - K  + H; 
	      d(1,0) =  2.00 * G / 3.00 - K + H;  d(1,1)  = -4.00 * G / 3.00 - K + H;  d(1,2)   =   2.00 * G / 3.00 - K  + H; 
              d(2,0) =  2.00 * G / 3.00 - K + H;  d(2,1)  =  2.00 * G / 3.00 - K + H;  d(2,2)   =  -4.00 * G / 3.00 - K  + H ;  
  
              singular =  SD_MathUtils<double>::InvertMatrix(d, d_inv);    


              while(iter++<=100 && norma>= toler) 
		  {
                      if(iter>=100){ KRATOS_ERROR(std::logic_error,  "RETURN MAPPING TO APEX  NOT CONVERGED IN RANKINE" , "")};
		      noalias(delta_lamda) =  delta_lamda - Vector(prod(d_inv, residual)); 
		
                        
		      delta_lamda_a = delta_lamda[0];
		      delta_lamda_b = delta_lamda[1];
                      delta_lamda_c = delta_lamda[2];
                      
		      /// normal
		      ///maccumulated_plastic_strain_current = maccumulated_plastic_strain_old + delta_lamda_a + delta_lamda_b +  delta_lamda_c;
		      
		      ///von mises
		      maccumulated_plastic_strain_current = maccumulated_plastic_strain_old + 0.81649658092773 * norm_2(delta_lamda); 
		      
		      
                       /// Updatinf mFt
                      mcurrent_Ft =  mFt - H * (maccumulated_plastic_strain_current - maccumulated_plastic_strain_old);  
                      
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
                                                
                  
	            residual[0] = PrincipalStress[0] - mcurrent_Ft;
                    residual[1] = PrincipalStress[1] - mcurrent_Ft;
	            residual[2] = PrincipalStress[2] - mcurrent_Ft;
                      
                    residual[0] =  residual[0] - delta_lamda_a*( 4.00 * G / 3.00 + K ) - delta_lamda_b*( -2.00 * G / 3.00 + K ) - delta_lamda_c*( -2.00 * G / 3.00 + K );   
                    residual[1] =  residual[1] - delta_lamda_a*(-2.00 * G / 3.00 + K ) - delta_lamda_b*( 4.00  * G / 3.00 + K ) - delta_lamda_c*( -2.00 * G / 3.00 + K );   
                    residual[2] =  residual[2] - delta_lamda_a*(-2.00 * G / 3.00 + K ) - delta_lamda_b*( -2.00 * G / 3.00 + K ) - delta_lamda_c*( 4.00  * G / 3.00 + K );         
                    norma      =  norm_2(residual); 

                        }
                      }


                     if(mcurrent_Ft <= 0.00) 
                        {
                           Sigma[0] = 0.00;
                           Sigma[1] = 0.00;
                           Sigma[2] = 0.00;
                        }
                     else 
                      {
		      Sigma[0] = PrincipalStress[0] - delta_lamda_a*( 4.00  * G / 3.00 + K ) - delta_lamda_b*( -2.00  * G / 3.00 + K ) - delta_lamda_c*( -2.00  * G / 3.00 + K );
		      Sigma[1] = PrincipalStress[1] - delta_lamda_a*(-2.00 * G  / 3.00 + K ) - delta_lamda_b*( 4.00  * G / 3.00 + K )  - delta_lamda_c*( -2.00  * G / 3.00 + K );
                      Sigma[2] = PrincipalStress[2] - delta_lamda_a*(-2.00 * G  / 3.00 + K ) - delta_lamda_b*( -2.00  * G / 3.00 + K )  - delta_lamda_c*( 4.00  * G / 3.00 + K ); 
                       }
                       
                     mPrincipalPlasticStrain_current[0] = mPrincipalPlasticStrain_old[0] +  delta_lamda_a;
	             mPrincipalPlasticStrain_current[1] = mPrincipalPlasticStrain_old[0] +  delta_lamda_b; 
		     mPrincipalPlasticStrain_current[2] = mPrincipalPlasticStrain_old[0] +  delta_lamda_c; 
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

void Isotropic_Rankine_Yield_Function::GetValue(const Variable<Matrix>& rVariable, Matrix& Result)
 {
      unsigned int size = 0; 
      
      if(this->mState==Plane_Strain)
         size = 4;
      else
         size = 6;
      
   if (rVariable==GREEN_LAGRANGE_PLASTIC_STRAIN_TENSOR)
   {
     Result.resize(1, size);
     for(unsigned int i = 0; i< size; i++ )
        Result(0,i) = mplastic_strain(i); 
   }
   
    return;
   }
   


}


