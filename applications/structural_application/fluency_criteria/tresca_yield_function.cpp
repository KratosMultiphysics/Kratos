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
#include "fluency_criteria/tresca_yield_function.h"
#include <cmath>



namespace Kratos
  {

  
	    typedef boost::numeric::ublas::vector<Vector> Second_Order_Tensor; // dos opciones: un tensor de segundo orden y/o un vector que almacena un vector		  
	    typedef boost::numeric::ublas::vector<Second_Order_Tensor> Third_Order_Tensor;
			  
            typedef boost::numeric::ublas::vector<boost::numeric::ublas::vector<Matrix> > Fourth_Order_Tensor;
			  
	    typedef matrix<Second_Order_Tensor> Matrix_Second_Tensor; // Acumulo un tensor de 2 orden en una matri
           
            Tresca_Yield_Function::Tresca_Yield_Function(myState State )
	    :FluencyCriteria()
	    {
              mState = State;
	    }

            Tresca_Yield_Function::~Tresca_Yield_Function() {}

    
//***********************************************************************
//***********************************************************************


      	void Tresca_Yield_Function::InitializeMaterial(const Properties& props)
	{
               mprops                 =  &props;
               msigma_y               = (*mprops)[YIELD_STRESS];
               miso_hardening_modulus = (*mprops)[ISOTROPIC_HARDENING_MODULUS];
               maccumulated_plastic_strain_current = 0.00;  
               maccumulated_plastic_strain_old     = 0.00;  
	}
  
 

      void Tresca_Yield_Function::UpdateVariables( const Vector& Variables)
         {
               mcurrent_sigma_y = msigma_y;   
         }

      void Tresca_Yield_Function::GetValue(Vector& Result)
	{}	    
	             
      void Tresca_Yield_Function::Finalize() 
          {
               
                msigma_y  = mcurrent_sigma_y;  
          }
		     

      void Tresca_Yield_Function::CalculateEquivalentUniaxialStress(
            const Vector& StressVector,double& Result)
            {

 	  array_1d<double,3> Trial_Stress_Vector = ZeroVector(3);
          CalculatePrincipalStressVector(StressVector, Trial_Stress_Vector);  
       
          mMultisurface_Platicity_Sigma       = ZeroVector(3);
          mMultisurface_Platicity_Yield       = ZeroVector(3);
 
	  ///* Multisurface Representation 
	  mMultisurface_Platicity_Sigma[0]    =   Trial_Stress_Vector[0] - Trial_Stress_Vector[2]; 
	  mMultisurface_Platicity_Yield[0]    =   mMultisurface_Platicity_Sigma[0] - mcurrent_sigma_y;

          mMultisurface_Platicity_Sigma[1]    =   Trial_Stress_Vector[1] - Trial_Stress_Vector[2]; 
	  mMultisurface_Platicity_Yield[1]    =   mMultisurface_Platicity_Sigma[1] - mcurrent_sigma_y; 


          mMultisurface_Platicity_Sigma[2]    =   Trial_Stress_Vector[0] - Trial_Stress_Vector[1]; 
	  mMultisurface_Platicity_Yield[2]    =   mMultisurface_Platicity_Sigma[2] - mcurrent_sigma_y; 
 

           Result = (*max_element(mMultisurface_Platicity_Yield.begin(), mMultisurface_Platicity_Yield.end()));                       

           }  

          void Tresca_Yield_Function::ReturnMapping(const Vector& StressVector,
                         const Vector& StrainVector,
                         Vector& delta_lamda,
                         array_1d<double,3>& Result)


           {


                double p_trial = 0.00;  
                Matrix Sigma_Tensor        = ZeroMatrix(3,3);
                State_Tensor(StressVector, Sigma_Tensor);
                p_trial = (Sigma_Tensor(0,0) + Sigma_Tensor(1,1) + Sigma_Tensor(2,2))/3.00 ;    
                array_1d<double,3> Trial_Stress_Vector     = ZeroVector(3); 
                array_1d<double,3> Aux_Trial_Stress_Vector = ZeroVector(3);   
                CalculatePrincipalStressVector(StressVector, Trial_Stress_Vector);  
                Trial_Stress_Vector[0] = Trial_Stress_Vector[0] - p_trial;
                Trial_Stress_Vector[1] = Trial_Stress_Vector[1] - p_trial;
                Trial_Stress_Vector[2] = Trial_Stress_Vector[2] - p_trial; 
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
                     double proof = Trial_Stress_Vector[0] + Trial_Stress_Vector[2] - 2.00 * Trial_Stress_Vector[1];
                     if(proof > 0.00 ) {mCases = right;} else {mCases = left;}
                     Two_Vector_Return_Mapping_To_Corner(StressVector,delta_lamda,  Trial_Stress_Vector);                        
                  }

                   Result[0] =  Trial_Stress_Vector[0] +  p_trial;    
                   Result[1] =  Trial_Stress_Vector[1] +  p_trial;
                   Result[2] =  Trial_Stress_Vector[2] +  p_trial;
 
           }   

      
        void Tresca_Yield_Function::One_Vector_Return_Mapping_To_Main_Plane
           (const Vector& StressVector, 
            Vector& delta_lamda ,array_1d<double,3>& Result) 
           {
                double E             = (*mprops)[YOUNG_MODULUS];
	        double NU            = (*mprops)[POISSON_RATIO];             
                double G             = 0.5*E / (1.00 + NU);
                double H             =  miso_hardening_modulus;           
		unsigned int iter = 0;     
                double  toler    = 1E-6; 
		double norma      = 1.00;  
                double d          = 0.00;   
                double residual   = 0.00;            
		delta_lamda       = ZeroVector(1);
		mcurrent_sigma_y  = msigma_y;
                

                while(iter++<=100 && norma>= toler) 
		  {

                      maccumulated_plastic_strain_current = maccumulated_plastic_strain_old; 
                      residual = this->mMultisurface_Platicity_Yield[0];
                      d = -4.00 * G - H;
                      delta_lamda[0] = delta_lamda[0] - residual/d;
                      maccumulated_plastic_strain_current+= delta_lamda[0]; 

                      ///* update teh current values
                      mcurrent_sigma_y+= H*delta_lamda[0]; //KRATOS_WATCH(mcurrent_sigma_y); 
                      //mcurrent_sigma_y= (*mprops)[YIELD_STRESS] + H * maccumulated_plastic_strain_current; //KRATOS_WATCH(mcurrent_sigma_y);
  
                     CalculateEquivalentUniaxialStress(StressVector, norma); 
                     residual =  norma - 4.00  * G * delta_lamda[0];
                     norma    = fabs(residual);                          
                  }

                 Result[0] = Result[0] - 2.00  * G * delta_lamda[0];
                 Result[2] = Result[2] + 2.00  * G * delta_lamda[0];
                    
         }
 

        void Tresca_Yield_Function::Two_Vector_Return_Mapping_To_Corner
           (const Vector& StressVector, 
            Vector& delta_lamda ,array_1d<double,3>& Result) 
           {
                double E          = (*mprops)[YOUNG_MODULUS];
	        double NU         = (*mprops)[POISSON_RATIO];             
                double G          = 0.5*E / (1.00 + NU);
                double H          =  miso_hardening_modulus;          
                int singular      = 0; 
		unsigned int iter = 0;     
                double  toler     = 1E-6; 
		double norma      = 1.00;   
                delta_lamda       = ZeroVector(2);
                mcurrent_sigma_y  = msigma_y; 
		Vector residual   = ZeroVector(2);
                Matrix d          = ZeroMatrix(2,2);  
                Matrix d_inv      = ZeroMatrix(2,2);   
               
                switch (mCases)
                      {
                         case right:
                            { 
                                residual[0] = this->mMultisurface_Platicity_Yield[0];   
                                residual[1] = this->mMultisurface_Platicity_Yield[2];   
                                d(0,0) = -4.00 * G - H; d(0,1) = -2.00 * G - H;  
                                d(1,0) = -2.00 * G - H; d(1,1) = -4.00 * G - H; 
                                singular =  SD_MathUtils<double>::InvertMatrix(d, d_inv); 
                                while(iter++<=100 && norma>= toler) 
		                 {
                                    maccumulated_plastic_strain_current = maccumulated_plastic_strain_old;
                                    noalias(delta_lamda) =  delta_lamda - Vector(prod(d_inv, residual));
                                    maccumulated_plastic_strain_current+= delta_lamda[0] + delta_lamda[1];
                                    mcurrent_sigma_y += H * (delta_lamda[0] + delta_lamda[1]);  
                                    CalculateEquivalentUniaxialStress(StressVector, norma); 
                                    residual[0] = (this->mMultisurface_Platicity_Yield[0]) - 2.00 * G * (2.00 * delta_lamda[0] +  delta_lamda[1] );     
                                    residual[1] = (this->mMultisurface_Platicity_Yield[2]) - 2.00 * G * (delta_lamda[0]  + 2.00 * delta_lamda[1] );
                                    norma    = norm_2(residual);    
                                  }
 
                              Result[0] = Result[0] - 2.00  * G * ( delta_lamda[0] + delta_lamda[1]) ;
                              Result[1] = Result[1] + 2.00  * G * delta_lamda[1]; 
                              Result[2] = Result[2] + 2.00  * G * delta_lamda[0];  
                              break;  
                          }

                         case left:
                            { 
                                residual[0] = this->mMultisurface_Platicity_Yield[0];   
                                residual[1] = this->mMultisurface_Platicity_Yield[1];   
                                d(0,0) = -4.00 * G - H; d(0,1) = -2.00 * G - H;  
                                d(1,0) = -2.00 * G - H; d(1,1) = -4.00 * G - H; 
                                singular =  SD_MathUtils<double>::InvertMatrix(d, d_inv); 
                                while(iter++<=100 && norma>= toler) 
		                 {
                                    maccumulated_plastic_strain_current = maccumulated_plastic_strain_old;
                                    noalias(delta_lamda) =  delta_lamda - Vector(prod(d_inv, residual));
                                    maccumulated_plastic_strain_current+= delta_lamda[0] + delta_lamda[1];
                                    mcurrent_sigma_y += H * (delta_lamda[0] + delta_lamda[1]);  
                                    CalculateEquivalentUniaxialStress(StressVector, norma); 
                                    residual[0] = (this->mMultisurface_Platicity_Yield[0]) - 2.00 * G * (2.00 * delta_lamda[0] +  delta_lamda[1] );     
                                    residual[1] = (this->mMultisurface_Platicity_Yield[1]) - 2.00 * G * (delta_lamda[0]  + 2.00 * delta_lamda[1] );
                                    norma    = norm_2(residual);    
                                  }
 
                              Result[0] = Result[0] - 2.00  * G * delta_lamda[0];
                              Result[1] = Result[1] - 2.00  * G * delta_lamda[1]; 
                              Result[2] = Result[2] + 2.00  * G * (delta_lamda[0] + delta_lamda[1]);  
                              break;  
                          } 
                }
           }


     bool Tresca_Yield_Function::CheckValidity( array_1d<double,3>& Sigma)
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



	void Tresca_Yield_Function::CalculateEquivalentUniaxialStressViaPrincipalStress(
 	const Vector& StressVector,double& Result)
	{

		      double crit      = 1E-10;
                      double zero      = 1E-10;  
                      double max       = 0.00;
                      unsigned int dim = 3;
		      
		      Matrix StressTensor     = ZeroMatrix(dim,dim);
                      Vector PrincipalStress  = ZeroVector(dim);
                      Vector Aux_Vector       = ZeroVector(dim);
   
                      this->State_Tensor(StressVector,StressTensor);
		      this->Comprobate_State_Tensor(StressTensor, StressVector); // funcion definida en clase base;
                      PrincipalStress         = SD_MathUtils<double>::EigenValues(StressTensor,crit, zero);
		      
                     Aux_Vector(0) =  fabs(PrincipalStress(0)-PrincipalStress(1));
		     Aux_Vector(1) =  fabs(PrincipalStress(0)-PrincipalStress(2));
		     Aux_Vector(2) =  fabs(PrincipalStress(1)-PrincipalStress(2));
		          
		      max = (*std::max_element(Aux_Vector.begin(),Aux_Vector.end()));
                     
                      Result = 0.50*max; // - msigma_max; 
                      //KRATOS_WATCH(Result)  
           }


	void Tresca_Yield_Function::CalculateEquivalentUniaxialStressViaInvariants(
	const Vector& StressVector,double& Result)

	{
                    
                      unsigned int dim  = 3;
                      double tetha_Lode = 0.00;
                      Vector I          = ZeroVector(3);
		      Vector J          = ZeroVector(3);
                      Vector J_des      = ZeroVector(3);		      

		      Matrix StressTensor     = ZeroMatrix(dim,dim);
		      Vector PrincipalStress  = ZeroVector(dim);

		      this->State_Tensor(StressVector,StressTensor);
		      this->Comprobate_State_Tensor(StressTensor, StressVector); // funcion definida en clase base;
		      Tensor_Utils<double>::TensorialInvariants(StressTensor, I, J, J_des);
		
		      if (J_des(1)==0.00 && J_des(2)==0.00) 
                        {
			    tetha_Lode = PI/2.00;                           	   
			}
                      else
		        {  
			tetha_Lode = (3.00*sqrt(3.00)*J_des(2))/(2.00*pow(J_des(1), 1.50));
			if(tetha_Lode > 1.00){tetha_Lode = 1.00; }
			tetha_Lode = asin(-tetha_Lode)/3.00;
		        }

		      Result = sqrt(J_des(1))*cos(tetha_Lode);// - msigma_max;
                      //KRATOS_WATCH("----------")
                      //KRATOS_WATCH(Result)

	 }




	void Tresca_Yield_Function::CalculateEquivalentUniaxialStressViaCilindricalCoordinate(
	const Vector& StressVector,double& Result) 


	{}

	  void Tresca_Yield_Function::CalculateDerivateFluencyCriteria(const Vector& StressVector, Vector& DerivateFluencyCriteria)
	{
	  		  Second_Order_Tensor a;  
                          Vector C = ZeroVector(3);
                          DerivateFluencyCriteria = ZeroVector(6);

			  double tetha_Lode = 0.00;  
			  Vector I          = ZeroVector(3);
			  Vector J          = ZeroVector(3);
			  Vector J_des      = ZeroVector(3);		      

			  Matrix StressTensor     = ZeroMatrix(3,3);
			  Vector PrincipalStress  = ZeroVector(3);
				  
			  
			  this->State_Tensor(StressVector,StressTensor);
			  this->Comprobate_State_Tensor(StressTensor, StressVector); // funcion definida en clase base;
			  Tensor_Utils<double>::TensorialInvariants(StressTensor, I, J, J_des);
		    
			  if (J_des(1)==0.00 && J_des(2)==0.00) 
			    {
				tetha_Lode = PI/2.00;                           	   
			    }
			  else
			    {  
			    tetha_Lode = -(3.00*sqrt(3.00)*J_des(2))/(2.00*pow(J_des(1), 1.50));
			    if(fabs(tetha_Lode) > 1.00){tetha_Lode = 1.00; }
			    tetha_Lode = asin(tetha_Lode)/3.00; 
			    }

			    if(fabs((fabs(tetha_Lode) - 0.523599)) < 1E-5)  // Angulo de  +-30.00
                            {
                            C(0) = 0.00; 
                            C(1) = sqrt(3.00);  
                            C(2) = 0.00;
                            }
                            else
                            { 
                            C(0) = 0.00; 
                            C(1) = 2.00*cos(tetha_Lode)*(1.00 + tan(tetha_Lode)*tan(3.00*tetha_Lode));  
                            C(2) = sqrt(3.00)*sin(tetha_Lode)/(J_des(1)*cos(3.00*tetha_Lode));
                            }
			    this->CalculateVectorFlowDerivate(StressVector, a);
			    for(unsigned  int i=0; i<3; i++)
                               {
                                 noalias(DerivateFluencyCriteria) = DerivateFluencyCriteria + a(i)*C(i); 
                               }

	}


    }


