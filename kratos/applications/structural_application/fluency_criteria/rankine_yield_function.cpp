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
#include "fluency_criteria/rankine_yield_function.h"
#include <cmath>



namespace Kratos
  {

  
	    typedef boost::numeric::ublas::vector<Vector> Second_Order_Tensor; // dos opciones: un tensor de segundo orden y/o un vector que almacena un vector		  
	    typedef boost::numeric::ublas::vector<Second_Order_Tensor> Third_Order_Tensor;
			  
            typedef boost::numeric::ublas::vector<boost::numeric::ublas::vector<Matrix> > Fourth_Order_Tensor;
			  
	    typedef matrix<Second_Order_Tensor> Matrix_Second_Tensor; // Acumulo un tensor de 2 orden en una matri
           
            Rankine_Yield_Function::Rankine_Yield_Function(myState State)
	    :FluencyCriteria()
	    {
              mState = State;
	    }

             Rankine_Yield_Function::~Rankine_Yield_Function() {}

    
//***********************************************************************
//***********************************************************************

	  void Rankine_Yield_Function::InitializeMaterial(const Properties& props) 
	  {   mprops = &props;
	  mFt[0] = (*mprops)[FT];
	  mFt[1] = (*mprops)[FT];
	  mFt[2] = (*mprops)[FT]; 
          minitialize = false; 
	  }

          
          void Rankine_Yield_Function::CalculateEquivalentUniaxialStressViaPrincipalStress(
	  const Vector& StressVector,double& Result){}


	  void Rankine_Yield_Function:: CalculateEquivalentUniaxialStress(
	  const Vector& StressVector,double& Result)
	  {
// 	  int    iter      = 50;
// 	  double zero      = 1.0E-9;
// 	  Matrix EigenVectors    = ZeroMatrix(3,3);
// 	  Matrix StressTensor    = ZeroMatrix(3,3);
 	  array_1d<double,3> Trial_Stress_Vector = ZeroVector(3);
// 	  this->State_Tensor(StressVector,StressTensor);
// 
// 	  SD_MathUtils<double>::EigenVectors(StressTensor, EigenVectors,Trial_Stress_Vector, zero, iter);
// 
// 	  ///*sigma_1 >  sigma_2 > sigma_3 
// 	  sort (Trial_Stress_Vector.begin(), Trial_Stress_Vector.end()); 
// 	  reverse(Trial_Stress_Vector.begin(),Trial_Stress_Vector.end()); 

          CalculatePrincipalStressVector(StressVector, Trial_Stress_Vector);  

          mMultisurface_Platicity_Sigma       = ZeroVector(3);
          mMultisurface_Platicity_Yield       = ZeroVector(3); 
	  ///* Multisurface Representation 
	  mMultisurface_Platicity_Sigma[0]    =   Trial_Stress_Vector[0]; 
	  mMultisurface_Platicity_Yield[0]    =   mMultisurface_Platicity_Sigma[0] - mFt[0];

	  mMultisurface_Platicity_Sigma[1]    =   Trial_Stress_Vector[1]; 
	  mMultisurface_Platicity_Yield[1]    =   mMultisurface_Platicity_Sigma[1] - mFt[1];

	  mMultisurface_Platicity_Sigma[2]    =   Trial_Stress_Vector[2]; 
	  mMultisurface_Platicity_Yield[2]    =   mMultisurface_Platicity_Sigma[2] - mFt[2];
  
          //KRATOS_WATCH(mRankine_Yield) 
          Result = (*max_element(mMultisurface_Platicity_Yield.begin(), mMultisurface_Platicity_Yield.end()));  
           
	  }

	  void Rankine_Yield_Function::UpdateVariables(const Vector& Variables)
	  {
	   mFt[0] = Variables[0];
           mFt[1] = Variables[1];
           mFt[2] = Variables[2];
 
           if(minitialize == false)
               { 
                 const double& ft     = (*mprops)[FT];
                 const double& gt     =  (*mprops)[FRACTURE_ENERGY]; 
                 const double& length = Variables[3];      
	         mH            = length * ft * ft / ( 2.00 * gt); 
                 minitialize = true;
               }
	  }



	  void Rankine_Yield_Function::CalculateEquivalentUniaxialStressViaInvariants(
	  const Vector& StressVector,double& Result)
	  {

	  }


	  void Rankine_Yield_Function::CalculateEquivalentUniaxialStressViaCilindricalCoordinate(
	  const Vector& StressVector,double& Result){}



	  void Rankine_Yield_Function::CalculateDerivateFluencyCriteria(const Vector& StressVector, Vector& DerivateFluencyCriteria)
	  {
	  }

          void Rankine_Yield_Function::ReturnMapping(const Vector& StressVector, 
            Vector& delta_lamda,
            array_1d<double,3>& Result)
            {              

              array_1d<double ,3> Trial_Stress_Vector; 
              CalculatePrincipalStressVector(StressVector, Trial_Stress_Vector);  
              double toler = 1E-6;  
              mcurrent_Ft = ZeroVector(3);            
              std::vector<int> active_surface;     
              active_surface.resize(0, false); 
	      active_surface.reserve(5); 

	      if(mMultisurface_Platicity_Yield[0] > toler ) {active_surface.push_back(0); }
	      if(mMultisurface_Platicity_Yield[1] > toler ) {active_surface.push_back(1); }      
	      if(mMultisurface_Platicity_Yield[2] > toler ) {active_surface.push_back(2); } 
                                            
              unsigned int iter    = 0;    
              double norma         = 1.00;   
	      double delta_gamma_a = 0.00; 
	      double delta_gamma_b = 0.00; 
	      double delta_gamma_c = 0.00;   
              
	      double E             = (*mprops)[YOUNG_MODULUS];
	      double NU            = (*mprops)[POISSON_RATIO];             
              double G             = 0.5*E / (1.00 + NU);
              double K             =  E / (3.00 * (1.00-2.00*NU));
              double H             =  mH;


	      Matrix d;  
	      Matrix d_inv;
	      Vector delta_gamma;
	      Vector residual;     
              noalias(mcurrent_Ft) = mFt;


	      ///* WARNING = Si el hablandamiento es no lineal usar newton Rapshon.  
	      ///* Una superficie activa
	      if(active_surface.size()==1)
	      {
                 iter  = 0;     
                 norma = 1.00;                
                 delta_gamma     = ZeroVector(1); 
                 residual        = ZeroVector(1); 
                                        
                 while(iter++<=100 && norma>= toler) 
	          {  
                    if(iter>=100){KRATOS_WATCH("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" )}
                    delta_gamma[0]  += (mMultisurface_Platicity_Yield[0]) / (4.00 * G /3.00 + K - H ); 
                    if(delta_gamma[0] < 0.00) {delta_gamma[0] = 0.00; }   
                    ///* Updatinf mFt
                    mcurrent_Ft[0] = mFt[0] - H * delta_gamma[0];   
                    ///* comprobando si mft se cumplio   
		    if(mcurrent_Ft[0] <= 0.00) {
		    mcurrent_Ft[0] = 0.00; } 


                    ///* update teh current values

                    UpdateVariables(mcurrent_Ft);  


                    delta_gamma_a   = delta_gamma[0];                        
                    CalculateEquivalentUniaxialStress(StressVector, norma); 
                    residual[0] =  norma - delta_gamma_a * (4.00  * G / 3.00 + K ) ;
                    norma    = fabs(residual[0]);   

                  }   
                     ///* Updating Stress  
		    if(mcurrent_Ft[0] == 0.0) 
                    {
		    Result[0] = 1E-14;
		    Result[1] = Trial_Stress_Vector[1] - delta_gamma_a*(-2.00 * G / 3.00 + K ); 
		    Result[2] = Trial_Stress_Vector[2] - delta_gamma_a*(-2.00 * G / 3.00 + K );  
		    } 
                    else
                    { 
	            Result[0] = Trial_Stress_Vector[0] - delta_gamma_a*(4.00  * G / 3.00 + K );    
		    Result[1] = Trial_Stress_Vector[1] - delta_gamma_a*(-2.00 * G / 3.00 + K ); 
		    Result[2] = Trial_Stress_Vector[2] - delta_gamma_a*(-2.00 * G / 3.00 + K );    
                    }      
              } 
              
	      ///*dos superficies activas    
	      if(active_surface.size()==2)
	      {
              int singular         =  0;    
	      delta_gamma = ZeroVector(2);
	      residual    = ZeroVector(2);   
              iter  = 0; 
              norma = 1.00;        
              residual[0] = mMultisurface_Platicity_Yield[0];
              residual[1] = mMultisurface_Platicity_Yield[1]; 
                        
              while(iter++<=100 && norma>= toler) 
		  {
                      if(iter>=100){KRATOS_WATCH("BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB" )} 
		      d.resize(2,2);
		      d_inv.resize(2,2);
		      d(0,0) = -( 4.00 * G / 3.00 + K ) + H;  d(0,1) = -(-2.00 * G / 3.00 + K );  
		      d(1,0) = -(-2.00 * G / 3.00 + K );      d(1,1) = -( 4.00 * G / 3.00 + K ) + H;  

		      singular             =  SD_MathUtils<double>::InvertMatrix(d, d_inv);
		      noalias(delta_gamma) =  delta_gamma - Vector(prod(d_inv, residual)); 
		      if(delta_gamma[0] < 0.00) {delta_gamma[0] = 0.00; }  
		      if(delta_gamma[1] < 0.00) {delta_gamma[1] = 0.00; }     

		      delta_gamma_a = delta_gamma[0];
		      delta_gamma_b = delta_gamma[1];
                     
                       ///* Updatinf mFt
                       mcurrent_Ft[0] = mFt[0] - H * delta_gamma[0];   
                       mcurrent_Ft[1] = mFt[1] - H * delta_gamma[1]; 
                       
                       ///* comprobando si mft se cumplio   
		       if(mcurrent_Ft[0] <= 0.00) {mcurrent_Ft[0] = 0.00; }
                       if(mcurrent_Ft[1] <= 0.00) {mcurrent_Ft[1] = 0.00; } 
                        

                       UpdateVariables(mcurrent_Ft);                      
                       CalculateEquivalentUniaxialStress(StressVector, norma); 
                               
                  
                       residual[0] = mMultisurface_Platicity_Yield[0];
                       residual[1] = mMultisurface_Platicity_Yield[1];
                       
                       residual[0]=  residual[0] - delta_gamma_a*( 4.00  * G / 3.00 + K ) - delta_gamma_b*( -2.00  * G / 3.00 + K );   
                       residual[1]=  residual[1] - delta_gamma_a*(-2.00 * G  / 3.00 + K ) - delta_gamma_b*( 4.00  * G / 3.00 + K );   
                       norma = norm_2(residual);

                       }

		      ///* Updating Stress
                     ///* Updating Stress  
                     Result[0] = Trial_Stress_Vector[0] - delta_gamma_a*( 4.00  * G / 3.00 + K ) - delta_gamma_b*( -2.00  * G / 3.00 + K );  
                     if(mcurrent_Ft[0] <= 0.00) {Result[0] = 1E-14;} 
                     Result[1] = Trial_Stress_Vector[1] - delta_gamma_a*(-2.00 * G  / 3.00 + K ) - delta_gamma_b*( 4.00  * G / 3.00 + K );   
		     if( mcurrent_Ft[1]<= 0.0)  {Result[0] = 1E-14;}  
		     Result[2] = Trial_Stress_Vector[2] - (delta_gamma_a + delta_gamma_b) * (-2.00 * G / 3.00 + K );    

                    
	      } 

	      ///* Muy poco probable en 2D
	      if(active_surface.size()==3)
	      {
              int singular         =  0;    
	      delta_gamma = ZeroVector(2);
	      residual    = ZeroVector(2);   
              iter  = 0; 
              norma = 1.00;  
	      delta_gamma = ZeroVector(3);
	      residual    = ZeroVector(3);   
	      residual[0] = mMultisurface_Platicity_Yield[0];
	      residual[1] = mMultisurface_Platicity_Yield[1];   
	      residual[2] = mMultisurface_Platicity_Yield[2];  
              while(iter++<=100 && norma>= toler) 
		  { 
                    if(iter>=100){KRATOS_WATCH("CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC" )} 
		    d.resize(3,3);
		    d_inv.resize(3,3);
		    d(0,0) = -( 4.00 * G / 3.00 + K ) + H;    d(0,1) = -(-2.00 * G / 3.00 + K );         d(0,2) = -(-2.00 * G / 3.00 + K ); 
		    d(1,0) = -(-2.00 * G / 3.00 + K );        d(1,1) = -( 4.00 * G / 3.00 + K ) + H;     d(1,2) = -(-2.00 * G / 3.00 + K );   
		    d(2,0) = -(-2.00 * G / 3.00 + K );        d(2,1) = -(-2.00 * G / 3.00 + K );         d(2,2) = -( 4.00 * G / 3.00 + K ) + H;   

		    singular             =  SD_MathUtils<double>::InvertMatrix(d, d_inv);
		    noalias(delta_gamma) =  delta_gamma - Vector(prod(d_inv, residual)); 

		    if(delta_gamma[0] < 0.00) {delta_gamma[0] = 0.00; }  
		    if(delta_gamma[1] < 0.00) {delta_gamma[1] = 0.00; }      
		    if(delta_gamma[2] < 0.00) {delta_gamma[2] = 0.00; }  
			

		    delta_gamma_a = delta_gamma[0];
		    delta_gamma_b = delta_gamma[1];
		    delta_gamma_c = delta_gamma[2];

		    ///* Updatinf mFt
		    mcurrent_Ft[0] = mFt[0] - H * delta_gamma[0];   
		    mcurrent_Ft[1] = mFt[1] - H * delta_gamma[1]; 
                    mcurrent_Ft[2] = mFt[2] - H * delta_gamma[2]; 
		    
		    ///* comprobando si mft se cumplio   
		    if(mcurrent_Ft[0] <= 0.00) {mcurrent_Ft[0] = 0.00; }
		    if(mcurrent_Ft[1] <= 0.00) {mcurrent_Ft[1] = 0.00; } 
		    if(mcurrent_Ft[2] <= 0.00) {mcurrent_Ft[2] = 0.00; } 
		    

		    UpdateVariables(mcurrent_Ft);                      
		    CalculateEquivalentUniaxialStress(StressVector, norma); 
                    residual[0] = mMultisurface_Platicity_Yield[0];
                    residual[1] = mMultisurface_Platicity_Yield[1];
                    residual[2] = mMultisurface_Platicity_Yield[1];
                       
		    residual[0] = residual[0] - delta_gamma_a*( 4.00  * G / 3.00 + K ) - (delta_gamma_b + delta_gamma_c) *( -2.00  * G / 3.00 + K );    
		    residual[1] = residual[1] - (delta_gamma_a + delta_gamma_c) * (-2.00 * G  / 3.00 + K ) - delta_gamma_b*( 4.00  * G / 3.00 + K ); 
		    residual[2] = residual[2] - (delta_gamma_a + delta_gamma_b) * (-2.00 * G / 3.00 + K )  - delta_gamma_c*( 4.00  * G / 3.00 + K );                    
                    norma = norm_2(residual);
                    }


		    ///* Updating Stress 
		    Result[0] = Trial_Stress_Vector[0] - delta_gamma_a*( 4.00  * G / 3.00 + K ) - (delta_gamma_b + delta_gamma_c) *( -2.00  * G / 3.00 + K );    
		    Result[1] = Trial_Stress_Vector[1] - (delta_gamma_a + delta_gamma_c) * (-2.00 * G  / 3.00 + K ) - delta_gamma_b*( 4.00  * G / 3.00 + K ); 
		    Result[2] = Trial_Stress_Vector[2] - (delta_gamma_a + delta_gamma_b) * (-2.00 * G / 3.00 + K )  - delta_gamma_c*( 4.00  * G / 3.00 + K );
                    if(mcurrent_Ft[0] <= 0.00) {Result[0] = 1E-14;} 
                    if(mcurrent_Ft[1] <= 0.00) {Result[1] = 1E-14;} 
                    if(mcurrent_Ft[2] <= 0.00) {Result[2] = 1E-14;} 
                 }
	        

	      
           }

            void Rankine_Yield_Function::GetValue(Vector& Result)
                 {
                    Result = mcurrent_Ft; 
                    //KRATOS_WATCH(mcurrent_Ft)   
                 } 

}


