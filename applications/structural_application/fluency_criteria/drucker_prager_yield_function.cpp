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
#include "fluency_criteria/drucker_prager_yield_function.h"
#include <cmath>



namespace Kratos
  {

  
	    typedef boost::numeric::ublas::vector<Vector> Second_Order_Tensor; // dos opciones: un tensor de segundo orden y/o un vector que almacena un vector		  
	    typedef boost::numeric::ublas::vector<Second_Order_Tensor> Third_Order_Tensor;
			  
            typedef boost::numeric::ublas::vector<boost::numeric::ublas::vector<Matrix> > Fourth_Order_Tensor;
			  
	    typedef matrix<Second_Order_Tensor> Matrix_Second_Tensor; // Acumulo un tensor de 2 orden en una matri
           
            Drucker_Prager_Yield_Function::Drucker_Prager_Yield_Function(myState State)
	    :FluencyCriteria()
	    {
	     
              mState = State;
	    }

             Drucker_Prager_Yield_Function::~Drucker_Prager_Yield_Function() {}

    
//***********************************************************************
//***********************************************************************

		    void Drucker_Prager_Yield_Function::InitializeMaterial(const Properties& props) {mprops = &props; }
		     

		    void Drucker_Prager_Yield_Function:: CalculateEquivalentUniaxialStress(
		    const Vector& StressVector,double& Result)  {}


		    void Drucker_Prager_Yield_Function::CalculateEquivalentUniaxialStressViaPrincipalStress(
		    const Vector& StressVector,double& Result)
                       {
                       }



		    void Drucker_Prager_Yield_Function::CalculateEquivalentUniaxialStressViaInvariants(
		    const Vector& StressVector,double& Result)
			   {
		
		      // Nota: la resistencia de comparacion de esta superficie de fluencia es FC.

                      
		      double frictional_internal = 0.00;
                      double Alfa   = 0.00;
                      double Beta  = 0.00;

			
                      Vector I          = ZeroVector(3);
		      Vector J          = ZeroVector(3);
                      Vector J_des      = ZeroVector(3);		      

		      Matrix StressTensor     = ZeroMatrix(3,3);
		      Vector PrincipalStress  = ZeroVector(3);

                      frictional_internal = (*mprops)[FRICTION_INTERNAL_ANGLE]*PI/180.00;
			    
		      double sin_phi = sin(frictional_internal);
                      double cos_phi = cos(frictional_internal);			
  

		      this->State_Tensor(StressVector,StressTensor);
		      this->Comprobate_State_Tensor(StressTensor, StressVector); // funcion definida en clase base;
		      Tensor_Utils<double>::TensorialInvariants(StressTensor, I, J, J_des);

		      Alfa = 2.00* sin_phi/(sqrt(3.00)*(3.00 - sin_phi));
                      Beta = 6.00*(*mprops)[COHESION]*cos_phi/(sqrt(3.00)*(3.00- sin_phi));
                      Result   = I(0)*Alfa + sqrt(J_des(1));  // Beta es la resitencia de comparacion
		      mSigma_e = Result; 
                      mSigma_y = Beta;
                      Result -= mSigma_y;  
                      //KRATOS_WATCH(Beta) 

                           }


		    void Drucker_Prager_Yield_Function::CalculateEquivalentUniaxialStressViaCilindricalCoordinate(
		    const Vector& StressVector,double& Result){}



		    void Drucker_Prager_Yield_Function::CalculateDerivateFluencyCriteria(const Vector& StressVector, Vector& DerivateFluencyCriteria)
			{
			    Second_Order_Tensor a;  
                            Vector C = ZeroVector(3);
                            DerivateFluencyCriteria = ZeroVector(6);
                            double frictional_internal = 0.00;
                            double Alfa   = 0.00;
        

                            frictional_internal = (*mprops)[FRICTION_INTERNAL_ANGLE]*PI/180.00;

			    
		            double sin_phi = sin(frictional_internal);
                            double cos_phi = cos(frictional_internal);

                            Alfa = 2.00*sin_phi/(sqrt(3.00)*(3.00 - sin_phi));
                            
			    C(0) =  Alfa; 
                            C(1) =  1.00;
                            C(2) =  0.00;
                            
			    this->CalculateVectorFlowDerivate(StressVector, a);
			    for(unsigned  int i=0; i<3; i++)
                               {
                                 noalias(DerivateFluencyCriteria) = DerivateFluencyCriteria + a(i)*C(i); 
                               }  
			      
			}


    }


