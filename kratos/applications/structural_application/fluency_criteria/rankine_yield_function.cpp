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

		    void Rankine_Yield_Function::InitializeMaterial(const Properties& props) { mprops = &props;}
		     

		    void Rankine_Yield_Function:: CalculateEquivalentUniaxialStress(
		    const Vector& StressVector,double& Result){}


		    void Rankine_Yield_Function::CalculateEquivalentUniaxialStressViaPrincipalStress(
		    const Vector& StressVector,double& Result)
                       {
		      double crit = 1E-15;
                      double zero = 1E-15;  
                      double max  = 0.00;
                      unsigned int dim  = 3;


		      Matrix StressTensor    = ZeroMatrix(dim,dim);
                      Vector PrincipalStress = ZeroVector(dim);
                      this->State_Tensor(StressVector,StressTensor);
		      this->Comprobate_State_Tensor(StressTensor, StressVector);
                      PrincipalStress = SD_MathUtils<double>::EigenValues(StressTensor,crit, zero);
		      max = (*std::max_element(PrincipalStress.begin(),PrincipalStress.end()));
                      //Result = max; // - msigma_max; 
                      //KRATOS_WATCH(Result) 
                   	}



		    void Rankine_Yield_Function::CalculateEquivalentUniaxialStressViaInvariants(
		    const Vector& StressVector,double& Result)
			   {
				      
	              unsigned int dim  = 3;
                      long double tetha_Lode = 0.00;
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

		      Result = 2.00*sqrt(3.00*J_des(1))*cos(tetha_Lode + PI/6.00) + I(0); // - msigma_max;
                      Result = Result/3.00;
		      //KRATOS_WATCH("----------")
                      //KRATOS_WATCH(Result)
				    

			    }


		    void Rankine_Yield_Function::CalculateEquivalentUniaxialStressViaCilindricalCoordinate(
		    const Vector& StressVector,double& Result){}



		    void Rankine_Yield_Function::CalculateDerivateFluencyCriteria(Vector DerivateFluencyCriteria){}


    }


