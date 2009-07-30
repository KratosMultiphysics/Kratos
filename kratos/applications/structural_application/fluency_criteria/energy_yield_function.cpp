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
      
		    void  Energy_Yield_Function::InitializeMaterial(const Properties& props) {mprops = &props;}
		     

		    void  Energy_Yield_Function::CalculateEquivalentUniaxialStress(
		    const Vector& StressVector, const Vector& StrainVector, const Matrix& Other, double& Result)

			    {
			    int    iter      = 50;
			    double zero      = 1.0E-9;
			    double teta_a    = 0.00;
			    double teta_b    = 0.00;
			    double teta      = 0.00;
			    //double norma     = 0.00; 
                            double n         = (*mprops)[FC]/(*mprops)[FT];
                       


			    Matrix StressTensor    = ZeroMatrix(3,3);
                            Matrix EigenVectors    = ZeroMatrix(3,3);
			    //Matrix InvertedMatrix  = ZeroMatrix(6,6);
			    Vector PrincipalStress = ZeroVector(3);
			    //Vector Aux_Vector      = ZeroVector(6);

// 			    if( mState == Plane_Stress || mState==Plane_Strain) 
// 				{
// 				  Aux_Vector      = ZeroVector(3);
//                                   InvertedMatrix  = ZeroMatrix(3,3);
// 				}  


		            this->State_Tensor(StressVector,StressTensor);
		            this->Comprobate_State_Tensor(StressTensor, StressVector); // funcion definida en clase base;
                           

			    //SD_MathUtils<double>::InvertMatrix(Other, InvertedMatrix);
			    
			    //norma = SD_MathUtils<double>::normTensor(StressTensor);

                            //if (norma>=1E-3)
			    //PrincipalStress  = SD_MathUtils<double>::EigenValues(StressTensor,crit, zero);
                            SD_MathUtils<double>::EigenVectors(StressTensor, EigenVectors,PrincipalStress, zero, iter);
			   
			    

			    teta_a =  Tensor_Utils<double>::Mc_aully(PrincipalStress);
			    teta_b = norm_1(PrincipalStress);  
			    
 
			    // Condicion para no tener divison de 0/0 
			    if (teta_a==0.00 && teta_b==0.00)
			    {teta = 1.00;}
			    else
			    {teta = teta_a/teta_b;}

			    
			    //noalias(Aux_Vector) = prod(trans(StressVector),InvertedMatrix);
			    //Result              = inner_prod(StressVector,Aux_Vector);
                            Result              = inner_prod(StressVector,StrainVector);
			    Result              = sqrt(Result);
			    Result              = (teta + (1.00-teta)/n)*Result;		    
			    }


		    void  Energy_Yield_Function::CalculateEquivalentUniaxialStressViaInvariants( 
		    const Vector& StressVector,double& Result){}


		    void  Energy_Yield_Function::CalculateEquivalentUniaxialStressViaCilindricalCoordinate( 
		    const Vector& StressVector,double& Result){}



		    void Energy_Yield_Function::CalculateDerivateFluencyCriteria(Vector DerivateFluencyCriteria){}

  
  
}



