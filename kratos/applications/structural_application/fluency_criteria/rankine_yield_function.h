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

#if !defined(RANKINE_FUNCTION_UTILS)
#define RANKINE_FUNCTION_UTILS


#include "custom_utilities/tensor_utils.h"
#include "fluency_criteria/fluency_criteria.h"
#include <cmath>



namespace Kratos
  {

    	
      class Rankine_Yield_Function: public FluencyCriteria    
      { 
    
        public:

	    typedef boost::numeric::ublas::vector<Vector> Second_Order_Tensor; // dos opciones: un tensor de segundo orden y/o un vector que almacena un vector		  
	    typedef boost::numeric::ublas::vector<Second_Order_Tensor> Third_Order_Tensor;
			  
            typedef boost::numeric::ublas::vector<boost::numeric::ublas::vector<Matrix> > Fourth_Order_Tensor;
			  
	    typedef matrix<Second_Order_Tensor> Matrix_Second_Tensor; // Acumulo un tensor de 2 orden en una matriz

            virtual boost::shared_ptr<FluencyCriteria> Clone() const
	        {
		      boost::shared_ptr<FluencyCriteria> p_clone(new Rankine_Yield_Function(mState));
		      return p_clone;
		}
  
            Rankine_Yield_Function(myState State);
	   
            ~Rankine_Yield_Function();

            KRATOS_CLASS_POINTER_DEFINITION( Rankine_Yield_Function );

//***********************************************************************
//***********************************************************************
// Energy_Criteria
// Diferent limits in traccion and compresion


		     void InitializeMaterial(const Properties& props);
		     

		    void CalculateEquivalentUniaxialStress(  
		    const Vector& StressVector,double& Result); 


		    void CalculateEquivalentUniaxialStressViaPrincipalStress(
		    const Vector& StressVector,double& Result);



		    void CalculateEquivalentUniaxialStressViaInvariants(
		    const Vector& StressVector,double& Result);


		    void CalculateEquivalentUniaxialStressViaCilindricalCoordinate(
		    const Vector& StressVector,double& Result);


		    void CalculateDerivateFluencyCriteria(const Vector& StressVector, Vector& DerivateFluencyCriteria);

                    void UpdateVariables(const Vector& Variables); 

                    void ReturnMapping(const Vector& StressVector, 
                    Vector& delta_lamda,
                    array_1d<double,3>& Result);

                    void GetValue(Vector& Result);

		    

      
	protected:

	  array_1d<double,3> mFt;
          array_1d<double,3> mcurrent_Ft;
          double mH;
          bool  minitialize;     
           
           

    };
}
#endif


