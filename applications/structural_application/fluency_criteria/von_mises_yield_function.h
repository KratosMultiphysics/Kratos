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

/* *********************************************************   
*          
*   Last Modified by:    $Author: Nelson $
*   Date:                $Date: 2011-11-9 
\*  Revision:           $Revision: 1.2 $
*
* ***********************************************************/

#if !defined(VON_MISES_FUNCTION_UTILS)
#define VON_MISES_FUNCTION_UTILS


#include "custom_utilities/tensor_utils.h"
#include "fluency_criteria/fluency_criteria.h"
#include "soft_hard_behavior/softening_hardening_criteria.h"
#include <cmath>



namespace Kratos
  {

    	
      class Von_Mises_Yield_Function: public FluencyCriteria    
      { 
    
        public:

	    typedef SofteningHardeningCriteria::Pointer SoftHardPointerType;   
	    typedef boost::numeric::ublas::vector<Vector> Second_Order_Tensor; // dos opciones: un tensor de segundo orden y/o un vector que almacena un vector		  
	    typedef boost::numeric::ublas::vector<Second_Order_Tensor> Third_Order_Tensor;
			  
            typedef boost::numeric::ublas::vector<boost::numeric::ublas::vector<Matrix> > Fourth_Order_Tensor;
			  
	    typedef matrix<Second_Order_Tensor> Matrix_Second_Tensor; // Acumulo un tensor de 2 orden en una matri    
  
            //Von_Mises_Yield_Function(const myState& State, myPotencialPlastic PotencialPlastic);
            Von_Mises_Yield_Function(const myState& State, const SoftHardPointerType& Soft);
	    
            virtual boost::shared_ptr<FluencyCriteria> Clone() const
	        {
		      //boost::shared_ptr<FluencyCriteria> p_clone(new Von_Mises_Yield_Function(mState,mPotencialPlastic));
		      boost::shared_ptr<FluencyCriteria> p_clone(new Von_Mises_Yield_Function(
		      mState, 
                      mpSigmaBehavior->Clone()));
		      return p_clone;
		}
	   
            ~Von_Mises_Yield_Function();

            KRATOS_CLASS_POINTER_DEFINITION( Von_Mises_Yield_Function );

//***********************************************************************
//***********************************************************************


		     void InitializeMaterial(const Properties& props);  
		     void ReturnMapping(const Vector& StrainVector, Vector& StressVector);
		     void FinalizeSolutionStep();
		     void UpdateMaterial();
		     void GetValue(const Variable<double>& rVariable, double& Result);
		     void CalculateElasticStrain(Vector& DesviatoricStress, double Pressure, Vector& Elastic_Strain);
		     double UniaxialTension(const Vector& Desviatoric_Trial_Stress);
		     void CalculatePlasticDamage(const Vector& Elastic_Stress, const Vector& Desviatoric_Trial_Stress);
		     //void GetValue(const Variable<Matrix>& rVariable, Matrix& Result);
		     //void GetValue(const Variable<Vector>& rVariable, Vector& Result);
		     void GetValue(double& Result);
		     //void GetValue(Matrix& Result);
		     
	private:
	double mplastic_damage_old;
	double mplastic_damage_new;
	double mhe;
        double mSigma_yold;
	double mSigma_ynew;
	double melastic_z_old;
	double melastic_z_new;
        bool   mcompute_tangent_matrix;
        double maccumulated_plastic_strain_current;   
        double maccumulated_plastic_strain_old;
        Vector  mold_plastic_strain;
        Vector  mcurrent_plastic_strain;
	SoftHardPointerType mpSigmaBehavior;

    };
}
#endif


