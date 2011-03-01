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

#if !defined(ISOTROPIC_RANKINE_FUNCTION)
#define ISOTROPIC_RANKINE_FUNCTION

#include "fluency_criteria/fluency_criteria.h"
#include <cmath>



namespace Kratos
  {

      
      /*! 
      STATE UPDATE PROCEDURE FOR ISOTROPIC RANKINE TYPE ELASTO-PLASTIC MATERIAL
      WITH ASSOCIATIVE FLOW RULE AND PIECE-WISE LINEAR
      ISOTROPIC SOFTENING
      IMPLICIT ELASTIC PREDICTOR/RETURN MAPPING ALGORITHM 
      PLANE STRAIN AND AXISYMMETRIC IMPLMENTATIONS
      */
      
      class Isotropic_Rankine_Yield_Function: public virtual FluencyCriteria    
      { 
    
        public:

            virtual boost::shared_ptr<FluencyCriteria> Clone() const
	        {
		      boost::shared_ptr<FluencyCriteria> p_clone(new Isotropic_Rankine_Yield_Function(mState));
		      return p_clone;
		}
  
            Isotropic_Rankine_Yield_Function();              
   
            Isotropic_Rankine_Yield_Function(myState State);
	   
            ~Isotropic_Rankine_Yield_Function();

            KRATOS_CLASS_POINTER_DEFINITION( Isotropic_Rankine_Yield_Function );

//***********************************************************************
//***********************************************************************



	void InitializeMaterial(const Properties& props);  
	bool CheckPlasticAdmisibility(const Vector& Stress); 
	void ReturnMapping(const Vector& StrainVector, Vector& StressVector);
	void FinalizeSolutionStep();
	void UpdateMaterial();
	bool CheckValidity( array_1d<double,3>&  Sigma); 
        void GetValue(const Variable<Matrix>& rVariable, Matrix& Result);

      
	public:

	double  mrankine_accumulated_plastic_strain_current;   
        double  mrankine_accumulated_plastic_strain_old;    
	double mFt;
	double mcurrent_Ft;
	double mH; 
	bool   minitialize;     
	array_1d<double, 3>    mPrincipalPlasticStrain_current;
	array_1d<double, 3>    mPrincipalPlasticStrain_old;  
	
        
        private:
        bool One_Vector_Return_Mapping_To_Main_Plane(const array_1d<double,3>& PrincipalStress, Vector& delta_lamda,    array_1d<double,3>& Sigma); 
        bool Two_Vector_Return_Mapping_To_Corner (   const array_1d<double,3>& PrincipalStress, Vector& delta_lamda ,    array_1d<double,3>& Sigma);        
        void Three_Vector_Return_Mapping_To_Apex (   const array_1d<double,3>& PrincipalStress, Vector& delta_lamda ,array_1d<double,3>& Sigma);
	
	enum   Cases {right, left};
	Cases  mCases; 
   
           

    };
}
#endif


