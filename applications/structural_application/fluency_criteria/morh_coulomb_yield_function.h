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

#if !defined(MORH_COULOMB_FUNCTION)
#define MORH_COULOMB_FUNCTION



#include "includes/ublas_interface.h"
#include "fluency_criteria/fluency_criteria.h"
#include "soft_hard_behavior/softening_hardening_criteria.h"
#include "utilities/math_utils.h"
#include "custom_utilities/sd_math_utils.h"
#include <cmath>



namespace Kratos
  {

      typedef SofteningHardeningCriteria::Pointer SoftHardPointerType;    
      
      /*! 
      STATE UPDATE PROCEDURE FOR MOHR-COULOMB TYPE ELASTO-PLASTIC MATERIAL
      WITH ASSOCIATIVE/NON-ASSOCIATIVE FLOW RULE AND PIECE-WISE LINEAR
      ISOTROPIC HARDENING:
      IMPLICIT ELASTIC PREDICTOR/RETURN MAPPING ALGORITHM (BOXES 8.4-7)
      PLANE STRAIN AND AXISYMMETRIC IMPLMENTATIONS
      */
      class Morh_Coulomb_Yield_Function: public FluencyCriteria    
      { 
    
        public:

              
            virtual boost::shared_ptr<FluencyCriteria> Clone() const
	        {
		      boost::shared_ptr<FluencyCriteria> p_clone(new Morh_Coulomb_Yield_Function(
		      mpSofteningBehavior,
		      mState, 
		      mPotencialPlastic));
		      return p_clone;
		}
  
            Morh_Coulomb_Yield_Function(); 
             
            Morh_Coulomb_Yield_Function(SoftHardPointerType SofteningBehavior,
                                        myState State, 
	                                myPotencialPlastic PotencialPlastic);
					
            ~Morh_Coulomb_Yield_Function();

            KRATOS_CLASS_POINTER_DEFINITION( Morh_Coulomb_Yield_Function );

	    
//***********************************************************************
//***********************************************************************



    bool CheckPlasticAdmisibility(const Vector& Stress);  
    void InitializeMaterial(const Properties& props);  
    void ReturnMapping(const Vector& StrainVector, Vector& StressVector);
    void FinalizeSolutionStep();
    void UpdateMaterial();
    void GetValue(const Variable<Matrix>& rVariable, Matrix& Result);
    void GetValue(const Variable<double>& rVariable, double& Result);
    
	
    private:
      
    bool ReturnMappingToMainPlane(const array_1d<double,3>& PrincipalStress, const array_1d<unsigned int,3>& order, array_1d<double,3>& Sigma);
    bool TwoVectorReturnToEdges  (const array_1d<double,3>& PrincipalStress, const array_1d<unsigned int,3>& order, const bool& edges, array_1d<double,3>& Sigma);
    void ReturnMappingToApex     (const array_1d<double,3>& PrincipalStress, array_1d<double, 3 >& Sigma);  
        
    bool CheckValidity(const array_1d<double,3>& Sigma);
    bool ReturnToEdges(const array_1d<double,3>& PrincipalStress);
  

      
    public:
    
    double mcohesion;
    double mcurrent_cohesion;
    double mdilatancy_angle;
    double mcurrent_dilatancy_angle;
    double minternal_friction_angle;
    double mcurrent_minternal_friction_angle;    
    double mmorh_coulomb_maccumulated_plastic_strain_old; 
    double mmorh_coulomb_maccumulated_plastic_strain_current;
    SoftHardPointerType mpSofteningBehavior;
    array_1d<double, 3>    mPrincipalPlasticStrain_current;
    array_1d<double, 3>    mPrincipalPlasticStrain_old;  
    
     
     ///@}
     ///@name Serialization
     ///@{
    
  

    };
}
#endif











