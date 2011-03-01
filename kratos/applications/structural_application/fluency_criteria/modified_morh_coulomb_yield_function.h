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

#if !defined(MODIFIED_MORH_COULOMB_FUNCTION)
#define MODIFIED_MORH_COULOMB_FUNCTION


#include "custom_utilities/tensor_utils.h"
#include "includes/ublas_interface.h"
#include "fluency_criteria/fluency_criteria.h"
#include "fluency_criteria/morh_coulomb_yield_function.h"
#include "fluency_criteria/isotropic_rankine_yield_function.h"
#include "fluency_criteria/fluency_criteria.h"
#include <cmath>
#include <math.h>



namespace Kratos
  {

      /*! Explicit elastic predictor/ return mapping algorithm for the associative and non-associative, isotropic
       softening Morh Coulomb with isotropic non-hardening Rankine Tension Cutt Off. */
            
      class Modified_Morh_Coulomb_Yield_Function:  public virtual FluencyCriteria 
          
      { 
    
	public:
	  
	    typedef  FluencyCriteria::Pointer                       FluencyPointerType;      
	    typedef  Morh_Coulomb_Yield_Function::Pointer           MorhCoulombPointerType;
	    typedef  Isotropic_Rankine_Yield_Function::Pointer      RankinePointerType;

            
	    Modified_Morh_Coulomb_Yield_Function();  
	    Modified_Morh_Coulomb_Yield_Function(
	    myState State, 
	    const MorhCoulombPointerType MorhCoulomb, 
	    const RankinePointerType RanKine);
	    
	    
	      
            virtual boost::shared_ptr<FluencyCriteria> Clone() const
	        {
		      boost::shared_ptr<FluencyCriteria> p_clone(new Modified_Morh_Coulomb_Yield_Function(
		      mState, 
		      mMorhCoulomb,
		      mRankine 
		      ));
		      return p_clone;
		}
   
            
	   
            ~Modified_Morh_Coulomb_Yield_Function();

            KRATOS_CLASS_POINTER_DEFINITION( Modified_Morh_Coulomb_Yield_Function );

//***********************************************************************
//***********************************************************************
   

      
      void InitializeMaterial(const Properties& props);  
      void ReturnMapping(const Vector& StrainVector, Vector& StressVector);
      void FinalizeSolutionStep();
      void UpdateMaterial();
      void GetValue(double Result);
      void GetValue(const Variable<Matrix>& rVariable, Matrix& Result);
      void GetValue(const Variable<double>& rVariable, double& Result);
   

      private:
	
	
      bool Return_Mapping_Intersection_Main_Plane_And_Sigma1_Tesile_Plane(const array_1d<double,3>& PrincipalStress, const array_1d<unsigned int,3>& order, array_1d<double,3>& Sigma);  
      bool Return_Mapping_Intersection_Main_Plane_Corner_And_Sigma1_Tesile_Plane(const array_1d<double,3>& PrincipalStress, const array_1d<unsigned int,3>& order, array_1d<double,3>& Sigma); 
      bool Return_Mapping_Intersection_Main_Plane_Corner_And_Sigma1_And_Sigma_2_Tesile_Plane(const array_1d<double,3>& PrincipalStress, const array_1d<unsigned int,3>& order, array_1d<double,3>& Sigma); 
      void CalculatePlasticDamage(const array_1d<double,3>& Sigma);
      
      MorhCoulombPointerType mMorhCoulomb;
      RankinePointerType     mRankine; 
      array_1d<double, 3>    mPrincipalPlasticStrain_current;
      array_1d<double, 3>    mPrincipalPlasticStrain_old;
      
      double m_modified_morh_coulomb_maccumulated_plastic_strain_old; 
      double m_modified_morh_coulomb_maccumulated_plastic_strain_current;
      //double mfracture_strain;
      double  mlength;
      double  mpastic_damage_old;
      double  mpastic_damage_current;
      
	  
	inline double cint(double x){
	  double intpart;
	if (modf(x,&intpart)>=0.5)
	   return  (x>=0)?  ceil(x):floor(x);
	else
	   return (x<0)?  ceil(x):floor(x);
	}

	inline double round(double r,unsigned places){
	double off=pow(10,places);
	return cint(r*off)/off;
        }
	

  

    };
}
#endif


