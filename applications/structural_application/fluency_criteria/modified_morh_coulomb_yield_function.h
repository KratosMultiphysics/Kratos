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
       //typedef Morh_Coulomb_Yield_Function MorhCoulomb;
       //typedef Isotropic_Rankine_Yield_Function Rankine;
      /*! Explicit elastic predictor/ return mapping algorithm for the associative and non-associative, isotropic
       softening Morh Coulomb with isotropic non-hardening Rankine Tension Cutt Off. */      
      class Modified_Morh_Coulomb_Yield_Function: public FluencyCriteria    //public Isotropic_Rankine_Yield_Function,  public Morh_Coulomb_Yield_Function
          
      { 
    
	public:
	  
	    typedef  FluencyCriteria::Pointer                       FluencyPointerType;      
	    typedef  Morh_Coulomb_Yield_Function::Pointer           MorhCoulombPointerType;
	    typedef  Isotropic_Rankine_Yield_Function::Pointer      RankinePointerType;

            
	    /// empty constructor
	    Modified_Morh_Coulomb_Yield_Function();  
	    
	    /// Two fucntion from python
	    Modified_Morh_Coulomb_Yield_Function(
	    myState State, 
	    const  MorhCoulombPointerType MorhCoulomb, 
	    const  RankinePointerType RanKine);
	    
	    /*
	    /// with derivate class
	    Modified_Morh_Coulomb_Yield_Function(
	    const SoftHardPointerType& SofteningBehavior,           /// For Ft
	    const SoftHardPointerType& SofteningBehaviorCohesion,   /// For Cohesion
	    const SoftHardPointerType& SofteningBehaviorFriction,   /// For Friction Angle
            const SoftHardPointerType& SofteningBehaviorDilatancy,  /// For Dilatancy
	    const myState& State, 
            const myPotencialPlastic& PotencialPlastic);
	    */
	    
            virtual boost::shared_ptr<FluencyCriteria> Clone() const
	        {
		      /*  
		      boost::shared_ptr<FluencyCriteria> p_clone(new Modified_Morh_Coulomb_Yield_Function(     
		      Rankine::mpSofteningBehaviorFt->Clone(),
		      MorhCoulomb::mpSofteningBehavior_Cohesion->Clone(),
		      MorhCoulomb::mpSofteningBehavior_Friction->Clone(),
		      MorhCoulomb::mpSofteningBehavior_Dilatancy->Clone(),
		      mState, 
		      MorhCoulomb::mPotencialPlastic));
		      return p_clone;
		      */
		      /*
		      MorhCoulombPointerType m_clone  = mMorhCoulomb->Clone();
		      RankinePointerType r_clone      = mRankine->Clone();
		      */
		      
		      MorhCoulombPointerType m_clone(new Morh_Coulomb_Yield_Function(                   
		      mMorhCoulomb->mpSofteningBehavior_Cohesion, 
		      mMorhCoulomb->mpSofteningBehavior_Friction,
		      mMorhCoulomb->mpSofteningBehavior_Dilatancy, 
		      mMorhCoulomb->mState, 
		      mMorhCoulomb->mPotencialPlastic));
		      
		      RankinePointerType r_clone(new Isotropic_Rankine_Yield_Function(
		      mRankine->mpSofteningBehaviorFt,
		      mState));
		      
		      
		      boost::shared_ptr<FluencyCriteria> p_clone(new Modified_Morh_Coulomb_Yield_Function(
		      mState, 
		      m_clone,
		      r_clone
		      ));
		      
 		      return p_clone;
		      
		}
   
            
	   
            ~Modified_Morh_Coulomb_Yield_Function();

            KRATOS_CLASS_POINTER_DEFINITION( Modified_Morh_Coulomb_Yield_Function );

//***********************************************************************
//***********************************************************************
   

      bool CheckPlasticAdmisibility(const Vector& Stress); 
      void InitializeMaterial(const Properties& props);  
      void ReturnMapping(const Vector& StrainVector, const Vector& TrialStress, Vector& StressVector);
      void FinalizeSolutionStep();
      void UpdateMaterial();
      void GetValue(const Variable<Matrix>& rVariable, Matrix& Result);
      void GetValue(const Variable<double>& rVariable, double& Result);
      void GetValue(const Variable<Vector>& rVariable, Vector& Result);
      void GetValue(double& Result);
      void GetValue(Matrix& Result);

      
      protected:
	
	
      bool Return_Mapping_Intersection_Main_Plane_And_Sigma1_Tesile_Plane(const array_1d<double,3>& PrincipalStress, const array_1d<unsigned int,3>& order, array_1d<double,3>& Sigma);  
      bool Return_Mapping_Intersection_Main_Plane_Corner_And_Sigma1_Tesile_Plane(const array_1d<double,3>& PrincipalStress, const array_1d<unsigned int,3>& order, array_1d<double,3>& Sigma); 
      bool Return_Mapping_Intersection_Main_Plane_Corner_And_Sigma1_And_Sigma_2_Tesile_Plane(const array_1d<double,3>& PrincipalStress, const array_1d<unsigned int,3>& order, array_1d<double,3>& Sigma); 
      void CalculatePlasticDamage(const array_1d<double,3>& Sigma);
      void ComputeActualStrees(const double& Ppvs, const array_1d<double,3>& Ppds, const array_1d<double,3>& PrincipalStress, array_1d<double,3>& Sigma);
      void ComputeActualStrain(const array_1d<double,3>& Pps);
      
      MorhCoulombPointerType mMorhCoulomb;
      RankinePointerType     mRankine; 
      array_1d<double, 3>    mPrincipalPlasticStrain_current;
      array_1d<double, 3>    mPrincipalPlasticStrain_old;
      
      double m_modified_morh_coulomb_maccumulated_plastic_strain_old; 
      double m_modified_morh_coulomb_maccumulated_plastic_strain_current;
      //double m_modified_morh_coulomb_maccumulated_plastic_strain_current_pos;
      //double m_modified_morh_coulomb_maccumulated_plastic_strain_current_neg;
      
      
      
      
      //double mfracture_strain;
      double  mlength;
      double  mpastic_damage_old;
      double  mpastic_damage_current;
      Vector mplastic_strain; 
      Vector mplastic_strain_old; 
      Vector mElastic_strain; 
      Vector mElastic_strain_old; 
      Matrix m_inv_DeltaF;
	  
	inline double cint(double x){
	  double intpart;
	if (modf(x,&intpart)>=0.5)
	   return  (x>=0)?  ceil(x):floor(x);
	else
	   return (x<0)?  ceil(x):floor(x);
	}

		inline double round(double r,unsigned int places){
	int off=1;
	for(unsigned int i=0; i<places; i++)
		off *= 10;



//		for(unsigned int i=0; i<places;i++
	//int off=pow(10,places);
	return cint(r*off)/off;
        }
	

  

    };
}
#endif


