//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Created by:          $Author:                    LHauser $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                    July 2018 $
//   Revision:            $Revision:                      0.0 $
//
//

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "custom_utilities/stress_invariants_utilities.hpp"
#include "custom_constitutive/custom_flow_rules/borja_casm_explicit_plastic_flow_rule.hpp"

#include "pfem_solid_mechanics_application_variables.h"

namespace Kratos
{


	// Version of Casm-like hiperelastic model. This is the Hyperelsatic model of Borja, that includes the Houlsby and the constant shear modulus as special cases.

	//************ CONSTRUCTOR ***********
	BorjaCasmExplicitFlowRule::BorjaCasmExplicitFlowRule()
		:NonAssociativeExplicitPlasticFlowRule()
	{
	}

	//*****************************INITIALIZATION CONSTRUCTOR*****************************
	//************************************************************************************

	BorjaCasmExplicitFlowRule::BorjaCasmExplicitFlowRule(YieldCriterionPointer pYieldCriterion)
		:NonAssociativeExplicitPlasticFlowRule(pYieldCriterion)
	{
		std::cout<<"   CASM FLOW RULE constructed"<<std::endl;
	}

  //********* ASSIGMENT OPERATOR
  BorjaCasmExplicitFlowRule& BorjaCasmExplicitFlowRule::operator=(BorjaCasmExplicitFlowRule const& rOther)
  {
		NonAssociativeExplicitPlasticFlowRule::operator=(rOther);
      return *this;
  }

  //********** COPY CONSTRUCTOR *********
  BorjaCasmExplicitFlowRule::BorjaCasmExplicitFlowRule(BorjaCasmExplicitFlowRule const& rOther)
		:NonAssociativeExplicitPlasticFlowRule(rOther)
  {
  }

  //*******   CLONE ********
  FlowRule::Pointer BorjaCasmExplicitFlowRule::Clone() const
  {
		FlowRule::Pointer p_clone(new BorjaCasmExplicitFlowRule(*this));
			return p_clone;
  }

  // ********** DESTRUCTOR **************
  BorjaCasmExplicitFlowRule::~BorjaCasmExplicitFlowRule()
  {
  }


  // *********** EVALUATE KIRCHHOFF STRESS VECTOR *****************************
  void BorjaCasmExplicitFlowRule::CalculateKirchhoffStressVector(const Vector& rHenckyStrainVector, Vector& rKirchhoffStressVector)
  {
//std::cout<<"BorjaCasmExplicitFlowRule::CalculateKirchhoffStressVector"<<std::endl;
		Vector IdentityVector = ZeroVector(6);

		for (unsigned int i = 0; i < 3; ++i)
			IdentityVector(i) = 1.0;

		double MeanStress; 
		double VolumetricStrain = MathUtils<double>::Dot( trans(rHenckyStrainVector), IdentityVector);

		Vector DeviatoricStrainVector; 
		DeviatoricStrainVector = rHenckyStrainVector -  (VolumetricStrain/3.0)*IdentityVector;   

		EvaluateMeanStress(VolumetricStrain, DeviatoricStrainVector,  MeanStress );

		EvaluateDeviatoricStress( VolumetricStrain, DeviatoricStrainVector, rKirchhoffStressVector);

		rKirchhoffStressVector += MeanStress*IdentityVector;
  }


  // ************ EVALUATE ONLY THE VOLUMETRIC PART OF THE HYPERELASTIC MODEL ****
  void BorjaCasmExplicitFlowRule::EvaluateMeanStress(const double& rVolumetricStrain, const Vector& rDeviatoricStrainVector, double& rMeanStress)
  {
		double SwellingSlope = mpYieldCriterion->GetHardeningLaw().GetProperties()[SWELLING_SLOPE];
    double AlphaShear = mpYieldCriterion->GetHardeningLaw().GetProperties()[ALPHA_SHEAR];

    double ReferencePressure = mpYieldCriterion->GetHardeningLaw().GetProperties()[PRE_CONSOLIDATION_STRESS];
    double OCR = mpYieldCriterion->GetHardeningLaw().GetProperties()[OVER_CONSOLIDATION_RATIO];
    ReferencePressure /= OCR;    

    double DeviatoricStrain2Norm = 0.0;
    for (unsigned int i = 0; i < 3; ++i)
      DeviatoricStrain2Norm += pow(rDeviatoricStrainVector(i), 2);

    for (unsigned int i = 3; i < 6; ++i)
      DeviatoricStrain2Norm += 2.0*pow(rDeviatoricStrainVector(i)/2.0, 2);

    rMeanStress = -ReferencePressure*std::exp( -rVolumetricStrain / SwellingSlope) * (1.0 + 1.0*AlphaShear*DeviatoricStrain2Norm/SwellingSlope);
  }

  void BorjaCasmExplicitFlowRule::EvaluateMeanStress(const Vector& rHenckyStrainVector, double& rMeanStress)
  {
    double VolumetricStrain = 0.0;
    for (unsigned int i = 0; i < 3; ++i)
			VolumetricStrain += rHenckyStrainVector(i);

    Vector DeviatoricStrain = rHenckyStrainVector;

    for (unsigned int i = 0; i < 3; ++i)
			DeviatoricStrain(i) -= VolumetricStrain/3.0;

    this->EvaluateMeanStress(VolumetricStrain, DeviatoricStrain, rMeanStress);
  }


	// ************* EVALUTE ONLY THE VOLUMETRIC PART OF THE HYPERELASTIC MODEL ***
  void BorjaCasmExplicitFlowRule::EvaluateDeviatoricStress(const double& rVolumetricStrain, const Vector & rDeviatoricStrainVector, Vector& rDeviatoricStress)
  {

    double ReferencePressure = mpYieldCriterion->GetHardeningLaw().GetProperties()[PRE_CONSOLIDATION_STRESS];
    double OCR = mpYieldCriterion->GetHardeningLaw().GetProperties()[OVER_CONSOLIDATION_RATIO];
    ReferencePressure /= OCR;    
    double SwellingSlope = mpYieldCriterion->GetHardeningLaw().GetProperties()[SWELLING_SLOPE];
    double AlphaShear = mpYieldCriterion->GetHardeningLaw().GetProperties()[ALPHA_SHEAR];
    double ConstantShearModulus = mpYieldCriterion->GetHardeningLaw().GetProperties()[INITIAL_SHEAR_MODULUS];


    rDeviatoricStress = rDeviatoricStrainVector;
    double ShearModulus = AlphaShear*ReferencePressure*std::exp( -rVolumetricStrain / SwellingSlope);
    rDeviatoricStress *= 2.0*( ShearModulus + ConstantShearModulus);

    for (unsigned int i = 3; i<6; ++i){
      rDeviatoricStress(i) /= 2.0;  // BECAUSE VOIGT NOTATION
    }
  }


	// *********** EVALUATE THE TANGENT ELASTIC MATRIX ************
	void BorjaCasmExplicitFlowRule::ComputeElasticMatrix(const Vector& rElasticStrainVector, Matrix& rElasticMatrix )
	{
//std::cout<<"BorjaCasmExplicitFlowRule::ComputeElasticMatrix"<<std::endl;
		Matrix FourthOrderIdentity = ZeroMatrix(6,6);
		for (unsigned int i = 0; i<3; ++i)
			 FourthOrderIdentity(i,i) = 1.0;

		for (unsigned int i = 3; i<6; ++i)
			 FourthOrderIdentity(i,i) = 0.50;
		// VOIGT NOTATION AND NOT KELVIN

		Matrix IdentityCross = ZeroMatrix(6,6);
		for (unsigned int i = 0; i<3; ++i) {
			 for (unsigned int j = 0; j<3; ++j) {
					IdentityCross(i,j) = 1.0;
			 }
		}


		Vector StressVector = ZeroVector(6);
		this->CalculateKirchhoffStressVector(rElasticStrainVector, StressVector);


		double MeanStress = 0.0;
		double VolumetricStrain = 0.0;
		for (unsigned int i = 0; i<3; i++) {
			 MeanStress += StressVector(i);
			 VolumetricStrain += rElasticStrainVector(i);
		}
		MeanStress /= 3.0;


		double ReferencePressure = mpYieldCriterion->GetHardeningLaw().GetProperties()[PRE_CONSOLIDATION_STRESS];
		double OCR = mpYieldCriterion->GetHardeningLaw().GetProperties()[OVER_CONSOLIDATION_RATIO];
		ReferencePressure /= OCR;    

		double SwellingSlope = mpYieldCriterion->GetHardeningLaw().GetProperties()[SWELLING_SLOPE];
		double AlphaShear = mpYieldCriterion->GetHardeningLaw().GetProperties()[ALPHA_SHEAR];

		double ConstantShearModulus = mpYieldCriterion->GetHardeningLaw().GetProperties()[INITIAL_SHEAR_MODULUS];


		rElasticMatrix  = (-1.0/SwellingSlope)*MeanStress*IdentityCross;
		rElasticMatrix += 2.0*AlphaShear*ReferencePressure*std::exp(-VolumetricStrain/SwellingSlope)*(FourthOrderIdentity - (1.0/3.0)*IdentityCross);



		Vector StrainVector = rElasticStrainVector; 
		for (unsigned int i = 0; i < 3; ++i)
			 StrainVector(i) -= VolumetricStrain / 3.0;

		double Modulus = 2.0 * ReferencePressure * exp( - VolumetricStrain/SwellingSlope) * ( AlphaShear / SwellingSlope);

		for (unsigned int i = 3; i < 6 ; ++i)
			 StrainVector(i) /= 2.0;


		// PARTE ASQUEROSA
		for (unsigned int i = 0; i<3; ++i) {
			 for (unsigned int j = 0; j<3; ++j) {
					rElasticMatrix(i,j) -= Modulus * (StrainVector(i) ); //-MeanStress);
					rElasticMatrix(i,j) -= Modulus * (StrainVector(j) ); //-MeanStress);
			 }
		}

		for (unsigned int i = 0; i<3; ++i) {
			 for (unsigned int j = 3; j < 6; ++j) {
					rElasticMatrix(i,j) -= Modulus*(StrainVector(j));///2.0;
			 }
		}

		for (unsigned int i = 3; i<6; ++i) {
			 for (unsigned int j = 0; j<3; ++j) {
					rElasticMatrix(i,j) -= Modulus*(StrainVector(i));///2.0;
			 }
		}

		// AND THE PART DUE TO THE INITIAL SHEAR MODULUS

		rElasticMatrix +=  2.0*ConstantShearModulus * ( FourthOrderIdentity - (1.0/3.0)*IdentityCross );

	}

		// COMPUTE THE PLASTIC HARDENING PARAMETER
		void BorjaCasmExplicitFlowRule::ComputePlasticHardeningParameter(const Vector& rHenckyStrainVector, const double& rAlpha, double& rH)
		{
//std::cout<<"BorjaCasmExplicitFlowRule::ComputePlasticHardeningParameter"<<std::endl;
			//
			double PreconsolidationStress = 0.0;
			PreconsolidationStress 				= mpYieldCriterion->GetHardeningLaw().CalculateHardening(PreconsolidationStress, rAlpha);
			const double SpacingR 				= mpYieldCriterion->GetHardeningLaw().GetProperties()[SPACING_RATIO];
			const double SwellingSlope		= mpYieldCriterion->GetHardeningLaw().GetProperties()[SWELLING_SLOPE];
      const double OtherSlope      	= mpYieldCriterion->GetHardeningLaw().GetProperties()[NORMAL_COMPRESSION_SLOPE];

			//calcualte Kirchhoff stress vector
			Vector StressVector = ZeroVector(6);
			this->CalculateKirchhoffStressVector(rHenckyStrainVector, StressVector);

			//calculate d_h/d_eps^p
			Vector PreconDerivativeEpsVol = ZeroVector(6);
			for (unsigned int i = 0; i<3; ++i)
				PreconDerivativeEpsVol(i) = 1.0;
			PreconDerivativeEpsVol *= -PreconsolidationStress/(OtherSlope - SwellingSlope);
			
			//calculate d_g/d_sig
			Vector PlasticPotentialD = ZeroVector(6);
			Matrix PlasticPotentialDD = ZeroMatrix(1,1);
			this->CalculatePlasticPotentialDerivatives(StressVector, PlasticPotentialD, PlasticPotentialDD, rAlpha);
			
			//compute H = - d_f/d_h * < d_h/d_eps^p, d_g/d_sig >
			rH = pow(log(SpacingR)*PreconsolidationStress,-1.0);
			rH *= MathUtils<double>::Dot( PreconDerivativeEpsVol, PlasticPotentialD );
/*

      double MeanStress;
      this->EvaluateMeanStress(rHenckyStrainVector, MeanStress);

      double PreconsolidationStress = mpYieldCriterion->GetHardeningLaw().CalculateHardening(PreconsolidationStress, rAlpha);
      double SwellingSlope          = mpYieldCriterion->GetHardeningLaw().GetProperties()[SWELLING_SLOPE];
      double OtherSlope             = mpYieldCriterion->GetHardeningLaw().GetProperties()[NORMAL_COMPRESSION_SLOPE];
      
      rH = (2.0*MeanStress-PreconsolidationStress) ;
      rH *= (-MeanStress);
      rH *= PreconsolidationStress/ ( OtherSlope - SwellingSlope);
*/
		}

	// EVALUATE THE PLASTIC POTENTIAL DERIVATIVES
	void BorjaCasmExplicitFlowRule::CalculatePlasticPotentialDerivatives( const Vector& rStressVector, Vector& rFirstDerivative, Matrix & rSecondDerivative, const double& rAlpha)
	{
		// set second derivatives to zero as not needed for explicit method
		//rFirstDerivative = ZeroVector(1);
		rSecondDerivative = ZeroMatrix(1,1);
		//return;
		
		// get material constants and preconsolidation pressure
		const double ShearM = mpYieldCriterion->GetHardeningLaw().GetProperties()[CRITICAL_STATE_LINE];
		
		// stress invariants & invariants derivatives
		double MeanStress, J2, LodeAngle;
		Vector V1, V2;
		StressInvariantsUtilities::CalculateStressInvariants( rStressVector, MeanStress, J2, LodeAngle);
		StressInvariantsUtilities::CalculateDerivativeVectors( rStressVector, V1, V2);
		
		//calculate third invariant effect 
		double ThirdInvariantEffect = 1.0;
		ThirdInvariantEffect = mpYieldCriterion->EvaluateThirdInvariantEffectMC( LodeAngle);
		
		// calculate d_g/d_sig = d_g/d_p * d_p/d_sig + d_g/d_J * d_J/d_sig + 0
		rFirstDerivative = 27.0*( ShearM/ThirdInvariantEffect*MeanStress + pow(3.0,0.5)*J2 ) / ( (3.0*MeanStress-2.0*pow(3.0,0.5)*J2) * (3.0*MeanStress+pow(3.0,0.5)*J2) ) * V1;
		rFirstDerivative += 3.0*pow(3.0,0.5)*(-9.0*MeanStress - ShearM/ThirdInvariantEffect*(3.0*MeanStress + 2.0*pow(3.0,0.5)*J2)) / ( (3.0*MeanStress-2.0*pow(3.0,0.5)*J2) * (3.0*MeanStress+pow(3.0,0.5)*J2) ) * V2;
	}

  void BorjaCasmExplicitFlowRule::save( Serializer& rSerializer) const 
  {
		KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, NonAssociativeExplicitPlasticFlowRule )
  }

  void BorjaCasmExplicitFlowRule::load( Serializer& rSerializer)
  {
		KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, NonAssociativeExplicitPlasticFlowRule )
  }

} //end namespace kratos
