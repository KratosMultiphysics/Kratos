//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Created by:          $Author:                    LHauser $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                    July 2018 $
//   Revision:            $Revision:                      0.0 $
//
//

// System includes
#include <iostream>

// External includes
#include<cmath>

// Project includes
#include "includes/properties.h"
#include "custom_constitutive/borja_hencky_casm_axisym_2D_law.hpp"
#include "pfem_solid_mechanics_application_variables.h"
#include "custom_utilities/solid_mechanics_math_utilities.hpp"

namespace Kratos
{

	//******************************CONSTRUCTOR*******************************************
	//************************************************************************************

	BorjaHenckyCasmPlasticAxisym2DLaw::BorjaHenckyCasmPlasticAxisym2DLaw()
		: NonLinearHenckyElasticPlasticAxisym2DLaw()
	{
		mpHardeningLaw   = HardeningLaw::Pointer( new CamClayHardeningLaw() );
		mpYieldCriterion = YieldCriterion::Pointer( new CasmYieldCriterion(mpHardeningLaw) );
		mpFlowRule       = FlowRule::Pointer( new BorjaCasmExplicitFlowRule(mpYieldCriterion) );
		std::cout<<"   CASM constructed"<<std::endl;
	}

	//******************************CONSTRUCTOR*******************************************
	//************************************************************************************

	BorjaHenckyCasmPlasticAxisym2DLaw::BorjaHenckyCasmPlasticAxisym2DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
	{
		mpHardeningLaw    =  pHardeningLaw;
		mpYieldCriterion  =  YieldCriterion::Pointer( new CasmYieldCriterion(mpHardeningLaw) );
		mpFlowRule        =  pFlowRule;
		std::cout<<"THE BAD, BAD, BAD CASM CONSTRUCTOR IS CALLED HERE"<<std::endl<<std::endl<<std::endl;
	}

	//******************************COPY CONSTRUCTOR**************************************
	//************************************************************************************

	BorjaHenckyCasmPlasticAxisym2DLaw::BorjaHenckyCasmPlasticAxisym2DLaw(const BorjaHenckyCasmPlasticAxisym2DLaw& rOther)
		: NonLinearHenckyElasticPlasticAxisym2DLaw(rOther)
	{
	}

	//********************************CLONE***********************************************
	//************************************************************************************

	ConstitutiveLaw::Pointer BorjaHenckyCasmPlasticAxisym2DLaw::Clone() const
	{
		BorjaHenckyCasmPlasticAxisym2DLaw::Pointer p_clone(new BorjaHenckyCasmPlasticAxisym2DLaw(*this));
		return p_clone;
	}

	//*******************************DESTRUCTOR*******************************************
	//************************************************************************************

	BorjaHenckyCasmPlasticAxisym2DLaw::~BorjaHenckyCasmPlasticAxisym2DLaw()
	{
	}


	//*************************************** GET VALUE *********************************
	double&  BorjaHenckyCasmPlasticAxisym2DLaw::GetValue(const Variable<double>& rThisVariable, double & rValue)
	{
		if ( rThisVariable == M_MODULUS ) {
			double Swelling = mpYieldCriterion->GetHardeningLaw().GetProperties()[ SWELLING_SLOPE ];
			double MeanStress ;
			MeanStress = this->GetValue( STRESS_INV_P, MeanStress);
			double K = MeanStress/Swelling;

			double Alpha = mpYieldCriterion->GetHardeningLaw().GetProperties()[ ALPHA_SHEAR ];
			double G = mpYieldCriterion->GetHardeningLaw().GetProperties()[ INITIAL_SHEAR_MODULUS ];
			G += Alpha*MeanStress;  // this modulus is approximated (not the real one)

			rValue = K + 4.0*G / 3.0;
		}
		else if ( ( rThisVariable == YOUNG_MODULUS) || ( rThisVariable == EQUIVALENT_YOUNG_MODULUS ) || ( rThisVariable == SHEAR_MODULUS ) || (rThisVariable == BULK_MODULUS ) )
		{
			double Swelling = mpYieldCriterion->GetHardeningLaw().GetProperties()[ SWELLING_SLOPE ];
			double MeanStress ;
			MeanStress = this->GetValue( STRESS_INV_P, MeanStress);
			double K = MeanStress/Swelling;
			if ( rThisVariable == BULK_MODULUS) 
			{
				rValue = K;
				return rValue;
			}

			double Alpha = mpYieldCriterion->GetHardeningLaw().GetProperties()[ ALPHA_SHEAR ];
			double G = mpYieldCriterion->GetHardeningLaw().GetProperties()[ INITIAL_SHEAR_MODULUS ];
			//G += Alpha*MeanStress;  // this modulus is approximated (not the real one)
			Matrix ELCG = mElasticLeftCauchyGreen;
			Vector Hencky;
			//Vector Hencky = ConvertCauchyGreenTensorToHenckyStrain( ELCG);
			// instead for calling the function, I just copy here the function. :'P --
			{
				Matrix & rCauchyGreenMatrix = ELCG;
				Matrix EigenVectors;
				Vector EigenValues;
				SolidMechanicsMathUtilities<double>::EigenVectors(rCauchyGreenMatrix, EigenVectors, EigenValues);

				Matrix Aux = ZeroMatrix(3,3);
				for (unsigned int i = 0; i < 3; ++i)
					Aux(i,i) = (std::log(EigenValues(i)))/2.0;

				Aux = prod(Aux, (EigenVectors));
				Aux = prod(trans(EigenVectors), Aux);
				Vector Result = MathUtils<double>::StrainTensorToVector(Aux, 6);
				Hencky = Result;
			}


			double volum = 0.0;
			for (unsigned int i = 0; i < 3; i++)
				volum += Hencky(i);

			double ReferencePressure = mpYieldCriterion->GetHardeningLaw().GetProperties()[PRE_CONSOLIDATION_STRESS];
			double OCR = mpYieldCriterion->GetHardeningLaw().GetProperties()[OVER_CONSOLIDATION_RATIO];
			ReferencePressure /= OCR;    
			G += Alpha* ReferencePressure * std::exp( - volum / Swelling) ; 

			if ( rThisVariable == SHEAR_MODULUS)
			{
				rValue = G;
				return rValue;
			}

			rValue = 9.0*K*G / ( 3.0*K + G);

		}
		else {
			rValue = NonLinearHenckyElasticPlasticAxisym2DLaw::GetValue( rThisVariable, rValue);
		}

		return rValue;
	}

	void BorjaHenckyCasmPlasticAxisym2DLaw::SetValue(const Variable<Vector>& rThisVariable, const Vector& rValue, const ProcessInfo& rCurrentProcessInfo)
	{

		if ( rThisVariable == ELASTIC_LEFT_CAUCHY_FROM_KIRCHHOFF_STRESS)
		{
			Vector rStressVector = rValue;
			// Sets the ElasticLeftCauchyGreen tensor based on a Kirchhoff Stress
			// for the moment only works in the principal directions to set the initial state

			double SwellingSlope = mpYieldCriterion->GetHardeningLaw().GetProperties()[SWELLING_SLOPE];
			double AlphaShear = mpYieldCriterion->GetHardeningLaw().GetProperties()[ALPHA_SHEAR];

			double ReferencePressure = mpYieldCriterion->GetHardeningLaw().GetProperties()[PRE_CONSOLIDATION_STRESS];
			double OCR = mpYieldCriterion->GetHardeningLaw().GetProperties()[OVER_CONSOLIDATION_RATIO];
			double ConstantShearModulus = mpYieldCriterion->GetHardeningLaw().GetProperties()[INITIAL_SHEAR_MODULUS];
			ReferencePressure /= OCR;    

			Vector Objective(3);
			Objective(0) = rStressVector(0) + rStressVector(1) + rStressVector(2);
			Objective(0) /= 3.0;                     // mean stress
			Objective(1) = rStressVector(1) - Objective(0); // large deviatoric stress
			Objective(2) = rStressVector(2) - Objective(0); // small deviatoric stress (two directions)

			Vector Guess = ZeroVector(3);

			bool NotConverged = true;

			double DeviatoricNorm2 = 0.0, ShearModulus;
			Vector Y = ZeroVector(3);

			Matrix TangentMatrix = ZeroMatrix(3,3);
			Matrix InverseTangent = ZeroMatrix(3,3);
			Vector Residual = ZeroVector(3);
			Vector dGuess = Residual;
			double error, detI; 
			unsigned nIter = 0;

			while (NotConverged) {

				//1 COMPUTE SOME ERROR
				DeviatoricNorm2 = pow( Guess(1), 2) + 2.0*pow( Guess(2), 2);
				ShearModulus = AlphaShear*ReferencePressure*std::exp( -Guess(0) / SwellingSlope) + ConstantShearModulus;
				Y(0) = -ReferencePressure *std::exp( - Guess(0) / SwellingSlope) * ( 1.0 + AlphaShear/SwellingSlope * DeviatoricNorm2) ;
				Y(1) = 2.0*ShearModulus*Guess(1);
				Y(2) = 2.0*ShearModulus*Guess(2);


				Residual = Y - Objective;
				error = 0.0;
				for (unsigned int i = 0; i < 3; ++i)
				   error += pow( Residual(i), 2);

				if (error < 1.0e-10) {
				   NotConverged = false;
				}
				
				//1.1 Compute the Tangent Matrix
				TangentMatrix(0,0) = ReferencePressure*std::exp( -Guess(0)/SwellingSlope) * (1.0 + AlphaShear/SwellingSlope*DeviatoricNorm2) / SwellingSlope;
				TangentMatrix(0,1) = -ReferencePressure*std::exp( -Guess(0)/SwellingSlope) * ( AlphaShear/SwellingSlope) * 2.0 *  Guess(1) ;
				TangentMatrix(0,2) = -ReferencePressure*std::exp( -Guess(0)/SwellingSlope) * ( AlphaShear/SwellingSlope) * 2.0 *  Guess(2) * 2.0 ;


				TangentMatrix(1,0) =  2.0* AlphaShear * ReferencePressure*std::exp(-Guess(0) / SwellingSlope)*(-1.0/SwellingSlope)*Guess(1);
				TangentMatrix(1,1) = 2.0*ShearModulus;
				TangentMatrix(1,2) = 0.0;

				TangentMatrix(2,0) =  2.0* AlphaShear * ReferencePressure*std::exp(-Guess(0) / SwellingSlope)*(-1.0/SwellingSlope)*Guess(2);
				TangentMatrix(2,1) = 0.0;
				TangentMatrix(2,2) = 2.0*ShearModulus;
				MathUtils<double>::InvertMatrix( TangentMatrix, InverseTangent, detI);

				dGuess = prod( InverseTangent, Residual);

				// (Try to solve some convergence problems)
				double dGNorm = 0.0;
				for (unsigned int i = 0; i <3 ; ++i)
				   dGNorm += pow( dGuess(i), 2);
				dGNorm = sqrt( dGNorm);
				if ( dGNorm > 0.01)
				   dGuess *= 0.1;


				Guess -= dGuess;
				nIter += 1;
				if (nIter > 1000) {
				   std::cout << " NONCONVERGING " << std::endl;
				   return;
				}
			}

			// 2. Set the ElasticLeftCauchy
			double Hencky1 = 0.0;
			double Hencky2 = 0.0; 

			Hencky1 = Guess(1) + Guess(0)/3.0;
			Hencky2 = Guess(2) + Guess(0)/3.0; // small strain 

			mElasticLeftCauchyGreen(0,0) = std::exp( 2.0*Hencky2);
			mElasticLeftCauchyGreen(1,1) = std::exp( 2.0*Hencky1);
			mElasticLeftCauchyGreen(2,2) = mElasticLeftCauchyGreen(0,0);

			// 3. Try to be inside the yield surface (this only make sense in the Initial State)

			const double ShearM = mpYieldCriterion->GetHardeningLaw().GetProperties()[CRITICAL_STATE_LINE];
			const double OtherSlope = mpYieldCriterion->GetHardeningLaw().GetProperties()[NORMAL_COMPRESSION_SLOPE];
			double PreconsolidationStress;

			double StressQ = 0.0;
			StressQ = pow( Objective(1), 2) + 2.0*pow( Objective(2), 2);
			StressQ = sqrt(1.5) * sqrt(StressQ);

			PreconsolidationStress = Objective(0) + pow ( StressQ / ShearM, 2) / Objective(0); 
			ReferencePressure *= OCR;

			PreconsolidationStress *= OCR;
			PreconsolidationStress = rStressVector(1)*OCR;
			if ( PreconsolidationStress > -5.0 ) {
				PreconsolidationStress = -5.0; // a treure;
			}
			else {
			}
			double VolumetricPlasticDef = - (OtherSlope - SwellingSlope) * std::log(-PreconsolidationStress / ReferencePressure);
			FlowRule::RadialReturnVariables ReturnMappingVariables;
			ReturnMappingVariables.DeltaGamma = VolumetricPlasticDef;
			ReturnMappingVariables.DeltaTime = 1.0;
			ReturnMappingVariables.DeltaBeta = 0.0;
			mpFlowRule->UpdateInternalVariables( ReturnMappingVariables);

			double ThisPrecon = 0.0;
			mpHardeningLaw->CalculateHardening( ThisPrecon, ReturnMappingVariables.DeltaGamma );
		}
		else
		{
			NonLinearHenckyElasticPlasticAxisym2DLaw::SetValue( rThisVariable, rValue, rCurrentProcessInfo );
		}
	}
   
	int BorjaHenckyCasmPlasticAxisym2DLaw::Check( const Properties& rMaterialProperties, const GeometryType& rElementGeometry, const ProcessInfo& rCurrentProcessInfo)
	{
		return 0;
	}


} // Namespace Kratos
