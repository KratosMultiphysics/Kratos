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
#include "includes/define.h"
#include "custom_utilities/stress_invariants_utilities.hpp"
#include "custom_constitutive/custom_yield_criteria/casm_cem_yield_criterion.hpp"

#include "pfem_solid_mechanics_application_variables.h"


namespace Kratos
{

	//*******************************CONSTRUCTOR******************************************
	//************************************************************************************
	CasmCemYieldCriterion::CasmCemYieldCriterion()
		:YieldCriterion()
	{
	}


	//*****************************INITIALIZATION CONSTRUCTOR*****************************
	//************************************************************************************
	CasmCemYieldCriterion::CasmCemYieldCriterion(HardeningLawPointer pHardeningLaw)
		:YieldCriterion(pHardeningLaw)
	{
		std::cout<<"   CASM CEM YIELD CRITERION constructed"<<std::endl;
	}


	//*******************************ASSIGMENT OPERATOR***********************************
	//************************************************************************************
	CasmCemYieldCriterion& CasmCemYieldCriterion::operator=(CasmCemYieldCriterion const& rOther)
	{
		YieldCriterion::operator=(rOther);
		return *this;
	}


	//*******************************COPY CONSTRUCTOR*************************************
	//************************************************************************************
	CasmCemYieldCriterion::CasmCemYieldCriterion(CasmCemYieldCriterion const& rOther)
		:YieldCriterion(rOther)
	{
	}


	//********************************DESTRUCTOR******************************************
	//************************************************************************************
	CasmCemYieldCriterion::~CasmCemYieldCriterion()
	{
	}


	//************************* CALCULATE YIELD FUNCTION  ******************
	//**********************************************************************
	double& CasmCemYieldCriterion::CalculateYieldCondition(double& rStateFunction, const Vector& rStressVector, const PlasticVariablesType& rPlasticVariables)
	{
		// calculate Kirchhoff invariants
		double MeanStress, LodeAngle;
		double DeviatoricQ; // == sqrt(3)*J2
		StressInvariantsUtilities::CalculateStressInvariants( rStressVector, MeanStress, DeviatoricQ, LodeAngle);
		DeviatoricQ *= sqrt(3.0);

		// slope of CS-line
		double ThirdInvariantEffect = EvaluateThirdInvariantEffectMC(LodeAngle);
		
		// get material constants
		const double ShearM 					= this->GetHardeningLaw().GetProperties()[CRITICAL_STATE_LINE];
		const double SpacingR 				= this->GetHardeningLaw().GetProperties()[SPACING_RATIO];
		const double ShapeN 					= this->GetHardeningLaw().GetProperties()[SHAPE_PARAMETER];
		const double AlphaTensile 		= this->GetHardeningLaw().GetProperties()[ALPHA_TENSILE];
		
		// calculate hardening parameters p0, b -> pt, pc
		double Pc, Pt;
		Pc = rPlasticVariables.PreconsolidationPressure*(1+rPlasticVariables.Bonding);
		Pt = rPlasticVariables.PreconsolidationPressure*(AlphaTensile*rPlasticVariables.Bonding);

		// evaluate yield function
		rStateFunction = pow(-DeviatoricQ/(ShearM/ThirdInvariantEffect*(MeanStress + Pt)), ShapeN );
		rStateFunction += 1/log(SpacingR)*log((MeanStress + Pt)/(Pc + Pt));
			/*
			std::cout<<"  CasmCemYieldCriterion::CalculateYieldCondition"<<std::endl;
			std::cout<<"   p:    "<<MeanStress<<"   J: "<<DeviatoricQ<<"   Lode: "<<LodeAngle<<std::endl;
			std::cout<<"   M:    "<<ShearM<<" Eff: "<<ThirdInvariantEffect<<std::endl;
			std::cout<<"   Sig:  "<<rStressVector<<std::endl;
			std::cout<<"   p0:   "<<HardeningVariables(0)<<std::endl;
			std::cout<<"   b:    "<<HardeningVariables(1)<<std::endl;
			std::cout<<"   Pc:   "<<Pc<<std::endl;
			std::cout<<"   Pt:   "<<Pt<<std::endl;
			std::cout<<"   Eff:  "<<ThirdInvariantEffect<<std::endl;
			std::cout<<"   f(sig): "<<rStateFunction<<std::endl<<std::endl;
			*/
		return rStateFunction; 
	}


	// ************************************************
	// ***** Calculate yield function derivatives *****
	// ************************************************
	void CasmCemYieldCriterion::CalculateYieldFunctionDerivative(const Vector& rStressVector, Vector& rYieldFunctionD, const PlasticVariablesType& rPlasticVariables)
	{
		// Kirchhoff stress invariants & invariants derivatives
		double MeanStress, J2, LodeAngle;
		Vector V1, V2;
		StressInvariantsUtilities::CalculateStressInvariants( rStressVector, MeanStress, J2, LodeAngle);
		StressInvariantsUtilities::CalculateDerivativeVectors( rStressVector, V1, V2);
		
		// slope CS-line
		const double ShearM = this->GetHardeningLaw().GetProperties()[CRITICAL_STATE_LINE];
		double ThirdInvariantEffect = EvaluateThirdInvariantEffectMC( LodeAngle);
		
		// get material constants
		const double SpacingR = this->GetHardeningLaw().GetProperties()[SPACING_RATIO];
		const double ShapeN = this->GetHardeningLaw().GetProperties()[SHAPE_PARAMETER];
		const double AlphaTensile 		= this->GetHardeningLaw().GetProperties()[ALPHA_TENSILE];
		
		// calculate hardening parameters p0, b -> pt, pc
		double Pt = rPlasticVariables.PreconsolidationPressure*(AlphaTensile*rPlasticVariables.Bonding);

		// calculate d_f/d_Sig = d_f/d_Inv * d_Inv/d_Sig 
		rYieldFunctionD = ( 1/( (MeanStress+Pt) * log(SpacingR) ) + ( ShapeN * pow( pow(3,1/2)*J2 , ShapeN) )/( pow(ShearM/ThirdInvariantEffect,ShapeN) * pow(-(MeanStress+Pt),ShapeN+1) ) ) * V1;
		rYieldFunctionD += ( ( ShapeN * pow(3,ShapeN/2) * pow(J2,ShapeN-1) )/( pow(ShearM/ThirdInvariantEffect,ShapeN) * pow(-(MeanStress+Pt),ShapeN) ) ) * V2;
		
		CalculateAndAddThirdInvDerivativeMC( rStressVector, rYieldFunctionD, Pt);
			/*
			std::cout<<"  CasmYieldCriterion::CalculateYieldFunctionDerivative"<<std::endl;
			std::cout<<"   p: "<<MeanStress<<" J: "<<J2<<" Lode: "<<LodeAngle<<std::endl;
			std::cout<<"   M: "<<ShearM<<" Eff: "<<ThirdInvariantEffect<<std::endl;
			std::cout<<"   Sig: "<<rStressVector<<std::endl;
			std::cout<<"   d_f/d_sig: "<<rYieldFunctionD<<std::endl;
			std::cout<<"   Pt: "<<Pt<<std::endl;
			*/
	}


	//******************************Evaluate Effect of Third Invariant *******************
	//************************************************************************************
	// NOT USED, CHECK BEFORE USING!!
	double CasmCemYieldCriterion::EvaluateThirdInvariantEffectSheng( const double& rLodeAngle)
	{
		// get phi_CS and catch phi of zero
		double Friction = this->GetHardeningLaw().GetProperties()[INTERNAL_FRICTION_ANGLE];
		if (Friction < 1.0E-3) {
			return 1.0;
		}
		Friction *= GetPI() / 180.0;
		
		// calculate 
		double auxAlpha = (3-std::sin(Friction))/(3+std::sin(Friction));
		if (auxAlpha > 1)
			return 1.0;
		else if (auxAlpha < 0.6)
			auxAlpha = 0.6;
		
		// influence of Lode angle according to Sheng et al. (2000)
		double Effect = 1.0;
		Effect = pow( 2*pow(auxAlpha,4) / ( 1+pow(auxAlpha,4)+(1-pow(auxAlpha,4))*std::sin(-3*rLodeAngle) ) , 1/4);

		return Effect;
	}


	double CasmCemYieldCriterion::EvaluateThirdInvariantEffectMC( const double& rLodeAngle)
	{

		double Effect = 1.0;
		double Friction = this->GetHardeningLaw().GetProperties()[INTERNAL_FRICTION_ANGLE];
		if (Friction < 1.0E-3) {
			return 1.0;
		}
		Friction *= GetPI() / 180.0;

		double LodeCut = GetSmoothingLodeAngle();

		if ( fabs( rLodeAngle)  < LodeCut) {
			Effect = std::cos( rLodeAngle) - 1.0/sqrt(3.0) * std::sin(Friction) * std::sin(rLodeAngle); 
		}
		else {

			double A, B;
			GetSmoothingConstants(A, B, rLodeAngle);
			Effect = A + B*std::sin(3.0*rLodeAngle);
		}

		// influence of Lode angle according to Abbo et al. (2011)
		Effect /= ( sqrt(3)/6) * (3.0 - std::sin(Friction) );
		
		//std::cout << "LODE: " << rLodeAngle << " " << rLodeAngle * 180.0 / GetPI() << "  " << "Effect" << " " << Effect << std::endl;
		
		return Effect;
		//return 1;
  }


	//*******************************Add  derivative of Effect Third Invariant ********************
	//*********************************************************************************************
	void CasmCemYieldCriterion::CalculateAndAddThirdInvDerivativeMC(const Vector& rStressVector, Vector& rYieldFunctionD, const double& rPt)
	{
		double KLode = 1.0;
		double Friction = this->GetHardeningLaw().GetProperties()[INTERNAL_FRICTION_ANGLE];
		if (Friction < 1.0E-3) {
			return ;
		}
		Friction *= GetPI() / 180.0;

		double MeanStress, J2, LodeAngle;
		Vector V1, V2, V3;

		StressInvariantsUtilities::CalculateStressInvariants( rStressVector, MeanStress, J2, LodeAngle);
		// since I will be dividing by J2
		if ( J2 < 1E-5)
			return;

		StressInvariantsUtilities::CalculateDerivativeVectors( rStressVector, V1, V2, V3);
		double C2, C3;
		double KLodeDeriv;
		const double ShearM = this->GetHardeningLaw().GetProperties()[CRITICAL_STATE_LINE];
		const double ShapeN = this->GetHardeningLaw().GetProperties()[SHAPE_PARAMETER];

		// calcualte K(Lode) and d_K/d_Lode
		double LodeCut = GetSmoothingLodeAngle();
		if ( fabs(LodeAngle)  < LodeCut) {

			KLode = std::cos(LodeAngle) - 1.0/sqrt(3.0) * std::sin(Friction) * std::sin(LodeAngle); 
			KLodeDeriv = -std::sin(LodeAngle) - 1.0/sqrt(3.0) * std::sin(Friction) * std::cos(LodeAngle);
		}
		else {

			double A, B;
			GetSmoothingConstants(A, B, LodeAngle);
			
			KLode = A + B * std::sin(3.0*LodeAngle);
			KLodeDeriv = 3.0 * B * std::cos(3.0*LodeAngle);
		}
		
		// calcualte the factors C2 and C3
		C2 = -std::tan(3.0*LodeAngle) * ShapeN * pow(6.0, ShapeN) * pow(J2, ShapeN-1) * pow(-(MeanStress+rPt) * ShearM * (3.0-std::sin(Friction)), -ShapeN);
		C2 *= pow(KLode, ShapeN-1) * KLodeDeriv;
			
		C3 = -ShapeN * pow(6.0, ShapeN) * sqrt(3.0) * pow(J2, ShapeN-3) * pow(-(MeanStress+rPt) * ShearM * (3.0-std::sin(Friction)), -ShapeN);
		C3 /= (2.0 * std::cos(3.0*LodeAngle));
		C3 *= pow(KLode, ShapeN-1) * KLodeDeriv;
		
		// d_f/d_Lode*d_Lode/d_sig = C2*d_J/d_sig + C3*d_J3/d_sig
		Vector ThisDerivative = C2 * V2 + C3 * V3;

      /*std::cout << " LODE " << LodeAngle <<" LODE " << LodeAngle * 180.0 / GetPI() << " EFFECT " << Effect <<  " DERIVATIVE " << EffectDeriv * std::cos( 3.0 * LodeAngle)  << std::endl;
        std::cout << " PREVIOUS ANAL DERIVATIVE " << rYieldFunctionD << std::endl;
        std::cout << " AND THIS NEW C2 " << C2 << " and C3 " << C3 << std::endl;
        std::cout << " V2 " << V2 << " V3 " << V3 << std::endl;*/
       
		rYieldFunctionD += ThisDerivative;


	}

	//************************* smoothingInvariants of something ***********
	//**********************************************************************

	void CasmCemYieldCriterion::GetSmoothingConstants(double& rA, double& rB, const double& rLodeAngle)
	{
		// based on C1-continuous smoothing of MC yield surface (Abbo et al 2011)
		double SmoothingAngle = this->GetSmoothingLodeAngle();
		double FrictionAngle = this->GetHardeningLaw().GetProperties()[INTERNAL_FRICTION_ANGLE];
		FrictionAngle *= GetPI() / 180.0;

		double Sign = 1.0;
		if ( rLodeAngle < 0.0)
			Sign = -1.0;

		rA = 3.0 +  std::tan(SmoothingAngle) * std::tan(3.0*SmoothingAngle) + Sign * (std::tan( 3.0*SmoothingAngle) - 3.0*std::tan(SmoothingAngle)) * std::sin( FrictionAngle) / sqrt(3.0);
		rA *= (1.0/3.0) * std::cos( SmoothingAngle );

		rB = -1.0 * ( Sign* std::sin(SmoothingAngle) + std::sin(FrictionAngle)*std::cos(SmoothingAngle) / sqrt(3.0) ) / ( 3.0*std::cos(3.0*SmoothingAngle) );
	}
	
	// ************************************************
	// ***** Calculate yield function derivatives *****
	// ************************************************
/*	void CasmCemYieldCriterion::CalculateYieldFunctionDerivative(const Vector& rStressVector, Vector& rYieldFunctionD, const double& rAlpha, const double& rBeta, const double& rAlphaCum, const double& rBetaCum)
	{
		// Kirchhoff stress invariants & invariants derivatives
		double MeanStress, J2, LodeAngle;
		Vector V1, V2;
		StressInvariantsUtilities::CalculateStressInvariants( rStressVector, MeanStress, J2, LodeAngle);
		StressInvariantsUtilities::CalculateDerivativeVectors( rStressVector, V1, V2);
		
		// slope CS-line
		const double ShearM = this->GetHardeningLaw().GetProperties()[CRITICAL_STATE_LINE];
		double ThirdInvariantEffect = EvaluateThirdInvariantEffectMC( LodeAngle);
		
		// get material constants
		const double SpacingR = this->GetHardeningLaw().GetProperties()[SPACING_RATIO];
		const double ShapeN = this->GetHardeningLaw().GetProperties()[SHAPE_PARAMETER];
		const double AlphaTensile 		= this->GetHardeningLaw().GetProperties()[ALPHA_TENSILE];
		
		// calculate hardening parameters p0, b -> pt, pc
		Vector HardeningVariables = ZeroVector(2);
		HardeningVariables = mpHardeningLaw->CalculateHardening(HardeningVariables, rAlpha, rBeta, rAlphaCum, rBetaCum);
		double Pt;//, Pc;
		//Pc = HardeningVariables(0)*(1+HardeningVariables(1));
		Pt = HardeningVariables(0)*(AlphaTensile*HardeningVariables(1));

		// calculate d_f/d_Sig = d_f/d_Inv * d_Inv/d_Sig 
		rYieldFunctionD = ( 1/( (MeanStress+Pt) * log(SpacingR) ) + ( ShapeN * pow( pow(3,1/2)*J2 , ShapeN) )/( pow(ShearM/ThirdInvariantEffect,ShapeN) * pow(-(MeanStress+Pt),ShapeN+1) ) ) * V1;
		rYieldFunctionD += ( ( ShapeN * pow(3,ShapeN/2) * pow(J2,ShapeN-1) )/( pow(ShearM/ThirdInvariantEffect,ShapeN) * pow(-(MeanStress+Pt),ShapeN) ) ) * V2;
		
		CalculateAndAddThirdInvDerivativeMC( rStressVector, rYieldFunctionD, Pt);
			*
			std::cout<<"  CasmYieldCriterion::CalculateYieldFunctionDerivative"<<std::endl;
			std::cout<<"   p: "<<MeanStress<<" J: "<<J2<<" Lode: "<<LodeAngle<<std::endl;
			std::cout<<"   M: "<<ShearM<<" Eff: "<<ThirdInvariantEffect<<std::endl;
			std::cout<<"   Sig: "<<rStressVector<<std::endl;
			std::cout<<"   d_f/d_sig: "<<rYieldFunctionD<<std::endl;
			std::cout<<"   Pt: "<<Pt<<std::endl;
			*
	}
*/

/*
	// ************************* CALCULATE YIELD FUNCTION  ******************
	// **********************************************************************
	double& CasmCemYieldCriterion::CalculateYieldCondition(double& rStateFunction, const Vector& rStressVector, const double& rAlpha, const double& rBeta, const double& rAlphaCum, const double& rBetaCum)
	{
		// calculate Kirchhoff invariants
		double MeanStress, LodeAngle;
		double DeviatoricQ; // == sqrt(3)*J2
		StressInvariantsUtilities::CalculateStressInvariants( rStressVector, MeanStress, DeviatoricQ, LodeAngle);
		DeviatoricQ *= sqrt(3.0);

		// slope of CS-line
		double ThirdInvariantEffect = EvaluateThirdInvariantEffectMC(LodeAngle);
		
		// get material constants
		const double ShearM 					= this->GetHardeningLaw().GetProperties()[CRITICAL_STATE_LINE];
		const double SpacingR 				= this->GetHardeningLaw().GetProperties()[SPACING_RATIO];
		const double ShapeN 					= this->GetHardeningLaw().GetProperties()[SHAPE_PARAMETER];
		const double AlphaTensile 		= this->GetHardeningLaw().GetProperties()[ALPHA_TENSILE];
		
		// calculate hardening parameters p0, b -> pt, pc
		Vector HardeningVariables = ZeroVector(2);
		HardeningVariables = mpHardeningLaw->CalculateHardening(HardeningVariables, rAlpha, rBeta, rAlphaCum, rBetaCum);
		double Pc, Pt;
		Pc = HardeningVariables(0)*(1+HardeningVariables(1));
		Pt = HardeningVariables(0)*(AlphaTensile*HardeningVariables(1));

		// evaluate yield function
		rStateFunction = pow(-DeviatoricQ/(ShearM/ThirdInvariantEffect*(MeanStress + Pt)), ShapeN );
		rStateFunction += 1/log(SpacingR)*log((MeanStress + Pt)/(Pc + Pt));
			*
			std::cout<<"  CasmCemYieldCriterion::CalculateYieldCondition"<<std::endl;
			std::cout<<"   p:    "<<MeanStress<<"   J: "<<DeviatoricQ<<"   Lode: "<<LodeAngle<<std::endl;
			std::cout<<"   M:    "<<ShearM<<" Eff: "<<ThirdInvariantEffect<<std::endl;
			std::cout<<"   Sig:  "<<rStressVector<<std::endl;
			std::cout<<"   p0:   "<<HardeningVariables(0)<<std::endl;
			std::cout<<"   b:    "<<HardeningVariables(1)<<std::endl;
			std::cout<<"   Pc:   "<<Pc<<std::endl;
			std::cout<<"   Pt:   "<<Pt<<std::endl;
			std::cout<<"   Eff:  "<<ThirdInvariantEffect<<std::endl;
			std::cout<<"   f(sig): "<<rStateFunction<<std::endl<<std::endl;
			*
		return rStateFunction; 
	}
*/	

	double CasmCemYieldCriterion::GetSmoothingLodeAngle()
	{
		return 27.0*GetPI()/180.0;
	}

	double CasmCemYieldCriterion::GetPI()
	{
		return 3.14159265359;
	}

	void CasmCemYieldCriterion::save( Serializer& rSerializer ) const
	{
		KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, YieldCriterion )
	}

	void CasmCemYieldCriterion::load( Serializer& rSerializer )
	{
		KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, YieldCriterion )
	}


}  // namespace Kratos.
