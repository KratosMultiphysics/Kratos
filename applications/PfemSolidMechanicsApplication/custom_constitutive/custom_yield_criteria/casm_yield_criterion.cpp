//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Created by:          $Author:                    LHauser $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                    June 2018 $
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
#include "custom_constitutive/custom_yield_criteria/casm_yield_criterion.hpp"

#include "pfem_solid_mechanics_application_variables.h"


namespace Kratos
{


	//*******************************CONSTRUCTOR******************************************
	//************************************************************************************
	CasmYieldCriterion::CasmYieldCriterion()
		:YieldCriterion()
	{

	}

	//*****************************INITIALIZATION CONSTRUCTOR*****************************
	//************************************************************************************

	CasmYieldCriterion::CasmYieldCriterion(HardeningLawPointer pHardeningLaw)
		:YieldCriterion(pHardeningLaw)
	{
		std::cout<<"   CASM YIELD CRITERION constructed"<<std::endl;
	}


	//*******************************ASSIGMENT OPERATOR***********************************
	//************************************************************************************

	CasmYieldCriterion& CasmYieldCriterion::operator=(CasmYieldCriterion const& rOther)
	{
		YieldCriterion::operator=(rOther);
		return *this;
	}

	//*******************************COPY CONSTRUCTOR*************************************
	//************************************************************************************

	CasmYieldCriterion::CasmYieldCriterion(CasmYieldCriterion const& rOther)
		:YieldCriterion(rOther)
	{

	}


	//********************************DESTRUCTOR******************************************
	//************************************************************************************

	CasmYieldCriterion::~CasmYieldCriterion()
	{
		
	}


	//************************* CALCULATE YIELD FUNCTION  ******************
	//**********************************************************************

	double& CasmYieldCriterion::CalculateYieldCondition(double& rStateFunction, const Vector& rStressVector, const double& rAlpha)
	{
		// calculate Kirchhoff invariants
		double MeanStress, LodeAngle;
		double DeviatoricQ; // == sqrt(3)*J2

		StressInvariantsUtilities::CalculateStressInvariants( rStressVector, MeanStress, DeviatoricQ, LodeAngle);
		DeviatoricQ *= sqrt(3.0);

		// slope of CS-line
		const double ShearM = this->GetHardeningLaw().GetProperties()[CRITICAL_STATE_LINE]; // for Lode = -30°
		double ThirdInvariantEffect = EvaluateThirdInvariantEffectMC(LodeAngle);
		
		// get spacing rario r & shape parameter n
		const double SpacingR = this->GetHardeningLaw().GetProperties()[SPACING_RATIO];
		const double ShapeN = this->GetHardeningLaw().GetProperties()[SHAPE_PARAMETER];

		// calculate hardening consolidation pressure (Kirchoff) with plastic volumetric strain increment rAlpha
		// h = h(F^p * F^pT)
		double PreconsolidationStress = 0.0;
		PreconsolidationStress = mpHardeningLaw->CalculateHardening(PreconsolidationStress, rAlpha);		
		// evaluate yield function
		rStateFunction = pow(-DeviatoricQ/(ShearM/ThirdInvariantEffect*(MeanStress)), ShapeN );
		rStateFunction += 1/log(SpacingR)*log(MeanStress/PreconsolidationStress);
			/*
			std::cout<<"  CasmYieldCriterion::CalculateYieldCondition"<<std::endl;
			std::cout<<"   p:    "<<MeanStress<<"   J: "<<DeviatoricQ<<"   Lode: "<<LodeAngle<<std::endl;
			std::cout<<"   M:    "<<ShearM<<" Eff: "<<ThirdInvariantEffect<<std::endl;
			std::cout<<"   Sig:  "<<rStressVector<<std::endl;
			std::cout<<"   p0:   "<<PreconsolidationStress<<std::endl;
			std::cout<<"   Eff:  "<<ThirdInvariantEffect<<std::endl;
			std::cout<<"   f(sig): "<<rStateFunction<<std::endl<<std::endl;
			*/
		return rStateFunction; 
	}


	//*******************************CALCULATE YIELD FUNCTION DERIVATIVE *****************
	//************************************************************************************
	void CasmYieldCriterion::CalculateYieldFunctionDerivative(const Vector& rStressVector, Vector& rYieldFunctionD, const double& rAlpha)
	{
		double PreconsolidationStress = 0.0;
		PreconsolidationStress = mpHardeningLaw->CalculateHardening(PreconsolidationStress, rAlpha);
		const double ShearM = this->GetHardeningLaw().GetProperties()[CRITICAL_STATE_LINE];
		const double SpacingR = this->GetHardeningLaw().GetProperties()[SPACING_RATIO];
		const double ShapeN = this->GetHardeningLaw().GetProperties()[SHAPE_PARAMETER];
		
		// stress invariants & invariants derivatives
		double MeanStress, J2, LodeAngle;
		Vector V1, V2;
		StressInvariantsUtilities::CalculateStressInvariants( rStressVector, MeanStress, J2, LodeAngle);
		StressInvariantsUtilities::CalculateDerivativeVectors( rStressVector, V1, V2);

		double ThirdInvariantEffect = EvaluateThirdInvariantEffectMC( LodeAngle);

		// calculate d_f/d_Sig = d_f/d_Inv * d_Inv/d_Sig 
		rYieldFunctionD = ( 1/( MeanStress * log(SpacingR) ) + ( ShapeN * pow( pow(3,1/2)*J2 , ShapeN) )/( pow(ShearM/ThirdInvariantEffect,ShapeN) * pow(-MeanStress,ShapeN+1) ) ) * V1;
		rYieldFunctionD += ( ( ShapeN * pow(3,ShapeN/2) * pow(J2,ShapeN-1) )/( pow(ShearM/ThirdInvariantEffect,ShapeN) * pow(-MeanStress,ShapeN) ) ) * V2;

		CalculateAndAddThirdInvDerivativeMC( rStressVector, rYieldFunctionD);
			/*
			std::cout<<"  CasmYieldCriterion::CalculateYieldFunctionDerivative"<<std::endl;
			std::cout<<"   p: "<<MeanStress<<" J: "<<J2<<" Lode: "<<LodeAngle<<std::endl;
			std::cout<<"   M: "<<ShearM<<" Eff: "<<ThirdInvariantEffect<<std::endl;
			std::cout<<"   Sig: "<<rStressVector<<std::endl;
			std::cout<<"   d_f/d_sig: "<<rYieldFunctionD<<std::endl;
			*/
	}

	//******************************Evaluate Effect of Third Invariant *******************
	//************************************************************************************
	double CasmYieldCriterion::EvaluateThirdInvariantEffectSheng( const double& rLodeAngle)
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
   
	double CasmYieldCriterion::EvaluateThirdInvariantEffectMC( const double& rLodeAngle)
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

  }


	//*******************************Add  derivative of Effect THird Invariant ********************
	//*********************************************************************************************
	void CasmYieldCriterion::CalculateAndAddThirdInvDerivativeMC(const Vector& rStressVector, Vector& rYieldFunctionD)
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
		C2 = -std::tan(3.0*LodeAngle) * ShapeN * pow(6.0, ShapeN) * pow(J2, ShapeN-1) * pow(-MeanStress * ShearM * (3.0-std::sin(Friction)), -ShapeN);
		C2 *= pow(KLode, ShapeN-1) * KLodeDeriv;
			
		C3 = -ShapeN * pow(6.0, ShapeN) * sqrt(3.0) * pow(J2, ShapeN-3) * pow(-MeanStress * ShearM * (3.0-std::sin(Friction)), -ShapeN);
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

	void CasmYieldCriterion::GetSmoothingConstants(double& rA, double& rB, const double& rLodeAngle)
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

	double CasmYieldCriterion::GetSmoothingLodeAngle()
	{
		return 27.0*GetPI()/180.0;
	}

	double CasmYieldCriterion::GetPI()
	{
		return 3.14159265359;
	}

	void CasmYieldCriterion::save( Serializer& rSerializer ) const
	{
		KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, YieldCriterion )
	}

	void CasmYieldCriterion::load( Serializer& rSerializer )
	{
		KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, YieldCriterion )
	}


}  // namespace Kratos.
