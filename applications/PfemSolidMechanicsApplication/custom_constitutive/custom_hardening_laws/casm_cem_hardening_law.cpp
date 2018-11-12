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
#include "includes/properties.h"
#include "custom_constitutive/custom_hardening_laws/casm_cem_hardening_law.hpp"
#include "pfem_solid_mechanics_application_variables.h"


namespace Kratos
{

	//*******************************CONSTRUCTOR******************************************
	//************************************************************************************

	CasmCemHardeningLaw::CasmCemHardeningLaw()
		:HardeningLaw()
	{
		std::cout<<"   CASM CEM HARDENING LAW constructed"<<std::endl;
	}


	//*******************************ASSIGMENT OPERATOR***********************************
	//************************************************************************************

	CasmCemHardeningLaw& CasmCemHardeningLaw::operator=(CasmCemHardeningLaw const& rOther)
	{
		HardeningLaw::operator=(rOther);
		return *this;
	}

	//*******************************COPY CONSTRUCTOR*************************************
	//************************************************************************************

	CasmCemHardeningLaw::CasmCemHardeningLaw(CasmCemHardeningLaw const& rOther)
		:HardeningLaw(rOther)
	{

	}


	//********************************DESTRUCTOR******************************************
	//************************************************************************************

	CasmCemHardeningLaw::~CasmCemHardeningLaw()
	{
		
	}

	/// Operations.

	//*******************************CALCULATE TOTAL HARDENING****************************
	//************************************************************************************
	void CasmCemHardeningLaw::CalculateHardening(PlasticVariablesType& rPlasticVariables, const double& rDeltaAlpha, const double& rDeltaBeta)
	{
		const double SwellingSlope 									= GetProperties()[SWELLING_SLOPE];
		const double OtherSlope 										= GetProperties()[NORMAL_COMPRESSION_SLOPE];
		const double Omega 													= GetProperties()[PLASTIC_DEVIATORIC_STRAIN_HARDENING];
		const double H1 														= GetProperties()[DEGRADATION_RATE_COMPRESSION];
		const double H2 														= GetProperties()[DEGRADATION_RATE_SHEAR];

		//update hardening parameter
		rPlasticVariables.PreconsolidationPressure += rPlasticVariables.PreconsolidationPressure/(OtherSlope-SwellingSlope)*(-rDeltaAlpha + Omega*rDeltaBeta );
		rPlasticVariables.Bonding += -rPlasticVariables.Bonding*( H1*fabs(rDeltaAlpha) + H2*fabs(rDeltaBeta) );
	}


	Vector& CasmCemHardeningLaw::CalculateHardening(Vector& rHardening, const double& rAlpha, const double& rBeta, const double& rAlphaCum, const double& rBetaCum, const double rTemperature)
	{
		//rAlpha ... inc vol strain, rBeta ... inc dev strain
		const double FirstPreconsolidationPressure 	= GetProperties()[PRE_CONSOLIDATION_STRESS];
		const double SwellingSlope 									= GetProperties()[SWELLING_SLOPE];
		const double OtherSlope 										= GetProperties()[NORMAL_COMPRESSION_SLOPE];
		const double Omega 													= GetProperties()[PLASTIC_DEVIATORIC_STRAIN_HARDENING];
		const double B0 														= GetProperties()[INITIAL_BONDING];
		const double H0 														= GetProperties()[DEGRADATION_THRESHOLD];
		const double H1 														= GetProperties()[DEGRADATION_RATE_COMPRESSION];
		const double H2 														= GetProperties()[DEGRADATION_RATE_SHEAR];

		
		rHardening(0) = -FirstPreconsolidationPressure*(std::exp ((-rAlpha + Omega*rBeta)/(OtherSlope-SwellingSlope)) );
		rHardening(1) = B0*(std::exp ( H0 - (H1*rAlphaCum + H2*rBetaCum) ) );
		
		return rHardening;
	}


	double& CasmCemHardeningLaw::CalculateHardening(double& rHardening, const double& rAlpha, const double& rBeta, const double& rAlphaCum, const double& rBetaCum, const double rTemperature)
	{
		Vector aux = ZeroVector(2);
		aux = CalculateHardening(aux, rAlpha, rBeta, rAlphaCum, rBetaCum);
		rHardening = aux(0);
		return rHardening;
	}
	
	double& CasmCemHardeningLaw::CalculateBonding(double& rHardening, const double& rAlpha, const double& rBeta, const double& rAlphaCum, const double& rBetaCum, const double rTemperature)
	{
		Vector aux = ZeroVector(2);
		aux = CalculateHardening(aux, rAlpha, rBeta, rAlphaCum, rBetaCum);
		rHardening = aux(1);
		return rHardening;
	}


	void CasmCemHardeningLaw::save( Serializer& rSerializer ) const
	{
		KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, HardeningLaw )
	}

	void CasmCemHardeningLaw::load( Serializer& rSerializer )
	{
		KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, HardeningLaw )
	}

}  // namespace Kratos.
