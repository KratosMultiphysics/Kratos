//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Created by:          $Author:                    LHauser $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                     Dec 2018 $
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
#include "custom_constitutive/custom_hardening_laws/casm_hardening_law.hpp"
#include "pfem_solid_mechanics_application_variables.h"


namespace Kratos
{

	//*******************************CONSTRUCTOR******************************************
	//************************************************************************************

	CasmHardeningLaw::CasmHardeningLaw()
		:HardeningLaw()
	{
		std::cout<<"   CASM HARDENING LAW constructed"<<std::endl;
	}

	//*******************************ASSIGMENT OPERATOR***********************************
	//************************************************************************************

	CasmHardeningLaw& CasmHardeningLaw::operator=(CasmHardeningLaw const& rOther)
	{
		HardeningLaw::operator=(rOther);
		return *this;
	}

	//*******************************COPY CONSTRUCTOR*************************************
	//************************************************************************************

	CasmHardeningLaw::CasmHardeningLaw(CasmHardeningLaw const& rOther)
		:HardeningLaw(rOther)
	{

	}

	//********************************DESTRUCTOR******************************************
	//************************************************************************************

	CasmHardeningLaw::~CasmHardeningLaw()
	{
		
	}

	/// Operations.

	//*******************************CALCULATE TOTAL HARDENING****************************
	//************************************************************************************
	void CasmHardeningLaw::CalculateHardening(PlasticVariablesType& rPlasticVariables, const double& rDeltaAlpha, const double& rDeltaBeta)
	{
		const double SwellingSlope 		= GetProperties()[SWELLING_SLOPE];
		const double OtherSlope 		= GetProperties()[NORMAL_COMPRESSION_SLOPE];

		//update hardening parameter
		rPlasticVariables.PreconsolidationPressure += rPlasticVariables.PreconsolidationPressure/(OtherSlope-SwellingSlope)*(-rDeltaAlpha );
	}

	Vector& CasmHardeningLaw::CalculateHardening(Vector& rHardening, const double& rAlpha, const double& rBeta, const double& rAlphaCum, const double& rBetaCum, const double rTemperature)
	{
		//rAlpha ... inc vol strain, rBeta ... inc dev strain
		const double FirstPreconsolidationPressure 	= GetProperties()[PRE_CONSOLIDATION_STRESS];
		const double SwellingSlope 					= GetProperties()[SWELLING_SLOPE];
		const double OtherSlope 					= GetProperties()[NORMAL_COMPRESSION_SLOPE];
		
		rHardening(0) = -FirstPreconsolidationPressure*(std::exp ((-rAlpha)/(OtherSlope-SwellingSlope)) );
		
		return rHardening;
	}

	double& CasmHardeningLaw::CalculateHardening(double& rHardening, const double& rAlpha, const double& rBeta, const double& rAlphaCum, const double& rBetaCum, const double rTemperature)
	{
		Vector aux = ZeroVector(1);
		aux = CalculateHardening(aux, rAlpha, rBeta, rAlphaCum, rBetaCum);
		rHardening = aux(0);
		return rHardening;
	}
	
	void CasmHardeningLaw::save( Serializer& rSerializer ) const
	{
		KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, HardeningLaw )
	}

	void CasmHardeningLaw::load( Serializer& rSerializer )
	{
		KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, HardeningLaw )
	}

}  // namespace Kratos.
