//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Created by:          $Author:                    LHauser $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                     Dec 2018 $
//   Revision:            $Revision:                      0.0 $
//
//

// System includes
#include <iostream>

// External includes
#include<cmath>

// Project includes
#include "includes/properties.h"
#include "custom_constitutive/borja_hencky_casm_cem_axisym_2D_law.hpp"
#include "pfem_solid_mechanics_application_variables.h"
#include "custom_utilities/solid_mechanics_math_utilities.hpp"

namespace Kratos
{

	//******************************CONSTRUCTOR*******************************************
	//************************************************************************************

	BorjaHenckyCasmCemPlasticAxisym2DLaw::BorjaHenckyCasmCemPlasticAxisym2DLaw()
		: BorjaHenckyCasmPlasticAxisym2DLaw()
	{
		mpHardeningLaw   = HardeningLaw::Pointer( new CasmCemHardeningLaw() );
		mpYieldCriterion = YieldCriterion::Pointer( new CasmCemYieldCriterion(mpHardeningLaw) );
		mpFlowRule       = FlowRule::Pointer( new BorjaCasmCemExplicitFlowRule(mpYieldCriterion) );
std::cout<<"   CASM cemented 2D axisym constructed"<<std::endl;
	}

	//******************************CONSTRUCTOR*******************************************
	//************************************************************************************

	BorjaHenckyCasmCemPlasticAxisym2DLaw::BorjaHenckyCasmCemPlasticAxisym2DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
	{
		mpHardeningLaw    =  pHardeningLaw;
		mpYieldCriterion  =  YieldCriterion::Pointer( new CasmCemYieldCriterion(mpHardeningLaw) );
		mpFlowRule        =  pFlowRule;
std::cout<<"   CASM cemented 2D axisym constructed"<<std::endl;
	}

	//******************************COPY CONSTRUCTOR**************************************
	//************************************************************************************

	BorjaHenckyCasmCemPlasticAxisym2DLaw::BorjaHenckyCasmCemPlasticAxisym2DLaw(const BorjaHenckyCasmCemPlasticAxisym2DLaw& rOther)
		: BorjaHenckyCasmPlasticAxisym2DLaw(rOther)
	{
	}

	//********************************CLONE***********************************************
	//************************************************************************************

	ConstitutiveLaw::Pointer BorjaHenckyCasmCemPlasticAxisym2DLaw::Clone() const
	{
		BorjaHenckyCasmCemPlasticAxisym2DLaw::Pointer p_clone(new BorjaHenckyCasmCemPlasticAxisym2DLaw(*this));
		return p_clone;
	}

	//*******************************DESTRUCTOR*******************************************
	//************************************************************************************

	BorjaHenckyCasmCemPlasticAxisym2DLaw::~BorjaHenckyCasmCemPlasticAxisym2DLaw()
	{
	}


	//*************************************** GET VALUE *********************************
	double&  BorjaHenckyCasmCemPlasticAxisym2DLaw::GetValue(const Variable<double>& rThisVariable, double & rValue)
	{
		if ( (rThisVariable == BONDING) )
		{
			const PlasticVariablesType& InternalVariables = mpFlowRule->GetPlasticVariables();
			rValue = InternalVariables.Bonding;
		}
		else {
			rValue = BorjaHenckyCasmPlasticAxisym2DLaw::GetValue( rThisVariable, rValue);
		}

		return rValue;
	}
	
	void BorjaHenckyCasmCemPlasticAxisym2DLaw::SetPlasticVariables ( const double& rInitialPreconPressure, const double& rInitialBonding)
	{
		mpFlowRule->SetPlasticVariables(rInitialPreconPressure, rInitialBonding);
	}
	
	const double BorjaHenckyCasmCemPlasticAxisym2DLaw::GetBonding()
	{
		return mpFlowRule->GetPlasticVariables().Bonding;
	}

	void BorjaHenckyCasmCemPlasticAxisym2DLaw::SetValue(const Variable<double>& rThisVariable, const double& rValue, const ProcessInfo& rCurrentProcessInfo)
	{
		if ( rThisVariable == BONDING)
		{
			mpFlowRule->SetBonding(rValue);
		}
		else
		{
			BorjaHenckyCasmPlasticAxisym2DLaw::SetValue( rThisVariable, rValue, rCurrentProcessInfo );
		}
	}

	int BorjaHenckyCasmCemPlasticAxisym2DLaw::Check( const Properties& rMaterialProperties, const GeometryType& rElementGeometry, const ProcessInfo& rCurrentProcessInfo)
	{
		return 0;
	}


} // Namespace Kratos
