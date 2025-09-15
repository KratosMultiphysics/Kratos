//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_constitutive/custom_flow_rules/linear_associative_plastic_flow_rule.hpp"

#include "solid_mechanics_application_variables.h"

namespace Kratos
{

//*******************************CONSTRUCTOR******************************************
//************************************************************************************

LinearAssociativePlasticFlowRule::LinearAssociativePlasticFlowRule()
	:NonLinearAssociativePlasticFlowRule()
{

}

//*****************************INITIALIZATION CONSTRUCTOR*****************************
//************************************************************************************

LinearAssociativePlasticFlowRule::LinearAssociativePlasticFlowRule(YieldCriterionPointer pYieldCriterion)
	:NonLinearAssociativePlasticFlowRule(pYieldCriterion)
{

}

//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

LinearAssociativePlasticFlowRule& LinearAssociativePlasticFlowRule::operator=(LinearAssociativePlasticFlowRule const& rOther)
{
   NonLinearAssociativePlasticFlowRule::operator=(rOther);
   return *this;
}

//*******************************COPY CONSTRUCTOR*************************************
//************************************************************************************

LinearAssociativePlasticFlowRule::LinearAssociativePlasticFlowRule(LinearAssociativePlasticFlowRule const& rOther)
	:NonLinearAssociativePlasticFlowRule(rOther)
{

}

//********************************CLONE***********************************************
//************************************************************************************

FlowRule::Pointer LinearAssociativePlasticFlowRule::Clone() const
{
  return Kratos::make_shared<LinearAssociativePlasticFlowRule>(*this);
}

//********************************DESTRUCTOR******************************************
//************************************************************************************

LinearAssociativePlasticFlowRule::~LinearAssociativePlasticFlowRule()
{
}

/// Operations.


//***************************CALCULATE LOCAL NEWTON PROCEDURE*************************
//************************************************************************************


bool LinearAssociativePlasticFlowRule::CalculateConsistencyCondition( RadialReturnVariables& rReturnMappingVariables, InternalVariables& rPlasticVariables, YieldCriterion::Parameters& rCriterionParameters )
{
	//start
	double DeltaStateFunction = 0;
	rReturnMappingVariables.DeltaGamma    = 0;

	double StateFunction                  = rReturnMappingVariables.TrialStateFunction;

	//Calculate Delta State Function:
	DeltaStateFunction = mpYieldCriterion->CalculateDeltaStateFunction( DeltaStateFunction, rCriterionParameters );

	//Calculate DeltaGamma:
	rReturnMappingVariables.DeltaGamma = StateFunction/DeltaStateFunction;

	//Update Equivalent Plastic Strain:
	rPlasticVariables.DeltaPlasticStrain       = sqrt(2.0/3.0) * rReturnMappingVariables.DeltaGamma;
	rPlasticVariables.EquivalentPlasticStrain += rPlasticVariables.DeltaPlasticStrain;

	//std::cout<<" Strain Rate "<<rPlasticVariables.DeltaPlasticStrain<<std::endl;

	return true;
}


//**************CALCULATE SCALING FACTORS FOR THE ELASTO PLASTIC MODULI***************
//************************************************************************************

// void LinearAssociativePlasticFlowRule::CalculateScalingFactors(const RadialReturnVariables& rReturnMappingVariables, PlasticFactors& rScalingFactors )
// {

// 	//1.-Identity build
// 	Matrix Identity  = identity_matrix<double> (3);

// 	//2.-Auxiliar matrices
// 	rScalingFactors.Normal      = rReturnMappingVariables.TrialIsoStressMatrix * ( 1.0 / rReturnMappingVariables.NormIsochoricStress );


// 	rScalingFactors.Dev_Normal  = zero_matrix<double> (3);


// 	//3.-Auxiliar constants
// 	double EquivalentPlasticStrain = mInternalVariables.EquivalentPlasticStrain + sqrt(2.0/3.0) * rReturnMappingVariables.DeltaGamma;
// 	double DeltaHardening = 0;
	// HardeningLaw::Parameters HardeningParameters;
	// HardeningParamters.SetTemperature(rReturnMappingVariables.Temperature);
	// HardeningParamters.SetEquivalentPlasticStrain(EquivalentPlasticStrain);

	// DeltaHardening = mpYieldCriterion->GetHardeningLaw().CalculateDeltaHardening( DeltaHardening, HardeningParameters );

// 	rScalingFactors.Beta0 = 1.0 + DeltaHardening/(3.0 * rReturnMappingVariables.LameMu_bar);

// 	rScalingFactors.Beta1 = 1.0 - ( 2.0 * rReturnMappingVariables.LameMu_bar * rReturnMappingVariables.DeltaGamma / rReturnMappingVariables.NormIsochoricStress );

// 	rScalingFactors.Beta2 = 0;

// 	rScalingFactors.Beta3 = ( ( 1.0 / rScalingFactors.Beta0 ) - ( 1 - rScalingFactors.Beta1 ) );

// 	rScalingFactors.Beta4 = 0;

// 	std::cout<<"LINEAR:: Beta0 "<<rScalingFactors.Beta0<<" Beta 1 "<<rScalingFactors.Beta1<<" Beta2 "<<rScalingFactors.Beta2<<" Beta 3 "<<rScalingFactors.Beta3<<" Beta4 "<<rScalingFactors.Beta4<<std::endl;

// };



void LinearAssociativePlasticFlowRule::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, NonLinearAssociativePlasticFlowRule )
}

void LinearAssociativePlasticFlowRule::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, NonLinearAssociativePlasticFlowRule )
}


}  // namespace Kratos.
