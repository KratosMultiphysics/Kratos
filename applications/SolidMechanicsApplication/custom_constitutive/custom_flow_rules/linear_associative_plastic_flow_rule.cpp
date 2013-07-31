//
//   Project Name:        KratosSolidMechanicsApplication $
//   Last modified by:    $Author:            JMCarbonell $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/properties.h"
#include "custom_constitutive/custom_flow_rules/linear_associative_plastic_flow_rule.hpp"

#include "solid_mechanics_application.h"

namespace Kratos
{

//*******************************CONSTRUCTOR******************************************
//************************************************************************************

LinearAssociativePlasticFlowRule::LinearAssociativePlasticFlowRule()
	:NonLinearAssociativePlasticFlowRule()
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


//********************************DESTRUCTOR******************************************
//************************************************************************************

LinearAssociativePlasticFlowRule::~LinearAssociativePlasticFlowRule()
{
}

/// Operations.


//***************************CALCULATE LOCAL NEWTON PROCEDURE*************************
//************************************************************************************


bool LinearAssociativePlasticFlowRule::CalculateConstistencyCondition( RadialReturnVariables& rReturnMappingVariables, InternalVariables& rPlasticVariables )
{
	//start
	double DeltaStateFunction = 0;
	rReturnMappingVariables.DeltaGamma    = 0;

	double StateFunction                  = rReturnMappingVariables.TrialStateFunction;
   
	//Calculate Delta State Function:
	DeltaStateFunction = mpYieldCriterion->CalculateDeltaStateFunction( DeltaStateFunction, rReturnMappingVariables.LameMu_bar, rPlasticVariables.EquivalentPlasticStrain );

	//Calculate DeltaGamma:
	rReturnMappingVariables.DeltaGamma = StateFunction/DeltaStateFunction;
	       
	//Update Equivalent Plastic Strain:
	rPlasticVariables.DeltaPlasticStrain       = sqrt(2.0/3.0) * rReturnMappingVariables.DeltaGamma;
	rPlasticVariables.EquivalentPlasticStrain += rPlasticVariables.DeltaPlasticStrain;
	       
	return true;	
}


//**************CALCULATE SCALING FACTORS FOR THE ELASTO PLASTIC MODULI***************
//************************************************************************************

void LinearAssociativePlasticFlowRule::CalculateScalingFactors( RadialReturnVariables& rReturnMappingVariables, PlasticFactors& rScalingFactors )
{
	
	//1.-Identity build
	Matrix IdentityMatrix  = identity_matrix<double> (3);

	//2.-Auxiliar matrices
	rScalingFactors.Normal      = rReturnMappingVariables.TrialIsoStressMatrix * ( 1.0 / rReturnMappingVariables.NormIsochoricStress );


	rScalingFactors.Dev_Normal  = zero_matrix<double> (3);


	//3.-Auxiliar constants
	double EquivalentPlasticStrain = mInternalVariables.EquivalentPlasticStrain + sqrt(2.0/3.0) * rReturnMappingVariables.DeltaGamma;
	double DeltaHardening = 0;
	DeltaHardening = mpHardeningLaw->CalculateDeltaHardening( DeltaHardening, EquivalentPlasticStrain );

	rScalingFactors.Beta0 = 1.0 + DeltaHardening/(3.0 * rReturnMappingVariables.LameMu_bar);
		
	rScalingFactors.Beta1 = 1.0 - ( 2.0 * LameMu_bar * DeltaGamma / rParameters.NormIsochoricStress );
		
	rScalingFactors.Beta2 = 0;
		
	rScalingFactors.Beta3 = ( ( 1.0 / Beta0 ) - ( 1 - Beta1 ) );
		
	rScalingFactors.Beta4 = 0;
	
};



void LinearAssociativePlasticFlowRule::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, NonLinearAssociativePlasticFlowRule );
}

void LinearAssociativePlasticFlowRule::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, NonLinearAssociativePlasticFlowRule );
}


}  // namespace Kratos.
