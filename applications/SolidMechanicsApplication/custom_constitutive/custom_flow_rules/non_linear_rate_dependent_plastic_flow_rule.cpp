//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2018 $
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
#include "custom_constitutive/custom_flow_rules/non_linear_rate_dependent_plastic_flow_rule.hpp"

#include "solid_mechanics_application_variables.h"

namespace Kratos
{

//*******************************CONSTRUCTOR******************************************
//************************************************************************************

NonLinearRateDependentPlasticFlowRule::NonLinearRateDependentPlasticFlowRule()
{

}

//*****************************INITIALIZATION CONSTRUCTOR*****************************
//************************************************************************************

NonLinearRateDependentPlasticFlowRule::NonLinearRateDependentPlasticFlowRule(YieldCriterionPointer pYieldCriterion)
	:NonLinearAssociativePlasticFlowRule(pYieldCriterion)
{

}

//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

NonLinearRateDependentPlasticFlowRule& NonLinearRateDependentPlasticFlowRule::operator=(NonLinearRateDependentPlasticFlowRule const& rOther)
{
   NonLinearAssociativePlasticFlowRule::operator=(rOther);
   return *this;
}

//*******************************COPY CONSTRUCTOR*************************************
//************************************************************************************

NonLinearRateDependentPlasticFlowRule::NonLinearRateDependentPlasticFlowRule(NonLinearRateDependentPlasticFlowRule const& rOther)
	:NonLinearAssociativePlasticFlowRule(rOther)
{

}

//********************************CLONE***********************************************
//************************************************************************************

FlowRule::Pointer NonLinearRateDependentPlasticFlowRule::Clone() const
{
  return Kratos::make_shared<NonLinearRateDependentPlasticFlowRule>(*this);
}


//********************************DESTRUCTOR******************************************
//************************************************************************************

NonLinearRateDependentPlasticFlowRule::~NonLinearRateDependentPlasticFlowRule()
{
}

/// Operations.


//***************************CALCULATE LOCAL NEWTON PROCEDURE*************************
//************************************************************************************

bool NonLinearRateDependentPlasticFlowRule::CalculateConsistencyCondition( RadialReturnVariables& rReturnMappingVariables, InternalVariables& rPlasticVariables, YieldCriterion::Parameters& rCriterionParameters )
{

	bool converged    = false;

	//Start 1rst Newton Raphson iteration
	rReturnMappingVariables.Options.Set(PLASTIC_RATE_REGION,true);
	rCriterionParameters.SetRateFactor(1); //plastic rate region on
	converged = this->CalculateRateDependentConsistency (rReturnMappingVariables, rPlasticVariables, rCriterionParameters);

	// if(!converged)
	//   std::cout<<" ConstitutiveLaw did not converge on the rate dependent return mapping"<<std::endl;

	const double& PlasticStrainRate = GetProperties()[PLASTIC_STRAIN_RATE];
	double MaterialDeltaPlasticStrain = PlasticStrainRate * rReturnMappingVariables.DeltaTime;

	//std::cout<<" DeltaPlasticStrain: "<<rPlasticVariables.DeltaPlasticStrain<<" MaterialDeltaPlasticStrain: "<<MaterialDeltaPlasticStrain<<std::endl;

	if( rPlasticVariables.DeltaPlasticStrain < MaterialDeltaPlasticStrain ){

	        //std::cout<<" DeltaPlasticStrain: "<<rPlasticVariables.DeltaPlasticStrain<<" MaterialDeltaPlasticStrain: "<<MaterialDeltaPlasticStrain<<std::endl;

	        //Start 2nd Newton Raphson iteration
		rReturnMappingVariables.Options.Set(PLASTIC_RATE_REGION,false);
		rCriterionParameters.SetRateFactor(0); //plastic rate region off

		converged = this->CalculateRateIndependentConsistency (rReturnMappingVariables, rPlasticVariables, rCriterionParameters);
		// if(!converged)
		//   std::cout<<" ConstitutiveLaw did not converge on the rate independent return mapping"<<std::endl;

	}

	return converged;

}


//***************************CALCULATE LOCAL NEWTON PROCEDURE (1)*********************
//************************************************************************************

bool NonLinearRateDependentPlasticFlowRule::CalculateRateDependentConsistency( RadialReturnVariables& rReturnMappingVariables, InternalVariables& rPlasticVariables, YieldCriterion::Parameters& rCriterionParameters)
{
	//Set convergence parameters
	unsigned int iter    = 0;
	double Tolerance     = 1e-5;
	double MaxIterations = 50;

	//start
	double DeltaDeltaGamma    = 0;
	double DeltaStateFunction = 0;
	rReturnMappingVariables.DeltaGamma  = 0;

	//initial guess
	//rReturnMappingVariables.DeltaGamma = sqrt(3.0*0.5) * ( rPlasticVariables.EquivalentPlasticStrain - rPlasticVariables.EquivalentPlasticStrainOld );


	const double& PlasticStrainRate = GetProperties()[PLASTIC_STRAIN_RATE];

	rReturnMappingVariables.DeltaGamma         = sqrt(3.0*0.5) * PlasticStrainRate * rReturnMappingVariables.DeltaTime;

	rPlasticVariables.DeltaPlasticStrain       = sqrt(2.0/3.0) * rReturnMappingVariables.DeltaGamma;

	rPlasticVariables.EquivalentPlasticStrain  = rPlasticVariables.EquivalentPlasticStrainOld + rPlasticVariables.DeltaPlasticStrain;


	//std::cout<<" Rate Dependent DeltaGamma "<<rReturnMappingVariables.DeltaGamma<<std::endl;

	double StateFunction = mpYieldCriterion->CalculateStateFunction( StateFunction, rCriterionParameters );

	double InitialStateFunction =  StateFunction;

	double alpha = 1;
	while ( fabs(StateFunction)>=Tolerance && iter<=MaxIterations )
	{
	       //Calculate Delta State Function:
		DeltaStateFunction = this->mpYieldCriterion->CalculateDeltaStateFunction( DeltaStateFunction, rCriterionParameters );

		//Calculate DeltaGamma:
		DeltaDeltaGamma = StateFunction/DeltaStateFunction;
		rReturnMappingVariables.DeltaGamma += DeltaDeltaGamma;

		//Update Equivalent Plastic Strain:
		rPlasticVariables.DeltaPlasticStrain       = sqrt(2.0/3.0) * rReturnMappingVariables.DeltaGamma;

		//alpha =  CalculateLineSearch( rReturnMappingVariables,rPlasticVariables,rCriterionParameters );

		rPlasticVariables.EquivalentPlasticStrain  = rPlasticVariables.EquivalentPlasticStrainOld + alpha * rPlasticVariables.DeltaPlasticStrain;


		//check negative equivalent plastic strain:
		// if( rPlasticVariables.EquivalentPlasticStrain < 0 )
		//   std::cout<<" EquivalentPlasticStrain NEGATIVE "<<rPlasticVariables.EquivalentPlasticStrain<<std::endl;

		//Calculate State Function:
		StateFunction = this->mpYieldCriterion->CalculateStateFunction( StateFunction, rCriterionParameters );

		iter++;
	}


	//std::cout<<" RateDependent Consistency [ PlasticStrain: "<<rPlasticVariables.EquivalentPlasticStrain<<" DeltaPlasticStrain "<<rPlasticVariables.DeltaPlasticStrain<<" DeltaGamma "<<rReturnMappingVariables.DeltaGamma<<" State Function "<<StateFunction<<"] (iters:"<<iter<<")"<<std::endl;

	// if( rReturnMappingVariables.DeltaGamma < 0 ){
	//   std::cout<<" DeltaGamma NEGATIVE "<<rReturnMappingVariables.DeltaGamma<<std::endl;
	// }


	if(iter>MaxIterations){
	  std::cout<<" IniStateFunction "<<InitialStateFunction<<" StateFunction "<<StateFunction<<" PlasticStrain "<<rPlasticVariables.EquivalentPlasticStrain<<" DeltaPlasticStrain "<<rPlasticVariables.DeltaPlasticStrain<<std::endl;
	  return false;
	}

	return true;
}



//***************************CALCULATE LOCAL NEWTON PROCEDURE (2)*********************
//************************************************************************************

bool NonLinearRateDependentPlasticFlowRule::CalculateRateIndependentConsistency( RadialReturnVariables& rReturnMappingVariables, InternalVariables& rPlasticVariables, YieldCriterion::Parameters& rCriterionParameters )
{
	//Set convergence parameters
	unsigned int iter    = 0;
	double Tolerance     = 1e-5;
	double MaxIterations = 50;

	//start
	double DeltaDeltaGamma    = 0;
	double DeltaStateFunction = 0;
	rReturnMappingVariables.DeltaGamma  = 1e-40;  //this can not be zero (zig-zag in the iterative loop if is zero)

	//initial guess
	//rReturnMappingVariables.DeltaGamma = sqrt(3.0*0.5) * ( rPlasticVariables.EquivalentPlasticStrain - rPlasticVariables.EquivalentPlasticStrainOld );


	double StateFunction = rReturnMappingVariables.TrialStateFunction;

	double InitialStateFunction =  StateFunction;

	//std::cout<<" StateFunction "<<StateFunction<<std::endl;

	rPlasticVariables.DeltaPlasticStrain       = sqrt(2.0/3.0) * rReturnMappingVariables.DeltaGamma;

	rPlasticVariables.EquivalentPlasticStrain  = rPlasticVariables.EquivalentPlasticStrainOld + rPlasticVariables.DeltaPlasticStrain;

	double alpha = 1;
	while ( fabs(StateFunction)>=Tolerance && iter<=MaxIterations)
	{
		//Calculate Delta State Function:
		DeltaStateFunction = this->mpYieldCriterion->CalculateDeltaStateFunction( DeltaStateFunction, rCriterionParameters );

		//Calculate DeltaGamma:
		DeltaDeltaGamma  = StateFunction/DeltaStateFunction;
		rReturnMappingVariables.DeltaGamma += DeltaDeltaGamma;

		//Update Equivalent Plastic Strain:
		rPlasticVariables.DeltaPlasticStrain       = sqrt(2.0/3.0) * rReturnMappingVariables.DeltaGamma;

		//alpha =  CalculateLineSearch( rReturnMappingVariables,rPlasticVariables,rCriterionParameters );

		rPlasticVariables.EquivalentPlasticStrain  = rPlasticVariables.EquivalentPlasticStrainOld + alpha * rPlasticVariables.DeltaPlasticStrain;

		//Calculate State Function:
		StateFunction = this->mpYieldCriterion->CalculateStateFunction( StateFunction, rCriterionParameters );

		//std::cout<<" StateFunction "<<StateFunction<<" PlasticStrain "<<rPlasticVariables.EquivalentPlasticStrain<<" DeltaPlasticStrain "<<rPlasticVariables.DeltaPlasticStrain<<" alpha "<<alpha<<" DeltaStateFunction "<<DeltaStateFunction<<std::endl;


		iter++;
	}

	// std::cout<<" RateDependent Independent Consistency [ PlasticStrain: "<<rPlasticVariables.EquivalentPlasticStrain<<" DeltaPlasticStrain "<<rPlasticVariables.DeltaPlasticStrain<<" DeltaGamma "<<rReturnMappingVariables.DeltaGamma<<" State Function "<<StateFunction<<"] (iters:"<<iter<<")"<<std::endl;

	// if( rReturnMappingVariables.DeltaGamma < 0 ){
	//   std::cout<<" ERROR: DeltaGamma NEGATIVE "<<rReturnMappingVariables.DeltaGamma<<std::endl;
	// }

	if(iter>MaxIterations){
	  std::cout<<" IniStateFunction "<<InitialStateFunction<<" Rate StateFunction "<<StateFunction<<" PlasticStrain "<<rPlasticVariables.EquivalentPlasticStrain<<std::endl;
	  return false;
	}


	return true;
}

//************************CALCULATE RETURN MAPPING LINE SEARCH************************
//************************************************************************************

double NonLinearRateDependentPlasticFlowRule::CalculateLineSearch( RadialReturnVariables& rReturnMappingVariables, InternalVariables& rPlasticVariables, YieldCriterion::Parameters& rCriterionParameters)
{
	//Set convergence parameters
	unsigned int iter    = 0;
	double MaxIterations = 10;

	//start preserve initial variables
	double DeltaPlasticStrain = rPlasticVariables.DeltaPlasticStrain;

	double StateFunction = this->mpYieldCriterion->CalculateStateFunction( StateFunction, rCriterionParameters );
        double R0 = rPlasticVariables.DeltaPlasticStrain * StateFunction;

	//double Residual0 = StateFunction;

	rPlasticVariables.EquivalentPlasticStrain  = rPlasticVariables.EquivalentPlasticStrainOld + rPlasticVariables.DeltaPlasticStrain;
	StateFunction = this->mpYieldCriterion->CalculateStateFunction( StateFunction, rCriterionParameters );

        double R1 = rPlasticVariables.DeltaPlasticStrain * StateFunction;

	double alpha = 1;

	if(R0*R1<0){

	  double R2 = R1;

	  if(fabs(R1)<fabs(R0))
	    R2=R0;
	  double R0start = R0;


	  double nabla = 0;
	  double delta = 1;

	  //if( Residual0 < StateFunction ){

	    while ( fabs(R2/R0start)>0.3 && iter<MaxIterations && (R1*R0)<0 && fabs(R1)>1e-7 && fabs(R0)>1e-7 )
	      {

		alpha = 0.5*(nabla+delta);

		rPlasticVariables.DeltaPlasticStrain *= alpha;

		rPlasticVariables.EquivalentPlasticStrain  = rPlasticVariables.EquivalentPlasticStrainOld + rPlasticVariables.DeltaPlasticStrain;
		StateFunction = this->mpYieldCriterion->CalculateStateFunction( StateFunction, rCriterionParameters );

		R2 = rPlasticVariables.DeltaPlasticStrain * StateFunction;

		rPlasticVariables.DeltaPlasticStrain /= alpha;


		if(R2*R1<0){
		  nabla = alpha;
		  R0 = R2;
		}
		else if(R2*R0<0){
		  delta = alpha;
		  R1 = R2;
		}
		else{
		  break;
		}

		iter++;
	      }
	    //}

	}

	rPlasticVariables.DeltaPlasticStrain = DeltaPlasticStrain;

	if( alpha != 1)
	  std::cout<<" [ LINE SEARCH: (Iterations: "<<iter<<", alpha: "<<alpha<<") ] "<<std::endl;


   	if(alpha>1 || alpha<=0)
   	  alpha=1;

	return alpha;
}


//***************************CALCULATE IMPLEX RETURN MAPPING**************************
//************************************************************************************

void NonLinearRateDependentPlasticFlowRule::CalculateImplexReturnMapping( RadialReturnVariables& rReturnMappingVariables, InternalVariables& rPlasticVariables, YieldCriterion::Parameters& rCriterionParameters, Matrix& rIsoStressMatrix )
{

        //1.-Computation of the plastic Multiplier
        rReturnMappingVariables.DeltaGamma = sqrt(3.0*0.5) * ( rPlasticVariables.EquivalentPlasticStrain - rPlasticVariables.EquivalentPlasticStrainOld );

	//2.- Update back stress, plastic strain and stress
	UpdateConfiguration( rReturnMappingVariables, rIsoStressMatrix );

	//3.- Calculate thermal dissipation and delta thermal dissipation
	if( rReturnMappingVariables.DeltaGamma > 0 ){

	  const double& PlasticStrainRate = GetProperties()[PLASTIC_STRAIN_RATE];
	  double MaterialDeltaPlasticStrain = PlasticStrainRate * rReturnMappingVariables.DeltaTime;

	  //plastic rate region on
	  rReturnMappingVariables.Options.Set(PLASTIC_RATE_REGION,true);
	  rCriterionParameters.SetRateFactor(1);

	  if( rPlasticVariables.DeltaPlasticStrain < MaterialDeltaPlasticStrain ){
	    //plastic rate region off
	    rReturnMappingVariables.Options.Set(PLASTIC_RATE_REGION,false);
	    rCriterionParameters.SetRateFactor(0);
	  }

	  this->CalculateImplexThermalDissipation( rCriterionParameters );

	  rReturnMappingVariables.Options.Set(PLASTIC_REGION,true);

	}
	else{

	  mThermalVariables.PlasticDissipation = 0;
	  mThermalVariables.DeltaPlasticDissipation = 0;

	}


}


void NonLinearRateDependentPlasticFlowRule::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, NonLinearAssociativePlasticFlowRule );
}

void NonLinearRateDependentPlasticFlowRule::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, NonLinearAssociativePlasticFlowRule );
}


}  // namespace Kratos.
