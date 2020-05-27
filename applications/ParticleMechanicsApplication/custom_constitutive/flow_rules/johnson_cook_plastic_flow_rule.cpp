//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    JMCarbonell
//					 (adapted to Particle Mechanics by Peter Wilson)
//

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/properties.h"
#include "custom_constitutive/flow_rules/johnson_cook_plastic_flow_rule.hpp"

#include "particle_mechanics_application_variables.h"

namespace Kratos
{

//*******************************CONSTRUCTOR******************************************
//************************************************************************************

	JohnsonCookPlasticFlowRule::JohnsonCookPlasticFlowRule()
{

}

//*****************************INITIALIZATION CONSTRUCTOR*****************************
//************************************************************************************

	JohnsonCookPlasticFlowRule::JohnsonCookPlasticFlowRule(YieldCriterionPointer pYieldCriterion)
	:ParticleFlowRule(pYieldCriterion)
{

}

//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

	JohnsonCookPlasticFlowRule& JohnsonCookPlasticFlowRule::operator=(JohnsonCookPlasticFlowRule const& rOther)
{
   ParticleFlowRule::operator=(rOther);
   return *this;
}

//*******************************COPY CONSTRUCTOR*************************************
//************************************************************************************

	JohnsonCookPlasticFlowRule::JohnsonCookPlasticFlowRule(JohnsonCookPlasticFlowRule const& rOther)
	:ParticleFlowRule(rOther)
{

}

//********************************CLONE***********************************************
//************************************************************************************

ParticleFlowRule::Pointer JohnsonCookPlasticFlowRule::Clone() const
{
  return Kratos::make_shared<JohnsonCookPlasticFlowRule>(*this);
}


//********************************DESTRUCTOR******************************************
//************************************************************************************

JohnsonCookPlasticFlowRule::~JohnsonCookPlasticFlowRule()
{
}

/// Operations.


//***************************CALCULATE LOCAL NEWTON PROCEDURE*************************
//************************************************************************************

bool JohnsonCookPlasticFlowRule::CalculateConsistencyCondition( RadialReturnVariables& rReturnMappingVariables, InternalVariables& rPlasticVariables, ParticleYieldCriterion::Parameters& rCriterionParameters )
{

	bool converged    = false;

	//Start 1rst Newton Raphson iteration
	rReturnMappingVariables.Options.Set(PLASTIC_RATE_REGION,true);
	rCriterionParameters.SetRateFactor(1); //plastic rate region on
	converged = this->CalculateRateDependentConsistency (rReturnMappingVariables, rPlasticVariables, rCriterionParameters);

	// if(!converged)
	//   std::cout<<" ConstitutiveLaw did not converge on the rate dependent return mapping"<<std::endl;

	const double& ReferenceStrainRate = GetProperties()[REFERENCE_STRAIN_RATE];
	double MaterialDeltaPlasticStrain = ReferenceStrainRate * rReturnMappingVariables.DeltaTime;

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

bool JohnsonCookPlasticFlowRule::CalculateRateDependentConsistency( RadialReturnVariables& rReturnMappingVariables, InternalVariables& rPlasticVariables, ParticleYieldCriterion::Parameters& rCriterionParameters)
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
	//rReturnMappingVariables.DeltaGamma = std::sqrt(3.0*0.5) * ( rPlasticVariables.EquivalentPlasticStrain - rPlasticVariables.EquivalentPlasticStrainOld );


	const double& ReferenceStrainRate = GetProperties()[REFERENCE_STRAIN_RATE];

	rReturnMappingVariables.DeltaGamma         = std::sqrt(3.0*0.5) * ReferenceStrainRate * rReturnMappingVariables.DeltaTime;

	rPlasticVariables.DeltaPlasticStrain       = std::sqrt(2.0/3.0) * rReturnMappingVariables.DeltaGamma;

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
		rPlasticVariables.DeltaPlasticStrain       = std::sqrt(2.0/3.0) * rReturnMappingVariables.DeltaGamma;

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

bool JohnsonCookPlasticFlowRule::CalculateRateIndependentConsistency( RadialReturnVariables& rReturnMappingVariables, InternalVariables& rPlasticVariables, ParticleYieldCriterion::Parameters& rCriterionParameters )
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
	//rReturnMappingVariables.DeltaGamma = std::sqrt(3.0*0.5) * ( rPlasticVariables.EquivalentPlasticStrain - rPlasticVariables.EquivalentPlasticStrainOld );


	double StateFunction = rReturnMappingVariables.TrialStateFunction;

	double InitialStateFunction =  StateFunction;

	//std::cout<<" StateFunction "<<StateFunction<<std::endl;

	rPlasticVariables.DeltaPlasticStrain       = std::sqrt(2.0/3.0) * rReturnMappingVariables.DeltaGamma;

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
		rPlasticVariables.DeltaPlasticStrain       = std::sqrt(2.0/3.0) * rReturnMappingVariables.DeltaGamma;

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

double JohnsonCookPlasticFlowRule::CalculateLineSearch( RadialReturnVariables& rReturnMappingVariables, InternalVariables& rPlasticVariables, ParticleYieldCriterion::Parameters& rCriterionParameters)
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

	if(R0*R1<0.0){

	  double R2 = R1;

	  if(fabs(R1)<fabs(R0))
	    R2=R0;
	  double R0start = R0;


	  double nabla = 0.0;
	  double delta = 1.0;

	  //if( Residual0 < StateFunction ){

	    while ( fabs(R2/R0start)>0.3 && iter<MaxIterations && (R1*R0)<0.0 && fabs(R1)>1e-7 && fabs(R0)>1e-7 )
	      {

		alpha = 0.5*(nabla+delta);

		rPlasticVariables.DeltaPlasticStrain *= alpha;

		rPlasticVariables.EquivalentPlasticStrain  = rPlasticVariables.EquivalentPlasticStrainOld + rPlasticVariables.DeltaPlasticStrain;
		StateFunction = this->mpYieldCriterion->CalculateStateFunction( StateFunction, rCriterionParameters );

		R2 = rPlasticVariables.DeltaPlasticStrain * StateFunction;

		rPlasticVariables.DeltaPlasticStrain /= alpha;


		if(R2*R1<0.0){
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

	if( alpha != 1.0)
	  std::cout<<" [ LINE SEARCH: (Iterations: "<<iter<<", alpha: "<<alpha<<") ] "<<std::endl;


   	if(alpha>1.0 || alpha<=0.0)
   	  alpha=1.0;

	return alpha;
}


//***************************CALCULATE IMPLEX RETURN MAPPING**************************
//************************************************************************************

void JohnsonCookPlasticFlowRule::CalculateImplexReturnMapping( RadialReturnVariables& rReturnMappingVariables, InternalVariables& rPlasticVariables, ParticleYieldCriterion::Parameters& rCriterionParameters, Matrix& rIsoStressMatrix )
{

        //1.-Computation of the plastic Multiplier
        rReturnMappingVariables.DeltaGamma = std::sqrt(3.0*0.5) * ( rPlasticVariables.EquivalentPlasticStrain - rPlasticVariables.EquivalentPlasticStrainOld );

	//2.- Update back stress, plastic strain and stress
	UpdateConfiguration( rReturnMappingVariables, rIsoStressMatrix );

	//3.- Calculate thermal dissipation and delta thermal dissipation
	if( rReturnMappingVariables.DeltaGamma > 0.0 ){

	  const double& PlasticStrainRate = GetProperties()[PLASTIC_STRAIN_RATE];
	  double MaterialDeltaPlasticStrain = PlasticStrainRate * rReturnMappingVariables.DeltaTime;

	  //plastic rate region on
	  rReturnMappingVariables.Options.Set(PLASTIC_RATE_REGION,true);
	  rCriterionParameters.SetRateFactor(1.0);

	  if( rPlasticVariables.DeltaPlasticStrain < MaterialDeltaPlasticStrain ){
	    //plastic rate region off
	    rReturnMappingVariables.Options.Set(PLASTIC_RATE_REGION,false);
	    rCriterionParameters.SetRateFactor(0.0);
	  }

	  this->CalculateImplexThermalDissipation( rCriterionParameters );

	  rReturnMappingVariables.Options.Set(PLASTIC_REGION,true);

	}
	else{

	  mThermalVariables.PlasticDissipation = 0.0;
	  mThermalVariables.DeltaPlasticDissipation = 0.0;

	}


}

//***************************CALCULATE THERMAL DISSIPATION****************************
//************************************************************************************

void JohnsonCookPlasticFlowRule::CalculateImplexThermalDissipation(ParticleYieldCriterion::Parameters& rCriterionParameters)
{

	//1.- Thermal Dissipation:

	mThermalVariables.PlasticDissipation = mpYieldCriterion->CalculateImplexPlasticDissipation(mThermalVariables.PlasticDissipation, rCriterionParameters);

	//2.- Thermal Dissipation Increment:

	mThermalVariables.DeltaPlasticDissipation = mpYieldCriterion->CalculateImplexDeltaPlasticDissipation(mThermalVariables.DeltaPlasticDissipation, rCriterionParameters);
}

//***************************UPDATE STRESS CONFIGURATION *****************************
//************************************************************************************

void JohnsonCookPlasticFlowRule::UpdateConfiguration(RadialReturnVariables& rReturnMappingVariables, Matrix& rIsoStressMatrix)
{
	//Back Stress update

	//Plastic Strain Update
	if (rReturnMappingVariables.NormIsochoricStress > 0) {

		//Stress Update:
		double Auxiliar = 2.0 * rReturnMappingVariables.LameMu_bar * rReturnMappingVariables.DeltaGamma;

		Matrix Normal = rIsoStressMatrix * (1.0 / rReturnMappingVariables.NormIsochoricStress);

		rIsoStressMatrix -= (Normal * Auxiliar);

	}
}

//***************************CALCULATE STRESS NORM ***********************************
//************************************************************************************

double& JohnsonCookPlasticFlowRule::CalculateStressNorm(Matrix& rStressMatrix, double& rStressNorm)
{
	rStressNorm = std::sqrt((rStressMatrix(0, 0) * rStressMatrix(0, 0)) + (rStressMatrix(1, 1) * rStressMatrix(1, 1)) + (rStressMatrix(2, 2) * rStressMatrix(2, 2)) +
	(rStressMatrix(0, 1) * rStressMatrix(0, 1)) + (rStressMatrix(0, 2) * rStressMatrix(0, 2)) + (rStressMatrix(1, 2) * rStressMatrix(1, 2)) +
	(rStressMatrix(1, 0) * rStressMatrix(1, 0)) + (rStressMatrix(2, 0) * rStressMatrix(2, 0)) + (rStressMatrix(2, 1) * rStressMatrix(2, 1)));


	return rStressNorm;
}


//***************************SET YIELD AND HARDENING VARIABLES************************
//************************************************************************************

void JohnsonCookPlasticFlowRule::SetCriterionParameters(RadialReturnVariables& rReturnMappingVariables, InternalVariables& rPlasticVariables, ParticleYieldCriterion::Parameters& rCriterionParameters)
{
	// constant variables during the return mapping
	rCriterionParameters.SetStressNorm(rReturnMappingVariables.NormIsochoricStress);

	rCriterionParameters.SetDeltaTime(rReturnMappingVariables.DeltaTime);

	rCriterionParameters.SetLameMu_bar(rReturnMappingVariables.LameMu_bar);

	rCriterionParameters.SetEquivalentPlasticStrainOld(rPlasticVariables.EquivalentPlasticStrainOld);

	rCriterionParameters.SetTemperature(rReturnMappingVariables.Temperature);


	// changing variables during the return mapping
	rReturnMappingVariables.DeltaGamma = 0;

	rCriterionParameters.SetDeltaGamma(rReturnMappingVariables.DeltaGamma);

	rCriterionParameters.SetEquivalentPlasticStrain(rPlasticVariables.EquivalentPlasticStrain);

	rCriterionParameters.SetRateFactor(0);

	// changing thermal variables during the return mapping
	rReturnMappingVariables.Thermal.clear();

}


//***************************CALCULATE THERMAL DISSIPATION****************************
//************************************************************************************

void JohnsonCookPlasticFlowRule::CalculateThermalDissipation(ParticleYieldCriterion::Parameters& rCriterionParameters, ThermalVariables& rThermalVariables)
{

	//1.- Thermal Dissipation:

	mThermalVariables.PlasticDissipation = mpYieldCriterion->CalculatePlasticDissipation(mThermalVariables.PlasticDissipation, rCriterionParameters);


	//2.- Thermal Dissipation Increment:

	mThermalVariables.DeltaPlasticDissipation = mpYieldCriterion->CalculateDeltaPlasticDissipation(mThermalVariables.DeltaPlasticDissipation, rCriterionParameters);

}

//***************************CALCULATE RADIAL RETURN MAPPING**************************
//************************************************************************************

bool JohnsonCookPlasticFlowRule::CalculateReturnMapping(RadialReturnVariables& rReturnMappingVariables, Matrix& rIsoStressMatrix)
{

	//0.- Initialize Variables
	bool PlasticityActive = false;
	rReturnMappingVariables.Options.Set(PLASTIC_REGION, false);

	InternalVariables PlasticVariables = mInternalVariables;
	ParticleYieldCriterion::Parameters CriterionParameters;
	this->SetCriterionParameters(rReturnMappingVariables, PlasticVariables, CriterionParameters);


	//1.- Isochoric stress norm
	rReturnMappingVariables.NormIsochoricStress = CalculateStressNorm(rIsoStressMatrix, rReturnMappingVariables.NormIsochoricStress);

	//2.- Check yield condition
	rReturnMappingVariables.TrialStateFunction = mpYieldCriterion->CalculateYieldCondition(rReturnMappingVariables.TrialStateFunction, CriterionParameters);


	//3.- Initialize PlasticDissipation
	mThermalVariables.PlasticDissipation = 0;
	mThermalVariables.DeltaPlasticDissipation = 0;



	if (rReturnMappingVariables.Options.Is(IMPLEX_ACTIVE))
	{

		this->CalculateImplexReturnMapping(rReturnMappingVariables, PlasticVariables, CriterionParameters, rIsoStressMatrix);

	}
	else {

		if (rReturnMappingVariables.TrialStateFunction <= 0)
		{

			PlasticityActive = false;
			PlasticVariables.DeltaPlasticStrain = 0;
			rReturnMappingVariables.Options.Set(PLASTIC_REGION, false);

		}
		else
		{

			//3.- Calculate the consistency condition
			bool converged = this->CalculateConsistencyCondition(rReturnMappingVariables, PlasticVariables, CriterionParameters);

			if (!converged)
				std::cout << " ConstitutiveLaw did not converge " << std::endl;


			//4.- Update back stress, plastic strain and stress
			UpdateConfiguration(rReturnMappingVariables, rIsoStressMatrix);


			//5.- Calculate thermal dissipation and delta thermal dissipation
			this->CalculateThermalDissipation(CriterionParameters, rReturnMappingVariables.Thermal);

			PlasticityActive = true;
			rReturnMappingVariables.Options.Set(PLASTIC_REGION, true);
		}

	}

	// std::cout<<" ReturnMapping "<<std::endl;
	// mInternalVariables.print();
	// mThermalVariables.print();
	// std::cout<<" rIsoStressMatrix "<<rIsoStressMatrix<<std::endl;

	rReturnMappingVariables.Options.Set(RETURN_MAPPING_COMPUTED, true);

	return 	PlasticityActive;
}


//**************CALCULATE SCALING FACTORS FOR THE ELASTO PLASTIC MODULI***************
//************************************************************************************

void JohnsonCookPlasticFlowRule::CalculateScalingFactors(const RadialReturnVariables& rReturnMappingVariables, PlasticFactors& rScalingFactors)
{

	//1.-Identity build
	Matrix Identity = identity_matrix<double>(3);

	//2.-Auxiliar matrices
	rScalingFactors.Normal = rReturnMappingVariables.TrialIsoStressMatrix * (1.0 / rReturnMappingVariables.NormIsochoricStress);

	Matrix Norm_Normal = prod(rScalingFactors.Normal, trans(rScalingFactors.Normal));

	double Trace_Norm_Normal = Norm_Normal(0, 0) + Norm_Normal(1, 1) + Norm_Normal(2, 2);

	rScalingFactors.Dev_Normal = Norm_Normal;
	rScalingFactors.Dev_Normal -= (1.0 / 3.0) * Trace_Norm_Normal * Identity;


	//3.-Auxiliar constants
	double EquivalentPlasticStrain = mInternalVariables.EquivalentPlasticStrain + std::sqrt(2.0 / 3.0) * rReturnMappingVariables.DeltaGamma;
	double DeltaHardening = 0;

	if (rReturnMappingVariables.Options.Is(IMPLEX_ACTIVE))
	{
		rScalingFactors.Beta0 = 0;

		rScalingFactors.Beta1 = 2.0 * rReturnMappingVariables.LameMu_bar * rReturnMappingVariables.DeltaGamma / rReturnMappingVariables.NormIsochoricStress;

		rScalingFactors.Beta2 = (2.0 / 3.0) * rReturnMappingVariables.NormIsochoricStress * rReturnMappingVariables.DeltaGamma / (rReturnMappingVariables.LameMu_bar);

		rScalingFactors.Beta3 = (-rScalingFactors.Beta1 + rScalingFactors.Beta2);

		rScalingFactors.Beta4 = (-rScalingFactors.Beta1) * rReturnMappingVariables.NormIsochoricStress / (rReturnMappingVariables.LameMu_bar);

	}
	else
	{

		ParticleHardeningLaw::Parameters HardeningParameters;
		HardeningParameters.SetTemperature(rReturnMappingVariables.Temperature);
		HardeningParameters.SetEquivalentPlasticStrain(EquivalentPlasticStrain);
		HardeningParameters.SetDeltaGamma(rReturnMappingVariables.DeltaGamma);
		HardeningParameters.SetDeltaTime(rReturnMappingVariables.DeltaTime);

		if (rReturnMappingVariables.Options.Is(PLASTIC_RATE_REGION))
			HardeningParameters.SetRateFactor(1);
		else if (rReturnMappingVariables.Options.IsNot(PLASTIC_RATE_REGION))
			HardeningParameters.SetRateFactor(0);

		DeltaHardening = mpYieldCriterion->GetHardeningLaw().CalculateDeltaHardening(DeltaHardening, HardeningParameters);

		rScalingFactors.Beta0 = 1.0 + DeltaHardening / (3.0 * rReturnMappingVariables.LameMu_bar);

		rScalingFactors.Beta1 = 2.0 * rReturnMappingVariables.LameMu_bar * rReturnMappingVariables.DeltaGamma / rReturnMappingVariables.NormIsochoricStress;

		rScalingFactors.Beta2 = ((1.0 - (1.0 / rScalingFactors.Beta0)) * (2.0 / 3.0) * rReturnMappingVariables.NormIsochoricStress * rReturnMappingVariables.DeltaGamma) / (rReturnMappingVariables.LameMu_bar);

		rScalingFactors.Beta3 = ((1.0 / rScalingFactors.Beta0) - rScalingFactors.Beta1 + rScalingFactors.Beta2);

		rScalingFactors.Beta4 = ((1.0 / rScalingFactors.Beta0) - rScalingFactors.Beta1) * rReturnMappingVariables.NormIsochoricStress / (rReturnMappingVariables.LameMu_bar);

	}

	//std::cout<<"FACTORS:: Beta0 "<<rScalingFactors.Beta0<<" Beta 1 "<<rScalingFactors.Beta1<<" Beta2 "<<rScalingFactors.Beta2<<" Beta 3 "<<rScalingFactors.Beta3<<" Beta4 "<<rScalingFactors.Beta4<<std::endl;
}


//***************************UPDATE INTERNAL VARIABLES********************************
//************************************************************************************

bool JohnsonCookPlasticFlowRule::UpdateInternalVariables(RadialReturnVariables& rReturnMappingVariables)
{

	mInternalVariables.EquivalentPlasticStrainOld = mInternalVariables.EquivalentPlasticStrain;

	mInternalVariables.DeltaPlasticStrain = std::sqrt(2.0 / 3.0) * rReturnMappingVariables.DeltaGamma;

	mInternalVariables.EquivalentPlasticStrain += mInternalVariables.DeltaPlasticStrain;

	mInternalVariables.DeltaPlasticStrain *= (1.0 / rReturnMappingVariables.DeltaTime);

	//update thermal variables
	// mThermalVariables = rReturnMappingVariables.Thermal;

	// mInternalVariables.print();
	// mThermalVariables.print();

	return true;
}


void JohnsonCookPlasticFlowRule::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, ParticleFlowRule );
}

void JohnsonCookPlasticFlowRule::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ParticleFlowRule );
}


}  // namespace Kratos.
