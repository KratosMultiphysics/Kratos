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
#include "custom_constitutive/custom_flow_rules/non_linear_associative_plastic_flow_rule.hpp"

#include "solid_mechanics_application_variables.h"

namespace Kratos
{

//*******************************CONSTRUCTOR******************************************
//************************************************************************************

NonLinearAssociativePlasticFlowRule::NonLinearAssociativePlasticFlowRule()
{

}

//*****************************INITIALIZATION CONSTRUCTOR*****************************
//************************************************************************************

NonLinearAssociativePlasticFlowRule::NonLinearAssociativePlasticFlowRule(YieldCriterionPointer pYieldCriterion)
	:FlowRule(pYieldCriterion)
{

}


//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

NonLinearAssociativePlasticFlowRule& NonLinearAssociativePlasticFlowRule::operator=(NonLinearAssociativePlasticFlowRule const& rOther)
{
   FlowRule::operator=(rOther);
   return *this;
}

//*******************************COPY CONSTRUCTOR*************************************
//************************************************************************************

NonLinearAssociativePlasticFlowRule::NonLinearAssociativePlasticFlowRule(NonLinearAssociativePlasticFlowRule const& rOther)
	:FlowRule(rOther)
{

}

//********************************CLONE***********************************************
//************************************************************************************

FlowRule::Pointer NonLinearAssociativePlasticFlowRule::Clone() const
{
  return Kratos::make_shared<NonLinearAssociativePlasticFlowRule>(*this);
}

//********************************DESTRUCTOR******************************************
//************************************************************************************

NonLinearAssociativePlasticFlowRule::~NonLinearAssociativePlasticFlowRule()
{
}

/// Operations.


//***************************CALCULATE STRESS NORM ***********************************
//************************************************************************************

double& NonLinearAssociativePlasticFlowRule::CalculateStressNorm ( Matrix & rStressMatrix, double& rStressNorm )
{

        rStressNorm =  sqrt((rStressMatrix(0,0)*rStressMatrix(0,0))+(rStressMatrix(1,1)*rStressMatrix(1,1))+(rStressMatrix(2,2)*rStressMatrix(2,2))+
		           (rStressMatrix(0,1)*rStressMatrix(0,1))+(rStressMatrix(0,2)*rStressMatrix(0,2))+(rStressMatrix(1,2)*rStressMatrix(1,2))+
		           (rStressMatrix(1,0)*rStressMatrix(1,0))+(rStressMatrix(2,0)*rStressMatrix(2,0))+(rStressMatrix(2,1)*rStressMatrix(2,1)));


	return rStressNorm;
}


//***************************SET YIELD AND HARDENING VARIABLES************************
//************************************************************************************

void NonLinearAssociativePlasticFlowRule::SetCriterionParameters( RadialReturnVariables& rReturnMappingVariables, InternalVariables& rPlasticVariables, YieldCriterion::Parameters& rCriterionParameters )
{
	// constant variables during the return mapping
	rCriterionParameters.SetStressNorm( rReturnMappingVariables.NormIsochoricStress );

	rCriterionParameters.SetDeltaTime( rReturnMappingVariables.DeltaTime );

	rCriterionParameters.SetLameMu_bar( rReturnMappingVariables.LameMu_bar );

	rCriterionParameters.SetEquivalentPlasticStrainOld( rPlasticVariables.EquivalentPlasticStrainOld );

	rCriterionParameters.SetTemperature( rReturnMappingVariables.Temperature );


	// changing variables during the return mapping
	rReturnMappingVariables.DeltaGamma = 0;

	rCriterionParameters.SetDeltaGamma( rReturnMappingVariables.DeltaGamma );

	rCriterionParameters.SetEquivalentPlasticStrain( rPlasticVariables.EquivalentPlasticStrain );

	rCriterionParameters.SetRateFactor(0);

	// changing thermal variables during the return mapping
	rReturnMappingVariables.Thermal.clear();

}


//***************************CALCULATE RADIAL RETURN MAPPING**************************
//************************************************************************************

bool NonLinearAssociativePlasticFlowRule::CalculateReturnMapping( RadialReturnVariables& rReturnMappingVariables, Matrix& rIsoStressMatrix )
{

	//0.- Initialize Variables
	bool PlasticityActive = false;
	rReturnMappingVariables.Options.Set(PLASTIC_REGION,false);

	InternalVariables PlasticVariables = mInternalVariables;
	YieldCriterion::Parameters CriterionParameters;
	this->SetCriterionParameters ( rReturnMappingVariables, PlasticVariables, CriterionParameters );


	//1.- Isochoric stress norm
	rReturnMappingVariables.NormIsochoricStress = CalculateStressNorm( rIsoStressMatrix, rReturnMappingVariables.NormIsochoricStress );

	//2.- Check yield condition
	rReturnMappingVariables.TrialStateFunction = mpYieldCriterion->CalculateYieldCondition( rReturnMappingVariables.TrialStateFunction, CriterionParameters );


	//3.- Initialize PlasticDissipation
	mThermalVariables.PlasticDissipation = 0;
	mThermalVariables.DeltaPlasticDissipation = 0;



	if( rReturnMappingVariables.Options.Is(IMPLEX_ACTIVE) )
	  {

   	     this->CalculateImplexReturnMapping ( rReturnMappingVariables, PlasticVariables, CriterionParameters, rIsoStressMatrix );

	  }
	else{

	  if( rReturnMappingVariables.TrialStateFunction <= 0 )
	    {

	      PlasticityActive = false;
	      PlasticVariables.DeltaPlasticStrain = 0;
	      rReturnMappingVariables.Options.Set(PLASTIC_REGION,false);

	    }
	  else
	    {

	      //3.- Calculate the consistency condition
  	      bool converged = this->CalculateConsistencyCondition( rReturnMappingVariables, PlasticVariables, CriterionParameters);

	      if(!converged)
	       	std::cout<<" ConstitutiveLaw did not converge "<<std::endl;


	      //4.- Update back stress, plastic strain and stress
	      UpdateConfiguration( rReturnMappingVariables, rIsoStressMatrix );


	      //5.- Calculate thermal dissipation and delta thermal dissipation
	      this->CalculateThermalDissipation( CriterionParameters, rReturnMappingVariables.Thermal );

	      PlasticityActive = true;
	      rReturnMappingVariables.Options.Set(PLASTIC_REGION,true);
	    }

	}

	// std::cout<<" ReturnMapping "<<std::endl;
	// mInternalVariables.print();
	// mThermalVariables.print();
	// std::cout<<" rIsoStressMatrix "<<rIsoStressMatrix<<std::endl;

	rReturnMappingVariables.Options.Set(RETURN_MAPPING_COMPUTED,true);

	return 	PlasticityActive;
}

//***************************CALCULATE LOCAL NEWTON PROCEDURE*************************
//************************************************************************************


bool NonLinearAssociativePlasticFlowRule::CalculateConsistencyCondition( RadialReturnVariables& rReturnMappingVariables, InternalVariables& rPlasticVariables, YieldCriterion::Parameters& rCriterionParameters )
{
	//Set convergence parameters
	unsigned int iter    = 0;
	double Tolerance     = 1e-5;
	double MaxIterations = 50;

	//start
	double DeltaDeltaGamma    = 0;
	double DeltaStateFunction = 0;
	rReturnMappingVariables.DeltaGamma    = 0;

	double StateFunction                  = rReturnMappingVariables.TrialStateFunction;

	while ( fabs(StateFunction)>=Tolerance && iter<=MaxIterations)
	{
		//Calculate Delta State Function:
		DeltaStateFunction = mpYieldCriterion->CalculateDeltaStateFunction( DeltaStateFunction, rCriterionParameters );

		//Calculate DeltaGamma:
		DeltaDeltaGamma  = StateFunction/DeltaStateFunction;
		rReturnMappingVariables.DeltaGamma += DeltaDeltaGamma;

		//Update Equivalent Plastic Strain:
		rPlasticVariables.DeltaPlasticStrain       = sqrt(2.0/3.0) * rReturnMappingVariables.DeltaGamma;
		rPlasticVariables.EquivalentPlasticStrain  = rPlasticVariables.EquivalentPlasticStrainOld + rPlasticVariables.DeltaPlasticStrain;

		//Calculate State Function:
		StateFunction = mpYieldCriterion->CalculateStateFunction( StateFunction, rCriterionParameters );


		iter++;
	}


	if(iter>MaxIterations)
	  return false;


	return true;
}



//***************************CALCULATE IMPLEX RETURN MAPPING**************************
//************************************************************************************

void NonLinearAssociativePlasticFlowRule::CalculateImplexReturnMapping( RadialReturnVariables& rReturnMappingVariables, InternalVariables& rPlasticVariables, YieldCriterion::Parameters& rCriterionParameters, Matrix& rIsoStressMatrix )
{


        //1.-Computation of the plastic Multiplier
        rReturnMappingVariables.DeltaGamma = sqrt(3.0/2.0) * ( rPlasticVariables.EquivalentPlasticStrain - rPlasticVariables.EquivalentPlasticStrainOld );

	//2.- Update back stress, plastic strain and stress
	UpdateConfiguration( rReturnMappingVariables, rIsoStressMatrix );

	//3.- Calculate thermal dissipation and delta thermal dissipation
	if( rReturnMappingVariables.DeltaGamma > 0 ){

  	  this->CalculateImplexThermalDissipation( rCriterionParameters );
	  rReturnMappingVariables.Options.Set(PLASTIC_REGION,true);
	}
	else{
	  mThermalVariables.PlasticDissipation = 0;
	  mThermalVariables.DeltaPlasticDissipation = 0;
	}


}


//***************************CALCULATE THERMAL DISSIPATION****************************
//************************************************************************************

void NonLinearAssociativePlasticFlowRule::CalculateThermalDissipation( YieldCriterion::Parameters& rCriterionParameters, ThermalVariables& rThermalVariables )
{

      //1.- Thermal Dissipation:

      mThermalVariables.PlasticDissipation = mpYieldCriterion->CalculatePlasticDissipation( mThermalVariables.PlasticDissipation, rCriterionParameters);


      //2.- Thermal Dissipation Increment:

      mThermalVariables.DeltaPlasticDissipation = mpYieldCriterion->CalculateDeltaPlasticDissipation( mThermalVariables.DeltaPlasticDissipation, rCriterionParameters );

}


//***************************CALCULATE THERMAL DISSIPATION****************************
//************************************************************************************

void NonLinearAssociativePlasticFlowRule::CalculateImplexThermalDissipation( YieldCriterion::Parameters& rCriterionParameters )
{

      //1.- Thermal Dissipation:

      mThermalVariables.PlasticDissipation = mpYieldCriterion->CalculateImplexPlasticDissipation( mThermalVariables.PlasticDissipation, rCriterionParameters );

      //2.- Thermal Dissipation Increment:

      mThermalVariables.DeltaPlasticDissipation = mpYieldCriterion->CalculateImplexDeltaPlasticDissipation( mThermalVariables.DeltaPlasticDissipation, rCriterionParameters );




}

//***************************UPDATE STRESS CONFIGURATION *****************************
//************************************************************************************

void NonLinearAssociativePlasticFlowRule::UpdateConfiguration( RadialReturnVariables& rReturnMappingVariables, Matrix & rIsoStressMatrix )
{
	//Back Stress update

	//Plastic Strain Update
        if( rReturnMappingVariables.NormIsochoricStress > 0 ){

	  //Stress Update:
	  double Auxiliar   = 2.0 * rReturnMappingVariables.LameMu_bar * rReturnMappingVariables.DeltaGamma;

	  Matrix Normal     = rIsoStressMatrix * ( 1.0 / rReturnMappingVariables.NormIsochoricStress );

	  rIsoStressMatrix -= ( Normal * Auxiliar );

	}

}

//***************************UPDATE INTERNAL VARIABLES********************************
//************************************************************************************

bool NonLinearAssociativePlasticFlowRule::UpdateInternalVariables( RadialReturnVariables& rReturnMappingVariables )
{

	mInternalVariables.EquivalentPlasticStrainOld  = mInternalVariables.EquivalentPlasticStrain;

	mInternalVariables.DeltaPlasticStrain          = sqrt(2.0/3.0) * rReturnMappingVariables.DeltaGamma;

	mInternalVariables.EquivalentPlasticStrain    += mInternalVariables.DeltaPlasticStrain;

	mInternalVariables.DeltaPlasticStrain         *= ( 1.0/rReturnMappingVariables.DeltaTime );

	//update thermal variables
	// mThermalVariables = rReturnMappingVariables.Thermal;

	// mInternalVariables.print();
	// mThermalVariables.print();

	return true;
}


//**************CALCULATE SCALING FACTORS FOR THE ELASTO PLASTIC MODULI***************
//************************************************************************************

void NonLinearAssociativePlasticFlowRule::CalculateScalingFactors(const RadialReturnVariables& rReturnMappingVariables, PlasticFactors& rScalingFactors )
{

 	//1.-Identity build
	Matrix Identity   = identity_matrix<double> (3);

	//2.-Auxiliar matrices
	rScalingFactors.Normal      = rReturnMappingVariables.TrialIsoStressMatrix * ( 1.0 / rReturnMappingVariables.NormIsochoricStress );

	Matrix Norm_Normal          = prod( rScalingFactors.Normal, trans(rScalingFactors.Normal) );

	double Trace_Norm_Normal    = Norm_Normal( 0, 0 ) + Norm_Normal( 1, 1 )	+ Norm_Normal( 2, 2 );

	rScalingFactors.Dev_Normal  = Norm_Normal;
	rScalingFactors.Dev_Normal -= (1.0/3.0) * Trace_Norm_Normal * Identity;


	//3.-Auxiliar constants
	double EquivalentPlasticStrain = mInternalVariables.EquivalentPlasticStrain + sqrt(2.0/3.0) * rReturnMappingVariables.DeltaGamma;
	double DeltaHardening = 0;

	if( rReturnMappingVariables.Options.Is(IMPLEX_ACTIVE) )
	  {
	    rScalingFactors.Beta0 = 0;

	    rScalingFactors.Beta1 = 2.0 * rReturnMappingVariables.LameMu_bar * rReturnMappingVariables.DeltaGamma / rReturnMappingVariables.NormIsochoricStress;

	    rScalingFactors.Beta2 = (2.0/3.0) * rReturnMappingVariables.NormIsochoricStress * rReturnMappingVariables.DeltaGamma / (rReturnMappingVariables.LameMu_bar) ;

	    rScalingFactors.Beta3 = ( - rScalingFactors.Beta1 + rScalingFactors.Beta2 );

	    rScalingFactors.Beta4 = ( - rScalingFactors.Beta1 ) * rReturnMappingVariables.NormIsochoricStress / ( rReturnMappingVariables.LameMu_bar ) ;

	  }
	else
	  {

	    HardeningLaw::Parameters HardeningParameters;
	    HardeningParameters.SetTemperature(rReturnMappingVariables.Temperature);
	    HardeningParameters.SetEquivalentPlasticStrain(EquivalentPlasticStrain);
	    HardeningParameters.SetDeltaGamma(rReturnMappingVariables.DeltaGamma);
	    HardeningParameters.SetDeltaTime(rReturnMappingVariables.DeltaTime);

	    if( rReturnMappingVariables.Options.Is(PLASTIC_RATE_REGION) )
	      HardeningParameters.SetRateFactor(1);
	    else if ( rReturnMappingVariables.Options.IsNot(PLASTIC_RATE_REGION) )
	      HardeningParameters.SetRateFactor(0);

	    DeltaHardening = mpYieldCriterion->GetHardeningLaw().CalculateDeltaHardening( DeltaHardening, HardeningParameters );

	    rScalingFactors.Beta0 = 1.0 + DeltaHardening/(3.0 * rReturnMappingVariables.LameMu_bar);

	    rScalingFactors.Beta1 = 2.0 * rReturnMappingVariables.LameMu_bar * rReturnMappingVariables.DeltaGamma / rReturnMappingVariables.NormIsochoricStress;

	    rScalingFactors.Beta2 = ( ( 1.0 - ( 1.0 / rScalingFactors.Beta0 ) ) * (2.0/3.0) * rReturnMappingVariables.NormIsochoricStress * rReturnMappingVariables.DeltaGamma )/(rReturnMappingVariables.LameMu_bar) ;

	    rScalingFactors.Beta3 = ( ( 1.0 / rScalingFactors.Beta0 ) - rScalingFactors.Beta1 + rScalingFactors.Beta2 );

	    rScalingFactors.Beta4 = ( ( 1.0 / rScalingFactors.Beta0 ) - rScalingFactors.Beta1 ) * rReturnMappingVariables.NormIsochoricStress / ( rReturnMappingVariables.LameMu_bar ) ;

	  }

	//std::cout<<"FACTORS:: Beta0 "<<rScalingFactors.Beta0<<" Beta 1 "<<rScalingFactors.Beta1<<" Beta2 "<<rScalingFactors.Beta2<<" Beta 3 "<<rScalingFactors.Beta3<<" Beta4 "<<rScalingFactors.Beta4<<std::endl;
}



void NonLinearAssociativePlasticFlowRule::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, FlowRule )
}

void NonLinearAssociativePlasticFlowRule::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, FlowRule )
}


}  // namespace Kratos.
