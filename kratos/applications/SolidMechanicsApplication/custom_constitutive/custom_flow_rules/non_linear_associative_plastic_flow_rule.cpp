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
#include "solid_mechanics_application.h"
#include "custom_constitutive/custom_flow_rules/non_linear_associative_plastic_flow_rule.hpp"

namespace Kratos
{

//*******************************CONSTRUCTOR******************************************
//************************************************************************************

NonLinearAssociativePlasticFlowRule::NonLinearAssociativePlasticFlowRule()
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


//********************************DESTRUCTOR******************************************
//************************************************************************************

NonLinearAssociativePlasticFlowRule::~NonLinearAssociativePlasticFlowRule()
{
}

/// Operations.


//***************************CALCULATE STRESS NORM ***********************************
//************************************************************************************

double& NonLinearAssociativePlasticFlowRule::CalculateNormStress ( Matrix & rStressMatrix, double& rNormStress )
{
	  
	rNormStress = sqrt(rStressMatrix( 0 , 0 )*rStressMatrix( 0 , 0 )+
			   rStressMatrix( 1 , 1 )*rStressMatrix( 1 , 1 )+
			   rStressMatrix( 2 , 2 )*rStressMatrix( 2 , 2 )+
			   2.0 * rStressMatrix( 0 , 1 )*rStressMatrix( 0 , 1 ) );

	return rNormStress;
}


//***************************CALCULATE RADIAL RETURN MAPPING**************************
//************************************************************************************


bool NonLinearAssociativePlasticFlowRule::CalculateReturnMapping( RadialReturnVariables& rReturnMappingVariables, Matrix& rIsoStressMatrix )
{
	  
	//0.- Initialize Variables
	rReturnMappingVariables.Control.PlasticRegion = false;
	InternalVariables PlasticVariables = mInternalVariables;
		
	//1.-Isochoric stress norm
	rReturnMappingVariables.NormIsochoricStress = CalculateNormStress( rIsoStressMatrix, rReturnMappingVariables.NormIsochoricStress );

	//2.- Check yield condition
	rReturnMappingVariables.TrialStateFunction = mpYieldCriterion->CalculateYieldCondition( rReturnMappingVariables.TrialStateFunction, rReturnMappingVariables.NormIsochoricStress, PlasticVariables.EquivalentPlasticStrain, rReturnMappingVariables.Temperature);

	if( rReturnMappingVariables.Control.ImplexActive == true ) 
	  {

	    CalculateImplexReturnMapping ( rReturnMappingVariables, PlasticVariables, rIsoStressMatrix );

	  }
	else{

	  if( rReturnMappingVariables.TrialStateFunction <= 0 )
	    {

	      PlasticVariables.DeltaPlasticStrain = 0;
	      rReturnMappingVariables.Control.PlasticRegion = false;

	    }
	  else
	    {

	      //3.- Calculate the consistency condition
	      bool converged = this->CalculateConsistencyCondition( rReturnMappingVariables, PlasticVariables );
	    
	      if(!converged)
		std::cout<<" ConstitutiveLaw did not converge "<<std::endl;


	      //4.- Update back stress, plastic strain and stress
	      UpdateConfiguration( rReturnMappingVariables, rIsoStressMatrix );
	    

	      //5.- Calculate thermal dissipation and delta thermal dissipation
	      this->CalculateThermalDissipation( rReturnMappingVariables, PlasticVariables );   
	    
	      rReturnMappingVariables.Control.PlasticRegion = true;

	    }

	}
	
	rReturnMappingVariables.Control.ReturnMappingComputed = true;

	return 	rReturnMappingVariables.Control.PlasticRegion;
}

//***************************CALCULATE LOCAL NEWTON PROCEDURE*************************
//************************************************************************************


bool NonLinearAssociativePlasticFlowRule::CalculateConsistencyCondition( RadialReturnVariables& rReturnMappingVariables, InternalVariables& rPlasticVariables )
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
	double InitialEquivalentPlasticStrain = rPlasticVariables.EquivalentPlasticStrain;
   
	while ( fabs(StateFunction)>=Tolerance && iter<=MaxIterations)
	{
		//Calculate Delta State Function:
		DeltaStateFunction = mpYieldCriterion->CalculateDeltaStateFunction( DeltaStateFunction, rReturnMappingVariables.LameMu_bar, rPlasticVariables.EquivalentPlasticStrain, rReturnMappingVariables.Temperature );

		//Calculate DeltaGamma:
		DeltaDeltaGamma  = StateFunction/DeltaStateFunction;
		rReturnMappingVariables.DeltaGamma += DeltaDeltaGamma;
	       
		//Update Equivalent Plastic Strain:
		rPlasticVariables.DeltaPlasticStrain       = sqrt(2.0/3.0) * rReturnMappingVariables.DeltaGamma;
		rPlasticVariables.EquivalentPlasticStrain  = InitialEquivalentPlasticStrain + rPlasticVariables.DeltaPlasticStrain;
	       
		//Calculate State Function:
		StateFunction = mpYieldCriterion->CalculateStateFunction( StateFunction, rReturnMappingVariables.NormIsochoricStress, rReturnMappingVariables.DeltaGamma, rReturnMappingVariables.LameMu_bar, rPlasticVariables.EquivalentPlasticStrain, InitialEquivalentPlasticStrain, rReturnMappingVariables.Temperature );
		

		iter++;
	}
	   

	if(iter>MaxIterations)
	  return false;


	return true;	
}


//***************************CALCULATE IMPLEX RETURN MAPPING**************************
//************************************************************************************

void NonLinearAssociativePlasticFlowRule::CalculateImplexReturnMapping( RadialReturnVariables& rReturnMappingVariables, InternalVariables& rPlasticVariables, Matrix& rIsoStressMatrix )
{

  
        //1.-Computation of the plastic Multiplier
        rReturnMappingVariables.DeltaGamma = sqrt(3.0/2.0) * ( rPlasticVariables.EquivalentPlasticStrain - rPlasticVariables.EquivalentPlasticStrainOld );
	
	//2.- Update back stress, plastic strain and stress
	UpdateConfiguration( rReturnMappingVariables, rIsoStressMatrix );

	//3.- Calculate thermal dissipation and delta thermal dissipation
	if( rReturnMappingVariables.DeltaGamma > 0 ){
	  
	  this->CalculateImplexThermalDissipation( rReturnMappingVariables, rPlasticVariables );
	  rReturnMappingVariables.Control.PlasticRegion = true;
	}
	else{
	  mThermalVariables.PlasticDissipation = 0;
	  mThermalVariables.DeltaPlasticDissipation = 0;
	} 
  

}


//***************************CALCULATE THERMAL DISSIPATION****************************
//************************************************************************************

void NonLinearAssociativePlasticFlowRule::CalculateThermalDissipation( RadialReturnVariables& rReturnMappingVariables, InternalVariables& rPlasticVariables )
{

      //1.- Thermal Dissipation:
 
      mThermalVariables.PlasticDissipation = mpYieldCriterion->CalculatePlasticDissipation( mThermalVariables.PlasticDissipation, rReturnMappingVariables.DeltaGamma, rReturnMappingVariables.TimeStep, rPlasticVariables.EquivalentPlasticStrain, rReturnMappingVariables.Temperature );
  

      //2.- Thermal Dissipation Increment:

      mThermalVariables.DeltaPlasticDissipation = mpYieldCriterion->CalculateDeltaPlasticDissipation( mThermalVariables.DeltaPlasticDissipation, rReturnMappingVariables.DeltaGamma, rReturnMappingVariables.TimeStep, rReturnMappingVariables.LameMu_bar, rPlasticVariables.EquivalentPlasticStrain, rReturnMappingVariables.Temperature );
		    		    
  
}


//***************************CALCULATE THERMAL DISSIPATION****************************
//************************************************************************************

void NonLinearAssociativePlasticFlowRule::CalculateImplexThermalDissipation( RadialReturnVariables& rReturnMappingVariables, InternalVariables& rPlasticVariables )
{
 
      //1.- Thermal Dissipation:
	
      mThermalVariables.PlasticDissipation = mpYieldCriterion->CalculateImplexPlasticDissipation( mThermalVariables.PlasticDissipation, rReturnMappingVariables.DeltaGamma, rReturnMappingVariables.TimeStep, rPlasticVariables.EquivalentPlasticStrain, rReturnMappingVariables.Temperature );
  
      //2.- Thermal Dissipation Increment:
      
      mThermalVariables.DeltaPlasticDissipation = mpYieldCriterion->CalculateImplexDeltaPlasticDissipation( mThermalVariables.DeltaPlasticDissipation, rReturnMappingVariables.DeltaGamma, rReturnMappingVariables.TimeStep, rReturnMappingVariables.LameMu_bar, rPlasticVariables.EquivalentPlasticStrain, rReturnMappingVariables.Temperature );


   
  
}

//***************************UPDATE STRESS CONFIGURATION *****************************
//************************************************************************************

void NonLinearAssociativePlasticFlowRule::UpdateConfiguration( RadialReturnVariables& rReturnMappingVariables, Matrix & rIsoStressMatrix )
{
	//Back Stress update
        
        //std::cout<< " ElasticIsoStress "<<rIsoStressMatrix<<std::endl;

	//Plastic Strain Update
        if( rReturnMappingVariables.NormIsochoricStress > 0 ){
    
	  //Stress Update: 
	  double Auxiliar   = 2.0 * rReturnMappingVariables.LameMu_bar * rReturnMappingVariables.DeltaGamma;

	  Matrix Normal     = rIsoStressMatrix * ( 1.0 / rReturnMappingVariables.NormIsochoricStress );
	  
	  rIsoStressMatrix -= ( Normal * Auxiliar );
	  
	}

	//std::cout<< " PlasticIsoStress "<<rIsoStressMatrix<<std::endl;
}

//***************************UPDATE INTERNAL VARIABLES********************************
//************************************************************************************

bool NonLinearAssociativePlasticFlowRule::UpdateInternalVariables( RadialReturnVariables& rReturnMappingVariables )
{
	
	mInternalVariables.EquivalentPlasticStrainOld  = mInternalVariables.EquivalentPlasticStrain;

	mInternalVariables.DeltaPlasticStrain          = sqrt(2.0/3.0) * rReturnMappingVariables.DeltaGamma;

	mInternalVariables.EquivalentPlasticStrain    += mInternalVariables.DeltaPlasticStrain;

	mInternalVariables.DeltaPlasticStrain         *= ( 1.0/rReturnMappingVariables.TimeStep );
 	
	return true;
}


//**************CALCULATE SCALING FACTORS FOR THE ELASTO PLASTIC MODULI***************
//************************************************************************************

void NonLinearAssociativePlasticFlowRule::CalculateScalingFactors(const RadialReturnVariables& rReturnMappingVariables, PlasticFactors& rScalingFactors )
{
  
 	//1.-Identity build
	Matrix IdentityMatrix  = identity_matrix<double> (3);

	//2.-Auxiliar matrices
	Matrix IsoStressMatrix      = MathUtils<double>::StressVectorToTensor( rReturnMappingVariables.TrialIsoStressVector );
	rScalingFactors.Normal      = IsoStressMatrix * ( 1.0 / rReturnMappingVariables.NormIsochoricStress );

	Matrix Norm_Normal          = prod( rScalingFactors.Normal, trans(rScalingFactors.Normal) );
	double Trace_Norm_Normal    = Norm_Normal( 0, 0 ) + Norm_Normal( 1, 1 ) + Norm_Normal( 2, 2 );

	rScalingFactors.Dev_Normal  = Norm_Normal;
	rScalingFactors.Dev_Normal -= (1.0/3.0) * Trace_Norm_Normal * IdentityMatrix;


	//3.-Auxiliar constants
	double EquivalentPlasticStrain = mInternalVariables.EquivalentPlasticStrain + sqrt(2.0/3.0) * rReturnMappingVariables.DeltaGamma;
	double DeltaHardening = 0;
	DeltaHardening = mpHardeningLaw->CalculateDeltaHardening( DeltaHardening, EquivalentPlasticStrain, rReturnMappingVariables.Temperature );

	rScalingFactors.Beta0 = 1.0 + DeltaHardening/(3.0 * rReturnMappingVariables.LameMu_bar);
		
	rScalingFactors.Beta1 = 2.0 * rReturnMappingVariables.LameMu_bar * rReturnMappingVariables.DeltaGamma / rReturnMappingVariables.NormIsochoricStress;
		
	rScalingFactors.Beta2 = ( ( 1.0 - ( 1.0 / rScalingFactors.Beta0 ) ) * (2.0/3.0) * rReturnMappingVariables.NormIsochoricStress * rReturnMappingVariables.DeltaGamma )/(rReturnMappingVariables.LameMu_bar) ;
		
	rScalingFactors.Beta3 = ( ( 1.0 / rScalingFactors.Beta0 ) - rScalingFactors.Beta1 + rScalingFactors.Beta2 );
		
	rScalingFactors.Beta4 = ( ( 1.0 / rScalingFactors.Beta0 ) - rScalingFactors.Beta1 ) * rReturnMappingVariables.NormIsochoricStress / ( rReturnMappingVariables.LameMu_bar ) ;
	
	//std::cout<<"FACTORS:: Beta0 "<<rScalingFactors.Beta0<<" Beta 1 "<<rScalingFactors.Beta1<<" Beta2 "<<rScalingFactors.Beta2<<" Beta 3 "<<rScalingFactors.Beta3<<" Beta4 "<<rScalingFactors.Beta4<<std::endl;
};



void NonLinearAssociativePlasticFlowRule::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, FlowRule );
}

void NonLinearAssociativePlasticFlowRule::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, FlowRule );
}


}  // namespace Kratos.
