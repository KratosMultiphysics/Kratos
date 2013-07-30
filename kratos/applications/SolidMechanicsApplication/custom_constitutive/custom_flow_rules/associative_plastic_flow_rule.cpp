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
#include "custom_constitutive/custom_flow_rules/associative_plastic_flow_rule.hpp"

#include "solid_mechanics_application.h"

namespace Kratos
{

//*******************************CONSTRUCTOR******************************************
//************************************************************************************

AssociativePlasticFlowRule::AssociativePlasticFlowRule()
	:FlowRule()
{
   
}


//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

AssociativePlasticFlowRule& AssociativePlasticFlowRule::operator=(AssociativePlasticFlowRule const& rOther)
{
   FlowRule::operator=(rOther);
   return *this;
}

//*******************************COPY CONSTRUCTOR*************************************
//************************************************************************************

AssociativePlasticFlowRule::AssociativePlasticFlowRule(AssociativePlasticFlowRule const& rOther)
	:HardeningLaw(rOther)
{

}


//********************************DESTRUCTOR******************************************
//************************************************************************************

AssociativePlasticFlowRule::~AssociativePlasticFlowRule()
{
}

/// Operations.

//*******************************CALCULATE TOTAL HARDENING****************************
//************************************************************************************

double& AssociativePlasticFlowRule::CalculateHardening(double &Hardening, double & alpha)
{
	
	//linear hardening properties
	const double& YieldStress                 =  GetProperties()[YIELD_STRESS];
	const double& KinematicHardeningConstant  =  GetProperties()[KINEMATIC_HARDENING];
	
	//exponential saturation properties
   	const double& K_reference           =  GetProperties()[REFERENCE_HARDENING];
	const double& K_infinity            =  GetProperties()[INFINITY_HARDENING];
	const double& Delta                 =  GetProperties()[HARDENING_EXPONENT];

	//Linear Hardening law:
	Hardening  = YieldStress + mTheta *  KinematicHardeningConstant;
	
	//Exponential Saturation:
	Hardening += (K_infinity-K_reference) * (1.0 - exp( (-1.0) * Delta * alpha ) );
	
	return Hardening;

}
  
//*******************************CALCULATE ISOTROPIC HARDENING************************
//************************************************************************************

double& AssociativePlasticFlowRule::CalculateIsotropicHardening(double &IsotropicHardening, double & alpha)
{

	//linear hardening properties
	const double& YieldStress                 =  GetProperties()[YIELD_STRESS];
	const double& KinematicHardeningConstant  =  GetProperties()[KINEMATIC_HARDENING];
	
	//exponential saturation properties
   	const double& K_reference           =  GetProperties()[REFERENCE_HARDENING];
	const double& K_infinity            =  GetProperties()[INFINITY_HARDENING];
	const double& Delta                 =  GetProperties()[HARDENING_EXPONENT];

	//Linear Hardening law: (mTheta = 1)
	IsotropicHardening  = YieldStress + KinematicHardeningConstant;
	
	//Exponential Saturation:
	IsotropicHardening += (K_infinity-K_reference) * (1.0 - exp( (-1.0) * Delta * alpha ) );
	
	return IsotropicHardening;	


}

//*******************************CALCULATE KINEMATIC HARDENING************************
//************************************************************************************

double& AssociativePlasticFlowRule::CalculateKinematicHardening(double &KinematicHardening, double & alpha)
{
	//linear hardening properties
	const double& KinematicHardeningConstant  =  GetProperties()[KINEMATIC_HARDENING];
	
	//Linear Hardening law:
	KinematicHardening  = (1.0 - mTheta) * KinematicHardeningConstant;
	
	return KinematicHardening;
}



//*******************************CALCULATE HARDENING DERIVATIVE***********************
//************************************************************************************

double& AssociativePlasticFlowRule::CalculateDeltaHardening(double &DeltaHardening, double & alpha)
{
      	//linear hardening properties
	const double& KinematicHardeningConstant  =  GetProperties()[KINEMATIC_HARDENING];
	
	//exponential saturation properties
   	const double& K_reference           =  GetProperties()[REFERENCE_HARDENING];
	const double& K_infinity            =  GetProperties()[INFINITY_HARDENING];
	const double& Delta                 =  GetProperties()[HARDENING_EXPONENT];

	//Linear Hardening law: (mTheta = 1)
	DeltaHardening  = mTheta * KinematicHardeningConstant;
	
	//Exponential Saturation:
	DeltaHardening += Delta * (K_infinity-K_reference) * ( exp( (-1.0) * Delta * alpha ) );
	
	return DeltaHardening;	
}

//***************************CALCULATE ISOTROPIC HARDENING DERIVATIVE*****************
//************************************************************************************

double& AssociativePlasticFlowRule::CalculateDeltaIsotropicHardening(double &DeltaIsotropicHardening, double & alpha)
{
       	//linear hardening properties
	const double& KinematicHardeningConstant  =  GetProperties()[KINEMATIC_HARDENING];
	
	//exponential saturation properties
   	const double& K_reference           =  GetProperties()[REFERENCE_HARDENING];
	const double& K_infinity            =  GetProperties()[INFINITY_HARDENING];
	const double& Delta                 =  GetProperties()[HARDENING_EXPONENT];

	//Linear Hardening law: (mTheta = 1)
	DeltaIsotropicHardening  = KinematicHardeningConstant;
	
	//Exponential Saturation:
	DeltaIsotropicHardening += Delta * (K_infinity-K_reference) * ( exp( (-1.0) * Delta * alpha ) );
	
	return DeltaIsotropicHardening;	

}


//***************************CALCULATE STRESS NORM ***********************************
//************************************************************************************

double& AssociativePlasticFlowRule::CalculateNormStress ( Matrix & rStressMatrix, double& rNormStress )
{
	  
	rNormStress = sqrt(rStressMatrix( 0 , 0 )*rStressMatrix( 0 , 0 )+
			   rStressMatrix( 1 , 1 )*rStressMatrix( 1 , 1 )+
			   rStressMatrix( 2 , 2 )*rStressMatrix( 2 , 2 )+
			   2.0 * rStressMatrix( 0 , 1 )*rStressMatrix( 0 , 1 ) );

	return rNormStress;

}


//***************************CALCULATE RADIAL RETURN MAPPING**************************
//************************************************************************************


bool AssociativePlasticFlowRule::CalculateReturnMapping( RadialReturnVariables& rReturnMappingVariables, Matrix& rIsoStressMatrix )
{
	  
	//0.- Initialize Variables
	bool Plasticity = false;
	InternalVariables PlasticVariables = mInternalVariables;
		
	//1.-Isochoric stress norm
	rReturnMappingVariables.NormIsochoricStress = CalculateNormStress( rIsoStressMatrix, ReturnMappingVariables.NormIsochoricStress );

	//2.- Check yield condition
	rReturnMappingVariables.TrialStateFunction = YieldCriterion->CheckYieldCondition( PlasticVariables, rReturnMappingVariables.NormIsochoricStress, rIsoStressMatrix )

		if( TrialStateFunction > 0 )
		{
			rPlasticVariables.DeltaPlasticStrain = 0;
			Plasticity = false;
		}
		else
		{
			bool converged = this->CalculateConsistencyCondition( PlasticVariables );

			if(!converged)
				std::cout<<" ConstitutiveLaw did not converge "<<std::endl;

			//3.- Update back stress, plastic strain and stress
			this->UpdateConfiguration( rReturnMappingVariables, rIsoStressMatrix );

	
			Plasticity = true;
		}


	return Plasticity;
}

//***************************CALCULATE LOCAL NEWTON PROCEDURE*************************
//************************************************************************************


bool AssociativePlasticFlowRule::CalculateConstistencyCondition( RadialReturnVariables& rReturnMappingVariables, InternalVariables& rPlasticVariables )
{
	//Set convergence parameters
	unsigned int iter    = 0;
	double Tolerance     = 1e-5;
	double MaxIterations = 50;

	//Initialize Parameters used in the determination of the Delta Plastic Strain
	double KinematicHardening = 0;
	double IsotropicHardening = 0;

	double DeltaKinematicHardening = 0;
	double DeltaIsotropicHardening = 0;

	//start
	double DeltaFunction = 0;
	rReturnMappingVariables.DeltaGamma    = 0;

	double StateFunction                  = rReturnMappingVariables.TrialStateFunction;
	double InitialEquivalentPlasticStrain = rPlasticVariables.EquivalentPlasticStrain ;
	double InitialKinematicHardening      = 0;
	InitialKinematicHardening = rHardeningLaw.CalculateKinematicHardening( InitialKinematicHardening, InitialEquivalentPlasticStrain );

   
	while ( fabs(StateFunction)>=Tolerance && iter<=MaxIterations)
	{
		//Calculate DeltaFunction:
		DeltaKinematicHardening = mpHardeningLaw->CalculateDeltaKinematicHardening( DeltaKinematicHardening, rPlasticVariables.EquivalentPlasticStrain );
		DeltaIsotropicHardening = mpHardeningLaw->CalculateDeltaIsotropicHardening( DeltaIsotropicHardening, rPlasticVariables.EquivalentPlasticStrain );

		DeltaFunction  = 2.0 * rReturnMappingVariables.LameMu_bar;
		DeltaFunction += (DeltaKinematicHardening + DeltaIsotropicHardening) * (2.0/3.0);

		//Calculate DeltaGamma:
		DeltaDeltaGamma  = StateFunction/DeltaFunction;
		rReturnMappingVariables.DeltaGamma += DeltaDeltaGamma;
	       
		//Update Equivalent Plastic Strain:
		rPlasticVariables.DeltaPlasticStrain       = sqrt(2.0/3.0) * rReturnMappingVariables.DeltaGamma;
		rPlasticVariables.EquivalentPlasticStrain  = InitialEquivalentPlasticStrain + rPlasticVariables.DeltaPlasticStrain;
	       
		//Calculate State Function:
		KinematicHardening = mpHardeningLaw->CalculateKinematicHardening( KinematicHardening, rPlasticVariables.EquivalentPlasticStrain );
		IsotropicHardening = mpHardeningLaw->CalculateIsotropicHardening( IsotropicHardening, rPlasticVariables.EquivalentPlasticStrain );

		StateFunction  = rReturnMappingVariables.NormIsochoricStress;
		StateFunction -= 2.0 * LameMu_bar * rReturnMappingVariables.DeltaGamma;
		StateFunction -= sqrt(2.0/3.0) * ( IsotropicHardening + ( KinematicHardening - InitialKinematicHardening ));

		iter++;
	}
	   

	if(iter>MaxIterations)
		return false;

	return true;	
}

//***************************UPDATE STRESS CONFIGURATION *****************************
//************************************************************************************

void AssociativePlasticFlowRule::UpdateConfiguration( RadialReturnVariables& rReturnMappingVariables, Matrix & rIsoStressMatrix )
{
	//Back Stress update
        
	//Plastic Strain Update

	//Stress Update: 
	double Auxiliar   = 2.0 * rReturnMappingVariables.LameMu_bar * rReturnMappingVariables.DeltaGamma;

	Matrix Normal     = rIsoStressMatrix * ( 1.0 / rReturnMappingVariables.NormIsochoricStress );

	rIsoStressMatrix -= ( Normal * Auxiliar );


}

//***************************UPDATE INTERNAL VARIABLES********************************
//************************************************************************************

bool AssociativePlasticFlowRule::UpdateInternalVariables( RadialReturnVariables& rReturnMappingVariables )
{
	
	//mInternalVariables.EquivalentPlasticStrainOld  = PlasticVariables.EquivalentPlasticStrain;

	mInternalVariables.DeltaPlasticStrain          = sqrt(2.0/3.0) * rReturnMappingVariables.DeltaGamma;

	mInternalVariables.EquivalentPlasticStrain    += mInternalVariables.DeltaPlasticStrain;

	mInternalVariables.DeltaPlasticStrain         *= ( 1.0/rReturnMappingVariables.TimeStep );
 	
}


//**************CALCULATE SCALING FACTORS FOR THE ELASTO PLASTIC MODULI***************
//************************************************************************************

void AssociativePlasticFlowRule::CalculateScalingFactors( RadialReturnVariables& rReturnMappingVariables, const Matrix & rIsoStressMatrix, PlasticFactors& rScalingFactors )
{
	
	//1.-Identity build
	Matrix IdentityMatrix  = identity_matrix<double> (3);

	//3.-Particular Parameters
	double YieldStress           =  GetProperties()[YIELD_STRESS];
	double KinematicHardening    =  GetProperties()[KINEMATIC_HARDENING];
	double Delta                 =  GetProperties()[HARDENING_EXPONENT];
	double K_reference           =  GetProperties()[REFERENCE_HARDENING];
	double K_infinity            =  GetProperties()[INFINITY_HARDENING];

	const double& YoungModulus       = GetProperties()[YOUNG_MODULUS];
	const double& PoissonCoefficient = GetProperties()[POISSON_RATIO];
		
	//1.-Auxiliar matrices
	rScalingFactors.Normal      = rIsoStressMatrix * ( 1.0 / rReturnMappingVariables.NormIsochoricStress );

	Matrix Norm_Normal          = prod( rScalingFactors.Normal, trans(rScalingFactors.Normal) );
	double Trace_Norm_Normal    = Norm_Normal( 0, 0 ) + Norm_Normal( 1, 1 ) + Norm_Normal( 2, 2 );

	rScalingFactors.Dev_Normal  = Norm_Normal;
	rScalingFactors.Dev_Normal -= (1.0/3.0) * Trace_Norm_Normal * IdentityMatrix;


	//2.-Auxiliar constants
	double EquivalentPlasticStrain = mInternalVariables.EquivalentPlasticStrain + sqrt(2.0/3.0) * rReturnMappingVariables.DeltaGamma;
	double DeltaHardening = 0;
	DeltaHardening = mpHardeningLaw->CalculateDeltaHardening( DeltaHardening, EquivalentPlasticStrain );

	rScalingFactors.Beta0 = 1.0 + DeltaHardening/(3.0 * rReturnMappingVariables.LameMu_bar);
		
	rScalingFactors.Beta1 = 2.0 * LameMu_bar * DeltaGamma / rParameters.NormIsochoricStress;
		
	rScalingFactors.Beta2 = ( ( 1.0 - ( 1.0 / Beta0 ) ) * (2.0/3.0) * rReturnMappingVariables.NormIsochoricStress * rReturnMappingVariables.DeltaGamma )/(rReturnMappingVariables.LameMu_bar) ;
		
	rScalingFactors.Beta3 = ( ( 1.0 / Beta0 ) - Beta1 + Beta2 );
		
	rScalingFactors.Beta4 = ( ( 1.0 / Beta0 ) - Beta1 ) * rReturnMappingVariables.NormIsochoricStress / ( rReturnMappingVariables.LameMu_bar ) ;
	
};



void AssociativePlasticFlowRule::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, FlowRule );
}

void AssociativePlasticFlowRule::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, FlowRule );
}


}  // namespace Kratos.
