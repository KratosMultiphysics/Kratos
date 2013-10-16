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
#include "solid_mechanics_application.h"
#include "custom_constitutive/custom_yield_criteria/mises_huber_yield_criterion.hpp"

namespace Kratos
{

//*******************************CONSTRUCTOR******************************************
//************************************************************************************

MisesHuberYieldCriterion::MisesHuberYieldCriterion()
	:YieldCriterion()
{
   
}


//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

MisesHuberYieldCriterion& MisesHuberYieldCriterion::operator=(MisesHuberYieldCriterion const& rOther)
{
   YieldCriterion::operator=(rOther);
   return *this;
}

//*******************************COPY CONSTRUCTOR*************************************
//************************************************************************************

MisesHuberYieldCriterion::MisesHuberYieldCriterion(MisesHuberYieldCriterion const& rOther)
	:YieldCriterion(rOther)
{

}


//********************************DESTRUCTOR******************************************
//************************************************************************************

MisesHuberYieldCriterion::~MisesHuberYieldCriterion()
{
}

/// Operations.


//***************************CALCULATE YIELD CONDITION********************************
//************************************************************************************

double& MisesHuberYieldCriterion::CalculateYieldCondition(double & rStateFunction, const double& rNormStress, const double& rAlpha, double rTemperature)
{
	double Hardening = 0;

	Hardening = mpHardeningLaw->CalculateHardening(Hardening, rAlpha, rTemperature);
		
	rStateFunction = rNormStress - sqrt(2.0/3.0) * Hardening;
		
	return rStateFunction;
};


//***************************CALCULATE YIELD CONDITION********************************
//************************************************************************************


double& MisesHuberYieldCriterion::CalculateYieldCondition(double & rStateFunction, const Matrix& rStressMatrix, const double& rAlpha, double rTemperature)
{
	double	NormStress = sqrt(rStressMatrix( 0 , 0 )*rStressMatrix( 0 , 0 )+
				  rStressMatrix( 1 , 1 )*rStressMatrix( 1 , 1 )+
				  rStressMatrix( 2 , 2 )*rStressMatrix( 2 , 2 )+
				  2.0 * rStressMatrix( 0 , 1 )*rStressMatrix( 0 , 1 ) );

	rStateFunction = this->CalculateYieldCondition( rStateFunction, NormStress, rAlpha, rTemperature );

	return rStateFunction;
};



//***************************CALCULATE STATE FUNCTION ********************************
//************************************************************************************

double& MisesHuberYieldCriterion::CalculateStateFunction(double & rStateFunction,const double& rNormStress, const double& rDeltaGamma, const double& rLameMu_bar, const double& rAlpha, const double& rAlphaOld, double rTemperature)
{
	double InitialKinematicHardening = 0;
	double KinematicHardening = 0;
	double IsotropicHardening = 0;

	InitialKinematicHardening = mpHardeningLaw->CalculateKinematicHardening(InitialKinematicHardening, rAlphaOld, rTemperature );
	
	KinematicHardening = mpHardeningLaw->CalculateKinematicHardening(KinematicHardening, rAlpha, rTemperature);

	IsotropicHardening = mpHardeningLaw->CalculateIsotropicHardening(IsotropicHardening, rAlpha, rTemperature);		

	rStateFunction = rNormStress - 2.0 * rLameMu_bar * rDeltaGamma - sqrt(2.0/3.0) * ( IsotropicHardening + ( KinematicHardening - InitialKinematicHardening ));
		
	return rStateFunction;
};


//***************************CALCULATE STATE FUNCTION ********************************
//************************************************************************************

double& MisesHuberYieldCriterion::CalculateDeltaStateFunction(double & rDeltaStateFunction, const double& rLameMu_bar, const double& rAlpha, double rTemperature)
{

	double DeltaKinematicHardening = 0;
	double DeltaIsotropicHardening = 0;


	DeltaKinematicHardening = mpHardeningLaw->CalculateDeltaKinematicHardening(DeltaKinematicHardening, rAlpha, rTemperature);

	DeltaIsotropicHardening = mpHardeningLaw->CalculateDeltaIsotropicHardening(DeltaIsotropicHardening, rAlpha, rTemperature);		

	rDeltaStateFunction = 2.0 * rLameMu_bar + (2.0/3.0) * (DeltaKinematicHardening + DeltaIsotropicHardening);
		
	return rDeltaStateFunction;
};


//***************************CALCULATE PLASTIC DISSIPATION****************************
//************************************************************************************

double& MisesHuberYieldCriterion::CalculatePlasticDissipation(double & rPlasticDissipation, const double& rDeltaGamma, const double& rDeltaTime, const double& rAlpha, const double& rTemperature)
{
	rPlasticDissipation = 0;
	return rPlasticDissipation;
};


//**********************CALCULATE DELTA PLASTIC DISSIPATION***************************
//************************************************************************************

double& MisesHuberYieldCriterion::CalculateDeltaPlasticDissipation(double & rDeltaPlasticDissipation, const double& rDeltaGamma, const double& rDeltaTime, const double& rLameMu_bar, const double& rAlpha, const double& rTemperature)
{
	rDeltaPlasticDissipation = 0;
	return rDeltaPlasticDissipation;
};


//***************************CALCULATE PLASTIC DISSIPATION****************************
//************************************************************************************

double& MisesHuberYieldCriterion::CalculateImplexPlasticDissipation(double & rPlasticDissipation, const double& rDeltaGamma, const double& rDeltaTime, const double& rAlpha, const double& rTemperature)
{
	rPlasticDissipation = 0;
	return rPlasticDissipation;
};


//**********************CALCULATE DELTA PLASTIC DISSIPATION***************************
//************************************************************************************

double& MisesHuberYieldCriterion::CalculateImplexDeltaPlasticDissipation(double & rDeltaPlasticDissipation, const double& rDeltaGamma, const double& rDeltaTime, const double& rLameMu_bar, const double& rAlpha, const double& rTemperature)
{
	rDeltaPlasticDissipation = 0;
	return rDeltaPlasticDissipation;
};

void MisesHuberYieldCriterion::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, YieldCriterion );
}

void MisesHuberYieldCriterion::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, YieldCriterion );
}


}  // namespace Kratos.
