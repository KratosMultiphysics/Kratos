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
#include "custom_constitutive/custom_yield_criteria/mises_huber_yield_criterion.hpp"

#include "solid_mechanics_application.h"

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

double& MisesHuberYieldCriterion::CalculateYieldCondition(double & rStateFunction, const double& rNormStress, const double& rAlpha)
{
	double Hardening = 0;

	Hardening = mpHardeningLaw->CalculateHardening(Hardening,rAlpha);
		
	rStateFunction = rNormStress - sqrt(2.0/3.0) * Hardening;
		
	return rStateFunction;
};


//***************************CALCULATE YIELD CONDITION********************************
//************************************************************************************


double& MisesHuberYieldCriterion::CalculateYieldCondition(double & rStateFunction, const Matrix& rStressMatrix, const double& rAlpha)
{
	double	NormStress = sqrt(rStressMatrix( 0 , 0 )*rStressMatrix( 0 , 0 )+
				  rStressMatrix( 1 , 1 )*rStressMatrix( 1 , 1 )+
				  rStressMatrix( 2 , 2 )*rStressMatrix( 2 , 2 )+
				  2.0 * rStressMatrix( 0 , 1 )*rStressMatrix( 0 , 1 ) );

	rStateFunction = this->CalculateYieldCondition( rStateFunction, NormStress, rAlpha );

	return rStateFunction;
};



//***************************CALCULATE STATE FUNCTION ********************************
//************************************************************************************

double& MisesHuberYieldCriterion::CalculateStateFunction(double & rStateFunction,const double& rNormStress, const &DeltaGamma, const double& LameMu_bar, const double& rAlpha, const double& rAlphaOld)
{
	double InitialKinematicHardening = 0;
	double KinematicHardening = 0;
	double IsotropicHardening = 0;

	InitialKinematicHardening = mpHardeningLaw->CalculateKinematicHardening(InitialKinematicHardening,rAlphaOld);
	
	KinematicHardening = mpHardeningLaw->CalculateKinematicHardening(KinematicHardening,rAlpha);

	IsotropicHardening = mpHardeningLaw->CalculateIsotropicHardening(IsotropicHardening,rAlpha);		

	rStateFunction = rNormStress - 2.0 * LameMu_bar * DeltaGamma - sqrt(2.0/3.0) * ( IsotropicHardening + ( KinematicHardening - InitialKinematicHardening ));
		
	return rStateFunction;
};


//***************************CALCULATE STATE FUNCTION ********************************
//************************************************************************************

double& MisesHuberYieldCriterion::CalculateDeltaStateFunction(double & rDeltaStateFunction, const double& LameMu_bar, const double& rAlpha)
{

	double DeltaKinematicHardening = 0;
	double DeltaIsotropicHardening = 0;


	DeltaKinematicHardening = mpHardeningLaw->CalculateDeltaKinematicHardening(KinematicHardening,rAlpha);

	DeltaIsotropicHardening = mpHardeningLaw->CalculateDeltaIsotropicHardening(IsotropicHardening,rAlpha);		

	rDeltaStateFunction = 2.0 * LameMu_bar + (2.0/3.0) * (DeltaKinematicHardening + DeltaIsotropicHardening);
		
	return rDeltaStateFunction;
};

void MisesHuberYieldCriterion::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, YiedlCriterion );
}

void MisesHuberYieldCriterion::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, YieldCriterion );
}


}  // namespace Kratos.
