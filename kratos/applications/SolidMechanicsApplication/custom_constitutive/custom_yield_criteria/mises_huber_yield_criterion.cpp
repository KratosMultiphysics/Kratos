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

double& MisesHuberYieldCriterion::CalculateYieldCondition(double & rStateFunction, double& rNormStress, double& rAlpha);
{
	double Hardening = 0;

	Hardening = mpHardeningLaw->CalculateHardening(Hardening,rAlpha);
		
	rStateFunction = rNormStress - sqrt(2.0/3.0) * Hardening;
		
	return rStateFunction;
};


//***************************CALCULATE YIELD CONDITION********************************
//************************************************************************************


double& MisesHuberYieldCriterion::CalculateYieldCondition(double & rStateFunction, Matrix& rStressMatrix, double& rAlpha);
{
	double	NormStress = sqrt(rStressMatrix( 0 , 0 )*rStressMatrix( 0 , 0 )+
				  rStressMatrix( 1 , 1 )*rStressMatrix( 1 , 1 )+
				  rStressMatrix( 2 , 2 )*rStressMatrix( 2 , 2 )+
				  2.0 * rStressMatrix( 0 , 1 )*rStressMatrix( 0 , 1 ) );

	rStateFunction = this->CalculateYieldCondition( rStateFunction, NormStress, rAlpha );

	return rStateFunction;
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
