//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
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
#include "solid_mechanics_application_variables.h"
#include "custom_constitutive/custom_hardening_laws/linear_isotropic_kinematic_hardening_law.hpp"


namespace Kratos
{

//*******************************CONSTRUCTOR******************************************
//************************************************************************************

LinearIsotropicKinematicHardeningLaw::LinearIsotropicKinematicHardeningLaw()
	:NonLinearIsotropicKinematicHardeningLaw()
{
   //Combined isotropic-kinematic 0<mTheta<1
   //Pure isotropic hardening mTheta=1;
   //Pure kinematic hardening mTheta=0;

   //Hardening law:
   this->mTheta = 0;

}


//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

LinearIsotropicKinematicHardeningLaw& LinearIsotropicKinematicHardeningLaw::operator=(LinearIsotropicKinematicHardeningLaw const& rOther)
{
   NonLinearIsotropicKinematicHardeningLaw::operator=(rOther);
   return *this;
}

//*******************************COPY CONSTRUCTOR*************************************
//************************************************************************************

LinearIsotropicKinematicHardeningLaw::LinearIsotropicKinematicHardeningLaw(LinearIsotropicKinematicHardeningLaw const& rOther)
	:NonLinearIsotropicKinematicHardeningLaw(rOther)
{

}


//********************************CLONE***********************************************
//************************************************************************************

HardeningLaw::Pointer LinearIsotropicKinematicHardeningLaw::Clone() const
{
  return Kratos::make_shared<LinearIsotropicKinematicHardeningLaw>(*this);
}


//********************************DESTRUCTOR******************************************
//************************************************************************************

LinearIsotropicKinematicHardeningLaw::~LinearIsotropicKinematicHardeningLaw()
{
}

/// Operations.


//*******************************CALCULATE ISOTROPIC HARDENING************************
//************************************************************************************

double& LinearIsotropicKinematicHardeningLaw::CalculateIsotropicHardening(double &rIsotropicHardening, const Parameters& rValues)
{

	//get values
	const double& rEquivalentPlasticStrain = rValues.GetEquivalentPlasticStrain();

	//linear hardening properties
	const double& YieldStress                 =  GetProperties()[YIELD_STRESS];
	const double& KinematicHardeningConstant  =  GetProperties()[KINEMATIC_HARDENING_MODULUS];


	//Linear Hardening law: (mTheta = 0)
	rIsotropicHardening  = YieldStress + (1.0 - mTheta) * KinematicHardeningConstant * rEquivalentPlasticStrain;


	return rIsotropicHardening;
}


//*******************************CALCULATE HARDENING DERIVATIVE***********************
//************************************************************************************

double& LinearIsotropicKinematicHardeningLaw::CalculateDeltaHardening(double &rDeltaHardening, const Parameters& rValues)
{
      	//linear hardening properties
	const double& KinematicHardeningConstant  =  GetProperties()[KINEMATIC_HARDENING_MODULUS];

	//Linear Hardening law: (mTheta = 0)
	rDeltaHardening  = (1.0 - mTheta) * KinematicHardeningConstant;

	return rDeltaHardening;
}

//***************************CALCULATE ISOTROPIC HARDENING DERIVATIVE*****************
//************************************************************************************

double& LinearIsotropicKinematicHardeningLaw::CalculateDeltaIsotropicHardening(double &rDeltaIsotropicHardening, const Parameters& rValues)
{
       	//linear hardening properties
	const double& KinematicHardeningConstant  =  GetProperties()[KINEMATIC_HARDENING_MODULUS];

	//Linear Hardening law: (mTheta = 0)
	rDeltaIsotropicHardening  = mTheta * KinematicHardeningConstant;

	return rDeltaIsotropicHardening;
}


void LinearIsotropicKinematicHardeningLaw::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, NonLinearIsotropicKinematicHardeningLaw )

}

void LinearIsotropicKinematicHardeningLaw::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, NonLinearIsotropicKinematicHardeningLaw )

}


}  // namespace Kratos.
