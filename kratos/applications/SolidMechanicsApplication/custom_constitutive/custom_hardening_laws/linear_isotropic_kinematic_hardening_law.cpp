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
   //mTheta = 1; 
   
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


//********************************DESTRUCTOR******************************************
//************************************************************************************

LinearIsotropicKinematicHardeningLaw::~LinearIsotropicKinematicHardeningLaw()
{
}

/// Operations.

//*******************************CALCULATE TOTAL HARDENING****************************
//************************************************************************************

double& LinearIsotropicKinematicHardeningLaw::CalculateHardening(double &rHardening, const double &rAlpha, double rTemperature)
{
	
	//linear hardening properties
	const double& YieldStress                 =  GetProperties()[YIELD_STRESS];
	const double& KinematicHardeningConstant  =  GetProperties()[KINEMATIC_HARDENING_MODULUS];
	

	//Linear Hardening law:
	rHardening  = YieldStress + mTheta *  KinematicHardeningConstant * rAlpha;
	
	return rHardening;

}
  
//*******************************CALCULATE ISOTROPIC HARDENING************************
//************************************************************************************

double& LinearIsotropicKinematicHardeningLaw::CalculateIsotropicHardening(double &rIsotropicHardening, const double &rAlpha, double rTemperature)
{

	//linear hardening properties
	const double& YieldStress                 =  GetProperties()[YIELD_STRESS];
	const double& KinematicHardeningConstant  =  GetProperties()[KINEMATIC_HARDENING_MODULUS];
	

	//Linear Hardening law: (mTheta = 1)
	rIsotropicHardening  = YieldStress + KinematicHardeningConstant * rAlpha;
	
	
	return rIsotropicHardening;	
}


//*******************************CALCULATE HARDENING DERIVATIVE***********************
//************************************************************************************

double& LinearIsotropicKinematicHardeningLaw::CalculateDeltaHardening(double &rDeltaHardening, const double &rAlpha, double rTemperature)
{
      	//linear hardening properties
	const double& KinematicHardeningConstant  =  GetProperties()[KINEMATIC_HARDENING_MODULUS];
	
	//Linear Hardening law: (mTheta = 1)
	rDeltaHardening  = mTheta * KinematicHardeningConstant;
		
	return rDeltaHardening;	
}

//***************************CALCULATE ISOTROPIC HARDENING DERIVATIVE*****************
//************************************************************************************

double& LinearIsotropicKinematicHardeningLaw::CalculateDeltaIsotropicHardening(double &rDeltaIsotropicHardening, const double &rAlpha, double rTemperature)
{
       	//linear hardening properties
	const double& KinematicHardeningConstant  =  GetProperties()[KINEMATIC_HARDENING_MODULUS];
	
	//Linear Hardening law: (mTheta = 1)
	rDeltaIsotropicHardening  = mTheta * KinematicHardeningConstant;
	
	return rDeltaIsotropicHardening;	
}


void LinearIsotropicKinematicHardeningLaw::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, NonLinearIsotropicKinematicHardeningLaw );

}

void LinearIsotropicKinematicHardeningLaw::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, NonLinearIsotropicKinematicHardeningLaw );

}


}  // namespace Kratos.
