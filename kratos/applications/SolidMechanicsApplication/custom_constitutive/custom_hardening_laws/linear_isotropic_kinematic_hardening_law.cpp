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
#include "custom_constitutive/custom_hardening_laws/linear_isotropic_kinematic_hardening_law.hpp"

#include "solid_mechanics_application.h"

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

double& LinearIsotropicKinematicHardeningLaw::CalculateHardening(double &Hardening,const double & rAlpha)
{
	
	//linear hardening properties
	const double& YieldStress                 =  GetProperties()[YIELD_STRESS];
	const double& KinematicHardeningConstant  =  GetProperties()[KINEMATIC_HARDENING];
	

	//Linear Hardening law:
	Hardening  = YieldStress + mTheta *  KinematicHardeningConstant;
	
	return Hardening;

}
  
//*******************************CALCULATE ISOTROPIC HARDENING************************
//************************************************************************************

double& LinearIsotropicKinematicHardeningLaw::CalculateIsotropicHardening(double &IsotropicHardening,const double & rAlpha)
{

	//linear hardening properties
	const double& YieldStress                 =  GetProperties()[YIELD_STRESS];
	const double& KinematicHardeningConstant  =  GetProperties()[KINEMATIC_HARDENING];
	

	//Linear Hardening law: (mTheta = 1)
	IsotropicHardening  = YieldStress + KinematicHardeningConstant;
	
	
	return IsotropicHardening;	
}


//*******************************CALCULATE HARDENING DERIVATIVE***********************
//************************************************************************************

double& LinearIsotropicKinematicHardeningLaw::CalculateDeltaHardening(double &DeltaHardening,const double & rAlpha)
{
      	//linear hardening properties
	const double& KinematicHardeningConstant  =  GetProperties()[KINEMATIC_HARDENING];
	
	//Linear Hardening law: (mTheta = 1)
	DeltaHardening  = mTheta * KinematicHardeningConstant;
		
	return DeltaHardening;	
}

//***************************CALCULATE ISOTROPIC HARDENING DERIVATIVE*****************
//************************************************************************************

double& LinearIsotropicKinematicHardeningLaw::CalculateDeltaIsotropicHardening(double &DeltaIsotropicHardening,const double & rAlpha)
{
       	//linear hardening properties
	const double& KinematicHardeningConstant  =  GetProperties()[KINEMATIC_HARDENING];
	
	//Linear Hardening law: (mTheta = 1)
	DeltaIsotropicHardening  = KinematicHardeningConstant;
	
	return DeltaIsotropicHardening;	
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
