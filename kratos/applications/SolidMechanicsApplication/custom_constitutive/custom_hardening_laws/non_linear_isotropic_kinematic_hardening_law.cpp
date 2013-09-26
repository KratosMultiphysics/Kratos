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
#include "custom_constitutive/custom_hardening_laws/non_linear_isotropic_kinematic_hardening_law.hpp"

namespace Kratos
{

//*******************************CONSTRUCTOR******************************************
//************************************************************************************

NonLinearIsotropicKinematicHardeningLaw::NonLinearIsotropicKinematicHardeningLaw()
	:HardeningLaw()
{
   //Combined isotropic-kinematic 0<mTheta<1
   //Pure isotropic hardening mTheta=1;  
   //Pure kinematic hardening mTheta=0; 

   //Hardening law:
   mTheta = 1; 
   
}


//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

NonLinearIsotropicKinematicHardeningLaw& NonLinearIsotropicKinematicHardeningLaw::operator=(NonLinearIsotropicKinematicHardeningLaw const& rOther)
{
   HardeningLaw::operator=(rOther);
   return *this;
}

//*******************************COPY CONSTRUCTOR*************************************
//************************************************************************************

NonLinearIsotropicKinematicHardeningLaw::NonLinearIsotropicKinematicHardeningLaw(NonLinearIsotropicKinematicHardeningLaw const& rOther)
	:HardeningLaw(rOther)
{

}


//********************************DESTRUCTOR******************************************
//************************************************************************************

NonLinearIsotropicKinematicHardeningLaw::~NonLinearIsotropicKinematicHardeningLaw()
{
}

/// Operations.

//*******************************CALCULATE TOTAL HARDENING****************************
//************************************************************************************

double& NonLinearIsotropicKinematicHardeningLaw::CalculateHardening(double &rHardening,const double & rAlpha)
{
	
	//linear hardening properties
	const double& YieldStress                 =  GetProperties()[YIELD_STRESS];
	const double& KinematicHardeningConstant  =  GetProperties()[KINEMATIC_HARDENING_MODULUS];
	
	//exponential saturation properties
   	const double& K_reference           =  GetProperties()[REFERENCE_HARDENING_MODULUS];
	const double& K_infinity            =  GetProperties()[INFINITY_HARDENING_MODULUS];
	const double& Delta                 =  GetProperties()[HARDENING_EXPONENT];

	//Linear Hardening law:
	rHardening  = YieldStress + mTheta *  KinematicHardeningConstant;
	
	//Exponential Saturation:
	rHardening += (K_infinity-K_reference) * (1.0 - exp( (-1.0) * Delta * rAlpha ) );
	
	return rHardening;

}
  
//*******************************CALCULATE ISOTROPIC HARDENING************************
//************************************************************************************

double& NonLinearIsotropicKinematicHardeningLaw::CalculateIsotropicHardening(double &rIsotropicHardening,const double & rAlpha)
{

	//linear hardening properties
	const double& YieldStress                 =  GetProperties()[YIELD_STRESS];
	const double& KinematicHardeningConstant  =  GetProperties()[KINEMATIC_HARDENING_MODULUS];
	
	//exponential saturation properties
   	const double& K_reference           =  GetProperties()[REFERENCE_HARDENING_MODULUS];
	const double& K_infinity            =  GetProperties()[INFINITY_HARDENING_MODULUS];
	const double& Delta                 =  GetProperties()[HARDENING_EXPONENT];

	//Linear Hardening law: (mTheta = 1)
	rIsotropicHardening  = YieldStress + KinematicHardeningConstant;
	
	//Exponential Saturation:
	rIsotropicHardening += (K_infinity-K_reference) * (1.0 - exp( (-1.0) * Delta * rAlpha ) );
	
	return rIsotropicHardening;	


}

//*******************************CALCULATE KINEMATIC HARDENING************************
//************************************************************************************

double& NonLinearIsotropicKinematicHardeningLaw::CalculateKinematicHardening(double &rKinematicHardening,const double & rAlpha)
{
	//linear hardening properties
	const double& KinematicHardeningConstant  =  GetProperties()[KINEMATIC_HARDENING_MODULUS];
	
	//Linear Hardening law:
	rKinematicHardening  = (1.0 - mTheta) * KinematicHardeningConstant;
	
	return rKinematicHardening;
}



//*******************************CALCULATE HARDENING DERIVATIVE***********************
//************************************************************************************

double& NonLinearIsotropicKinematicHardeningLaw::CalculateDeltaHardening(double &rDeltaHardening,const double & rAlpha)
{
      	//linear hardening properties
	const double& KinematicHardeningConstant  =  GetProperties()[KINEMATIC_HARDENING_MODULUS];
	
	//exponential saturation properties
   	const double& K_reference           =  GetProperties()[REFERENCE_HARDENING_MODULUS];
	const double& K_infinity            =  GetProperties()[INFINITY_HARDENING_MODULUS];
	const double& Delta                 =  GetProperties()[HARDENING_EXPONENT];

	//Linear Hardening law: (mTheta = 1)
	rDeltaHardening  = mTheta * KinematicHardeningConstant;
	
	//Exponential Saturation:
	rDeltaHardening += Delta * (K_infinity-K_reference) * ( exp( (-1.0) * Delta * rAlpha ) );
	
	return rDeltaHardening;	
}

//***************************CALCULATE ISOTROPIC HARDENING DERIVATIVE*****************
//************************************************************************************

double& NonLinearIsotropicKinematicHardeningLaw::CalculateDeltaIsotropicHardening(double &rDeltaIsotropicHardening,const double & rAlpha)
{
       	//linear hardening properties
	const double& KinematicHardeningConstant  =  GetProperties()[KINEMATIC_HARDENING_MODULUS];
	
	//exponential saturation properties
   	const double& K_reference           =  GetProperties()[REFERENCE_HARDENING_MODULUS];
	const double& K_infinity            =  GetProperties()[INFINITY_HARDENING_MODULUS];
	const double& Delta                 =  GetProperties()[HARDENING_EXPONENT];

	//Linear Hardening law: (mTheta = 1)
	rDeltaIsotropicHardening  = KinematicHardeningConstant;
	
	//Exponential Saturation:
	rDeltaIsotropicHardening += Delta * (K_infinity-K_reference) * ( exp( (-1.0) * Delta * rAlpha ) );
	
	return rDeltaIsotropicHardening;	

}


//***************************CALCULATE KINEMATIC HARDENING DERIVATIVE*****************
//************************************************************************************

double& NonLinearIsotropicKinematicHardeningLaw::CalculateDeltaKinematicHardening(double &rDeltaKinematicHardening,const double & rAlpha)
{
	rDeltaKinematicHardening = 0;

	return rDeltaKinematicHardening;
}


void NonLinearIsotropicKinematicHardeningLaw::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, HardeningLaw );
    rSerializer.save("Theta",mTheta);
}

void NonLinearIsotropicKinematicHardeningLaw::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, HardeningLaw );
    rSerializer.load("Theta",mTheta);
}


}  // namespace Kratos.
