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
#include "custom_constitutive/custom_hardening_laws/non_linear_isotropic_kinematic_hardening_law.hpp"

#include "solid_mechanics_application.h"

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

double& NonLinearIsotropicKinematicHardeningLaw::CalculateHardening(double &Hardening, double & rAlpha)
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
	Hardening += (K_infinity-K_reference) * (1.0 - exp( (-1.0) * Delta * rAlpha ) );
	
	return Hardening;

}
  
//*******************************CALCULATE ISOTROPIC HARDENING************************
//************************************************************************************

double& NonLinearIsotropicKinematicHardeningLaw::CalculateIsotropicHardening(double &IsotropicHardening, double & rAlpha)
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
	IsotropicHardening += (K_infinity-K_reference) * (1.0 - exp( (-1.0) * Delta * rAlpha ) );
	
	return IsotropicHardening;	


}

//*******************************CALCULATE KINEMATIC HARDENING************************
//************************************************************************************

double& NonLinearIsotropicKinematicHardeningLaw::CalculateKinematicHardening(double &KinematicHardening, double & rAlpha)
{
	//linear hardening properties
	const double& KinematicHardeningConstant  =  GetProperties()[KINEMATIC_HARDENING];
	
	//Linear Hardening law:
	KinematicHardening  = (1.0 - mTheta) * KinematicHardeningConstant;
	
	return KinematicHardening;
}



//*******************************CALCULATE HARDENING DERIVATIVE***********************
//************************************************************************************

double& NonLinearIsotropicKinematicHardeningLaw::CalculateDeltaHardening(double &DeltaHardening, double & rAlpha)
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
	DeltaHardening += Delta * (K_infinity-K_reference) * ( exp( (-1.0) * Delta * rAlpha ) );
	
	return DeltaHardening;	
}

//***************************CALCULATE ISOTROPIC HARDENING DERIVATIVE*****************
//************************************************************************************

double& NonLinearIsotropicKinematicHardeningLaw::CalculateDeltaIsotropicHardening(double &DeltaIsotropicHardening, double & rAlpha)
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
	DeltaIsotropicHardening += Delta * (K_infinity-K_reference) * ( exp( (-1.0) * Delta * rAlpha ) );
	
	return DeltaIsotropicHardening;	

}


//***************************CALCULATE KINEMATIC HARDENING DERIVATIVE*****************
//************************************************************************************

double& NonLinearIsotropicKinematicHardeningLaw::CalculateDeltaKinematicHardening(double &DeltaKinematicHardening, double & rAlpha)
{
	DeltaKinematicHardening = 0;

	return DeltakinematicHardening;
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
