//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_constitutive/custom_hardening_laws/non_linear_isotropic_kinematic_hardening_law.hpp"

#include "solid_mechanics_application_variables.h"

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
	,mTheta(rOther.mTheta)
{

}


//********************************CLONE***********************************************
//************************************************************************************

HardeningLaw::Pointer NonLinearIsotropicKinematicHardeningLaw::Clone() const
{
  return Kratos::make_shared<NonLinearIsotropicKinematicHardeningLaw>(*this);
}


//********************************DESTRUCTOR******************************************
//************************************************************************************

NonLinearIsotropicKinematicHardeningLaw::~NonLinearIsotropicKinematicHardeningLaw()
{
}

/// Operations.

//*******************************CALCULATE TOTAL HARDENING****************************
//************************************************************************************

double& NonLinearIsotropicKinematicHardeningLaw::CalculateHardening(double &rHardening, const Parameters& rValues)
{

	double IsotropicHardening = this->CalculateIsotropicHardening( IsotropicHardening, rValues );

	double KinematicHardening = this->CalculateKinematicHardening( KinematicHardening, rValues );

	rHardening = IsotropicHardening + KinematicHardening;



	return rHardening;

}

//*******************************CALCULATE ISOTROPIC HARDENING************************
//************************************************************************************

double& NonLinearIsotropicKinematicHardeningLaw::CalculateIsotropicHardening(double &rIsotropicHardening, const Parameters& rValues)
{
	//get values
	const double& rEquivalentPlasticStrain = rValues.GetEquivalentPlasticStrain();
	const double& rTemperature             = rValues.GetTemperature();

        //linear hardening properties
	double  YieldStress                 =  GetProperties()[YIELD_STRESS];
	double  KinematicHardeningConstant  =  GetProperties()[KINEMATIC_HARDENING_MODULUS];

	//exponential saturation properties
   	double  K_reference           =  GetProperties()[REFERENCE_HARDENING_MODULUS];
	double  K_infinity            =  GetProperties()[INFINITY_HARDENING_MODULUS];
	const double& Delta           =  GetProperties()[HARDENING_EXPONENT];


	YieldStress                *= this->CalculateThermalReferenceEffect(rTemperature);
	K_reference                *= this->CalculateThermalReferenceEffect(rTemperature);

	K_infinity                 *= this->CalculateThermalCurrentEffect(rTemperature);
	KinematicHardeningConstant *= this->CalculateThermalCurrentEffect(rTemperature);


	//Linear Hardening law: (mTheta = 1)
	rIsotropicHardening  = YieldStress + mTheta * KinematicHardeningConstant * rEquivalentPlasticStrain;

	//Exponential Saturation:
	rIsotropicHardening += (K_infinity-K_reference) * (1.0 - exp( (-1.0) * Delta * rEquivalentPlasticStrain ) );

	return rIsotropicHardening;


}

//*******************************CALCULATE KINEMATIC HARDENING************************
//************************************************************************************

double& NonLinearIsotropicKinematicHardeningLaw::CalculateKinematicHardening(double &rKinematicHardening, const Parameters& rValues)
{
	//get values
	const double& rTemperature          =  rValues.GetTemperature();

      	//linear hardening properties
	double  KinematicHardeningConstant  =  GetProperties()[KINEMATIC_HARDENING_MODULUS];

	KinematicHardeningConstant *= this->CalculateThermalCurrentEffect(rTemperature);

	//Linear Hardening law:
	rKinematicHardening  = (1.0 - mTheta) * KinematicHardeningConstant;

	return rKinematicHardening;
}



//*******************************CALCULATE HARDENING DERIVATIVE***********************
//************************************************************************************

double& NonLinearIsotropicKinematicHardeningLaw::CalculateDeltaHardening(double &rDeltaHardening, const Parameters& rValues)
{

	double DeltaIsotropicHardening = this->CalculateDeltaIsotropicHardening( DeltaIsotropicHardening, rValues );

	double DeltaKinematicHardening = this->CalculateDeltaKinematicHardening( DeltaKinematicHardening, rValues );

	rDeltaHardening = DeltaIsotropicHardening + DeltaKinematicHardening;


	return rDeltaHardening;
}

//***************************CALCULATE ISOTROPIC HARDENING DERIVATIVE*****************
//************************************************************************************

double& NonLinearIsotropicKinematicHardeningLaw::CalculateDeltaIsotropicHardening(double &rDeltaIsotropicHardening, const Parameters& rValues)
{
	//get values
	const double& rEquivalentPlasticStrain = rValues.GetEquivalentPlasticStrain();
	const double& rTemperature             = rValues.GetTemperature();

       	//linear hardening properties
	double  KinematicHardeningConstant  =  GetProperties()[KINEMATIC_HARDENING_MODULUS];

	//exponential saturation properties
   	double  K_reference           =  GetProperties()[REFERENCE_HARDENING_MODULUS];
	double  K_infinity            =  GetProperties()[INFINITY_HARDENING_MODULUS];
	const double& Delta           =  GetProperties()[HARDENING_EXPONENT];


	K_reference                *= this->CalculateThermalReferenceEffect(rTemperature);
	K_infinity                 *= this->CalculateThermalCurrentEffect(rTemperature);
	KinematicHardeningConstant *= this->CalculateThermalCurrentEffect(rTemperature);


	//Linear Hardening law: (mTheta = 1)
	rDeltaIsotropicHardening  = mTheta * KinematicHardeningConstant;

	//Exponential Saturation:
	rDeltaIsotropicHardening += Delta * (K_infinity-K_reference) * ( exp( (-1.0) * Delta * rEquivalentPlasticStrain ) );

	return rDeltaIsotropicHardening;

}


//***************************CALCULATE KINEMATIC HARDENING DERIVATIVE*****************
//************************************************************************************

double& NonLinearIsotropicKinematicHardeningLaw::CalculateDeltaKinematicHardening(double &rDeltaKinematicHardening, const Parameters& rValues)
{
	rDeltaKinematicHardening = 0;

	return rDeltaKinematicHardening;
}


//***************************CALCULATE TEMPERATURE EVOLUTION PROPERTIES***************
//************************************************************************************


double NonLinearIsotropicKinematicHardeningLaw::CalculateThermalReferenceEffect(const double &rTemperature)
{
  return 1;
}

//***************************CALCULATE TEMPERATURE EVOLUTION PROPERTIES***************
//************************************************************************************

double NonLinearIsotropicKinematicHardeningLaw::CalculateThermalCurrentEffect(const double &rTemperature)
{
  return 1;
}


void NonLinearIsotropicKinematicHardeningLaw::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, HardeningLaw )
    rSerializer.save("Theta",mTheta);
}

void NonLinearIsotropicKinematicHardeningLaw::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, HardeningLaw )
    rSerializer.load("Theta",mTheta);
}


}  // namespace Kratos.
