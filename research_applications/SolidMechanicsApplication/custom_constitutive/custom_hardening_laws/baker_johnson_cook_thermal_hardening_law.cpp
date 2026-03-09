//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2018 $
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
#include "custom_constitutive/custom_hardening_laws/baker_johnson_cook_thermal_hardening_law.hpp"

#include "solid_mechanics_application_variables.h"

namespace Kratos
{

//*******************************CONSTRUCTOR******************************************
//************************************************************************************

BakerJohnsonCookThermalHardeningLaw::BakerJohnsonCookThermalHardeningLaw()
	:HardeningLaw()
{

}


//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

BakerJohnsonCookThermalHardeningLaw& BakerJohnsonCookThermalHardeningLaw::operator=(BakerJohnsonCookThermalHardeningLaw const& rOther)
{
   HardeningLaw::operator=(rOther);
   return *this;
}

//*******************************COPY CONSTRUCTOR*************************************
//************************************************************************************

BakerJohnsonCookThermalHardeningLaw::BakerJohnsonCookThermalHardeningLaw(BakerJohnsonCookThermalHardeningLaw const& rOther)
	:HardeningLaw(rOther)
{

}

//********************************CLONE***********************************************
//************************************************************************************

HardeningLaw::Pointer BakerJohnsonCookThermalHardeningLaw::Clone() const
{
  return Kratos::make_shared<BakerJohnsonCookThermalHardeningLaw>(*this);
}

//********************************DESTRUCTOR******************************************
//************************************************************************************

BakerJohnsonCookThermalHardeningLaw::~BakerJohnsonCookThermalHardeningLaw()
{
}

/// Operations.

//*******************************CALCULATE TOTAL HARDENING****************************
//************************************************************************************

double& BakerJohnsonCookThermalHardeningLaw::CalculateHardening(double &rHardening, const Parameters& rValues)
{
	//get values
	const double& rRateFactor              = rValues.GetRateFactor();
	const double& rEquivalentPlasticStrain = rValues.GetEquivalentPlasticStrain();
	const double& rTemperature             = rValues.GetTemperature();
	const double& rDeltaGamma              = rValues.GetDeltaGamma();
	const double& rDeltaTime               = rValues.GetDeltaTime();

	//Constant Parameters of the -- Johnson and Cook --:
	double K = GetProperties()[JC_PARAMETER_K];
	double C = GetProperties()[JC_PARAMETER_C];

	double n = GetProperties()[JC_PARAMETER_n];
	double m = GetProperties()[JC_PARAMETER_m];

	double ReferenceTemperature = GetProperties()[REFERENCE_TEMPERATURE];
	double MeldTemperature      = GetProperties()[MELD_TEMPERATURE];
	double PlasticStrainRate    = GetProperties()[PLASTIC_STRAIN_RATE];

	if(rTemperature - ReferenceTemperature < 0){
		std::cout<<" Initial Temperature conditions not defined properly ("<<rTemperature<<" < "<<ReferenceTemperature<<")"<<std::endl;

	}

	double NormalizedTemperature =  pow( (rTemperature)/(MeldTemperature), m );

	double Kv = K * exp((-1)*NormalizedTemperature);
	double nv = n * exp((-1)*NormalizedTemperature);


	rHardening   = ( Kv * pow(rEquivalentPlasticStrain, nv) );

	if( rRateFactor != 0 ){

	  if( rDeltaGamma == 0 )
	    std::cout<<" Something is wrong in the Johnson_Cook_hardening variables supplied "<<std::endl;

	  rHardening  *= (1.0 + rRateFactor * C * std::log( (rDeltaGamma * sqrt(2.0/3.0))/(PlasticStrainRate * rDeltaTime) ) );

	}

	return rHardening;

}

//*******************************CALCULATE ISOTROPIC HARDENING************************
//************************************************************************************

double& BakerJohnsonCookThermalHardeningLaw::CalculateIsotropicHardening(double &rIsotropicHardening, const Parameters& rValues)
{
	KRATOS_THROW_ERROR(std::logic_error, "Call to IsotropicHardening instead of Hardening","");

	return rIsotropicHardening;
}

//*******************************CALCULATE KINEMATIC HARDENING************************
//************************************************************************************

double& BakerJohnsonCookThermalHardeningLaw::CalculateKinematicHardening(double &rKinematicHardening, const Parameters& rValues)
{
	KRATOS_THROW_ERROR(std::logic_error, "Call to KinematicHardening instead of Hardening","");

	return rKinematicHardening;
}



//*******************************CALCULATE HARDENING DERIVATIVE***********************
//************************************************************************************

double& BakerJohnsonCookThermalHardeningLaw::CalculateDeltaHardening(double &rDeltaHardening, const Parameters& rValues)
{

	//get values
	const double& rRateFactor                 = rValues.GetRateFactor();
	const double& rEquivalentPlasticStrain    = rValues.GetEquivalentPlasticStrain();
	const double& rTemperature                = rValues.GetTemperature();
	const double& rDeltaGamma                 = rValues.GetDeltaGamma();
	const double& rDeltaTime                  = rValues.GetDeltaTime();

	//Constant Parameters of the -- Johnson and Cook --:
	double K = GetProperties()[JC_PARAMETER_K];
	double C = GetProperties()[JC_PARAMETER_C];

	double n = GetProperties()[JC_PARAMETER_n];
	double m = GetProperties()[JC_PARAMETER_m];

	double ReferenceTemperature = GetProperties()[REFERENCE_TEMPERATURE];
	double MeldTemperature      = GetProperties()[MELD_TEMPERATURE];
	double PlasticStrainRate    = GetProperties()[PLASTIC_STRAIN_RATE];

	if(rTemperature - ReferenceTemperature < 0){
		std::cout<<" Initial Temperature conditions not defined properly ("<<rTemperature<<" < "<<ReferenceTemperature<<")"<<std::endl;

	}

	double NormalizedTemperature =  pow( (rTemperature)/(MeldTemperature), m );

	double Kv = K * exp((-1)*NormalizedTemperature);
	double nv = n * exp((-1)*NormalizedTemperature);


	rDeltaHardening  = ( nv * Kv * pow(rEquivalentPlasticStrain, nv-1) );

	if( rRateFactor != 0 ){

	  if( rDeltaGamma == 0 )
	    std::cout<<" Something is wrong in the Johnson_Cook_hardening variables supplied "<<std::endl;

	  rDeltaHardening *= (1.0 + rRateFactor * C * std::log( (rDeltaGamma * sqrt(2.0/3.0))/(PlasticStrainRate * rDeltaTime) ) );

	  rDeltaHardening += rRateFactor * ( sqrt(3.0/2.0) * ( Kv * pow( rEquivalentPlasticStrain, nv ) ) * C / rDeltaGamma );
	}

	return rDeltaHardening;
}

//***************************CALCULATE ISOTROPIC HARDENING DERIVATIVE*****************
//************************************************************************************

double& BakerJohnsonCookThermalHardeningLaw::CalculateDeltaIsotropicHardening(double &rDeltaIsotropicHardening, const Parameters& rValues)
{
	KRATOS_THROW_ERROR(std::logic_error, "Call to DeltaIsotropicHardening instead of Hardening","");

	return rDeltaIsotropicHardening;

}


//***************************CALCULATE KINEMATIC HARDENING DERIVATIVE*****************
//************************************************************************************

double& BakerJohnsonCookThermalHardeningLaw::CalculateDeltaKinematicHardening(double &rDeltaKinematicHardening, const Parameters& rValues)
{

	KRATOS_THROW_ERROR(std::logic_error, "Call to DeltaKinematicHardening instead of Hardening","");

	return rDeltaKinematicHardening;
}


//***************************CALCULATE HARDENING DERIVATIVE TEMPERATURE***************
//************************************************************************************


double& BakerJohnsonCookThermalHardeningLaw::CalculateDeltaThermalHardening(double &rDeltaThermalHardening, const Parameters& rValues)
{
	//get values
	const double& rRateFactor                 = rValues.GetRateFactor();
	const double& rEquivalentPlasticStrain    = rValues.GetEquivalentPlasticStrain();
	const double& rTemperature                = rValues.GetTemperature();
	const double& rDeltaGamma                 = rValues.GetDeltaGamma();
	const double& rDeltaTime                  = rValues.GetDeltaTime();

	//Constant Parameters of the -- Johnson and Cook --:
	double K = GetProperties()[JC_PARAMETER_K];
	double C = GetProperties()[JC_PARAMETER_C];

	double n = GetProperties()[JC_PARAMETER_n];
	double m = GetProperties()[JC_PARAMETER_m];

	double ReferenceTemperature = GetProperties()[REFERENCE_TEMPERATURE];
	double MeldTemperature      = GetProperties()[MELD_TEMPERATURE];
	double PlasticStrainRate    = GetProperties()[PLASTIC_STRAIN_RATE];

	if(rTemperature - ReferenceTemperature < 0){
		std::cout<<" Initial Temperature conditions not defined properly ("<<rTemperature<<" < "<<ReferenceTemperature<<")"<<std::endl;

	}

	double NormalizedTemperature =  pow( (rTemperature)/(MeldTemperature), m );

	double Kv = K * exp((-1)*NormalizedTemperature);
	double nv = n * exp((-1)*NormalizedTemperature);


	double DeltaNormalizedTemperature = (m / MeldTemperature) * pow( (rTemperature/MeldTemperature), m-1 );

        rDeltaThermalHardening  = (1 + nv *std::log( rEquivalentPlasticStrain ) );

	rDeltaThermalHardening *= ( Kv * pow ( rEquivalentPlasticStrain, nv) * DeltaNormalizedTemperature );

	if( rRateFactor != 0 ){

	  if( rDeltaGamma == 0 )
	    std::cout<<" Something is wrong in the Johnson_Cook_hardening variables supplied "<<std::endl;
	  rDeltaThermalHardening *= ( 1.0 + rRateFactor * C * std::log( (rDeltaGamma * sqrt(2.0/3.0))/(PlasticStrainRate * rDeltaTime) ) );
	}

	return rDeltaThermalHardening;
}


void BakerJohnsonCookThermalHardeningLaw::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, HardeningLaw );
}

void BakerJohnsonCookThermalHardeningLaw::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, HardeningLaw );
}


}  // namespace Kratos.
