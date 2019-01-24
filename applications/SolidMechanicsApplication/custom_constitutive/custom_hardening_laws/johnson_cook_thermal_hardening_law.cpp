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
#include "custom_constitutive/custom_hardening_laws/johnson_cook_thermal_hardening_law.hpp"

#include "solid_mechanics_application_variables.h"

namespace Kratos
{

//*******************************CONSTRUCTOR******************************************
//************************************************************************************

JohnsonCookThermalHardeningLaw::JohnsonCookThermalHardeningLaw()
	:HardeningLaw()
{

}


//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

JohnsonCookThermalHardeningLaw& JohnsonCookThermalHardeningLaw::operator=(JohnsonCookThermalHardeningLaw const& rOther)
{
   HardeningLaw::operator=(rOther);
   return *this;
}

//*******************************COPY CONSTRUCTOR*************************************
//************************************************************************************

JohnsonCookThermalHardeningLaw::JohnsonCookThermalHardeningLaw(JohnsonCookThermalHardeningLaw const& rOther)
	:HardeningLaw(rOther)
{

}

//********************************CLONE***********************************************
//************************************************************************************

HardeningLaw::Pointer JohnsonCookThermalHardeningLaw::Clone() const
{
  return Kratos::make_shared<JohnsonCookThermalHardeningLaw>(*this);
}

//********************************DESTRUCTOR******************************************
//************************************************************************************

JohnsonCookThermalHardeningLaw::~JohnsonCookThermalHardeningLaw()
{
}

/// Operations.

//*******************************CALCULATE TOTAL HARDENING****************************
//************************************************************************************

double& JohnsonCookThermalHardeningLaw::CalculateHardening(double &rHardening, const Parameters& rValues)
{
	//get values
	const double& rRateFactor              = rValues.GetRateFactor();
	const double& rEquivalentPlasticStrain = rValues.GetEquivalentPlasticStrain();
	const double& rTemperature             = rValues.GetTemperature();
	const double& rDeltaGamma              = rValues.GetDeltaGamma();
	const double& rDeltaTime               = rValues.GetDeltaTime();

	//Constant Parameters of the -- Johnson and Cook --:
	double A = GetProperties()[JC_PARAMETER_A];
	double B = GetProperties()[JC_PARAMETER_B];
	double C = GetProperties()[JC_PARAMETER_C];

	double n = GetProperties()[JC_PARAMETER_n];
	double m = GetProperties()[JC_PARAMETER_m];

	double ReferenceTemperature = GetProperties()[REFERENCE_TEMPERATURE];
	double MeldTemperature      = GetProperties()[MELD_TEMPERATURE];
	double PlasticStrainRate    = GetProperties()[PLASTIC_STRAIN_RATE];

	double DeltaTemperature = rTemperature - ReferenceTemperature;
	if( DeltaTemperature < 0){
	  if( DeltaTemperature < -1.0 )
	    std::cout<<" Initial Temperature conditions not defined properly ("<<rTemperature<<" < "<<ReferenceTemperature<<") :"<<(rTemperature - ReferenceTemperature)<<std::endl;
	  DeltaTemperature = 0;
	}

	double NormalizedTemperature = (1.0 - pow( (DeltaTemperature/(MeldTemperature - ReferenceTemperature)), m) );

	// if( NormalizedTemperature < 0 )
	//   NormalizedTemperature = 0;

	if( rEquivalentPlasticStrain <= 0 )
	  rHardening   =  A * NormalizedTemperature;
	else
	  rHardening   = ( A + B * pow(rEquivalentPlasticStrain, n) ) * NormalizedTemperature;

	if( rRateFactor != 0 ){

	  if( rDeltaGamma == 0 )
	    std::cout<<" H Something is wrong in the Johnson_Cook_hardening variables supplied :: DeltaGamma= "<<rDeltaGamma<<" RateFactor= "<<rRateFactor<<std::endl;

	  double RateComparisson = (rDeltaGamma * sqrt(2.0/3.0))/(PlasticStrainRate * rDeltaTime);

	  if( RateComparisson <= 0 )
	    RateComparisson = 1e-40;

	  rHardening  *= ( 1 + rRateFactor * C * std::log( RateComparisson ) );

	}

	//std::cout<< " rHardening "<<rHardening<<std::endl;

	return rHardening;

}

//*******************************CALCULATE ISOTROPIC HARDENING************************
//************************************************************************************

double& JohnsonCookThermalHardeningLaw::CalculateIsotropicHardening(double &rIsotropicHardening, const Parameters& rValues)
{
	KRATOS_THROW_ERROR(std::logic_error, "Call to IsotropicHardening instead of Hardening","");

	return rIsotropicHardening;
}

//*******************************CALCULATE KINEMATIC HARDENING************************
//************************************************************************************

double& JohnsonCookThermalHardeningLaw::CalculateKinematicHardening(double &rKinematicHardening, const Parameters& rValues)
{
	KRATOS_THROW_ERROR(std::logic_error, "Call to KinematicHardening instead of Hardening","");

	return rKinematicHardening;
}



//*******************************CALCULATE HARDENING DERIVATIVE***********************
//************************************************************************************

double& JohnsonCookThermalHardeningLaw::CalculateDeltaHardening(double &rDeltaHardening, const Parameters& rValues)
{

	//get values
	const double& rRateFactor                 = rValues.GetRateFactor();
	const double& rEquivalentPlasticStrain    = rValues.GetEquivalentPlasticStrain();
	const double& rTemperature                = rValues.GetTemperature();
	const double& rDeltaGamma                 = rValues.GetDeltaGamma();
	const double& rDeltaTime                  = rValues.GetDeltaTime();

	//Constant Parameters of the -- Johnson and Cook --:
	double A = GetProperties()[JC_PARAMETER_A];
	double B = GetProperties()[JC_PARAMETER_B];
	double C = GetProperties()[JC_PARAMETER_C];

	double n = GetProperties()[JC_PARAMETER_n];
	double m = GetProperties()[JC_PARAMETER_m];

	double ReferenceTemperature = GetProperties()[REFERENCE_TEMPERATURE];
	double MeldTemperature      = GetProperties()[MELD_TEMPERATURE];
	double PlasticStrainRate    = GetProperties()[PLASTIC_STRAIN_RATE];

	double DeltaTemperature = rTemperature - ReferenceTemperature;
	if( DeltaTemperature < 0){
	  if( DeltaTemperature < -1.0 )
	    std::cout<<" Initial Temperature conditions not defined properly ("<<rTemperature<<" < "<<ReferenceTemperature<<") :"<<(rTemperature - ReferenceTemperature)<<std::endl;
	  DeltaTemperature = 0;
	}

	double NormalizedTemperature = (1.0 - pow( (DeltaTemperature/(MeldTemperature - ReferenceTemperature)), m) );

	// if( NormalizedTemperature < 0 )
	//   NormalizedTemperature = 0;

	if( rEquivalentPlasticStrain <= 0 )
	  rDeltaHardening = 0;
	else
	  rDeltaHardening = ( B * n * pow( rEquivalentPlasticStrain, n-1 ) ) * NormalizedTemperature;


	if( rRateFactor != 0 ){

	  if( rDeltaGamma == 0 )
	    std::cout<<" DH Something is wrong in the Johnson_Cook_hardening variables supplied :: DeltaGamma= "<<rDeltaGamma<<" RateFactor= "<<rRateFactor<<std::endl;

	  double RateComparisson = (rDeltaGamma * sqrt(2.0/3.0))/(PlasticStrainRate * rDeltaTime);

          if( RateComparisson <= 0 )
	    RateComparisson = 1e-40;

	  rDeltaHardening *= ( 1.0 + rRateFactor * C * std::log( RateComparisson ) );

	  if( rEquivalentPlasticStrain <= 0 )
	    rDeltaHardening += rRateFactor * ( sqrt(3.0/2.0) * ( A ) * NormalizedTemperature * C / rDeltaGamma );
	  else
	    rDeltaHardening += rRateFactor * ( sqrt(3.0/2.0) * ( A + B * pow( rEquivalentPlasticStrain, n ) ) * NormalizedTemperature * C / rDeltaGamma );


	}

	//std::cout<< " DeltaHardening "<<rDeltaHardening<<std::endl;

	return rDeltaHardening;
}

//***************************CALCULATE ISOTROPIC HARDENING DERIVATIVE*****************
//************************************************************************************

double& JohnsonCookThermalHardeningLaw::CalculateDeltaIsotropicHardening(double &rDeltaIsotropicHardening, const Parameters& rValues)
{
	KRATOS_THROW_ERROR(std::logic_error, "Call to DeltaIsotropicHardening instead of Hardening","");

	return rDeltaIsotropicHardening;

}


//***************************CALCULATE KINEMATIC HARDENING DERIVATIVE*****************
//************************************************************************************

double& JohnsonCookThermalHardeningLaw::CalculateDeltaKinematicHardening(double &rDeltaKinematicHardening, const Parameters& rValues)
{

	KRATOS_THROW_ERROR(std::logic_error, "Call to DeltaKinematicHardening instead of Hardening","");

	return rDeltaKinematicHardening;
}


//***************************CALCULATE HARDENING DERIVATIVE TEMPERATURE***************
//************************************************************************************


double& JohnsonCookThermalHardeningLaw::CalculateDeltaThermalHardening(double &rDeltaThermalHardening, const Parameters& rValues)
{
	//get values
	const double& rRateFactor                 = rValues.GetRateFactor();
	const double& rEquivalentPlasticStrain    = rValues.GetEquivalentPlasticStrain();
	const double& rTemperature                = rValues.GetTemperature();
	const double& rDeltaGamma                 = rValues.GetDeltaGamma();
	const double& rDeltaTime                  = rValues.GetDeltaTime();

	//Constant Parameters of the -- Johnson and Cook --:
	double A = GetProperties()[JC_PARAMETER_A];
	double B = GetProperties()[JC_PARAMETER_B];
	double C = GetProperties()[JC_PARAMETER_C];

	double n = GetProperties()[JC_PARAMETER_n];
	double m = GetProperties()[JC_PARAMETER_m];

	double ReferenceTemperature = GetProperties()[REFERENCE_TEMPERATURE];
	double MeldTemperature      = GetProperties()[MELD_TEMPERATURE];
	double PlasticStrainRate    = GetProperties()[PLASTIC_STRAIN_RATE];

	double DeltaTemperature = rTemperature - ReferenceTemperature;
	double DeltaNormalizedTemperature = 0;

	if( DeltaTemperature <= 0 ){
	  if( DeltaTemperature < -1.0 )
	    std::cout<<" Initial Temperature conditions not defined properly ("<<rTemperature<<" < "<<ReferenceTemperature<<") :"<<(rTemperature - ReferenceTemperature)<<std::endl;
	  DeltaTemperature = 0;
	}
	else{
	  DeltaNormalizedTemperature = ( pow( (DeltaTemperature/(MeldTemperature - ReferenceTemperature)), m-1)/(MeldTemperature - ReferenceTemperature) );
	}


	// std::cout<<" DeltaNormalizedTemperature "<<DeltaNormalizedTemperature<<" DeltaTemperature "<<DeltaTemperature<<std::endl;
	// std::cout<<" MeldTemperature "<<MeldTemperature<<" Reference Temperature "<<ReferenceTemperature<<std::endl;
	// std::cout<<" T/(Tm-Tr) "<<(DeltaTemperature/(MeldTemperature - ReferenceTemperature))<<" (T/(Tm-Tr))**m-1 " << pow( (DeltaTemperature/(MeldTemperature - ReferenceTemperature)), m-1)<<" m-1 "<<m-1<<std::endl;

	if( rEquivalentPlasticStrain < 0 )
	  rDeltaThermalHardening  =  m *( A ) * DeltaNormalizedTemperature;
	else
	  rDeltaThermalHardening  =  m *( A + B * pow ( rEquivalentPlasticStrain, n) ) * DeltaNormalizedTemperature;

	// std::cout<<" DeltaThermalHardening "<<rDeltaThermalHardening<<" EquivalentPlasticStrain "<<rEquivalentPlasticStrain<<std::endl;

	if( rRateFactor != 0 ){

	  if( rDeltaGamma == 0 )
	    std::cout<<" DTH Something is wrong in the Johnson_Cook_hardening variables supplied "<<std::endl;

	  double RateComparisson = (rDeltaGamma * sqrt(2.0/3.0))/(PlasticStrainRate * rDeltaTime);

	  // std::cout<<" RateComparison "<<RateComparisson<<" PlasticStrainRate "<<PlasticStrainRate<<std::endl;

	  if( RateComparisson <= 0 )
	    RateComparisson =  1e-40;

	  rDeltaThermalHardening *=  ( 1.0 + rRateFactor * C * std::log( RateComparisson ) );


	}

	return rDeltaThermalHardening;
}


void JohnsonCookThermalHardeningLaw::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, HardeningLaw );
}

void JohnsonCookThermalHardeningLaw::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, HardeningLaw );
}


}  // namespace Kratos.
