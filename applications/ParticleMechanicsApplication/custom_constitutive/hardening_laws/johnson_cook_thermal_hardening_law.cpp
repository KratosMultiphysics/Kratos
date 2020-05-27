//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    JMCarbonell
//					 (adapted to Particle Mechanics by Peter Wilson)
//

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/properties.h"
#include "custom_constitutive/hardening_laws/johnson_cook_thermal_hardening_law.hpp"

#include "particle_mechanics_application_variables.h"

namespace Kratos
{

//*******************************CONSTRUCTOR******************************************
//************************************************************************************

JohnsonCookThermalHardeningLaw::JohnsonCookThermalHardeningLaw()
	:ParticleHardeningLaw()
{

}


//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

JohnsonCookThermalHardeningLaw& JohnsonCookThermalHardeningLaw::operator=(JohnsonCookThermalHardeningLaw const& rOther)
{
	ParticleHardeningLaw::operator=(rOther);
   return *this;
}

//*******************************COPY CONSTRUCTOR*************************************
//************************************************************************************

JohnsonCookThermalHardeningLaw::JohnsonCookThermalHardeningLaw(JohnsonCookThermalHardeningLaw const& rOther)
	:ParticleHardeningLaw(rOther)
{

}

//********************************CLONE***********************************************
//************************************************************************************

ParticleHardeningLaw::Pointer JohnsonCookThermalHardeningLaw::Clone() const
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
	const double rPlasticStrainRate	      = rValues.GetPlasticStrainRate(); // plastic strain rate
	const double rEquivalentPlasticStrain = rValues.GetEquivalentPlasticStrain();
	const double rTemperature             = rValues.GetTemperature();

	//Constant Parameters of the -- Johnson and Cook --:
	const double A = GetProperties()[JC_PARAMETER_A];
	const double B = GetProperties()[JC_PARAMETER_B];
	const double C = GetProperties()[JC_PARAMETER_C];

	const double n = GetProperties()[JC_PARAMETER_n];
	const double m = GetProperties()[JC_PARAMETER_m];

	const double ReferenceTemperature = GetProperties()[REFERENCE_TEMPERATURE];
	const double MeldTemperature      = GetProperties()[MELD_TEMPERATURE];
	const double ReferenceStrainRate    = GetProperties()[REFERENCE_STRAIN_RATE];

	// Hardening formula is: = (A + B* ep^n) * strain_rate_hardening_factor * thermal_hardening_factor

	// Calculate thermal hardening factor
	double thermal_hardening_factor;
	if (rTemperature < ReferenceTemperature)
	{
		thermal_hardening_factor = 1.0;
	}
	else if (rTemperature >= MeldTemperature)
	{
		thermal_hardening_factor = 0.0;
	}
	else
	{
		thermal_hardening_factor = 1.0 - std::pow((rTemperature-ReferenceTemperature) / (MeldTemperature-ReferenceTemperature), m);
	}

	// Calculate strain rate hardening factor
	double strain_rate_hardening_factor = 1.0;
	if (rPlasticStrainRate > ReferenceStrainRate)
	{
		strain_rate_hardening_factor += C * std::log(rPlasticStrainRate / ReferenceStrainRate);
	}

	// Store results into hardening
	rHardening = thermal_hardening_factor * strain_rate_hardening_factor;

	return rHardening;
}

void JohnsonCookThermalHardeningLaw::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, ParticleHardeningLaw );
}

void JohnsonCookThermalHardeningLaw::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ParticleHardeningLaw );
}


}  // namespace Kratos.
