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
#include "custom_constitutive/custom_hardening_laws/non_linear_isotropic_kinematic_thermal_hardening_law.hpp"

#include "solid_mechanics_application_variables.h"

namespace Kratos
{

//*******************************CONSTRUCTOR******************************************
//************************************************************************************

NonLinearIsotropicKinematicThermalHardeningLaw::NonLinearIsotropicKinematicThermalHardeningLaw()
	:NonLinearIsotropicKinematicHardeningLaw()
{

}


//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

NonLinearIsotropicKinematicThermalHardeningLaw& NonLinearIsotropicKinematicThermalHardeningLaw::operator=(NonLinearIsotropicKinematicThermalHardeningLaw const& rOther)
{
   NonLinearIsotropicKinematicHardeningLaw::operator=(rOther);
   return *this;
}

//*******************************COPY CONSTRUCTOR*************************************
//************************************************************************************

NonLinearIsotropicKinematicThermalHardeningLaw::NonLinearIsotropicKinematicThermalHardeningLaw(NonLinearIsotropicKinematicThermalHardeningLaw const& rOther)
	:NonLinearIsotropicKinematicHardeningLaw(rOther)
{

}

//********************************CLONE***********************************************
//************************************************************************************

HardeningLaw::Pointer NonLinearIsotropicKinematicThermalHardeningLaw::Clone() const
{
  return Kratos::make_shared<NonLinearIsotropicKinematicThermalHardeningLaw>(*this);
}

//********************************DESTRUCTOR******************************************
//************************************************************************************

NonLinearIsotropicKinematicThermalHardeningLaw::~NonLinearIsotropicKinematicThermalHardeningLaw()
{
}

/// Operations.


//***************************CALCULATE TEMPERATURE EVOLUTION PROPERTIES***************
//************************************************************************************


double NonLinearIsotropicKinematicThermalHardeningLaw::CalculateThermalReferenceEffect(const double &rTemperature)
{

	//parameters for the thermal solution
	const double& ReferenceConductivity = GetProperties()[REFERENCE_CONDUCTIVITY];
	const double& ReferenceTemperature  = GetProperties()[REFERENCE_TEMPERATURE];

	//thermal effect in the initial parameters
	double reference_temp_effect = ( 1.0 - ReferenceConductivity * (rTemperature - ReferenceTemperature) );

	return reference_temp_effect;
}

//***************************CALCULATE TEMPERATURE EVOLUTION PROPERTIES***************
//************************************************************************************

double NonLinearIsotropicKinematicThermalHardeningLaw::CalculateThermalCurrentEffect(const double &rTemperature)
{

	//parameters for the thermal solution
	const double& HardnessConductivity  = GetProperties()[HARDNESS_CONDUCTIVITY];
	const double& ReferenceTemperature  = GetProperties()[REFERENCE_TEMPERATURE];

	//thermal effect in the final parameters
	double current_temp_effect   = ( 1.0 - HardnessConductivity * (rTemperature - ReferenceTemperature) );

	return current_temp_effect;
}

//***************************CALCULATE HARDENING DERIVATIVE TEMPERATURE***************
//************************************************************************************


double& NonLinearIsotropicKinematicThermalHardeningLaw::CalculateDeltaThermalHardening(double &rDeltaThermalHardening, const Parameters& rValues)
{
	//get values
	const double& rEquivalentPlasticStrain = rValues.GetEquivalentPlasticStrain();

        //linear hardening properties
	double  YieldStress                 =  GetProperties()[YIELD_STRESS];
	double  KinematicHardeningConstant  =  GetProperties()[KINEMATIC_HARDENING_MODULUS];

	//exponential saturation properties
   	double  K_reference           =  GetProperties()[REFERENCE_HARDENING_MODULUS];
	double  K_infinity            =  GetProperties()[INFINITY_HARDENING_MODULUS];
	const double& Delta           =  GetProperties()[HARDENING_EXPONENT];


	//parameters for the thermal solution
	const double& ReferenceConductivity = GetProperties()[REFERENCE_CONDUCTIVITY];
	const double& HardnessConductivity  = GetProperties()[HARDNESS_CONDUCTIVITY];


	//Linear Hardening law: (mTheta = 1)
	rDeltaThermalHardening  = ( YieldStress * ReferenceConductivity + this->mTheta * KinematicHardeningConstant * HardnessConductivity * rEquivalentPlasticStrain );

	//Exponential Saturation:
	rDeltaThermalHardening += ( K_infinity * HardnessConductivity -K_reference * ReferenceConductivity ) * (1.0 - exp( (-1.0) * Delta * rEquivalentPlasticStrain ) );

	return rDeltaThermalHardening;
}


void NonLinearIsotropicKinematicThermalHardeningLaw::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, NonLinearIsotropicKinematicHardeningLaw );
}

void NonLinearIsotropicKinematicThermalHardeningLaw::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, NonLinearIsotropicKinematicHardeningLaw );
}


}  // namespace Kratos.
