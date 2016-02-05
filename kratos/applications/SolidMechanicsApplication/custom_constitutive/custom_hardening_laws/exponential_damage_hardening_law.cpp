//
//   Project Name:        KratosSolidMechanicsApplication $
//   Last modified by:    $Author:              IPouplana $
//   Date:                $Date:                July 2015 $
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
#include "custom_constitutive/custom_hardening_laws/exponential_damage_hardening_law.hpp"

namespace Kratos
{

//*******************************CONSTRUCTOR******************************************
//************************************************************************************

ExponentialDamageHardeningLaw::ExponentialDamageHardeningLaw()
	:HardeningLaw()
{
       
}


//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

ExponentialDamageHardeningLaw& ExponentialDamageHardeningLaw::operator=(ExponentialDamageHardeningLaw const& rOther)
{
   HardeningLaw::operator=(rOther);
   return *this;
}

//*******************************COPY CONSTRUCTOR*************************************
//************************************************************************************

ExponentialDamageHardeningLaw::ExponentialDamageHardeningLaw(ExponentialDamageHardeningLaw const& rOther)
	:HardeningLaw(rOther)
{

}


//********************************CLONE***********************************************
//************************************************************************************

HardeningLaw::Pointer ExponentialDamageHardeningLaw::Clone() const
{
  HardeningLaw::Pointer p_clone(new ExponentialDamageHardeningLaw(*this));
  return p_clone;
}


//********************************DESTRUCTOR******************************************
//************************************************************************************

ExponentialDamageHardeningLaw::~ExponentialDamageHardeningLaw()
{
}

/// Operations.

//****************************** CALCULATE DAMAGE PARAMETER **************************
//************************************************************************************

double& ExponentialDamageHardeningLaw::CalculateHardening(double &rHardening, const Parameters& rValues)
{
	const Properties& MaterialProperties = this->GetProperties();
    const double& FractureEnergy  = MaterialProperties[FRACTURE_ENERGY];
    const double& DamageThreshold = MaterialProperties[DAMAGE_THRESHOLD];

    const double& CharacteristicLength = rValues.GetCharacteristicSize();
    const double& InternalVariable = rValues.GetDeltaGamma();

    double A = 1.0/(FractureEnergy/(CharacteristicLength*DamageThreshold*DamageThreshold)-0.5);

    if(A < 0.0) A = 0.0;
    
    //Compute Damage variable from the internal historical variable
    rHardening = 1.0-DamageThreshold/InternalVariable*exp(A*(1.0-InternalVariable/DamageThreshold));

    if(rHardening < 0.0) rHardening = 0.0;
    if(rHardening > 1.0) rHardening = 1.0;
    
	return rHardening;
}


//***************************** CALCULATE DAMAGE DERIVATIVE **************************
//************************************************************************************

double& ExponentialDamageHardeningLaw::CalculateDeltaHardening(double &rDeltaHardening, const Parameters& rValues)
{
	const Properties& MaterialProperties = this->GetProperties();
    const double& FractureEnergy  = MaterialProperties[FRACTURE_ENERGY];
    const double& DamageThreshold = MaterialProperties[DAMAGE_THRESHOLD];

    const double& CharacteristicLength = rValues.GetCharacteristicSize();
    const double& InternalVariable = rValues.GetDeltaGamma();

    double A = 1.0/(FractureEnergy/(CharacteristicLength*DamageThreshold*DamageThreshold)-0.5);

    if(A < 0.0) A = 0.0;
    
    //Damage derivative with respect to the internal historical variable
    rDeltaHardening = (DamageThreshold + A*InternalVariable)/(InternalVariable*InternalVariable)*exp(A*(1.0 - InternalVariable/DamageThreshold));

    if(rDeltaHardening < 0.0) rDeltaHardening = 0.0;
    
	return rDeltaHardening;
}


//************************************************************************************
//************************************************************************************


void ExponentialDamageHardeningLaw::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, HardeningLaw )
}

void ExponentialDamageHardeningLaw::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, HardeningLaw )
}


}  // namespace Kratos.
