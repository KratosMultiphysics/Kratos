//
//   Project Name:        KratosPoromechanicsApplication $
//   Created by:          $Author:              IPouplana $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2015 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_constitutive/continuum_laws/custom_hardening_laws/exponential_damage_hardening_law.hpp"

#include "poromechanics_application_variables.h"

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
  return Kratos::make_shared<ExponentialDamageHardeningLaw>(*this);
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

    const double& CharacteristicSize = rValues.GetCharacteristicSize();
    const double& StateVariable = rValues.GetDeltaGamma();

    double A = 1.0/(FractureEnergy/(CharacteristicSize*DamageThreshold*DamageThreshold)-0.5);

    if(A < 0.0) A = 0.0;

    //Compute Damage variable from the internal historical variable
    rHardening = 1.0-DamageThreshold/StateVariable*exp(A*(1.0-StateVariable/DamageThreshold));

    if(rHardening < 0.0)
    {
        rHardening = 0.0;
    }
    else if(rHardening > 1.0)
    {
        rHardening = 1.0;
    }

	return rHardening;
}


//***************************** CALCULATE DAMAGE DERIVATIVE **************************
//************************************************************************************

double& ExponentialDamageHardeningLaw::CalculateDeltaHardening(double &rDeltaHardening, const Parameters& rValues)
{
	const Properties& MaterialProperties = this->GetProperties();
    const double& FractureEnergy  = MaterialProperties[FRACTURE_ENERGY];
    const double& DamageThreshold = MaterialProperties[DAMAGE_THRESHOLD];

    const double& CharacteristicSize = rValues.GetCharacteristicSize();
    const double& StateVariable = rValues.GetDeltaGamma();

    double A = 1.0/(FractureEnergy/(CharacteristicSize*DamageThreshold*DamageThreshold)-0.5);

    if(A < 0.0) A = 0.0;

    //Damage derivative with respect to the internal historical variable
    rDeltaHardening = (DamageThreshold + A*StateVariable)/(StateVariable*StateVariable)*exp(A*(1.0-StateVariable/DamageThreshold));

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
