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
#include "custom_constitutive/continuum_laws/custom_hardening_laws/modified_exponential_damage_hardening_law.hpp"

#include "poromechanics_application_variables.h"

namespace Kratos
{

//*******************************CONSTRUCTOR******************************************
//************************************************************************************

ModifiedExponentialDamageHardeningLaw::ModifiedExponentialDamageHardeningLaw()
	:HardeningLaw()
{

}


//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

ModifiedExponentialDamageHardeningLaw& ModifiedExponentialDamageHardeningLaw::operator=(ModifiedExponentialDamageHardeningLaw const& rOther)
{
   HardeningLaw::operator=(rOther);
   return *this;
}

//*******************************COPY CONSTRUCTOR*************************************
//************************************************************************************

ModifiedExponentialDamageHardeningLaw::ModifiedExponentialDamageHardeningLaw(ModifiedExponentialDamageHardeningLaw const& rOther)
	:HardeningLaw(rOther)
{

}


//********************************CLONE***********************************************
//************************************************************************************

HardeningLaw::Pointer ModifiedExponentialDamageHardeningLaw::Clone() const
{
  return Kratos::make_shared<ModifiedExponentialDamageHardeningLaw>(*this);
}


//********************************DESTRUCTOR******************************************
//************************************************************************************

ModifiedExponentialDamageHardeningLaw::~ModifiedExponentialDamageHardeningLaw()
{
}

/// Operations.

//****************************** CALCULATE DAMAGE PARAMETER **************************
//************************************************************************************

double& ModifiedExponentialDamageHardeningLaw::CalculateHardening(double &rHardening, const Parameters& rValues)
{
	const Properties& MaterialProperties = this->GetProperties();
    const double& DamageThreshold = MaterialProperties[DAMAGE_THRESHOLD];
    const double& ResidualStrength = MaterialProperties[RESIDUAL_STRENGTH];
    const double& SofteningSlope = MaterialProperties[SOFTENING_SLOPE];

    const double& StateVariable = rValues.GetDeltaGamma();

    //Compute Damage variable from the internal historical variable
    rHardening = 1.0-DamageThreshold*(1.0-ResidualStrength)/StateVariable-ResidualStrength*exp(-SofteningSlope*(StateVariable-DamageThreshold));

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

double& ModifiedExponentialDamageHardeningLaw::CalculateDeltaHardening(double &rDeltaHardening, const Parameters& rValues)
{
	const Properties& MaterialProperties = this->GetProperties();
    const double& DamageThreshold = MaterialProperties[DAMAGE_THRESHOLD];
    const double& ResidualStrength = MaterialProperties[RESIDUAL_STRENGTH];
    const double& SofteningSlope = MaterialProperties[SOFTENING_SLOPE];

    const double& StateVariable = rValues.GetDeltaGamma();

    //Damage derivative with respect to the internal historical variable
    rDeltaHardening = DamageThreshold*(1.0-ResidualStrength)/(StateVariable*StateVariable)+ResidualStrength*SofteningSlope*
                        exp(-SofteningSlope*(StateVariable-DamageThreshold));

    if(rDeltaHardening < 0.0) rDeltaHardening = 0.0;

	return rDeltaHardening;
}


//************************************************************************************
//************************************************************************************


void ModifiedExponentialDamageHardeningLaw::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, HardeningLaw )
}

void ModifiedExponentialDamageHardeningLaw::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, HardeningLaw )
}


}  // namespace Kratos.
