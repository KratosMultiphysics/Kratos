//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    Bodhinanda Chandra
//

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "custom_constitutive/hardening_laws/exponential_strain_softening_law.hpp"
#include "particle_mechanics_application.h"
#include "utilities/math_utils.h"


namespace Kratos
{

//*******************************CONSTRUCTOR******************************************
//************************************************************************************

ExponentialStrainSofteningLaw::ExponentialStrainSofteningLaw()
	:HardeningLaw()
{
   
}


//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

ExponentialStrainSofteningLaw& ExponentialStrainSofteningLaw::operator=(ExponentialStrainSofteningLaw const& rOther)
{
   HardeningLaw::operator=(rOther);
   return *this;
}

//*******************************COPY CONSTRUCTOR*************************************
//************************************************************************************

ExponentialStrainSofteningLaw::ExponentialStrainSofteningLaw(ExponentialStrainSofteningLaw const& rOther)
	:HardeningLaw(rOther)
{

}


//********************************DESTRUCTOR******************************************
//************************************************************************************

ExponentialStrainSofteningLaw::~ExponentialStrainSofteningLaw()
{
}

/// Operations.

void ExponentialStrainSofteningLaw::InitializeMaterial (const Properties& rMaterialProperties)
{
    HardeningLaw::InitializeMaterial(rMaterialProperties);

    this->InitializeHardeningParameters(mHardeningParameters);

}

void ExponentialStrainSofteningLaw::InitializeHardeningParameters(HardeningParameters& rHardeningParameters){
    rHardeningParameters.DeltaCohesion        = 0.0;
    rHardeningParameters.DeltaFrictionAngle   = 0.0;
    rHardeningParameters.DeltaDilatancyAngle  = 0.0;
}

//*******************************CALCULATE TOTAL HARDENING****************************
//************************************************************************************

double& ExponentialStrainSofteningLaw::CalculateHardening(double &rHardening, const double &rAlpha, const Variable<double>& rThisVariable)
{
    if (rThisVariable==INTERNAL_FRICTION_ANGLE)
    {
        rHardening = mHardeningParameters.DeltaFrictionAngle;
    }
    else if (rThisVariable==INTERNAL_DILATANCY_ANGLE)
    {
        rHardening = mHardeningParameters.DeltaDilatancyAngle;
    }
    else if (rThisVariable==COHESION)
    {
        rHardening = mHardeningParameters.DeltaCohesion;
    }

    return rHardening;
}
  


void ExponentialStrainSofteningLaw::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, HardeningLaw )

}

void ExponentialStrainSofteningLaw::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, HardeningLaw )

}


}  // namespace Kratos.
