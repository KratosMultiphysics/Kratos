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

//*******************************CALCULATE TOTAL HARDENING****************************
//************************************************************************************

double& ExponentialStrainSofteningLaw::CalculateHardening(double &rHardening, const double &rAlpha, const Variable<double>& rThisVariable)
{
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
