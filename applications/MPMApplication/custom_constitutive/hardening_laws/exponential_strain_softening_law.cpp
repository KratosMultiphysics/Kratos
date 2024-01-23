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
#include "mpm_application_variables.h"
#include "utilities/math_utils.h"


namespace Kratos
{

//*******************************CONSTRUCTOR******************************************
//************************************************************************************

ExponentialStrainSofteningLaw::ExponentialStrainSofteningLaw()
	:MPMHardeningLaw()
{

}


//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

ExponentialStrainSofteningLaw& ExponentialStrainSofteningLaw::operator=(ExponentialStrainSofteningLaw const& rOther)
{
   MPMHardeningLaw::operator=(rOther);
   return *this;
}

//*******************************COPY CONSTRUCTOR*************************************
//************************************************************************************

ExponentialStrainSofteningLaw::ExponentialStrainSofteningLaw(ExponentialStrainSofteningLaw const& rOther)
	:MPMHardeningLaw(rOther)
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

double& ExponentialStrainSofteningLaw::CalculateHardening(double &rHardening, const double &rAlpha, const Variable<double>& rThisVariable, const Properties& rProperties)
{
    rHardening = 0.0;

    const double beta = rProperties[SHAPE_FUNCTION_BETA];

    if (rThisVariable == COHESION)
    {
        const double peak_cohesion = rProperties[COHESION];
        const double residual_cohesion = rProperties[COHESION_RESIDUAL];
        rHardening = -1.0 * beta * (peak_cohesion - residual_cohesion) * std::exp( -1.0 * beta * rAlpha );
    }

    else if (rThisVariable == INTERNAL_FRICTION_ANGLE)
    {
        const double peak_friction_angle = rProperties[INTERNAL_FRICTION_ANGLE];
        const double residual_friction_angle = rProperties[INTERNAL_FRICTION_ANGLE_RESIDUAL];
        rHardening = -1.0 * beta * (peak_friction_angle - residual_friction_angle) * std::exp( -1.0 * beta * rAlpha );
    }

    else if (rThisVariable == INTERNAL_DILATANCY_ANGLE)
    {
        const double peak_dilatancy_angle = rProperties[INTERNAL_DILATANCY_ANGLE];
        const double residual_dilatancy_angle = rProperties[INTERNAL_DILATANCY_ANGLE_RESIDUAL];
        rHardening = -1.0 * beta * (peak_dilatancy_angle - residual_dilatancy_angle) * std::exp( -1.0 * beta * rAlpha );
    }

    return rHardening;
}


void ExponentialStrainSofteningLaw::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, MPMHardeningLaw )

}

void ExponentialStrainSofteningLaw::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, MPMHardeningLaw )

}


}  // namespace Kratos.
