//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    Ilaria Iaconeta, Bodhinanda Chandra
//


// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "custom_constitutive/yield_criteria/mc_yield_criterion.hpp"
#include "mpm_application_variables.h"
#include "utilities/math_utils.h"

namespace Kratos
{


//*******************************CONSTRUCTOR******************************************
//************************************************************************************
MCYieldCriterion::MCYieldCriterion()
    :MPMYieldCriterion()
{

}

//*****************************INITIALIZATION CONSTRUCTOR*****************************
//************************************************************************************

MCYieldCriterion::MCYieldCriterion(HardeningLawPointer pHardeningLaw)
    :MPMYieldCriterion(pHardeningLaw)
{

}

//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

MCYieldCriterion& MCYieldCriterion::operator=(MCYieldCriterion const& rOther)
{
    MPMYieldCriterion::operator=(rOther);
    return *this;
}

//*******************************COPY CONSTRUCTOR*************************************
//************************************************************************************

MCYieldCriterion::MCYieldCriterion(MCYieldCriterion const& rOther)
    :MPMYieldCriterion(rOther)
{

}


//********************************DESTRUCTOR******************************************
//************************************************************************************

MCYieldCriterion::~MCYieldCriterion()
{
}



//************************* CALCULATE YIELD FUNCTION  ******************
//**********************************************************************

double& MCYieldCriterion::CalculateYieldCondition(double& rStateFunction, const Vector& rStressVector, const double& rCohesion, const double& rFrictionAngle, const Properties& rProp)
{
    const double friction_coefficient = (1 + std::sin(rFrictionAngle))/(1 - std::sin(rFrictionAngle));
    const double cohesion_coefficient = 2 * rCohesion * std::sqrt(friction_coefficient);

    // f = k*s1 -s3 - comp
    rStateFunction = friction_coefficient * rStressVector[0] - rStressVector[2] - cohesion_coefficient;

    return rStateFunction;
}

void MCYieldCriterion::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, MPMYieldCriterion )
}

void MCYieldCriterion::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, MPMYieldCriterion )
}


}  // namespace Kratos.
