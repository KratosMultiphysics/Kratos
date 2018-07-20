//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    Ilaria Iaconeta
//

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "custom_constitutive/yield_criteria/mc_yield_criterion.hpp"
#include "particle_mechanics_application.h"
#include "utilities/math_utils.h"

// ROUNDED MohrCoulomb YIELD CRITERION (Sloan & Booker, 1986 ) (in the octahedral)
// Hyperbolic in the other plane... (Abbo & Lyamin, Sloan, Hambleton )
namespace Kratos
{


//*******************************CONSTRUCTOR******************************************
//************************************************************************************
MCYieldCriterion::MCYieldCriterion()
    :YieldCriterion()
{

}

//*****************************INITIALIZATION CONSTRUCTOR*****************************
//************************************************************************************

MCYieldCriterion::MCYieldCriterion(HardeningLawPointer pHardeningLaw)
    :YieldCriterion(pHardeningLaw)
{

}

//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

MCYieldCriterion& MCYieldCriterion::operator=(MCYieldCriterion const& rOther)
{
    YieldCriterion::operator=(rOther);
    return *this;
}

//*******************************COPY CONSTRUCTOR*************************************
//************************************************************************************

MCYieldCriterion::MCYieldCriterion(MCYieldCriterion const& rOther)
    :YieldCriterion(rOther)
{

}


//********************************DESTRUCTOR******************************************
//************************************************************************************

MCYieldCriterion::~MCYieldCriterion()
{
}



//************************* CALCULATE YIELD FUNCTION  ******************
//**********************************************************************

double& MCYieldCriterion::CalculateYieldCondition(double& rStateFunction, const Vector& rStressVector, const double& rAlpha)
{
    const double Cohesion = this->GetHardeningLaw().GetProperties()[COHESION];
    const double FrictionAngle = this->GetHardeningLaw().GetProperties()[INTERNAL_FRICTION_ANGLE];

    const double FrictionCoefficient = (1 + std::sin(FrictionAngle))/(1 - std::sin(FrictionAngle));
    const double CohesionCoefficient = 2 * Cohesion * sqrt(FrictionCoefficient);

    // f = k*s1 -s3 - comp
    rStateFunction = FrictionCoefficient * rStressVector(0) - rStressVector(2) - CohesionCoefficient;

    return rStateFunction;
}

double MCYieldCriterion::GetSmoothingLodeAngle()
{
    return 29.0*GetPI()/180.0;
}

double MCYieldCriterion::GetPI()
{
    return std::atan(1.0)*4.0;
}

double MCYieldCriterion::GetSmoothingHiperbolic()
{
    return 2.0;
}

void MCYieldCriterion::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, YieldCriterion )
}

void MCYieldCriterion::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, YieldCriterion )
}


}  // namespace Kratos.
