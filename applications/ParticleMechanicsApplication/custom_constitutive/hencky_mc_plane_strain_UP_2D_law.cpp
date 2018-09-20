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
#include <iostream>
#include <cmath>

// External includes

// Project includes
#include "includes/properties.h"
#include "custom_constitutive/hencky_mc_plane_strain_UP_2D_law.hpp"
#include "particle_mechanics_application.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

HenckyMCPlasticPlaneStrainUP2DLaw::HenckyMCPlasticPlaneStrainUP2DLaw()
    : HenckyElasticPlasticPlaneStrainUP2DLaw()
{
    mpHardeningLaw   = HardeningLaw::Pointer( new HardeningLaw() );
    mpYieldCriterion = YieldCriterion::Pointer( new MCYieldCriterion(mpHardeningLaw) );
    mpMPMFlowRule    = MPMFlowRule::Pointer( new MCPlasticFlowRule(mpYieldCriterion) );
}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

HenckyMCPlasticPlaneStrainUP2DLaw::HenckyMCPlasticPlaneStrainUP2DLaw(MPMFlowRulePointer pMPMFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
{
    mpHardeningLaw    =  pHardeningLaw;
    mpYieldCriterion  =  YieldCriterion::Pointer( new MCYieldCriterion(mpHardeningLaw) );
    mpMPMFlowRule     =  pMPMFlowRule;
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

HenckyMCPlasticPlaneStrainUP2DLaw::HenckyMCPlasticPlaneStrainUP2DLaw(const HenckyMCPlasticPlaneStrainUP2DLaw& rOther)
    : HenckyElasticPlasticPlaneStrainUP2DLaw(rOther)
{

}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer HenckyMCPlasticPlaneStrainUP2DLaw::Clone() const
{
    HenckyMCPlasticPlaneStrainUP2DLaw::Pointer p_clone(new HenckyMCPlasticPlaneStrainUP2DLaw(*this));
    return p_clone;
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

HenckyMCPlasticPlaneStrainUP2DLaw::~HenckyMCPlasticPlaneStrainUP2DLaw()
{
}


} // Namespace Kratos
