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

// External includes
#include<cmath>

// Project includes
#include "includes/properties.h"
#include "custom_constitutive/hencky_mc_plane_strain_2D_law.hpp"
#include "particle_mechanics_application.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

HenckyMCPlasticPlaneStrain2DLaw::HenckyMCPlasticPlaneStrain2DLaw()
    : HenckyElasticPlasticPlaneStrain2DLaw()
{
    mpHardeningLaw   = HardeningLaw::Pointer( new HardeningLaw() );
    mpYieldCriterion = YieldCriterion::Pointer( new MCYieldCriterion(mpHardeningLaw) );
    mpMPMFlowRule       = MPMFlowRule::Pointer( new MCPlasticFlowRule(mpYieldCriterion) );
}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

HenckyMCPlasticPlaneStrain2DLaw::HenckyMCPlasticPlaneStrain2DLaw(FlowRulePointer pMPMFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
{
    mpHardeningLaw    =  pHardeningLaw;
    mpYieldCriterion  =  YieldCriterion::Pointer( new MCYieldCriterion(mpHardeningLaw) );
    mpMPMFlowRule        =  pMPMFlowRule;
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

HenckyMCPlasticPlaneStrain2DLaw::HenckyMCPlasticPlaneStrain2DLaw(const HenckyMCPlasticPlaneStrain2DLaw& rOther)
    : HenckyElasticPlasticPlaneStrain2DLaw(rOther)
{

}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer HenckyMCPlasticPlaneStrain2DLaw::Clone() const
{
    HenckyMCPlasticPlaneStrain2DLaw::Pointer p_clone(new HenckyMCPlasticPlaneStrain2DLaw(*this));
    return p_clone;
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

HenckyMCPlasticPlaneStrain2DLaw::~HenckyMCPlasticPlaneStrain2DLaw()
{
}


} // Namespace Kratos
