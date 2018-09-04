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
#include <iostream>

// External includes
#include<cmath>

// Project includes
#include "includes/properties.h"
#include "custom_constitutive/hencky_mc_strain_softening_3D_law.hpp"
#include "particle_mechanics_application.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

HenckyMCStrainSofteningPlastic3DLaw::HenckyMCStrainSofteningPlastic3DLaw()
    : HenckyElasticPlastic3DLaw()
{
    mpHardeningLaw   = HardeningLaw::Pointer( new ExponentialStrainSofteningLaw() );
    mpYieldCriterion = YieldCriterion::Pointer( new MCYieldCriterion(mpHardeningLaw) );
    mpMPMFlowRule    = MPMFlowRule::Pointer( new MCStrainSofteningPlasticFlowRule(mpYieldCriterion) );
}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

HenckyMCStrainSofteningPlastic3DLaw::HenckyMCStrainSofteningPlastic3DLaw(FlowRulePointer pMPMFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
{
    mpHardeningLaw    =  pHardeningLaw;
    mpYieldCriterion  =  YieldCriterion::Pointer( new MCYieldCriterion(mpHardeningLaw) );
    mpMPMFlowRule        =  pMPMFlowRule;
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

HenckyMCStrainSofteningPlastic3DLaw::HenckyMCStrainSofteningPlastic3DLaw(const HenckyMCStrainSofteningPlastic3DLaw& rOther)
    : HenckyElasticPlastic3DLaw(rOther)
{

}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer HenckyMCStrainSofteningPlastic3DLaw::Clone() const
{
    HenckyMCStrainSofteningPlastic3DLaw::Pointer p_clone(new HenckyMCStrainSofteningPlastic3DLaw(*this));
    return p_clone;
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

HenckyMCStrainSofteningPlastic3DLaw::~HenckyMCStrainSofteningPlastic3DLaw()
{
}


} // Namespace Kratos
