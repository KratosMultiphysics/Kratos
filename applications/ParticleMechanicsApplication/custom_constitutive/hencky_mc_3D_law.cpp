//
//   Project Name:        KratosParticleMechanicsApplication $
//   Created by:          $Author:                 IIaconeta $
//   Date:                $Date:               February 2017 $
//   Revision:            $Revision:                     0.0 $
//
//

// System includes
#include <iostream>

// External includes
#include<cmath>

// Project includes
#include "includes/properties.h"
#include "custom_constitutive/hencky_mc_3D_law.hpp"
#include "particle_mechanics_application.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

HenckyMCPlastic3DLaw::HenckyMCPlastic3DLaw()
    : HenckyElasticPlastic3DLaw()
{
  mpHardeningLaw   = HardeningLaw::Pointer( new HardeningLaw() );
  mpYieldCriterion = YieldCriterion::Pointer( new MCYieldCriterion(mpHardeningLaw) );
  mpMPMFlowRule       = MPMFlowRule::Pointer( new MCPlasticFlowRule(mpYieldCriterion) );
}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

HenckyMCPlastic3DLaw::HenckyMCPlastic3DLaw(FlowRulePointer pMPMFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
{
  mpHardeningLaw    =  pHardeningLaw;
  mpYieldCriterion  =  YieldCriterion::Pointer( new MCYieldCriterion(mpHardeningLaw) );
  mpMPMFlowRule        =  pMPMFlowRule;
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

HenckyMCPlastic3DLaw::HenckyMCPlastic3DLaw(const HenckyMCPlastic3DLaw& rOther)
    : HenckyElasticPlastic3DLaw(rOther)
{

}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer HenckyMCPlastic3DLaw::Clone() const
{
    HenckyMCPlastic3DLaw::Pointer p_clone(new HenckyMCPlastic3DLaw(*this));
    return p_clone;
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

HenckyMCPlastic3DLaw::~HenckyMCPlastic3DLaw()
{
}


} // Namespace Kratos
