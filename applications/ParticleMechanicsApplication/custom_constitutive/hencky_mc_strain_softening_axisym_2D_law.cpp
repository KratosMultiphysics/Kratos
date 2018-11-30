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
#include "custom_constitutive/hencky_mc_strain_softening_axisym_2D_law.hpp"
#include "particle_mechanics_application.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

HenckyMCStrainSofteningPlasticAxisym2DLaw::HenckyMCStrainSofteningPlasticAxisym2DLaw()
    : HenckyElasticPlasticAxisym2DLaw()
{
  mpHardeningLaw      = MPMHardeningLaw::Pointer( new ExponentialStrainSofteningLaw() );
  mpYieldCriterion    = MPMYieldCriterion::Pointer( new MCYieldCriterion(mpHardeningLaw) );
  mpMPMFlowRule       = MPMFlowRule::Pointer( new MCPlasticFlowRule(mpYieldCriterion) );
}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

HenckyMCStrainSofteningPlasticAxisym2DLaw::HenckyMCStrainSofteningPlasticAxisym2DLaw(FlowRulePointer pMPMFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
{
  mpHardeningLaw       =  pHardeningLaw;
  mpYieldCriterion     =  MPMYieldCriterion::Pointer( new MCYieldCriterion(mpHardeningLaw) );
  mpMPMFlowRule        =  pMPMFlowRule;
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

HenckyMCStrainSofteningPlasticAxisym2DLaw::HenckyMCStrainSofteningPlasticAxisym2DLaw(const HenckyMCStrainSofteningPlasticAxisym2DLaw& rOther)
    : HenckyElasticPlasticAxisym2DLaw(rOther)
{

}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer HenckyMCStrainSofteningPlasticAxisym2DLaw::Clone() const
{
    HenckyMCStrainSofteningPlasticAxisym2DLaw::Pointer p_clone(new HenckyMCStrainSofteningPlasticAxisym2DLaw(*this));
    return p_clone;
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

HenckyMCStrainSofteningPlasticAxisym2DLaw::~HenckyMCStrainSofteningPlasticAxisym2DLaw()
{
}


//*********************************CHECK**********************************************
//************************************************************************************

int HenckyMCStrainSofteningPlasticAxisym2DLaw::Check(const Properties& rProperties, const GeometryType& rGeometry, const ProcessInfo& rCurrentProcessInfo)
{
    HenckyElasticPlasticAxisym2DLaw::Check(rProperties, rGeometry, rCurrentProcessInfo);

    KRATOS_ERROR_IF(YOUNG_MODULUS.Key() == 0 || rProperties[YOUNG_MODULUS]<= 0.00) << "YOUNG_MODULUS has Key zero or invalid value " << std::endl;

    const double& nu = rProperties[POISSON_RATIO];
    const double tolerance = 10.e-7;
    const bool check = bool( (nu > 0.5-tolerance ) || (nu < (-1.0 + tolerance)) );

    KRATOS_ERROR_IF(POISSON_RATIO.Key() == 0 || check==true) << "POISSON_RATIO has Key zero invalid value " << std::endl;

    KRATOS_ERROR_IF(COHESION.Key() == 0 || rProperties[COHESION]< 0.00) << "COHESION has Key zero or invalid value " << std::endl;
    KRATOS_ERROR_IF(INTERNAL_FRICTION_ANGLE.Key() == 0 || rProperties[INTERNAL_FRICTION_ANGLE]< 0.00) << "INTERNAL_FRICTION_ANGLE has Key zero or invalid value " << std::endl;

    return 0;
}


} // Namespace Kratos
