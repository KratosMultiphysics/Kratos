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
    mpMPMFlowRule    = MPMFlowRule::Pointer( new MCPlasticFlowRule(mpYieldCriterion) );
}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

HenckyMCPlasticPlaneStrain2DLaw::HenckyMCPlasticPlaneStrain2DLaw(FlowRulePointer pMPMFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
{
    mpHardeningLaw    =  pHardeningLaw;
    mpYieldCriterion  =  YieldCriterion::Pointer( new MCYieldCriterion(mpHardeningLaw) );
    mpMPMFlowRule     =  pMPMFlowRule;
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

//*********************************CHECK**********************************************
//************************************************************************************

int HenckyMCPlasticPlaneStrain2DLaw::Check(const Properties& rProperties, const GeometryType& rGeometry, const ProcessInfo& rCurrentProcessInfo)
{
    HenckyElasticPlasticPlaneStrain2DLaw::Check(rProperties, rGeometry, rCurrentProcessInfo);
    
    KRATOS_ERROR_IF(YOUNG_MODULUS.Key() == 0 || rProperties[YOUNG_MODULUS]<= 0.00) << "YOUNG_MODULUS has Key zero or invalid value " << std::endl;

    const double& nu = rProperties[POISSON_RATIO];
    const bool check = bool( (nu >0.499 && nu<0.501 ) || (nu < -0.999 && nu > -1.01 ) );

    KRATOS_ERROR_IF(POISSON_RATIO.Key() == 0 || check==true) << "POISSON_RATIO has Key zero invalid value " << std::endl;

    KRATOS_ERROR_IF(COHESION.Key() == 0 || rProperties[COHESION]< 0.00) << "COHESION has Key zero or invalid value " << std::endl;
    KRATOS_ERROR_IF(INTERNAL_FRICTION_ANGLE.Key() == 0 || rProperties[INTERNAL_FRICTION_ANGLE]< 0.00) << "INTERNAL_FRICTION_ANGLE has Key zero or invalid value " << std::endl;
    
    return 0;
}

} // Namespace Kratos
