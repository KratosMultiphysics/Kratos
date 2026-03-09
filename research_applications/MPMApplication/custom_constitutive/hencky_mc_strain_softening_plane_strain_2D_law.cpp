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
#include <cmath>

// External includes

// Project includes
#include "includes/properties.h"
#include "custom_constitutive/hencky_mc_strain_softening_plane_strain_2D_law.hpp"
#include "mpm_application_variables.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

HenckyMCStrainSofteningPlasticPlaneStrain2DLaw::HenckyMCStrainSofteningPlasticPlaneStrain2DLaw()
    : HenckyElasticPlasticPlaneStrain2DLaw()
{
    mpHardeningLaw   = MPMHardeningLaw::Pointer( new ExponentialStrainSofteningLaw() );
    mpYieldCriterion = MPMYieldCriterion::Pointer( new MCYieldCriterion(mpHardeningLaw) );
    mpMPMFlowRule    = MPMFlowRule::Pointer( new MCStrainSofteningPlasticFlowRule(mpYieldCriterion) );
}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

HenckyMCStrainSofteningPlasticPlaneStrain2DLaw::HenckyMCStrainSofteningPlasticPlaneStrain2DLaw(FlowRulePointer pMPMFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
{
    mpHardeningLaw    =  pHardeningLaw;
    mpYieldCriterion  =  MPMYieldCriterion::Pointer( new MCYieldCriterion(mpHardeningLaw) );
    mpMPMFlowRule     =  pMPMFlowRule;
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

HenckyMCStrainSofteningPlasticPlaneStrain2DLaw::HenckyMCStrainSofteningPlasticPlaneStrain2DLaw(const HenckyMCStrainSofteningPlasticPlaneStrain2DLaw& rOther)
    : HenckyElasticPlasticPlaneStrain2DLaw(rOther)
{

}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer HenckyMCStrainSofteningPlasticPlaneStrain2DLaw::Clone() const
{
    HenckyMCStrainSofteningPlasticPlaneStrain2DLaw::Pointer p_clone(new HenckyMCStrainSofteningPlasticPlaneStrain2DLaw(*this));
    return p_clone;
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

HenckyMCStrainSofteningPlasticPlaneStrain2DLaw::~HenckyMCStrainSofteningPlasticPlaneStrain2DLaw()
{
}

//*********************************CHECK**********************************************
//************************************************************************************

int HenckyMCStrainSofteningPlasticPlaneStrain2DLaw::Check(const Properties& rProperties, const GeometryType& rGeometry, const ProcessInfo& rCurrentProcessInfo) const
{
    HenckyElasticPlasticPlaneStrain2DLaw::Check(rProperties, rGeometry, rCurrentProcessInfo);

    KRATOS_ERROR_IF(YOUNG_MODULUS.Key() == 0 || rProperties[YOUNG_MODULUS]<= 0.00) << "YOUNG_MODULUS has Key zero or invalid value " << std::endl;

    const double& nu = rProperties[POISSON_RATIO];
    const double tolerance = 10.e-7;
    const bool check = bool( (nu > 0.5-tolerance ) || (nu < (-1.0 + tolerance)) );

    KRATOS_ERROR_IF(POISSON_RATIO.Key() == 0 || check==true) << "POISSON_RATIO has Key zero invalid value " << std::endl;

    KRATOS_ERROR_IF(COHESION.Key() == 0 || rProperties[COHESION]< 0.00) << "COHESION has Key zero or invalid value " << std::endl;
    KRATOS_ERROR_IF(INTERNAL_FRICTION_ANGLE.Key() == 0 || rProperties[INTERNAL_FRICTION_ANGLE]< 0.00) << "INTERNAL_FRICTION_ANGLE has Key zero or invalid value " << std::endl;

    KRATOS_ERROR_IF(COHESION_RESIDUAL.Key() == 0 || rProperties[COHESION_RESIDUAL]< 0.00) << "COHESION_RESIDUAL has Key zero or invalid value " << std::endl;
    KRATOS_ERROR_IF(INTERNAL_FRICTION_ANGLE_RESIDUAL.Key() == 0 || rProperties[INTERNAL_FRICTION_ANGLE_RESIDUAL]< 0.00) << "INTERNAL_FRICTION_ANGLE_RESIDUAL has Key zero or invalid value " << std::endl;

    KRATOS_ERROR_IF(INTERNAL_DILATANCY_ANGLE_RESIDUAL.Key() == 0 || rProperties[INTERNAL_DILATANCY_ANGLE_RESIDUAL] != 0.00) << "INTERNAL_DILATANCY_ANGLE_RESIDUAL has Key zero or invalid value " << std::endl;

    KRATOS_ERROR_IF(SHAPE_FUNCTION_BETA.Key() == 0 || rProperties[SHAPE_FUNCTION_BETA] < 0.00) << "SHAPE_FUNCTION_BETA has Key zero or invalid value " << std::endl;

    return 0;
}


} // Namespace Kratos
