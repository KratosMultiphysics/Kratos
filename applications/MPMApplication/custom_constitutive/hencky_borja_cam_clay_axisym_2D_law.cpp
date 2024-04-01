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
#include "custom_constitutive/hencky_borja_cam_clay_axisym_2D_law.hpp"
#include "mpm_application_variables.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

HenckyBorjaCamClayPlasticAxisym2DLaw::HenckyBorjaCamClayPlasticAxisym2DLaw()
    : HenckyElasticPlasticAxisym2DLaw()
{
    mpHardeningLaw   = MPMHardeningLaw::Pointer( new CamClayHardeningLaw() );
    mpYieldCriterion = MPMYieldCriterion::Pointer( new ModifiedCamClayYieldCriterion(mpHardeningLaw) );
    mpMPMFlowRule    = MPMFlowRule::Pointer( new BorjaCamClayPlasticFlowRule(mpYieldCriterion) );
}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

HenckyBorjaCamClayPlasticAxisym2DLaw::HenckyBorjaCamClayPlasticAxisym2DLaw(FlowRulePointer pMPMFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
{
    mpHardeningLaw    =  pHardeningLaw;
    mpYieldCriterion  =  MPMYieldCriterion::Pointer( new ModifiedCamClayYieldCriterion(mpHardeningLaw) );
    mpMPMFlowRule     =  pMPMFlowRule;
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

HenckyBorjaCamClayPlasticAxisym2DLaw::HenckyBorjaCamClayPlasticAxisym2DLaw(const HenckyBorjaCamClayPlasticAxisym2DLaw& rOther)
    : HenckyElasticPlasticAxisym2DLaw(rOther)
{

}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer HenckyBorjaCamClayPlasticAxisym2DLaw::Clone() const
{
    HenckyBorjaCamClayPlasticAxisym2DLaw::Pointer p_clone(new HenckyBorjaCamClayPlasticAxisym2DLaw(*this));
    return p_clone;
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

HenckyBorjaCamClayPlasticAxisym2DLaw::~HenckyBorjaCamClayPlasticAxisym2DLaw()
{
}


//*********************************CHECK**********************************************
//************************************************************************************

int HenckyBorjaCamClayPlasticAxisym2DLaw::Check(const Properties& rProperties, const GeometryType& rGeometry, const ProcessInfo& rCurrentProcessInfo) const
{
    HenckyElasticPlasticAxisym2DLaw::Check(rProperties, rGeometry, rCurrentProcessInfo);

    KRATOS_ERROR_IF(PRE_CONSOLIDATION_STRESS.Key() == 0 || rProperties[PRE_CONSOLIDATION_STRESS] >= 0.00) << "PRE_CONSOLIDATION_STRESS has Key zero or invalid value (Expected negative value) " << std::endl;
    KRATOS_ERROR_IF(OVER_CONSOLIDATION_RATIO.Key() == 0 || rProperties[OVER_CONSOLIDATION_RATIO] <= 0.00) << "OVER_CONSOLIDATION_RATIO has Key zero invalid value " << std::endl;

    KRATOS_ERROR_IF(SWELLING_SLOPE.Key() == 0 || rProperties[SWELLING_SLOPE] <= 0.00) << "SWELLING_SLOPE has Key zero or invalid value " << std::endl;
    KRATOS_ERROR_IF(NORMAL_COMPRESSION_SLOPE.Key() == 0 || rProperties[NORMAL_COMPRESSION_SLOPE] <= 0.00) << "NORMAL_COMPRESSION_SLOPE has Key zero or invalid value " << std::endl;
    KRATOS_ERROR_IF(CRITICAL_STATE_LINE.Key() == 0 || rProperties[CRITICAL_STATE_LINE] <= 0.00) << "CRITICAL_STATE_LINE has Key zero or invalid value " << std::endl;
    KRATOS_ERROR_IF(INITIAL_SHEAR_MODULUS.Key() == 0 || rProperties[INITIAL_SHEAR_MODULUS] <= 0.00) << "INITIAL_SHEAR_MODULUS has Key zero or invalid value " << std::endl;

    KRATOS_ERROR_IF(ALPHA_SHEAR.Key() == 0 ) << "ALPHA_SHEAR has Key zero " << std::endl;

    return 0;
}


} // Namespace Kratos
