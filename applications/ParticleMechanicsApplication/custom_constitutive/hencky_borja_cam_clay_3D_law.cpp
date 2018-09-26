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
#include "custom_constitutive/hencky_borja_cam_clay_3D_law.hpp"
#include "particle_mechanics_application.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

HenckyBorjaCamClayPlastic3DLaw::HenckyBorjaCamClayPlastic3DLaw()
    : HenckyElasticPlastic3DLaw()
{
    mpHardeningLaw   = HardeningLaw::Pointer( new CamClayHardeningLaw() );
    mpYieldCriterion = YieldCriterion::Pointer( new ModifiedCamClayYieldCriterion(mpHardeningLaw) );
    mpMPMFlowRule    = MPMFlowRule::Pointer( new BorjaCamClayPlasticFlowRule(mpYieldCriterion) );
}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

HenckyBorjaCamClayPlastic3DLaw::HenckyBorjaCamClayPlastic3DLaw(FlowRulePointer pMPMFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
{
    mpHardeningLaw    =  pHardeningLaw;
    mpYieldCriterion  =  YieldCriterion::Pointer( new ModifiedCamClayYieldCriterion(mpHardeningLaw) );
    mpMPMFlowRule     =  pMPMFlowRule;
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

HenckyBorjaCamClayPlastic3DLaw::HenckyBorjaCamClayPlastic3DLaw(const HenckyBorjaCamClayPlastic3DLaw& rOther)
    : HenckyElasticPlastic3DLaw(rOther)
{

}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer HenckyBorjaCamClayPlastic3DLaw::Clone() const
{
    HenckyBorjaCamClayPlastic3DLaw::Pointer p_clone(new HenckyBorjaCamClayPlastic3DLaw(*this));
    return p_clone;
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

HenckyBorjaCamClayPlastic3DLaw::~HenckyBorjaCamClayPlastic3DLaw()
{
}

//*********************************CHECK**********************************************
//************************************************************************************

int HenckyBorjaCamClayPlastic3DLaw::Check(const Properties& rProperties, const GeometryType& rGeometry, const ProcessInfo& rCurrentProcessInfo)
{
    HenckyElasticPlastic3DLaw::Check(rProperties, rGeometry, rCurrentProcessInfo);

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
