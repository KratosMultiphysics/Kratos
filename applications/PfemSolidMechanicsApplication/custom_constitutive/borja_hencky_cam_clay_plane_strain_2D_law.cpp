//
//   Project Name:        KratosSolidMechanicsApplication $
//   Last modified by:    $Author:            JMCarbonell $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes
#include <iostream>

// External includes
#include<cmath>

// Project includes
#include "includes/properties.h"
#include "custom_constitutive/borja_hencky_cam_clay_plane_strain_2D_law.hpp"
#include "pfem_solid_mechanics_application.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

BorjaHenckyCamClayPlasticPlaneStrain2DLaw::BorjaHenckyCamClayPlasticPlaneStrain2DLaw()
    : NonLinearHenckyElasticPlasticPlaneStrain2DLaw()
{
  mpHardeningLaw   = HardeningLaw::Pointer( new CamClayKinematicHardeningLaw() );
  mpYieldCriterion = YieldCriterion::Pointer( new CamClayYieldCriterion(mpHardeningLaw) );
  mpFlowRule       = FlowRule::Pointer( new BorjaCamClayExplicitFlowRule(mpYieldCriterion) );
}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

BorjaHenckyCamClayPlasticPlaneStrain2DLaw::BorjaHenckyCamClayPlasticPlaneStrain2DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
{
  mpHardeningLaw    =  pHardeningLaw;
  mpYieldCriterion  =  YieldCriterion::Pointer( new CamClayYieldCriterion(mpHardeningLaw) );
  mpFlowRule        =  pFlowRule;
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

BorjaHenckyCamClayPlasticPlaneStrain2DLaw::BorjaHenckyCamClayPlasticPlaneStrain2DLaw(const BorjaHenckyCamClayPlasticPlaneStrain2DLaw& rOther)
    : NonLinearHenckyElasticPlasticPlaneStrain2DLaw(rOther)
{

}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer BorjaHenckyCamClayPlasticPlaneStrain2DLaw::Clone() const
{
    BorjaHenckyCamClayPlasticPlaneStrain2DLaw::Pointer p_clone(new BorjaHenckyCamClayPlasticPlaneStrain2DLaw(*this));
    return p_clone;
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

BorjaHenckyCamClayPlasticPlaneStrain2DLaw::~BorjaHenckyCamClayPlasticPlaneStrain2DLaw()
{
}


} // Namespace Kratos
