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
#include "custom_constitutive/hencky_cam_clay_plane_strain_2D_law.hpp"
#include "pfem_solid_mechanics_application.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

NonLinearHenckyCamClayPlasticPlaneStrain2DLaw::NonLinearHenckyCamClayPlasticPlaneStrain2DLaw()
    : NonLinearHenckyElasticPlasticPlaneStrain2DLaw()
{
  mpFlowRule       = FlowRule::Pointer( new CamClayExplicitFlowRule() );
  mpYieldCriterion = YieldCriterion::Pointer( new CamClayYieldCriterion() );
  mpHardeningLaw   = HardeningLaw::Pointer( new CamClayKinematicHardeningLaw() );
}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

NonLinearHenckyCamClayPlasticPlaneStrain2DLaw::NonLinearHenckyCamClayPlasticPlaneStrain2DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
{
  mpFlowRule        =  pFlowRule;
  mpYieldCriterion  =  YieldCriterion::Pointer( new CamClayYieldCriterion() );
  mpHardeningLaw    =  pHardeningLaw;
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

NonLinearHenckyCamClayPlasticPlaneStrain2DLaw::NonLinearHenckyCamClayPlasticPlaneStrain2DLaw(const NonLinearHenckyCamClayPlasticPlaneStrain2DLaw& rOther)
    : NonLinearHenckyElasticPlasticPlaneStrain2DLaw(rOther)
{

}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer NonLinearHenckyCamClayPlasticPlaneStrain2DLaw::Clone() const
{
    NonLinearHenckyCamClayPlasticPlaneStrain2DLaw::Pointer p_clone(new NonLinearHenckyCamClayPlasticPlaneStrain2DLaw(*this));
    return p_clone;
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

NonLinearHenckyCamClayPlasticPlaneStrain2DLaw::~NonLinearHenckyCamClayPlasticPlaneStrain2DLaw()
{
}


} // Namespace Kratos
