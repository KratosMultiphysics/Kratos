//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Created by:          $Author:                  LMonforte $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                    July 2015 $
//   Revision:            $Revision:                      0.0 $
//
//

// System includes
#include <iostream>

// External includes
#include<cmath>

// Project includes
#include "includes/properties.h"
#include "custom_constitutive/linear_hencky_cam_clay_axisym_2D_law.hpp"
#include "pfem_solid_mechanics_application_variables.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

LinearHenckyCamClayPlasticAxisym2DLaw::LinearHenckyCamClayPlasticAxisym2DLaw()
    : NonLinearHenckyElasticPlasticAxisym2DLaw()
{
  mpHardeningLaw   = HardeningLaw::Pointer( new CamClayKinematicHardeningLaw() );
  mpYieldCriterion = YieldCriterion::Pointer( new CamClayYieldCriterion(mpHardeningLaw) );
  mpFlowRule       = FlowRule::Pointer( new LinearCamClayExplicitFlowRule(mpYieldCriterion) );
}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

LinearHenckyCamClayPlasticAxisym2DLaw::LinearHenckyCamClayPlasticAxisym2DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
{
  mpHardeningLaw    =  pHardeningLaw;
  mpYieldCriterion  =  YieldCriterion::Pointer( new CamClayYieldCriterion(mpHardeningLaw) );
  mpFlowRule        =  pFlowRule;
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

LinearHenckyCamClayPlasticAxisym2DLaw::LinearHenckyCamClayPlasticAxisym2DLaw(const LinearHenckyCamClayPlasticAxisym2DLaw& rOther)
    : NonLinearHenckyElasticPlasticAxisym2DLaw(rOther)
{

}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer LinearHenckyCamClayPlasticAxisym2DLaw::Clone() const
{
    LinearHenckyCamClayPlasticAxisym2DLaw::Pointer p_clone(new LinearHenckyCamClayPlasticAxisym2DLaw(*this));
    return p_clone;
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

LinearHenckyCamClayPlasticAxisym2DLaw::~LinearHenckyCamClayPlasticAxisym2DLaw()
{
}


} // Namespace Kratos
