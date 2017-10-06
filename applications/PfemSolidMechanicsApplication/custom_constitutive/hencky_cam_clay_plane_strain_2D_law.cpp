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
#include "custom_constitutive/hencky_cam_clay_plane_strain_2D_law.hpp"
#include "pfem_solid_mechanics_application_variables.h"

// TO DO: Remove. This constitutive model is a special case of Borja

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

NonLinearHenckyCamClayPlasticPlaneStrain2DLaw::NonLinearHenckyCamClayPlasticPlaneStrain2DLaw()
    : NonLinearHenckyElasticPlasticPlaneStrain2DLaw()
{
  mpHardeningLaw   = HardeningLaw::Pointer( new CamClayKinematicHardeningLaw() );
  mpYieldCriterion = YieldCriterion::Pointer( new CamClayYieldCriterion(mpHardeningLaw) );
  mpFlowRule       = FlowRule::Pointer( new CamClayExplicitFlowRule(mpYieldCriterion) );
}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

NonLinearHenckyCamClayPlasticPlaneStrain2DLaw::NonLinearHenckyCamClayPlasticPlaneStrain2DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
{
  mpHardeningLaw    =  pHardeningLaw;
  mpYieldCriterion  =  YieldCriterion::Pointer( new CamClayYieldCriterion(mpHardeningLaw) );
  mpFlowRule        =  pFlowRule;
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
