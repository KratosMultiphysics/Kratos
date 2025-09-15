//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Created by:          $Author:                  LMonforte $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                    July 2015 $
//   Revision:            $Revision:                      0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_constitutive/hencky_J2_plane_strain_2D_law.hpp"

#include "pfem_solid_mechanics_application_variables.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

HenckyJ2PlasticPlaneStrain2DLaw::HenckyJ2PlasticPlaneStrain2DLaw()
    : NonLinearHenckyElasticPlasticPlaneStrain2DLaw()
{
  mpHardeningLaw   = HardeningLaw::Pointer( new HardeningLaw() );
  mpYieldCriterion = YieldCriterion::Pointer( new J2YieldCriterion(mpHardeningLaw) );
  mpFlowRule       = FlowRule::Pointer( new J2ExplicitFlowRule(mpYieldCriterion) );
}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

HenckyJ2PlasticPlaneStrain2DLaw::HenckyJ2PlasticPlaneStrain2DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
{
  mpHardeningLaw    =  pHardeningLaw;
  mpYieldCriterion  =  YieldCriterion::Pointer( new J2YieldCriterion(mpHardeningLaw) );
  mpFlowRule        =  pFlowRule;
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

HenckyJ2PlasticPlaneStrain2DLaw::HenckyJ2PlasticPlaneStrain2DLaw(const HenckyJ2PlasticPlaneStrain2DLaw& rOther)
    : NonLinearHenckyElasticPlasticPlaneStrain2DLaw(rOther)
{

}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer HenckyJ2PlasticPlaneStrain2DLaw::Clone() const
{
    HenckyJ2PlasticPlaneStrain2DLaw::Pointer p_clone(new HenckyJ2PlasticPlaneStrain2DLaw(*this));
    return p_clone;
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

HenckyJ2PlasticPlaneStrain2DLaw::~HenckyJ2PlasticPlaneStrain2DLaw()
{
}


} // Namespace Kratos
