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
#include "custom_constitutive/hencky_J2_plane_strain_2D_law.hpp"
#include "pfem_solid_mechanics_application.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

HenckyJ2PlasticPlaneStrain2DLaw::HenckyJ2PlasticPlaneStrain2DLaw()
    : NonLinearHenckyElasticPlasticPlaneStrain2DLaw()
{
  mpFlowRule       = FlowRule::Pointer( new J2ExplicitFlowRule() );
  mpYieldCriterion = YieldCriterion::Pointer( new J2YieldCriterion() );
  mpHardeningLaw   = HardeningLaw::Pointer( new HardeningLaw() );
}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

HenckyJ2PlasticPlaneStrain2DLaw::HenckyJ2PlasticPlaneStrain2DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
{
  mpFlowRule        =  pFlowRule;
  mpYieldCriterion  =  YieldCriterion::Pointer( new J2YieldCriterion() );
  mpHardeningLaw    =  pHardeningLaw;
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
