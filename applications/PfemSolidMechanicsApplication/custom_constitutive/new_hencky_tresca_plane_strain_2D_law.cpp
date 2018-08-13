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
#include "custom_constitutive/new_hencky_tresca_plane_strain_2D_law.hpp"
#include "pfem_solid_mechanics_application_variables.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

NewHenckyTrescaPlasticPlaneStrain2DLaw::NewHenckyTrescaPlasticPlaneStrain2DLaw()
    : HenckyTrescaPlasticPlaneStrain2DLaw()
{
  mpFlowRule       = FlowRule::Pointer( new TrescaExplicitFlowRule() );  //not relevant since it is associated
  mpYieldCriterion = YieldCriterion::Pointer( new NewTrescaYieldCriterion() );
  mpHardeningLaw   = HardeningLaw::Pointer( new HardeningLaw() );
}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

NewHenckyTrescaPlasticPlaneStrain2DLaw::NewHenckyTrescaPlasticPlaneStrain2DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
{
  mpHardeningLaw    =  pHardeningLaw;
  mpYieldCriterion  =  YieldCriterion::Pointer( new NewTrescaYieldCriterion(mpHardeningLaw) );
  mpFlowRule        =  pFlowRule;
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

NewHenckyTrescaPlasticPlaneStrain2DLaw::NewHenckyTrescaPlasticPlaneStrain2DLaw(const NewHenckyTrescaPlasticPlaneStrain2DLaw& rOther)
    : HenckyTrescaPlasticPlaneStrain2DLaw(rOther)
{

}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer NewHenckyTrescaPlasticPlaneStrain2DLaw::Clone() const
{
    NewHenckyTrescaPlasticPlaneStrain2DLaw::Pointer p_clone(new NewHenckyTrescaPlasticPlaneStrain2DLaw(*this));
    return p_clone;
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

NewHenckyTrescaPlasticPlaneStrain2DLaw::~NewHenckyTrescaPlasticPlaneStrain2DLaw()
{
}


} // Namespace Kratos
