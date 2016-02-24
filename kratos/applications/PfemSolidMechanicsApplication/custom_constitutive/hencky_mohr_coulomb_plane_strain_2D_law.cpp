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
#include "custom_constitutive/hencky_mohr_coulomb_plane_strain_2D_law.hpp"
#include "pfem_solid_mechanics_application_variables.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

HenckyMohrCoulombPlasticPlaneStrain2DLaw::HenckyMohrCoulombPlasticPlaneStrain2DLaw()
    : NonLinearHenckyElasticPlasticPlaneStrain2DLaw()
{
  mpHardeningLaw   = HardeningLaw::Pointer( new HardeningLaw() );
  mpYieldCriterion = YieldCriterion::Pointer( new MohrCoulombYieldCriterion(mpHardeningLaw) );
  mpFlowRule       = FlowRule::Pointer( new MohrCoulombExplicitFlowRule(mpYieldCriterion) );
}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

HenckyMohrCoulombPlasticPlaneStrain2DLaw::HenckyMohrCoulombPlasticPlaneStrain2DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
{
  mpHardeningLaw    =  pHardeningLaw;
  mpYieldCriterion  =  YieldCriterion::Pointer( new MohrCoulombYieldCriterion(mpHardeningLaw) );
  mpFlowRule        =  pFlowRule;
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

HenckyMohrCoulombPlasticPlaneStrain2DLaw::HenckyMohrCoulombPlasticPlaneStrain2DLaw(const HenckyMohrCoulombPlasticPlaneStrain2DLaw& rOther)
    : NonLinearHenckyElasticPlasticPlaneStrain2DLaw(rOther)
{

}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer HenckyMohrCoulombPlasticPlaneStrain2DLaw::Clone() const
{
    HenckyMohrCoulombPlasticPlaneStrain2DLaw::Pointer p_clone(new HenckyMohrCoulombPlasticPlaneStrain2DLaw(*this));
    return p_clone;
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

HenckyMohrCoulombPlasticPlaneStrain2DLaw::~HenckyMohrCoulombPlasticPlaneStrain2DLaw()
{
}


} // Namespace Kratos
