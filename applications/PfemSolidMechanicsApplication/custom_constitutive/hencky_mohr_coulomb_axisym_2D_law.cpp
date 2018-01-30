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
#include "custom_constitutive/hencky_mohr_coulomb_axisym_2D_law.hpp"
#include "pfem_solid_mechanics_application_variables.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

HenckyMohrCoulombPlasticAxisym2DLaw::HenckyMohrCoulombPlasticAxisym2DLaw()
    : NonLinearHenckyElasticPlasticAxisym2DLaw()
{
  mpFlowRule       = FlowRule::Pointer( new MohrCoulombExplicitFlowRule() );
  mpYieldCriterion = YieldCriterion::Pointer( new MohrCoulombYieldCriterion() );
  mpHardeningLaw   = HardeningLaw::Pointer( new HardeningLaw() );
}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

HenckyMohrCoulombPlasticAxisym2DLaw::HenckyMohrCoulombPlasticAxisym2DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
{
  mpFlowRule        =  pFlowRule;
  mpYieldCriterion  =  YieldCriterion::Pointer( new MohrCoulombYieldCriterion() );
  mpHardeningLaw    =  pHardeningLaw;
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

HenckyMohrCoulombPlasticAxisym2DLaw::HenckyMohrCoulombPlasticAxisym2DLaw(const HenckyMohrCoulombPlasticAxisym2DLaw& rOther)
    : NonLinearHenckyElasticPlasticAxisym2DLaw(rOther)
{

}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer HenckyMohrCoulombPlasticAxisym2DLaw::Clone() const
{
    HenckyMohrCoulombPlasticAxisym2DLaw::Pointer p_clone(new HenckyMohrCoulombPlasticAxisym2DLaw(*this));
    return p_clone;
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

HenckyMohrCoulombPlasticAxisym2DLaw::~HenckyMohrCoulombPlasticAxisym2DLaw()
{
}


} // Namespace Kratos
