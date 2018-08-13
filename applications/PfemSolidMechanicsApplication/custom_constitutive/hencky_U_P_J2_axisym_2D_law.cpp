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
#include "custom_constitutive/hencky_U_P_J2_axisym_2D_law.hpp"

#include "pfem_solid_mechanics_application_variables.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

HenckyPlasticUPJ2Axisym2DLaw::HenckyPlasticUPJ2Axisym2DLaw()
    : NonLinearHenckyElasticPlasticUPAxisym2DLaw()
{
  mpFlowRule       = FlowRule::Pointer( new J2ExplicitFlowRule() );
  mpYieldCriterion = YieldCriterion::Pointer( new J2YieldCriterion() );
  mpHardeningLaw   = HardeningLaw::Pointer( new HardeningLaw() );
}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

HenckyPlasticUPJ2Axisym2DLaw::HenckyPlasticUPJ2Axisym2DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
{
  mpFlowRule        =  pFlowRule;
  mpYieldCriterion  =  YieldCriterion::Pointer( new J2YieldCriterion() );
  mpHardeningLaw    =  pHardeningLaw;
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

HenckyPlasticUPJ2Axisym2DLaw::HenckyPlasticUPJ2Axisym2DLaw(const HenckyPlasticUPJ2Axisym2DLaw& rOther)
    : NonLinearHenckyElasticPlasticUPAxisym2DLaw(rOther)
{

}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer HenckyPlasticUPJ2Axisym2DLaw::Clone() const
{
    HenckyPlasticUPJ2Axisym2DLaw::Pointer p_clone(new HenckyPlasticUPJ2Axisym2DLaw(*this));
    return p_clone;
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

HenckyPlasticUPJ2Axisym2DLaw::~HenckyPlasticUPJ2Axisym2DLaw()
{
}


} // Namespace Kratos
