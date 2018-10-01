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
#include "custom_constitutive/hencky_U_P_Tresca_plane_strain_2D_law.hpp"

#include "pfem_solid_mechanics_application_variables.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

HenckyPlasticUPTrescaPlaneStrain2DLaw::HenckyPlasticUPTrescaPlaneStrain2DLaw()
    : NonLinearHenckyElasticPlasticUPPlaneStrain2DLaw()
{
  mpFlowRule       = FlowRule::Pointer( new TrescaExplicitFlowRule() );
  mpYieldCriterion = YieldCriterion::Pointer( new TrescaYieldCriterion() );
  mpHardeningLaw   = HardeningLaw::Pointer( new HardeningLaw() );
}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

HenckyPlasticUPTrescaPlaneStrain2DLaw::HenckyPlasticUPTrescaPlaneStrain2DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
{
  mpFlowRule        =  pFlowRule;
  mpYieldCriterion  =  YieldCriterion::Pointer( new TrescaYieldCriterion() );
  mpHardeningLaw    =  pHardeningLaw;
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

HenckyPlasticUPTrescaPlaneStrain2DLaw::HenckyPlasticUPTrescaPlaneStrain2DLaw(const HenckyPlasticUPTrescaPlaneStrain2DLaw& rOther)
    : NonLinearHenckyElasticPlasticUPPlaneStrain2DLaw(rOther)
{

}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer HenckyPlasticUPTrescaPlaneStrain2DLaw::Clone() const
{
    HenckyPlasticUPTrescaPlaneStrain2DLaw::Pointer p_clone(new HenckyPlasticUPTrescaPlaneStrain2DLaw(*this));
    return p_clone;
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

HenckyPlasticUPTrescaPlaneStrain2DLaw::~HenckyPlasticUPTrescaPlaneStrain2DLaw()
{
}


} // Namespace Kratos
