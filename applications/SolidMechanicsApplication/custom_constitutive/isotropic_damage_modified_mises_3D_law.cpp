//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:              IPouplana $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2015 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes

// Project includes
#include "custom_constitutive/custom_flow_rules/isotropic_damage_flow_rule.hpp"
#include "custom_constitutive/custom_yield_criteria/modified_mises_yield_criterion.hpp"
#include "custom_constitutive/custom_hardening_laws/modified_exponential_damage_hardening_law.hpp"
#include "custom_constitutive/isotropic_damage_modified_mises_3D_law.hpp"

#include "solid_mechanics_application_variables.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

IsotropicDamageModifiedMises3DLaw::IsotropicDamageModifiedMises3DLaw()
    : LinearElasticPlastic3DLaw()
{
  mpHardeningLaw   = HardeningLaw::Pointer( new ModifiedExponentialDamageHardeningLaw() );
  mpYieldCriterion = YieldCriterion::Pointer( new ModifiedMisesYieldCriterion(mpHardeningLaw) );
  mpFlowRule       = FlowRule::Pointer( new IsotropicDamageFlowRule(mpYieldCriterion) );
}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

IsotropicDamageModifiedMises3DLaw::IsotropicDamageModifiedMises3DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
    : LinearElasticPlastic3DLaw(pFlowRule, pYieldCriterion, pHardeningLaw)
{

}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

IsotropicDamageModifiedMises3DLaw::IsotropicDamageModifiedMises3DLaw(const IsotropicDamageModifiedMises3DLaw& rOther)
    : LinearElasticPlastic3DLaw(rOther)
{

}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer IsotropicDamageModifiedMises3DLaw::Clone() const
{
    return Kratos::make_shared<IsotropicDamageModifiedMises3DLaw>(*this);
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

IsotropicDamageModifiedMises3DLaw::~IsotropicDamageModifiedMises3DLaw()
{
}

//************* COMPUTING  METHODS
//************************************************************************************
//************************************************************************************

//******************CHECK CONSISTENCY IN THE CONSTITUTIVE LAW*************************
//************************************************************************************

int IsotropicDamageModifiedMises3DLaw::Check(const Properties& rMaterialProperties,
                             const GeometryType& rElementGeometry,
                             const ProcessInfo& rCurrentProcessInfo)
{
    int ierr = HyperElasticPlastic3DLaw::Check(rMaterialProperties,rElementGeometry,rCurrentProcessInfo);
    if(ierr != 0) return ierr;

    if(DAMAGE_THRESHOLD.Key() == 0 || rMaterialProperties.Has( DAMAGE_THRESHOLD ) == false || rMaterialProperties[DAMAGE_THRESHOLD] <= 0.0)
        KRATOS_THROW_ERROR( std::invalid_argument,"DAMAGE_THRESHOLD has Key zero, is not defined or has an invalid value for property", rMaterialProperties.Id() )
    if(STRENGTH_RATIO.Key() == 0 || rMaterialProperties.Has( STRENGTH_RATIO ) == false || rMaterialProperties[STRENGTH_RATIO] <= 0.0)
        KRATOS_THROW_ERROR( std::invalid_argument,"STRENGTH_RATIO has Key zero, is not defined or has an invalid value for property", rMaterialProperties.Id() )
    if(RESIDUAL_STRENGTH.Key() == 0 || rMaterialProperties.Has( RESIDUAL_STRENGTH ) == false || rMaterialProperties[RESIDUAL_STRENGTH] < 0.0)
        KRATOS_THROW_ERROR( std::invalid_argument,"RESIDUAL_STRENGTH has Key zero, is not defined or has an invalid value for property", rMaterialProperties.Id() )
    if(SOFTENING_SLOPE.Key() == 0 || rMaterialProperties.Has( SOFTENING_SLOPE ) == false || rMaterialProperties[SOFTENING_SLOPE] < 0.0)
        KRATOS_THROW_ERROR( std::invalid_argument,"SOFTENING_SLOPE has Key zero, is not defined or has an invalid value for property", rMaterialProperties.Id() )

    return ierr;
}

} // Namespace Kratos
