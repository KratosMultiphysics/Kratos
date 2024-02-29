
//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ignasi de Pouplana
//


#if !defined(KRATOS_POROMECHANICS_APPLICATION_H_INCLUDED )
#define  KRATOS_POROMECHANICS_APPLICATION_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"

// Application includes
#include "poromechanics_application_variables.h"

#include "custom_conditions/one-phase_flow/U_Pl_force_condition.hpp"
#include "custom_conditions/one-phase_flow/U_Pl_face_load_condition.hpp"
#include "custom_conditions/one-phase_flow/U_Pl_normal_face_load_condition.hpp"
#include "custom_conditions/one-phase_flow/U_Pl_liquid_discharge_condition.hpp"
#include "custom_conditions/one-phase_flow/U_Pl_normal_liquid_flux_condition.hpp"
#include "custom_conditions/one-phase_flow/U_Pl_normal_liquid_flux_FIC_condition.hpp"
#include "custom_conditions/one-phase_flow/U_Pl_face_load_interface_condition.hpp"
#include "custom_conditions/one-phase_flow/U_Pl_normal_liquid_flux_interface_condition.hpp"
#include "custom_conditions/one-phase_flow/line_load_2D_diff_order_condition.hpp"
#include "custom_conditions/one-phase_flow/line_normal_load_2D_diff_order_condition.hpp"
#include "custom_conditions/one-phase_flow/line_normal_liquid_flux_2D_diff_order_condition.hpp"
#include "custom_conditions/one-phase_flow/surface_load_3D_diff_order_condition.hpp"
#include "custom_conditions/one-phase_flow/surface_normal_load_3D_diff_order_condition.hpp"
#include "custom_conditions/one-phase_flow/surface_normal_liquid_flux_3D_diff_order_condition.hpp"

#include "custom_elements/one-phase_flow/U_Pl_small_strain_element.hpp"
#include "custom_elements/one-phase_flow/U_Pl_small_strain_interface_element.hpp"
#include "custom_elements/one-phase_flow/U_Pl_small_strain_link_interface_element.hpp"
#include "custom_elements/one-phase_flow/U_Pl_small_strain_FIC_element.hpp"
#include "custom_elements/one-phase_flow/small_strain_U_Pl_diff_order_element.hpp"

#include "custom_constitutive/interface_element_laws/bilinear_cohesive_3D_law.hpp"
#include "custom_constitutive/interface_element_laws/bilinear_cohesive_2D_law.hpp"
#include "custom_constitutive/interface_element_laws/elastic_cohesive_3D_law.hpp"
#include "custom_constitutive/interface_element_laws/elastic_cohesive_2D_law.hpp"
#include "custom_constitutive/interface_element_laws/isotropic_damage_cohesive_3D_law.hpp"
#include "custom_constitutive/interface_element_laws/isotropic_damage_cohesive_2D_law.hpp"
#include "custom_constitutive/interface_element_laws/elastoplastic_mohr_coulomb_cohesive_3D_law.hpp"
#include "custom_constitutive/interface_element_laws/elastoplastic_mohr_coulomb_cohesive_2D_law.hpp"
#include "custom_constitutive/interface_element_laws/elastoplastic_mod_mohr_coulomb_cohesive_3D_law.hpp"
#include "custom_constitutive/interface_element_laws/elastoplastic_mod_mohr_coulomb_cohesive_2D_law.hpp"
#include "custom_constitutive/interface_element_laws/exponential_cohesive_3D_law.hpp"
#include "custom_constitutive/interface_element_laws/exponential_cohesive_2D_law.hpp"

#include "custom_constitutive/continuum_laws/custom_flow_rules/local_damage_flow_rule.hpp"
#include "custom_constitutive/continuum_laws/custom_flow_rules/nonlocal_damage_flow_rule.hpp"

#include "custom_constitutive/continuum_laws/simo_ju_local_damage_3D_law.hpp"
#include "custom_constitutive/continuum_laws/simo_ju_local_damage_plane_strain_2D_law.hpp"
#include "custom_constitutive/continuum_laws/simo_ju_local_damage_plane_stress_2D_law.hpp"

#include "custom_constitutive/continuum_laws/simo_ju_nonlocal_damage_3D_law.hpp"
#include "custom_constitutive/continuum_laws/simo_ju_nonlocal_damage_plane_strain_2D_law.hpp"
#include "custom_constitutive/continuum_laws/simo_ju_nonlocal_damage_plane_stress_2D_law.hpp"

#include "custom_constitutive/continuum_laws/modified_mises_nonlocal_damage_3D_law.hpp"
#include "custom_constitutive/continuum_laws/modified_mises_nonlocal_damage_plane_strain_2D_law.hpp"
#include "custom_constitutive/continuum_laws/modified_mises_nonlocal_damage_plane_stress_2D_law.hpp"

#include "custom_constitutive/continuum_laws/history_linear_elastic_3D_law.hpp"
#include "custom_constitutive/continuum_laws/history_linear_elastic_plane_strain_2D_law.hpp"
#include "custom_constitutive/continuum_laws/history_linear_elastic_plane_stress_2D_law.hpp"

#include "custom_constitutive/continuum_laws/custom_flow_rules/isotropic_damage_flow_rule.hpp"
#include "custom_constitutive/continuum_laws/custom_yield_criteria/simo_ju_yield_criterion.hpp"
#include "custom_constitutive/continuum_laws/custom_yield_criteria/modified_mises_yield_criterion.hpp"
#include "custom_constitutive/continuum_laws/custom_hardening_laws/exponential_damage_hardening_law.hpp"
#include "custom_constitutive/continuum_laws/custom_hardening_laws/modified_exponential_damage_hardening_law.hpp"
#include "custom_constitutive/continuum_laws/hyperelastic_3D_law.hpp"
#include "custom_constitutive/continuum_laws/linear_elastic_3D_law.hpp"
#include "custom_constitutive/continuum_laws/linear_elastic_plane_strain_2D_law.hpp"
#include "custom_constitutive/continuum_laws/linear_elastic_plane_stress_2D_law.hpp"

namespace Kratos
{

class KRATOS_API(POROMECHANICS_APPLICATION) KratosPoromechanicsApplication : public KratosApplication
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(KratosPoromechanicsApplication);

    // Default constructor
    KratosPoromechanicsApplication();

    // Destructor
    ~KratosPoromechanicsApplication() override {}


    void Register() override;

    // Turn back information as a string
    std::string Info() const override
    {
        return "KratosPoromechanicsApplication";
    }

    // Print information about this object
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
        PrintData(rOStream);
    }

    // Print object's data
    void PrintData(std::ostream& rOStream) const override
    {
        KRATOS_WATCH("in my application");
        KRATOS_WATCH(KratosComponents<VariableData>::GetComponents().size() );
        rOStream << "Variables:" << std::endl;
        KratosComponents<VariableData>().PrintData(rOStream);
        rOStream << std::endl;
        rOStream << "Elements:" << std::endl;
        KratosComponents<Element>().PrintData(rOStream);
        rOStream << std::endl;
        rOStream << "Conditions:" << std::endl;
        KratosComponents<Condition>().PrintData(rOStream);
    }

private:

// Member Variables

const UPlSmallStrainElement<2,3> mUPlSmallStrainElement2D3N;
const UPlSmallStrainElement<2,4> mUPlSmallStrainElement2D4N;
const UPlSmallStrainElement<3,4> mUPlSmallStrainElement3D4N;
const UPlSmallStrainElement<3,8> mUPlSmallStrainElement3D8N;

const UPlSmallStrainInterfaceElement<2,4> mUPlSmallStrainInterfaceElement2D4N;
const UPlSmallStrainInterfaceElement<3,6> mUPlSmallStrainInterfaceElement3D6N;
const UPlSmallStrainInterfaceElement<3,8> mUPlSmallStrainInterfaceElement3D8N;

const UPlSmallStrainLinkInterfaceElement<2,4> mUPlSmallStrainLinkInterfaceElement2D4N;
const UPlSmallStrainLinkInterfaceElement<3,6> mUPlSmallStrainLinkInterfaceElement3D6N;
const UPlSmallStrainLinkInterfaceElement<3,8> mUPlSmallStrainLinkInterfaceElement3D8N;

const UPlSmallStrainFICElement<2,3> mUPlSmallStrainFICElement2D3N;
const UPlSmallStrainFICElement<2,4> mUPlSmallStrainFICElement2D4N;
const UPlSmallStrainFICElement<3,4> mUPlSmallStrainFICElement3D4N;
const UPlSmallStrainFICElement<3,8> mUPlSmallStrainFICElement3D8N;

const SmallStrainUPlDiffOrderElement mSmallStrainUPlDiffOrderElement2D6N;
const SmallStrainUPlDiffOrderElement mSmallStrainUPlDiffOrderElement2D8N;
const SmallStrainUPlDiffOrderElement mSmallStrainUPlDiffOrderElement2D9N;
const SmallStrainUPlDiffOrderElement mSmallStrainUPlDiffOrderElement3D10N;
const SmallStrainUPlDiffOrderElement mSmallStrainUPlDiffOrderElement3D20N;
const SmallStrainUPlDiffOrderElement mSmallStrainUPlDiffOrderElement3D27N;

const UPlForceCondition<2,1> mUPlForceCondition2D1N;
const UPlForceCondition<3,1> mUPlForceCondition3D1N;
const UPlFaceLoadCondition<2,2> mUPlFaceLoadCondition2D2N;
const UPlFaceLoadCondition<3,3> mUPlFaceLoadCondition3D3N;
const UPlFaceLoadCondition<3,4> mUPlFaceLoadCondition3D4N;
const UPlNormalFaceLoadCondition<2,2> mUPlNormalFaceLoadCondition2D2N;
const UPlNormalFaceLoadCondition<3,3> mUPlNormalFaceLoadCondition3D3N;
const UPlNormalFaceLoadCondition<3,4> mUPlNormalFaceLoadCondition3D4N;
const UPlLiquidDischargeCondition<2,1> mUPlLiquidDischargeCondition2D1N;
const UPlLiquidDischargeCondition<3,1> mUPlLiquidDischargeCondition3D1N;
const UPlNormalLiquidFluxCondition<2,2> mUPlNormalLiquidFluxCondition2D2N;
const UPlNormalLiquidFluxCondition<3,3> mUPlNormalLiquidFluxCondition3D3N;
const UPlNormalLiquidFluxCondition<3,4> mUPlNormalLiquidFluxCondition3D4N;

const UPlFaceLoadInterfaceCondition<2,2> mUPlFaceLoadInterfaceCondition2D2N;
const UPlFaceLoadInterfaceCondition<3,4> mUPlFaceLoadInterfaceCondition3D4N;
const UPlNormalLiquidFluxInterfaceCondition<2,2> mUPlNormalLiquidFluxInterfaceCondition2D2N;
const UPlNormalLiquidFluxInterfaceCondition<3,4> mUPlNormalLiquidFluxInterfaceCondition3D4N;

const UPlNormalLiquidFluxFICCondition<2,2> mUPlNormalLiquidFluxFICCondition2D2N;
const UPlNormalLiquidFluxFICCondition<3,3> mUPlNormalLiquidFluxFICCondition3D3N;
const UPlNormalLiquidFluxFICCondition<3,4> mUPlNormalLiquidFluxFICCondition3D4N;

const LineLoad2DDiffOrderCondition mLineLoadDiffOrderCondition2D3N;
const LineNormalLoad2DDiffOrderCondition mLineNormalLoadDiffOrderCondition2D3N;
const LineNormalLiquidFlux2DDiffOrderCondition mLineNormalLiquidFluxDiffOrderCondition2D3N;
const SurfaceLoad3DDiffOrderCondition mSurfaceLoadDiffOrderCondition3D6N;
const SurfaceLoad3DDiffOrderCondition mSurfaceLoadDiffOrderCondition3D8N;
const SurfaceLoad3DDiffOrderCondition mSurfaceLoadDiffOrderCondition3D9N;
const SurfaceNormalLoad3DDiffOrderCondition mSurfaceNormalLoadDiffOrderCondition3D6N;
const SurfaceNormalLoad3DDiffOrderCondition mSurfaceNormalLoadDiffOrderCondition3D8N;
const SurfaceNormalLoad3DDiffOrderCondition mSurfaceNormalLoadDiffOrderCondition3D9N;
const SurfaceNormalLiquidFlux3DDiffOrderCondition mSurfaceNormalLiquidFluxDiffOrderCondition3D6N;
const SurfaceNormalLiquidFlux3DDiffOrderCondition mSurfaceNormalLiquidFluxDiffOrderCondition3D8N;
const SurfaceNormalLiquidFlux3DDiffOrderCondition mSurfaceNormalLiquidFluxDiffOrderCondition3D9N;

const ElastoPlasticMohrCoulombCohesive3DLaw mElastoPlasticMohrCoulombCohesive3DLaw;
const ElastoPlasticMohrCoulombCohesive2DLaw mElastoPlasticMohrCoulombCohesive2DLaw;
const ElastoPlasticModMohrCoulombCohesive3DLaw mElastoPlasticModMohrCoulombCohesive3DLaw;
const ElastoPlasticModMohrCoulombCohesive2DLaw mElastoPlasticModMohrCoulombCohesive2DLaw;
const IsotropicDamageCohesive3DLaw mIsotropicDamageCohesive3DLaw;
const IsotropicDamageCohesive2DLaw mIsotropicDamageCohesive2DLaw;
const BilinearCohesive3DLaw mBilinearCohesive3DLaw;
const BilinearCohesive2DLaw mBilinearCohesive2DLaw;
const ElasticCohesive3DLaw mElasticCohesive3DLaw;
const ElasticCohesive2DLaw mElasticCohesive2DLaw;
const ExponentialCohesive3DLaw mExponentialCohesive3DLaw;
const ExponentialCohesive2DLaw mExponentialCohesive2DLaw;

const LocalDamageFlowRule mLocalDamageFlowRule;
const NonlocalDamageFlowRule mNonlocalDamageFlowRule;

const SimoJuLocalDamage3DLaw mSimoJuLocalDamage3DLaw;
const SimoJuLocalDamagePlaneStrain2DLaw mSimoJuLocalDamagePlaneStrain2DLaw;
const SimoJuLocalDamagePlaneStress2DLaw mSimoJuLocalDamagePlaneStress2DLaw;

const SimoJuNonlocalDamage3DLaw mSimoJuNonlocalDamage3DLaw;
const SimoJuNonlocalDamagePlaneStrain2DLaw mSimoJuNonlocalDamagePlaneStrain2DLaw;
const SimoJuNonlocalDamagePlaneStress2DLaw mSimoJuNonlocalDamagePlaneStress2DLaw;

const ModifiedMisesNonlocalDamage3DLaw mModifiedMisesNonlocalDamage3DLaw;
const ModifiedMisesNonlocalDamagePlaneStrain2DLaw mModifiedMisesNonlocalDamagePlaneStrain2DLaw;
const ModifiedMisesNonlocalDamagePlaneStress2DLaw mModifiedMisesNonlocalDamagePlaneStress2DLaw;

const HistoryLinearElastic3DLaw mHistoryLinearElastic3DLaw;
const HistoryLinearElasticPlaneStrain2DLaw mHistoryLinearElasticPlaneStrain2DLaw;
const HistoryLinearElasticPlaneStress2DLaw mHistoryLinearElasticPlaneStress2DLaw;

const HyperElastic3DLaw                       mHyperElastic3DLaw;
const LinearElastic3DLaw                      mLinearElastic3DLaw;
const LinearElasticPlaneStrain2DLaw           mLinearElasticPlaneStrain2DLaw;
const LinearElasticPlaneStress2DLaw           mLinearElasticPlaneStress2DLaw;
const IsotropicDamageFlowRule                 mIsotropicDamageFlowRule;
const SimoJuYieldCriterion                    mSimoJuYieldCriterion;
const ModifiedMisesYieldCriterion             mModifiedMisesYieldCriterion;
const ExponentialDamageHardeningLaw           mExponentialDamageHardeningLaw;
const ModifiedExponentialDamageHardeningLaw   mModifiedExponentialDamageHardeningLaw;

// Assignment operator.
KratosPoromechanicsApplication& operator=(KratosPoromechanicsApplication const& rOther);

// Copy constructor.
KratosPoromechanicsApplication(KratosPoromechanicsApplication const& rOther);

}; // Class KratosPoromechanicsApplication
}  // namespace Kratos.

#endif // KRATOS_POROMECHANICS_APPLICATION_H_INCLUDED  defined


