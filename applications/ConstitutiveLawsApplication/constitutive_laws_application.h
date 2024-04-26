// KRATOS ___                _   _ _         _   _             __                       _
//       / __\___  _ __  ___| |_(_) |_ _   _| |_(_)_   _____  / /  __ ___      _____   /_\  _ __  _ __
//      / /  / _ \| '_ \/ __| __| | __| | | | __| \ \ / / _ \/ /  / _` \ \ /\ / / __| //_\\| '_ \| '_  |
//     / /__| (_) | | | \__ \ |_| | |_| |_| | |_| |\ V /  __/ /__| (_| |\ V  V /\__ \/  _  \ |_) | |_) |
//     \____/\___/|_| |_|___/\__|_|\__|\__,_|\__|_| \_/ \___\____/\__,_| \_/\_/ |___/\_/ \_/ .__/| .__/
//                                                                                         |_|   |_|
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Alejandro Cornejo Velazquez
//                   Riccardo Rossi
//

#pragma once


// System includes


// External includes


// Project includes
#include "includes/kratos_application.h"

// Constitutive laws
#include "custom_constitutive/small_strains/linear/wrinkling_linear_2d_law.h"
#include "custom_constitutive/small_strains/plasticity/truss_plasticity_constitutive_law.h"
#include "custom_constitutive/finite_strains/hyperelasticity/hyper_elastic_isotropic_ogden_1d.h"
#include "custom_constitutive/finite_strains/hyperelasticity/hyper_elastic_isotropic_henky_1d.h"
#include "custom_constitutive/small_strains/linear/elastic_isotropic_plane_stress_uncoupled_shear.h"
#include "custom_constitutive/finite_strains/hyperelasticity/hyper_elastic_isotropic_kirchhoff_3d.h"
#include "custom_constitutive/finite_strains/hyperelasticity/hyper_elastic_isotropic_kirchhoff_plane_stress_2d.h"
#include "custom_constitutive/finite_strains/hyperelasticity/hyper_elastic_isotropic_kirchhoff_plane_strain_2d.h"
#include "custom_constitutive/finite_strains/hyperelasticity/hyper_elastic_isotropic_neo_hookean_3d.h"
#include "custom_constitutive/finite_strains/hyperelasticity/hyper_elastic_isotropic_q_incomp_isoch_neo_hook_3d.h"
#include "custom_constitutive/finite_strains/hyperelasticity/hyper_elastic_isotropic_neo_hookean_plane_strain_2d.h"
#include "custom_constitutive/finite_strains/hyperelasticity/hyper_elastic_simo_taylor_neo_hookean_3d.h"
#include "custom_constitutive/finite_strains/hyperelasticity/hyper_elastic_simo_taylor_neo_hookean_plane_strain_2d.h"
#include "custom_constitutive/small_strains/linear/linear_elastic_orthotropic_2D_law.h"
#include "custom_constitutive/small_strains/plasticity/small_strain_j2_plasticity_plane_strain_2d.h"
#include "custom_constitutive/small_strains/plasticity/small_strain_j2_plasticity_3d.h"
#include "custom_constitutive/small_strains/damage/small_strain_isotropic_damage_3d.h"
#include "custom_constitutive/small_strains/damage/small_strain_isotropic_damage_implex_3d.h"
#include "custom_constitutive/small_strains/damage/small_strain_isotropic_damage_plane_strain_2d.h"
#include "custom_constitutive/small_strains/damage/small_strain_isotropic_damage_traction_only_3d.h"
#include "custom_constitutive/small_strains/damage/small_strain_isotropic_damage_traction_only_implex_3d.h"
#include "custom_constitutive/small_strains/damage/plane_stress_d_plus_d_minus_damage_masonry_2d.h"
#include "custom_constitutive/small_strains/damage/d_plus_d_minus_damage_masonry_3d.h"
#include "custom_constitutive/small_strains/plasticity/small_strain_isotropic_plasticity_factory.h"
#include "custom_constitutive/small_strains/plasticity/small_strain_kinematic_plasticity_factory.h"
#include "custom_constitutive/small_strains/plasticity/generic_small_strain_isotropic_plasticity.h"
#include "custom_constitutive/small_strains/plasticity/generic_small_strain_kinematic_plasticity.h"
#include "custom_constitutive/finite_strains/plasticity/finite_strain_isotropic_plasticity_factory.h"
#include "custom_constitutive/finite_strains/plasticity/finite_strain_kinematic_plasticity_factory.h"
#include "custom_constitutive/finite_strains/plasticity/generic_finite_strain_isotropic_plasticity.h"
#include "custom_constitutive/finite_strains/plasticity/generic_finite_strain_kinematic_plasticity.h"
#include "custom_constitutive/small_strains/damage/generic_small_strain_isotropic_damage.h"
#include "custom_constitutive/small_strains/damage/small_strain_isotropic_damage_factory.h"
#include "custom_constitutive/small_strains/viscous/viscous_generalized_kelvin.h"
#include "custom_constitutive/small_strains/viscous/viscous_generalized_maxwell.h"
#include "custom_constitutive/small_strains/viscoplasticity/generic_small_strain_viscoplasticity_3d.h"
#include "custom_constitutive/small_strains/damage/generic_small_strain_d_plus_d_minus_damage.h"
#include "custom_constitutive/small_strains/plasticity/plasticity_isotropic_kinematic_j2.h"
#include "custom_constitutive/small_strains/plastic_damage/generic_small_strain_plastic_damage_model.h"
#include "custom_constitutive/small_strains/damage/generic_small_strain_orthotropic_damage.h"
#include "custom_constitutive/composites/serial_parallel_rule_of_mixtures_law.h"
#include "custom_constitutive/small_strains/anisotropy_orthotropy/generic_anisotropic_3d_law.h"
#include "custom_constitutive/small_strains/linear/multi_linear_elastic_1d_law.h"
#include "custom_constitutive/small_strains/linear/multi_linear_isotropic_plane_stress_2d.h"
#include "custom_constitutive/small_strains/damage/generic_small_strain_isotropic_damage_plane_stress.h"

// Integrators
#include "custom_constitutive/auxiliary_files/cl_integrators/generic_cl_integrator_damage.h"
#include "custom_constitutive/auxiliary_files/cl_integrators/generic_cl_integrator_plasticity.h"
#include "custom_constitutive/auxiliary_files/cl_integrators/generic_cl_integrator_kinematic_plasticity.h"
#include "custom_constitutive/auxiliary_files/cl_integrators/d+d-cl_integrators/generic_compression_cl_integrator.h"
#include "custom_constitutive/auxiliary_files/cl_integrators/d+d-cl_integrators/generic_tension_cl_integrator.h"
#include "custom_constitutive/small_strains/fatigue/generic_small_strain_high_cycle_fatigue_law.h"

// Yield surfaces
#include "custom_constitutive/auxiliary_files/yield_surfaces/generic_yield_surface.h"
#include "custom_constitutive/auxiliary_files/yield_surfaces/von_mises_yield_surface.h"
#include "custom_constitutive/auxiliary_files/yield_surfaces/modified_mohr_coulomb_yield_surface.h"
#include "custom_constitutive/auxiliary_files/yield_surfaces/mohr_coulomb_yield_surface.h"
#include "custom_constitutive/auxiliary_files/yield_surfaces/rankine_yield_surface.h"
#include "custom_constitutive/auxiliary_files/yield_surfaces/simo_ju_yield_surface.h"
#include "custom_constitutive/auxiliary_files/yield_surfaces/drucker_prager_yield_surface.h"
#include "custom_constitutive/auxiliary_files/yield_surfaces/tresca_yield_surface.h"

// Thermal yield surfaces
#include "custom_constitutive/thermal/auxiliary_files/thermal_yield_surfaces/thermal_von_mises_yield_surface.h"
#include "custom_constitutive/thermal/auxiliary_files/thermal_yield_surfaces/thermal_tresca_yield_surface.h"
#include "custom_constitutive/thermal/auxiliary_files/thermal_yield_surfaces/thermal_rankine_yield_surface.h"
#include "custom_constitutive/thermal/auxiliary_files/thermal_yield_surfaces/thermal_simo_ju_yield_surface.h"
#include "custom_constitutive/thermal/auxiliary_files/thermal_yield_surfaces/thermal_mohr_coulomb_yield_surface.h"
#include "custom_constitutive/thermal/auxiliary_files/thermal_yield_surfaces/thermal_modified_mohr_coulomb_yield_surface.h"
#include "custom_constitutive/thermal/auxiliary_files/thermal_yield_surfaces/thermal_drucker_prager_yield_surface.h"


// Plastic potentials
#include "custom_constitutive/auxiliary_files/plastic_potentials/generic_plastic_potential.h"
#include "custom_constitutive/auxiliary_files/plastic_potentials/von_mises_plastic_potential.h"
#include "custom_constitutive/auxiliary_files/plastic_potentials/tresca_plastic_potential.h"
#include "custom_constitutive/auxiliary_files/plastic_potentials/modified_mohr_coulomb_plastic_potential.h"
#include "custom_constitutive/auxiliary_files/plastic_potentials/mohr_coulomb_plastic_potential.h"
#include "custom_constitutive/auxiliary_files/plastic_potentials/drucker_prager_plastic_potential.h"
#include "custom_constitutive/auxiliary_files/plastic_potentials/rankine_plastic_potential.h"

// Rules of mixtures
#include "custom_constitutive/composites/rule_of_mixtures_law.h"
#include "custom_constitutive/composites/traction_separation_law.h"

#include "custom_constitutive/small_strains/plastic_damage/associative_plastic_damage_model.h"

// Thermal CL
#include "custom_constitutive/thermal/small_strains/elastic/thermal_elastic_isotropic_3d.h"
#include "custom_constitutive/thermal/small_strains/elastic/thermal_linear_plane_strain.h"
#include "custom_constitutive/thermal/small_strains/elastic/thermal_linear_plane_stress.h"

#include "custom_constitutive/thermal/small_strains/damage/generic_small_strain_thermal_isotropic_damage.h"
#include "custom_constitutive/thermal/small_strains/damage/generic_small_strain_thermal_isotropic_damage_plane_stress.h"

namespace Kratos {

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
*/
class KRATOS_API(CONSTITUTIVE_LAWS_APPLICATION) KratosConstitutiveLawsApplication : public KratosApplication {
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of KratosConstitutiveLawsApplication
    KRATOS_CLASS_POINTER_DEFINITION(KratosConstitutiveLawsApplication);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    KratosConstitutiveLawsApplication();

    /// Destructor.
    ~KratosConstitutiveLawsApplication() override {}

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    void Register() override;

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "KratosConstitutiveLawsApplication";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
        PrintData(rOStream);
    }

    ///// Print object's data.
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

    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{


    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{

    // static const ApplicationCondition  msApplicationCondition;

    ///@}
    ///@name Member Variables
    ///@{


    // Damage and plasticity laws
    const ElasticIsotropicPlaneStressUncoupledShear  mElasticIsotropicPlaneStressUncoupledShear;
    const HyperElasticIsotropicKirchhoff3D  mHyperElasticIsotropicKirchhoff3D;
    const HyperElasticIsotropicKirchhoffPlaneStress2D  mHyperElasticIsotropicKirchhoffPlaneStress2D;
    const HyperElasticIsotropicKirchhoffPlaneStrain2D  mHyperElasticIsotropicKirchhoffPlaneStrain2D;
    const HyperElasticIsotropicNeoHookean3D  mHyperElasticIsotropicNeoHookean3D;
    const HyperElasticIsotropicQuasiIncompressibleIshochoricNeoHookean3D  mHyperElasticIsotropicQuasiIncompressibleIshochoricNeoHookean3D;
    const HyperElasticIsotropicNeoHookeanPlaneStrain2D  mHyperElasticIsotropicNeoHookeanPlaneStrain2D;
    const HyperElasticSimoTaylorNeoHookean3D mHyperElasticSimoTaylorNeoHookean3D;
    const HyperElasticSimoTaylorNeoHookeanPlaneStrain2D mHyperElasticSimoTaylorNeoHookeanPlaneStrain2D;
    const LinearElasticOrthotropic2DLaw mLinearElasticOrthotropic2DLaw;

    const SmallStrainJ2Plasticity3D mSmallStrainJ2Plasticity3D;
    const SmallStrainJ2PlasticityPlaneStrain2D mSmallStrainJ2PlasticityPlaneStrain2D;
    const SmallStrainIsotropicDamage3D mSmallStrainIsotropicDamage3D;
    const SmallStrainIsotropicDamageImplex3D mSmallStrainIsotropicDamageImplex3D;
    const SmallStrainIsotropicDamagePlaneStrain2D mSmallStrainIsotropicDamagePlaneStrain2D;
    const SmallStrainIsotropicDamageTractionOnly3D mSmallStrainIsotropicDamageTractionOnly3D;
    const SmallStrainIsotropicDamageTractionOnlyImplex3D mSmallStrainIsotropicDamageTractionOnlyImplex3D;
    const TrussPlasticityConstitutiveLaw mTrussPlasticityConstitutiveLaw;
    const HyperElasticIsotropicOgden1D mHyperElasticIsotropicOgden1D;
    const HyperElasticIsotropicHenky1D mHyperElasticIsotropicHenky1D;
    const WrinklingLinear2DLaw mWrinklingLinear2DLaw;
    const MultiLinearElastic1DLaw mMultiLinearElastic1DLaw;
    const MultiLinearIsotropicPlaneStress2D mMultiLinearIsotropicPlaneStress2D;

    // Damage and plasticity laws
    const SerialParallelRuleOfMixturesLaw mSerialParallelRuleOfMixturesLaw;
    const SmallStrainIsotropicPlasticityFactory mSmallStrainIsotropicPlasticityFactory;
    const SmallStrainKinematicPlasticityFactory mSmallStrainKinematicPlasticityFactory;
    const FiniteStrainIsotropicPlasticityFactory mFiniteStrainIsotropicPlasticityFactory;
    const SmallStrainIsotropicDamageFactory mSmallStrainIsotropicDamageFactory;
    const ViscousGeneralizedKelvin<ElasticIsotropic3D> mViscousGeneralizedKelvin3D;
    const ViscousGeneralizedMaxwell<ElasticIsotropic3D> mViscousGeneralizedMaxwell3D;
    const GenericSmallStrainViscoplasticity3D mGenericSmallStrainViscoplasticity3D;
    const PlasticityIsotropicKinematicJ2 mPlasticityIsotropicKinematicJ2;

    /// Plasticity
    /* Small strain */
    const GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainIsotropicPlasticityPlaneStrainVonMisesVonMises;
    const GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<ModifiedMohrCoulombPlasticPotential<3>>>> mSmallStrainIsotropicPlasticityPlaneStrainVonMisesModifiedMohrCoulomb;
    const GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<DruckerPragerPlasticPotential<3>>>> mSmallStrainIsotropicPlasticityPlaneStrainVonMisesDruckerPrager;
    const GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<TrescaPlasticPotential<3>>>> mSmallStrainIsotropicPlasticityPlaneStrainVonMisesTresca;
    const GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainIsotropicPlasticityPlaneStrainModifiedMohrCoulombVonMises;
    const GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<3>>>> mSmallStrainIsotropicPlasticityPlaneStrainModifiedMohrCoulombModifiedMohrCoulomb;
    const GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<DruckerPragerPlasticPotential<3>>>> mSmallStrainIsotropicPlasticityPlaneStrainModifiedMohrCoulombDruckerPrager;
    const GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<TrescaPlasticPotential<3>>>> mSmallStrainIsotropicPlasticityPlaneStrainModifiedMohrCoulombTresca;
    const GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainIsotropicPlasticityPlaneStrainTrescaVonMises;
    const GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<ModifiedMohrCoulombPlasticPotential<3>>>> mSmallStrainIsotropicPlasticityPlaneStrainTrescaModifiedMohrCoulomb;
    const GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<DruckerPragerPlasticPotential<3>>>> mSmallStrainIsotropicPlasticityPlaneStrainTrescaDruckerPrager;
    const GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<TrescaPlasticPotential<3>>>> mSmallStrainIsotropicPlasticityPlaneStrainTrescaTresca;
    const GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainIsotropicPlasticityPlaneStrainDruckerPragerVonMises;
    const GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<ModifiedMohrCoulombPlasticPotential<3>>>> mSmallStrainIsotropicPlasticityPlaneStrainDruckerPragerModifiedMohrCoulomb;
    const GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<3>>>> mSmallStrainIsotropicPlasticityPlaneStrainDruckerPragerDruckerPrager;
    const GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<TrescaPlasticPotential<3>>>> mSmallStrainIsotropicPlasticityPlaneStrainDruckerPragerTresca;
    const GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<MohrCoulombPlasticPotential<3>>>> mSmallStrainIsotropicPlasticityPlaneStrainVonMisesMohrCoulomb;
    const GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainIsotropicPlasticityPlaneStrainMohrCoulombVonMises;
    const GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<MohrCoulombPlasticPotential<3>>>> mSmallStrainIsotropicPlasticityPlaneStrainMohrCoulombMohrCoulomb;
    const GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<DruckerPragerPlasticPotential<3>>>> mSmallStrainIsotropicPlasticityPlaneStrainMohrCoulombDruckerPrager;
    const GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<TrescaPlasticPotential<3>>>> mSmallStrainIsotropicPlasticityPlaneStrainMohrCoulombTresca;
    const GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<MohrCoulombPlasticPotential<3>>>> mSmallStrainIsotropicPlasticityPlaneStrainTrescaMohrCoulomb;
    const GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<MohrCoulombPlasticPotential<3>>>> mSmallStrainIsotropicPlasticityPlaneStrainDruckerPragerMohrCoulomb;

    const GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainIsotropicPlasticity3DVonMisesVonMises;
    const GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>> mSmallStrainIsotropicPlasticity3DVonMisesModifiedMohrCoulomb;
    const GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<DruckerPragerPlasticPotential<6>>>> mSmallStrainIsotropicPlasticity3DVonMisesDruckerPrager;
    const GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<TrescaPlasticPotential<6>>>> mSmallStrainIsotropicPlasticity3DVonMisesTresca;
    const GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainIsotropicPlasticity3DModifiedMohrCoulombVonMises;
    const GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>> mSmallStrainIsotropicPlasticity3DModifiedMohrCoulombModifiedMohrCoulomb;
    const GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<DruckerPragerPlasticPotential<6>>>> mSmallStrainIsotropicPlasticity3DModifiedMohrCoulombDruckerPrager;
    const GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<TrescaPlasticPotential<6>>>> mSmallStrainIsotropicPlasticity3DModifiedMohrCoulombTresca;
    const GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainIsotropicPlasticity3DTrescaVonMises;
    const GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>> mSmallStrainIsotropicPlasticity3DTrescaModifiedMohrCoulomb;
    const GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<DruckerPragerPlasticPotential<6>>>> mSmallStrainIsotropicPlasticity3DTrescaDruckerPrager;
    const GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<TrescaPlasticPotential<6>>>> mSmallStrainIsotropicPlasticity3DTrescaTresca;
    const GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainIsotropicPlasticity3DDruckerPragerVonMises;
    const GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>> mSmallStrainIsotropicPlasticity3DDruckerPragerModifiedMohrCoulomb;
    const GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<6>>>> mSmallStrainIsotropicPlasticity3DDruckerPragerDruckerPrager;
    const GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<TrescaPlasticPotential<6>>>> mSmallStrainIsotropicPlasticity3DDruckerPragerTresca;
    const GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<MohrCoulombPlasticPotential<6>>>> mSmallStrainIsotropicPlasticity3DVonMisesMohrCoulomb;
    const GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainIsotropicPlasticity3DMohrCoulombVonMises;
    const GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<MohrCoulombPlasticPotential<6>>>> mSmallStrainIsotropicPlasticity3DMohrCoulombMohrCoulomb;
    const GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<DruckerPragerPlasticPotential<6>>>> mSmallStrainIsotropicPlasticity3DMohrCoulombDruckerPrager;
    const GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<TrescaPlasticPotential<6>>>> mSmallStrainIsotropicPlasticity3DMohrCoulombTresca;
    const GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<MohrCoulombPlasticPotential<6>>>> mSmallStrainIsotropicPlasticity3DTrescaMohrCoulomb;
    const GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<MohrCoulombPlasticPotential<6>>>> mSmallStrainIsotropicPlasticity3DDruckerPragerMohrCoulomb;

    const GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainKinematicPlasticityPlaneStrainVonMisesVonMises;
    const GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<ModifiedMohrCoulombPlasticPotential<3>>>> mSmallStrainKinematicPlasticityPlaneStrainVonMisesModifiedMohrCoulomb;
    const GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<DruckerPragerPlasticPotential<3>>>> mSmallStrainKinematicPlasticityPlaneStrainVonMisesDruckerPrager;
    const GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<TrescaPlasticPotential<3>>>> mSmallStrainKinematicPlasticityPlaneStrainVonMisesTresca;
    const GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainKinematicPlasticityPlaneStrainModifiedMohrCoulombVonMises;
    const GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<3>>>> mSmallStrainKinematicPlasticityPlaneStrainModifiedMohrCoulombModifiedMohrCoulomb;
    const GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<ModifiedMohrCoulombYieldSurface<DruckerPragerPlasticPotential<3>>>> mSmallStrainKinematicPlasticityPlaneStrainModifiedMohrCoulombDruckerPrager;
    const GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<ModifiedMohrCoulombYieldSurface<TrescaPlasticPotential<3>>>> mSmallStrainKinematicPlasticityPlaneStrainModifiedMohrCoulombTresca;
    const GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainKinematicPlasticityPlaneStrainTrescaVonMises;
    const GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<ModifiedMohrCoulombPlasticPotential<3>>>> mSmallStrainKinematicPlasticityPlaneStrainTrescaModifiedMohrCoulomb;
    const GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<DruckerPragerPlasticPotential<3>>>> mSmallStrainKinematicPlasticityPlaneStrainTrescaDruckerPrager;
    const GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<TrescaPlasticPotential<3>>>> mSmallStrainKinematicPlasticityPlaneStrainTrescaTresca;
    const GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainKinematicPlasticityPlaneStrainDruckerPragerVonMises;
    const GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<ModifiedMohrCoulombPlasticPotential<3>>>> mSmallStrainKinematicPlasticityPlaneStrainDruckerPragerModifiedMohrCoulomb;
    const GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<3>>>> mSmallStrainKinematicPlasticityPlaneStrainDruckerPragerDruckerPrager;
    const GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<TrescaPlasticPotential<3>>>> mSmallStrainKinematicPlasticityPlaneStrainDruckerPragerTresca;
    const GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<MohrCoulombPlasticPotential<3>>>> mSmallStrainKinematicPlasticityPlaneStrainVonMisesMohrCoulomb;
    const GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<MohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainKinematicPlasticityPlaneStrainMohrCoulombVonMises;
    const GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<MohrCoulombYieldSurface<MohrCoulombPlasticPotential<3>>>> mSmallStrainKinematicPlasticityPlaneStrainMohrCoulombMohrCoulomb;
    const GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<MohrCoulombYieldSurface<DruckerPragerPlasticPotential<3>>>> mSmallStrainKinematicPlasticityPlaneStrainMohrCoulombDruckerPrager;
    const GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<MohrCoulombYieldSurface<TrescaPlasticPotential<3>>>> mSmallStrainKinematicPlasticityPlaneStrainMohrCoulombTresca;
    const GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<MohrCoulombPlasticPotential<3>>>> mSmallStrainKinematicPlasticityPlaneStrainTrescaMohrCoulomb;
    const GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<MohrCoulombPlasticPotential<3>>>> mSmallStrainKinematicPlasticityPlaneStrainDruckerPragerMohrCoulomb;

    const GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainKinematicPlasticity3DVonMisesVonMises;
    const GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>> mSmallStrainKinematicPlasticity3DVonMisesModifiedMohrCoulomb;
    const GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<DruckerPragerPlasticPotential<6>>>> mSmallStrainKinematicPlasticity3DVonMisesDruckerPrager;
    const GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<TrescaPlasticPotential<6>>>> mSmallStrainKinematicPlasticity3DVonMisesTresca;
    const GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainKinematicPlasticity3DModifiedMohrCoulombVonMises;
    const GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>> mSmallStrainKinematicPlasticity3DModifiedMohrCoulombModifiedMohrCoulomb;
    const GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<ModifiedMohrCoulombYieldSurface<DruckerPragerPlasticPotential<6>>>> mSmallStrainKinematicPlasticity3DModifiedMohrCoulombDruckerPrager;
    const GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<ModifiedMohrCoulombYieldSurface<TrescaPlasticPotential<6>>>> mSmallStrainKinematicPlasticity3DModifiedMohrCoulombTresca;
    const GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainKinematicPlasticity3DTrescaVonMises;
    const GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>> mSmallStrainKinematicPlasticity3DTrescaModifiedMohrCoulomb;
    const GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<DruckerPragerPlasticPotential<6>>>> mSmallStrainKinematicPlasticity3DTrescaDruckerPrager;
    const GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<TrescaPlasticPotential<6>>>> mSmallStrainKinematicPlasticity3DTrescaTresca;
    const GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainKinematicPlasticity3DDruckerPragerVonMises;
    const GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>> mSmallStrainKinematicPlasticity3DDruckerPragerModifiedMohrCoulomb;
    const GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<6>>>> mSmallStrainKinematicPlasticity3DDruckerPragerDruckerPrager;
    const GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<TrescaPlasticPotential<6>>>> mSmallStrainKinematicPlasticity3DDruckerPragerTresca;
    const GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<MohrCoulombPlasticPotential<6>>>> mSmallStrainKinematicPlasticity3DVonMisesMohrCoulomb;
    const GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<MohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainKinematicPlasticity3DMohrCoulombVonMises;
    const GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<MohrCoulombYieldSurface<MohrCoulombPlasticPotential<6>>>> mSmallStrainKinematicPlasticity3DMohrCoulombMohrCoulomb;
    const GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<MohrCoulombYieldSurface<DruckerPragerPlasticPotential<6>>>> mSmallStrainKinematicPlasticity3DMohrCoulombDruckerPrager;
    const GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<MohrCoulombYieldSurface<TrescaPlasticPotential<6>>>> mSmallStrainKinematicPlasticity3DMohrCoulombTresca;
    const GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<MohrCoulombPlasticPotential<6>>>> mSmallStrainKinematicPlasticity3DTrescaMohrCoulomb;
    const GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<MohrCoulombPlasticPotential<6>>>> mSmallStrainKinematicPlasticity3DDruckerPragerMohrCoulomb;

    // Plastic Damage Model
    const GenericSmallStrainPlasticDamageModel <GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>, GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainPlasticDamageModel3DVonMisesVonMisesVonMises;
    const GenericSmallStrainPlasticDamageModel <GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>, GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainPlasticDamageModel3DVonMisesVonMisesDruckerPrager;


    /* Finite strain Isotropic plasticity*/
    const GenericFiniteStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>> mFiniteStrainIsotropicPlasticity3DVonMisesVonMises;
    const GenericFiniteStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>> mFiniteStrainIsotropicPlasticity3DVonMisesModifiedMohrCoulomb;
    const GenericFiniteStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<DruckerPragerPlasticPotential<6>>>> mFiniteStrainIsotropicPlasticity3DVonMisesDruckerPrager;
    const GenericFiniteStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<TrescaPlasticPotential<6>>>> mFiniteStrainIsotropicPlasticity3DVonMisesTresca;
    const GenericFiniteStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>> mFiniteStrainIsotropicPlasticity3DModifiedMohrCoulombVonMises;
    const GenericFiniteStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>> mFiniteStrainIsotropicPlasticity3DModifiedMohrCoulombModifiedMohrCoulomb;
    const GenericFiniteStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<DruckerPragerPlasticPotential<6>>>> mFiniteStrainIsotropicPlasticity3DModifiedMohrCoulombDruckerPrager;
    const GenericFiniteStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<TrescaPlasticPotential<6>>>> mFiniteStrainIsotropicPlasticity3DModifiedMohrCoulombTresca;
    const GenericFiniteStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<VonMisesPlasticPotential<6>>>> mFiniteStrainIsotropicPlasticity3DTrescaVonMises;
    const GenericFiniteStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>> mFiniteStrainIsotropicPlasticity3DTrescaModifiedMohrCoulomb;
    const GenericFiniteStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<DruckerPragerPlasticPotential<6>>>> mFiniteStrainIsotropicPlasticity3DTrescaDruckerPrager;
    const GenericFiniteStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<TrescaPlasticPotential<6>>>> mFiniteStrainIsotropicPlasticity3DTrescaTresca;
    const GenericFiniteStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>> mFiniteStrainIsotropicPlasticity3DDruckerPragerVonMises;
    const GenericFiniteStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>> mFiniteStrainIsotropicPlasticity3DDruckerPragerModifiedMohrCoulomb;
    const GenericFiniteStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<6>>>> mFiniteStrainIsotropicPlasticity3DDruckerPragerDruckerPrager;
    const GenericFiniteStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<TrescaPlasticPotential<6>>>> mFiniteStrainIsotropicPlasticity3DDruckerPragerTresca;
    const GenericFiniteStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<MohrCoulombPlasticPotential<6>>>> mFiniteStrainIsotropicPlasticity3DVonMisesMohrCoulomb;
    const GenericFiniteStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>> mFiniteStrainIsotropicPlasticity3DMohrCoulombVonMises;
    const GenericFiniteStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<MohrCoulombPlasticPotential<6>>>> mFiniteStrainIsotropicPlasticity3DMohrCoulombMohrCoulomb;
    const GenericFiniteStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<DruckerPragerPlasticPotential<6>>>> mFiniteStrainIsotropicPlasticity3DMohrCoulombDruckerPrager;
    const GenericFiniteStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<TrescaPlasticPotential<6>>>> mFiniteStrainIsotropicPlasticity3DMohrCoulombTresca;
    const GenericFiniteStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<MohrCoulombPlasticPotential<6>>>> mFiniteStrainIsotropicPlasticity3DTrescaMohrCoulomb;
    const GenericFiniteStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<MohrCoulombPlasticPotential<6>>>> mFiniteStrainIsotropicPlasticity3DDruckerPragerMohrCoulomb;

    /* Finite strain Kinematic plasticity*/
    const GenericFiniteStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>> mFiniteStrainKinematicPlasticity3DVonMisesVonMises;
    const GenericFiniteStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>> mFiniteStrainKinematicPlasticity3DVonMisesModifiedMohrCoulomb;
    const GenericFiniteStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<DruckerPragerPlasticPotential<6>>>> mFiniteStrainKinematicPlasticity3DVonMisesDruckerPrager;
    const GenericFiniteStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<TrescaPlasticPotential<6>>>> mFiniteStrainKinematicPlasticity3DVonMisesTresca;
    const GenericFiniteStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>> mFiniteStrainKinematicPlasticity3DModifiedMohrCoulombVonMises;
    const GenericFiniteStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>> mFiniteStrainKinematicPlasticity3DModifiedMohrCoulombModifiedMohrCoulomb;
    const GenericFiniteStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<ModifiedMohrCoulombYieldSurface<DruckerPragerPlasticPotential<6>>>> mFiniteStrainKinematicPlasticity3DModifiedMohrCoulombDruckerPrager;
    const GenericFiniteStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<ModifiedMohrCoulombYieldSurface<TrescaPlasticPotential<6>>>> mFiniteStrainKinematicPlasticity3DModifiedMohrCoulombTresca;
    const GenericFiniteStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<VonMisesPlasticPotential<6>>>> mFiniteStrainKinematicPlasticity3DTrescaVonMises;
    const GenericFiniteStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>> mFiniteStrainKinematicPlasticity3DTrescaModifiedMohrCoulomb;
    const GenericFiniteStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<DruckerPragerPlasticPotential<6>>>> mFiniteStrainKinematicPlasticity3DTrescaDruckerPrager;
    const GenericFiniteStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<TrescaPlasticPotential<6>>>> mFiniteStrainKinematicPlasticity3DTrescaTresca;
    const GenericFiniteStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>> mFiniteStrainKinematicPlasticity3DDruckerPragerVonMises;
    const GenericFiniteStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>> mFiniteStrainKinematicPlasticity3DDruckerPragerModifiedMohrCoulomb;
    const GenericFiniteStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<6>>>> mFiniteStrainKinematicPlasticity3DDruckerPragerDruckerPrager;
    const GenericFiniteStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<TrescaPlasticPotential<6>>>> mFiniteStrainKinematicPlasticity3DDruckerPragerTresca;
    const GenericFiniteStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<MohrCoulombPlasticPotential<6>>>> mFiniteStrainKinematicPlasticity3DVonMisesMohrCoulomb;
    const GenericFiniteStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<MohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>> mFiniteStrainKinematicPlasticity3DMohrCoulombVonMises;
    const GenericFiniteStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<MohrCoulombYieldSurface<MohrCoulombPlasticPotential<6>>>> mFiniteStrainKinematicPlasticity3DMohrCoulombMohrCoulomb;
    const GenericFiniteStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<MohrCoulombYieldSurface<DruckerPragerPlasticPotential<6>>>> mFiniteStrainKinematicPlasticity3DMohrCoulombDruckerPrager;
    const GenericFiniteStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<MohrCoulombYieldSurface<TrescaPlasticPotential<6>>>> mFiniteStrainKinematicPlasticity3DMohrCoulombTresca;
    const GenericFiniteStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<MohrCoulombPlasticPotential<6>>>> mFiniteStrainKinematicPlasticity3DTrescaMohrCoulomb;
    const GenericFiniteStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<MohrCoulombPlasticPotential<6>>>> mFiniteStrainKinematicPlasticity3DDruckerPragerMohrCoulomb;

    /// Damage
    /* Small strain 3D */
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainIsotropicDamage3DVonMises;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>> mSmallStrainIsotropicDamage3DModifiedMohrCoulomb;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<TrescaPlasticPotential<6>>>> mSmallStrainIsotropicDamage3DTresca;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<6>>>> mSmallStrainIsotropicDamage3DDruckerPrager;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<RankinePlasticPotential<6>>>> mSmallStrainIsotropicDamage3DRankine;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainIsotropicDamage3DSimoJu;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<MohrCoulombYieldSurface<MohrCoulombPlasticPotential<6>>>> mSmallStrainIsotropicDamage3DMohrCoulomb;

    /* Small strain 2D */
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainIsotropicDamagePlaneStrainVonMises;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<3>>>> mSmallStrainIsotropicDamagePlaneStrainModifiedMohrCoulomb;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<TrescaPlasticPotential<3>>>> mSmallStrainIsotropicDamagePlaneStrainTresca;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<3>>>> mSmallStrainIsotropicDamagePlaneStrainDruckerPrager;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<RankinePlasticPotential<3>>>> mSmallStrainIsotropicDamagePlaneStrainRankine;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainIsotropicDamagePlaneStrainSimoJu;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<MohrCoulombYieldSurface<MohrCoulombPlasticPotential<3>>>> mSmallStrainIsotropicDamagePlaneStrainMohrCoulomb;

    /* Small strain  plane stress 2D */
    const GenericSmallStrainIsotropicDamagePlaneStress <GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainIsotropicDamagePlaneStressVonMises;
    const GenericSmallStrainIsotropicDamagePlaneStress <GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<3>>>> mSmallStrainIsotropicDamagePlaneStressModifiedMohrCoulomb;
    const GenericSmallStrainIsotropicDamagePlaneStress <GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<TrescaPlasticPotential<3>>>> mSmallStrainIsotropicDamagePlaneStressTresca;
    const GenericSmallStrainIsotropicDamagePlaneStress <GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<3>>>> mSmallStrainIsotropicDamagePlaneStressDruckerPrager;
    const GenericSmallStrainIsotropicDamagePlaneStress <GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<RankinePlasticPotential<3>>>> mSmallStrainIsotropicDamagePlaneStressRankine;
    const GenericSmallStrainIsotropicDamagePlaneStress <GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainIsotropicDamagePlaneStressSimoJu;
    const GenericSmallStrainIsotropicDamagePlaneStress <GenericConstitutiveLawIntegratorDamage<MohrCoulombYieldSurface<MohrCoulombPlasticPotential<3>>>> mSmallStrainIsotropicDamagePlaneStressMohrCoulomb;

    // HCF (High Cycle Fatigue)
    const GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainHighCycleFatigue3DLawVonMises;
    const GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>> mSmallStrainHighCycleFatigue3DLawModifiedMohrCoulomb;
    const GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<TrescaPlasticPotential<6>>>> mSmallStrainHighCycleFatigue3DLawTresca;
    const GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<6>>>> mSmallStrainHighCycleFatigue3DLawDruckerPrager;
    const GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<RankinePlasticPotential<6>>>> mSmallStrainHighCycleFatigue3DLawRankine;
    const GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainHighCycleFatigue3DLawSimoJu;
    const GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<MohrCoulombYieldSurface<MohrCoulombPlasticPotential<6>>>> mSmallStrainHighCycleFatigue3DLawMohrCoulomb;


    // d+d- laws (3D)
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainDplusDminusDamageModifiedMohrCoulombModifiedMohrCoulomb3D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<RankineYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainDplusDminusDamageModifiedMohrCoulombRankine3D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainDplusDminusDamageModifiedMohrCoulombSimoJu3D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainDplusDminusDamageModifiedMohrCoulombVonMises3D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<TrescaYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainDplusDminusDamageModifiedMohrCoulombTresca3D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainDplusDminusDamageModifiedMohrCoulombDruckerPrager3D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<RankineYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainDplusDminusDamageRankineModifiedMohrCoulomb3D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<RankineYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<RankineYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainDplusDminusDamageRankineRankine3D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<RankineYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainDplusDminusDamageRankineSimoJu3D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<RankineYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainDplusDminusDamageRankineVonMises3D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<RankineYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<TrescaYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainDplusDminusDamageRankineTresca3D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<RankineYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainDplusDminusDamageRankineDruckerPrager3D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainDplusDminusDamageSimoJuModifiedMohrCoulomb3D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<RankineYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainDplusDminusDamageSimoJuRankine3D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainDplusDminusDamageSimoJuSimoJu3D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainDplusDminusDamageSimoJuVonMises3D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<TrescaYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainDplusDminusDamageSimoJuTresca3D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainDplusDminusDamageSimoJuDruckerPrager3D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainDplusDminusDamageVonMisesModifiedMohrCoulomb3D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<RankineYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainDplusDminusDamageVonMisesRankine3D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainDplusDminusDamageVonMisesSimoJu3D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainDplusDminusDamageVonMisesVonMises3D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<TrescaYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainDplusDminusDamageVonMisesTresca3D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainDplusDminusDamageVonMisesDruckerPrager3D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<TrescaYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainDplusDminusDamageTrescaModifiedMohrCoulomb3D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<TrescaYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<RankineYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainDplusDminusDamageTrescaRankine3D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<TrescaYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainDplusDminusDamageTrescaSimoJu3D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<TrescaYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainDplusDminusDamageTrescaVonMises3D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<TrescaYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<TrescaYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainDplusDminusDamageTrescaTresca3D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<TrescaYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainDplusDminusDamageTrescaDruckerPrager3D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainDplusDminusDamageDruckerPragerModifiedMohrCoulomb3D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<RankineYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainDplusDminusDamageDruckerPragerRankine3D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainDplusDminusDamageDruckerPragerSimoJu3D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainDplusDminusDamageDruckerPragerVonMises3D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<TrescaYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainDplusDminusDamageDruckerPragerTresca3D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainDplusDminusDamageDruckerPragerDruckerPrager3D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainDplusDminusDamageMohrCoulombMohrCoulomb3D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<RankineYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainDplusDminusDamageMohrCoulombRankine3D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainDplusDminusDamageMohrCoulombSimoJu3D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainDplusDminusDamageMohrCoulombVonMises3D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<TrescaYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainDplusDminusDamageMohrCoulombTresca3D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainDplusDminusDamageMohrCoulombDruckerPrager3D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<RankineYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainDplusDminusDamageRankineMohrCoulomb3D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainDplusDminusDamageSimoJuMohrCoulomb3D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainDplusDminusDamageVonMisesMohrCoulomb3D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<TrescaYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainDplusDminusDamageTrescaMohrCoulomb3D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainDplusDminusDamageDruckerPragerMohrCoulomb3D;

    // d+d- laws (2D)
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainDplusDminusDamageModifiedMohrCoulombModifiedMohrCoulomb2D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<RankineYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainDplusDminusDamageModifiedMohrCoulombRankine2D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<SimoJuYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainDplusDminusDamageModifiedMohrCoulombSimoJu2D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<VonMisesYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainDplusDminusDamageModifiedMohrCoulombVonMises2D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<TrescaYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainDplusDminusDamageModifiedMohrCoulombTresca2D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainDplusDminusDamageModifiedMohrCoulombDruckerPrager2D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<RankineYieldSurface<VonMisesPlasticPotential<3>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainDplusDminusDamageRankineModifiedMohrCoulomb2D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<RankineYieldSurface<VonMisesPlasticPotential<3>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<RankineYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainDplusDminusDamageRankineRankine2D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<RankineYieldSurface<VonMisesPlasticPotential<3>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<SimoJuYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainDplusDminusDamageRankineSimoJu2D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<RankineYieldSurface<VonMisesPlasticPotential<3>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<VonMisesYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainDplusDminusDamageRankineVonMises2D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<RankineYieldSurface<VonMisesPlasticPotential<3>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<TrescaYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainDplusDminusDamageRankineTresca2D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<RankineYieldSurface<VonMisesPlasticPotential<3>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainDplusDminusDamageRankineDruckerPrager2D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<SimoJuYieldSurface<VonMisesPlasticPotential<3>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainDplusDminusDamageSimoJuModifiedMohrCoulomb2D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<SimoJuYieldSurface<VonMisesPlasticPotential<3>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<RankineYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainDplusDminusDamageSimoJuRankine2D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<SimoJuYieldSurface<VonMisesPlasticPotential<3>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<SimoJuYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainDplusDminusDamageSimoJuSimoJu2D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<SimoJuYieldSurface<VonMisesPlasticPotential<3>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<VonMisesYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainDplusDminusDamageSimoJuVonMises2D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<SimoJuYieldSurface<VonMisesPlasticPotential<3>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<TrescaYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainDplusDminusDamageSimoJuTresca2D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<SimoJuYieldSurface<VonMisesPlasticPotential<3>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainDplusDminusDamageSimoJuDruckerPrager2D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<VonMisesYieldSurface<VonMisesPlasticPotential<3>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainDplusDminusDamageVonMisesModifiedMohrCoulomb2D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<VonMisesYieldSurface<VonMisesPlasticPotential<3>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<RankineYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainDplusDminusDamageVonMisesRankine2D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<VonMisesYieldSurface<VonMisesPlasticPotential<3>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<SimoJuYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainDplusDminusDamageVonMisesSimoJu2D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<VonMisesYieldSurface<VonMisesPlasticPotential<3>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<VonMisesYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainDplusDminusDamageVonMisesVonMises2D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<VonMisesYieldSurface<VonMisesPlasticPotential<3>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<TrescaYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainDplusDminusDamageVonMisesTresca2D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<VonMisesYieldSurface<VonMisesPlasticPotential<3>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainDplusDminusDamageVonMisesDruckerPrager2D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<TrescaYieldSurface<VonMisesPlasticPotential<3>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainDplusDminusDamageTrescaModifiedMohrCoulomb2D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<TrescaYieldSurface<VonMisesPlasticPotential<3>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<RankineYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainDplusDminusDamageTrescaRankine2D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<TrescaYieldSurface<VonMisesPlasticPotential<3>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<SimoJuYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainDplusDminusDamageTrescaSimoJu2D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<TrescaYieldSurface<VonMisesPlasticPotential<3>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<VonMisesYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainDplusDminusDamageTrescaVonMises2D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<TrescaYieldSurface<VonMisesPlasticPotential<3>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<TrescaYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainDplusDminusDamageTrescaTresca2D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<TrescaYieldSurface<VonMisesPlasticPotential<3>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainDplusDminusDamageTrescaDruckerPrager2D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<3>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainDplusDminusDamageDruckerPragerModifiedMohrCoulomb2D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<3>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<RankineYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainDplusDminusDamageDruckerPragerRankine2D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<3>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<SimoJuYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainDplusDminusDamageDruckerPragerSimoJu2D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<3>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<VonMisesYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainDplusDminusDamageDruckerPragerVonMises2D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<3>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<TrescaYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainDplusDminusDamageDruckerPragerTresca2D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<3>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainDplusDminusDamageDruckerPragerDruckerPrager2D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainDplusDminusDamageMohrCoulombMohrCoulomb2D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<RankineYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainDplusDminusDamageMohrCoulombRankine2D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<SimoJuYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainDplusDminusDamageMohrCoulombSimoJu2D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<VonMisesYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainDplusDminusDamageMohrCoulombVonMises2D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<TrescaYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainDplusDminusDamageMohrCoulombTresca2D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainDplusDminusDamageMohrCoulombDruckerPrager2D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<RankineYieldSurface<VonMisesPlasticPotential<3>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainDplusDminusDamageRankineMohrCoulomb2D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<SimoJuYieldSurface<VonMisesPlasticPotential<3>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainDplusDminusDamageSimoJuMohrCoulomb2D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<VonMisesYieldSurface<VonMisesPlasticPotential<3>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainDplusDminusDamageVonMisesMohrCoulomb2D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<TrescaYieldSurface<VonMisesPlasticPotential<3>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainDplusDminusDamageTrescaMohrCoulomb2D;
    const GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<3>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainDplusDminusDamageDruckerPragerMohrCoulomb2D;
    const DamageDPlusDMinusMasonry2DLaw mDamageDPlusDMinusPlaneStressMasonry2DLaw;
    const DamageDPlusDMinusMasonry3DLaw mDamageDPlusDMinusMasonry3DLaw;

    // Orthotropic Damage
    const GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainOrthotropicDamageRankine3D;
    const GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainOrthotropicDamageVonMises3D;
    const GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainOrthotropicDamageDruckerPrager3D;
    const GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainOrthotropicDamageTresca3D;
    const GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainOrthotropicDamageMohrCoulomb3D;
    const GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainOrthotropicDamageModifiedMohrCoulomb3D;
    const GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainOrthotropicDamageSimoJu3D;
    const GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainOrthotropicDamageRankine2D;
    const GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainOrthotropicDamageVonMises2D;
    const GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainOrthotropicDamageDruckerPrager2D;
    const GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainOrthotropicDamageTresca2D;
    const GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainOrthotropicDamageMohrCoulomb2D;
    const GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainOrthotropicDamageModifiedMohrCoulomb2D;
    const GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainOrthotropicDamageSimoJu2D;

    // Rules of mixtures
    const ParallelRuleOfMixturesLaw<3> mParallelRuleOfMixturesLaw3D;
    const TractionSeparationLaw3D<3> mTractionSeparationLaw3D;
	const ParallelRuleOfMixturesLaw<2> mParallelRuleOfMixturesLaw2D;

    // Anisotropic law
    const GenericAnisotropic3DLaw mGenericAnisotropic3DLaw;

    const AssociativePlasticDamageModel <VonMisesYieldSurface<VonMisesPlasticPotential<6>>> mAssociativePlasticDamageModel3DVonMises;
    const AssociativePlasticDamageModel <DruckerPragerYieldSurface<DruckerPragerPlasticPotential<6>>> mAssociativePlasticDamageModel3DDruckerPrager;
    const AssociativePlasticDamageModel <ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>> mAssociativePlasticDamageModel3DModifiedMohrCoulomb;
    const AssociativePlasticDamageModel <RankineYieldSurface<RankinePlasticPotential<6>>> mAssociativePlasticDamageModel3DRankine;

    // Thermal CL
    const ThermalElasticIsotropic3D mThermalElasticIsotropic3D;
    const ThermalLinearPlaneStrain mThermalLinearPlaneStrain;
    const ThermalLinearPlaneStress mThermalLinearPlaneStress;

    /// Thermal Damage
    /* Small strain 3D */
    const GenericSmallStrainThermalIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ThermalVonMisesYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainThermalIsotropicDamage3DVonMises;
    const GenericSmallStrainThermalIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ThermalModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainThermalIsotropicDamage3DModifiedMohrCoulomb;
    const GenericSmallStrainThermalIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ThermalTrescaYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainThermalIsotropicDamage3DTresca;
    const GenericSmallStrainThermalIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ThermalDruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainThermalIsotropicDamage3DDruckerPrager;
    const GenericSmallStrainThermalIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ThermalRankineYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainThermalIsotropicDamage3DRankine;
    const GenericSmallStrainThermalIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ThermalSimoJuYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainThermalIsotropicDamage3DSimoJu;
    const GenericSmallStrainThermalIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ThermalMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainThermalIsotropicDamage3DMohrCoulomb;

    /* Small strain 2D */
    const GenericSmallStrainThermalIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ThermalVonMisesYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainThermalIsotropicDamagePlaneStrainVonMises;
    const GenericSmallStrainThermalIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ThermalModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainThermalIsotropicDamagePlaneStrainModifiedMohrCoulomb;
    const GenericSmallStrainThermalIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ThermalTrescaYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainThermalIsotropicDamagePlaneStrainTresca;
    const GenericSmallStrainThermalIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ThermalDruckerPragerYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainThermalIsotropicDamagePlaneStrainDruckerPrager;
    const GenericSmallStrainThermalIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ThermalRankineYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainThermalIsotropicDamagePlaneStrainRankine;
    const GenericSmallStrainThermalIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ThermalSimoJuYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainThermalIsotropicDamagePlaneStrainSimoJu;
    const GenericSmallStrainThermalIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ThermalMohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainThermalIsotropicDamagePlaneStrainMohrCoulomb;

    // Thermal plane stress damage
    const GenericSmallStrainThermalIsotropicDamagePlaneStress <GenericConstitutiveLawIntegratorDamage<ThermalVonMisesYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainThermalIsotropicDamagePlaneStressVonMises;
    const GenericSmallStrainThermalIsotropicDamagePlaneStress <GenericConstitutiveLawIntegratorDamage<ThermalModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainThermalIsotropicDamagePlaneStressModifiedMohrCoulomb;
    const GenericSmallStrainThermalIsotropicDamagePlaneStress <GenericConstitutiveLawIntegratorDamage<ThermalTrescaYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainThermalIsotropicDamagePlaneStressTresca;
    const GenericSmallStrainThermalIsotropicDamagePlaneStress <GenericConstitutiveLawIntegratorDamage<ThermalDruckerPragerYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainThermalIsotropicDamagePlaneStressDruckerPrager;
    const GenericSmallStrainThermalIsotropicDamagePlaneStress <GenericConstitutiveLawIntegratorDamage<ThermalRankineYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainThermalIsotropicDamagePlaneStressRankine;
    const GenericSmallStrainThermalIsotropicDamagePlaneStress <GenericConstitutiveLawIntegratorDamage<ThermalSimoJuYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainThermalIsotropicDamagePlaneStressSimoJu;
    const GenericSmallStrainThermalIsotropicDamagePlaneStress <GenericConstitutiveLawIntegratorDamage<ThermalMohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainThermalIsotropicDamagePlaneStressMohrCoulomb;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    KratosConstitutiveLawsApplication& operator=(KratosConstitutiveLawsApplication const& rOther);

    /// Copy constructor.
    KratosConstitutiveLawsApplication(KratosConstitutiveLawsApplication const& rOther);


    ///@}

}; // Class KratosConstitutiveLawsApplication

///@}


///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.
