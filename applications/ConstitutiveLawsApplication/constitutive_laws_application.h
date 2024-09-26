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

#if !defined(KRATOS_CONSTITUTIVE_LAWS_APPLICATION_H_INCLUDED )
#define  KRATOS_CONSTITUTIVE_LAWS_APPLICATION_H_INCLUDED


// System includes


// External includes


// Project includes
#include "includes/kratos_application.h"

// Constitutive laws
#include "custom_constitutive/wrinkling_linear_2d_law.h"
#include "custom_constitutive/truss_plasticity_constitutive_law.h"
#include "custom_constitutive/hyper_elastic_isotropic_ogden_1d.h"
#include "custom_constitutive/hyper_elastic_isotropic_henky_1d.h"
#include "custom_constitutive/elastic_isotropic_plane_stress_uncoupled_shear.h"
#include "custom_constitutive/hyper_elastic_isotropic_kirchhoff_3d.h"
#include "custom_constitutive/hyper_elastic_isotropic_kirchhoff_plane_stress_2d.h"
#include "custom_constitutive/hyper_elastic_isotropic_kirchhoff_plane_strain_2d.h"
#include "custom_constitutive/hyper_elastic_isotropic_neo_hookean_3d.h"
#include "custom_constitutive/hyper_elastic_isotropic_neo_hookean_plane_strain_2d.h"
#include "custom_constitutive/linear_elastic_orthotropic_2D_law.h"
#include "custom_constitutive/small_strain_j2_plasticity_plane_strain_2d.h"
#include "custom_constitutive/small_strain_j2_plasticity_3d.h"
#include "custom_constitutive/small_strain_plasticity_3d_ath.h"
#include "custom_constitutive/small_strain_isotropic_damage_3d.h"
#include "custom_constitutive/small_strain_isotropic_damage_3d_athira.h"
#include "custom_constitutive/small_strain_isotropic_damage_3d_mazars.h"
#include "custom_constitutive/small_strain_isotropic_damage_implex_3d.h"
#include "custom_constitutive/small_strain_isotropic_damage_plane_strain_2d.h"
#include "custom_constitutive/small_strain_isotropic_damage_traction_only_3d.h"
#include "custom_constitutive/small_strain_isotropic_damage_traction_only_implex_3d.h"
#include "custom_constitutive/plane_stress_d_plus_d_minus_damage_masonry_2d.h"
#include "custom_constitutive/d_plus_d_minus_damage_masonry_3d.h"
#include "custom_constitutive/small_strain_isotropic_plasticity_factory.h"
#include "custom_constitutive/small_strain_kinematic_plasticity_factory.h"
#include "custom_constitutive/generic_small_strain_isotropic_plasticity.h"
#include "custom_constitutive/generic_small_strain_kinematic_plasticity.h"
#include "custom_constitutive/finite_strain_isotropic_plasticity_factory.h"
#include "custom_constitutive/finite_strain_kinematic_plasticity_factory.h"
#include "custom_constitutive/generic_finite_strain_isotropic_plasticity.h"
#include "custom_constitutive/generic_finite_strain_kinematic_plasticity.h"
#include "custom_constitutive/generic_small_strain_isotropic_damage.h"
#include "custom_constitutive/small_strain_isotropic_damage_factory.h"
#include "custom_constitutive/viscous_generalized_kelvin.h"
#include "custom_constitutive/viscous_generalized_maxwell.h"
#include "custom_constitutive/generic_small_strain_viscoplasticity_3d.h"
#include "custom_constitutive/generic_small_strain_d_plus_d_minus_damage.h"
#include "custom_constitutive/plasticity_isotropic_kinematic_j2.h"
#include "custom_constitutive/generic_small_strain_plastic_damage_model.h"
#include "custom_constitutive/generic_small_strain_orthotropic_damage.h"
#include "custom_constitutive/serial_parallel_rule_of_mixtures_law.h"
#include "custom_constitutive/generic_anisotropic_3d_law.h"
#include "custom_constitutive/multi_linear_elastic_1d_law.h"
#include "custom_constitutive/multi_linear_isotropic_plane_stress_2d.h"
#include "custom_constitutive/non_local_elastic_isotropic_damage.h"
#include "custom_constitutive/elastic_local_anisotropic_damage_3d.h"
#include "custom_constitutive/elastic_isotropic_damage_3d_nonlocal_equivalent_strain.h"
#include "custom_constitutive/elastic_anisotropic_damage_3d_nonlocal_equivalent_strain.h"
#include "custom_constitutive/elastic_anisotropic_damage_3d_two_nl_equivalent_strains.h"
#include "custom_constitutive/elastic_anisotropic_damage_3d_two_nl_equivalent_strains_incremental.h"
#include "custom_constitutive/elastic_anisotropic_damage_3d_nonlocal_equivalent_strains_tension_compression.h"
#include "custom_constitutive/elastic_local_anisotropic_damage_3d_tension_compression.h"



// Integrators
#include "custom_constitutive/constitutive_laws_integrators/generic_constitutive_law_integrator_damage.h"
#include "custom_constitutive/constitutive_laws_integrators/generic_constitutive_law_integrator_plasticity.h"
#include "custom_constitutive/constitutive_laws_integrators/generic_constitutive_law_integrator_kinematic_plasticity.h"
#include "custom_constitutive/constitutive_laws_integrators/d+d-constitutive_law_integrators/generic_compression_constitutive_law_integrator.h"
#include "custom_constitutive/constitutive_laws_integrators/d+d-constitutive_law_integrators/generic_tension_constitutive_law_integrator.h"
#include "custom_constitutive/generic_small_strain_high_cycle_fatigue_law.h"

// Yield surfaces
#include "custom_constitutive/yield_surfaces/generic_yield_surface.h"
#include "custom_constitutive/yield_surfaces/von_mises_yield_surface.h"
#include "custom_constitutive/yield_surfaces/modified_mohr_coulomb_yield_surface.h"
#include "custom_constitutive/yield_surfaces/mohr_coulomb_yield_surface.h"
#include "custom_constitutive/yield_surfaces/rankine_yield_surface.h"
#include "custom_constitutive/yield_surfaces/simo_ju_yield_surface.h"
#include "custom_constitutive/yield_surfaces/drucker_prager_yield_surface.h"
#include "custom_constitutive/yield_surfaces/tresca_yield_surface.h"

// Plastic potentials
#include "custom_constitutive/plastic_potentials/generic_plastic_potential.h"
#include "custom_constitutive/plastic_potentials/von_mises_plastic_potential.h"
#include "custom_constitutive/plastic_potentials/tresca_plastic_potential.h"
#include "custom_constitutive/plastic_potentials/modified_mohr_coulomb_plastic_potential.h"
#include "custom_constitutive/plastic_potentials/mohr_coulomb_plastic_potential.h"
#include "custom_constitutive/plastic_potentials/drucker_prager_plastic_potential.h"
#include "custom_constitutive/plastic_potentials/rankine_plastic_potential.h"

// Rules of mixtures
#include "custom_constitutive/rule_of_mixtures_law.h"

#include "custom_constitutive/associative_plastic_damage_model.h"

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
    const HyperElasticIsotropicNeoHookeanPlaneStrain2D  mHyperElasticIsotropicNeoHookeanPlaneStrain2D;
    const LinearElasticOrthotropic2DLaw mLinearElasticOrthotropic2DLaw;

    const SmallStrainJ2Plasticity3D mSmallStrainJ2Plasticity3D;
    const SmallStrainPlasticity3DAth mSmallStrainPlasticity3DAth;
    const SmallStrainJ2PlasticityPlaneStrain2D mSmallStrainJ2PlasticityPlaneStrain2D;
    const SmallStrainIsotropicDamage3D mSmallStrainIsotropicDamage3D;
    const SmallStrainIsotropicDamageAthira3D mSmallStrainIsotropicDamageAthira3D;
    const SmallStrainIsotropicDamage3DMazars mSmallStrainIsotropicDamage3DMazars;
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
    const NonLocalElasticIsotropicDamage mNonLocalElasticIsotropicDamage;
    const ElasticAnisotropicDamage mElasticAnisotropicDamage;
    const ElasticIsotropicDamage3DNonLocalEquivalentStrain mElasticIsotropicDamage3DNonLocalEquivalentStrain;
    const ElasticAnisotropicDamage3DNonLocalEquivalentStrain mElasticAnisotropicDamage3DNonLocalEquivalentStrain;
    const ElasticAnisotropicDamage3DTwoNLEquivalentStrains mElasticAnisotropicDamage3DTwoNLEquivalentStrains;
    const ElasticAnisotropicDamage3DTwoNLEquivalentStrainsIncremental mElasticAnisotropicDamage3DTwoNLEquivalentStrainsIncremental;
    const ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC mElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC;
    const ElasticLocalAnisotropicDamage3DTensionCompression mElasticLocalAnisotropicDamage3DTensionCompression;

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
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainIsotropicDamage3DVonMisesVonMises;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>> mSmallStrainIsotropicDamage3DVonMisesModifiedMohrCoulomb;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<DruckerPragerPlasticPotential<6>>>> mSmallStrainIsotropicDamage3DVonMisesDruckerPrager;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<TrescaPlasticPotential<6>>>> mSmallStrainIsotropicDamage3DVonMisesTresca;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainIsotropicDamage3DModifiedMohrCoulombVonMises;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>> mSmallStrainIsotropicDamage3DModifiedMohrCoulombModifiedMohrCoulomb;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<DruckerPragerPlasticPotential<6>>>> mSmallStrainIsotropicDamage3DModifiedMohrCoulombDruckerPrager;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<TrescaPlasticPotential<6>>>> mSmallStrainIsotropicDamage3DModifiedMohrCoulombTresca;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainIsotropicDamage3DTrescaVonMises;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>> mSmallStrainIsotropicDamage3DTrescaModifiedMohrCoulomb;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<DruckerPragerPlasticPotential<6>>>> mSmallStrainIsotropicDamage3DTrescaDruckerPrager;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<TrescaPlasticPotential<6>>>> mSmallStrainIsotropicDamage3DTrescaTresca;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainIsotropicDamage3DDruckerPragerVonMises;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>> mSmallStrainIsotropicDamage3DDruckerPragerModifiedMohrCoulomb;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<6>>>> mSmallStrainIsotropicDamage3DDruckerPragerDruckerPrager;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<TrescaPlasticPotential<6>>>> mSmallStrainIsotropicDamage3DDruckerPragerTresca;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainIsotropicDamage3DRankineVonMises;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>> mSmallStrainIsotropicDamage3DRankineModifiedMohrCoulomb;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<DruckerPragerPlasticPotential<6>>>> mSmallStrainIsotropicDamage3DRankineDruckerPrager;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<TrescaPlasticPotential<6>>>> mSmallStrainIsotropicDamage3DRankineTresca;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainIsotropicDamage3DSimoJuVonMises;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>> mSmallStrainIsotropicDamage3DSimoJuModifiedMohrCoulomb;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<DruckerPragerPlasticPotential<6>>>> mSmallStrainIsotropicDamage3DSimoJuDruckerPrager;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<TrescaPlasticPotential<6>>>> mSmallStrainIsotropicDamage3DSimoJuTresca;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<MohrCoulombPlasticPotential<6>>>> mSmallStrainIsotropicDamage3DVonMisesMohrCoulomb;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainIsotropicDamage3DMohrCoulombVonMises;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<MohrCoulombYieldSurface<MohrCoulombPlasticPotential<6>>>> mSmallStrainIsotropicDamage3DMohrCoulombMohrCoulomb;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<MohrCoulombYieldSurface<DruckerPragerPlasticPotential<6>>>> mSmallStrainIsotropicDamage3DMohrCoulombDruckerPrager;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<MohrCoulombYieldSurface<TrescaPlasticPotential<6>>>> mSmallStrainIsotropicDamage3DMohrCoulombTresca;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<MohrCoulombPlasticPotential<6>>>> mSmallStrainIsotropicDamage3DTrescaMohrCoulomb;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<MohrCoulombPlasticPotential<6>>>> mSmallStrainIsotropicDamage3DDruckerPragerMohrCoulomb;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<MohrCoulombPlasticPotential<6>>>> mSmallStrainIsotropicDamage3DRankineMohrCoulomb;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<MohrCoulombPlasticPotential<6>>>> mSmallStrainIsotropicDamage3DSimoJuMohrCoulomb;

    /* Small strain 2D */
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainIsotropicDamage2DVonMisesVonMises;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<ModifiedMohrCoulombPlasticPotential<3>>>> mSmallStrainIsotropicDamage2DVonMisesModifiedMohrCoulomb;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<DruckerPragerPlasticPotential<3>>>> mSmallStrainIsotropicDamage2DVonMisesDruckerPrager;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<TrescaPlasticPotential<3>>>> mSmallStrainIsotropicDamage2DVonMisesTresca;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainIsotropicDamage2DModifiedMohrCoulombVonMises;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<3>>>> mSmallStrainIsotropicDamage2DModifiedMohrCoulombModifiedMohrCoulomb;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<DruckerPragerPlasticPotential<3>>>> mSmallStrainIsotropicDamage2DModifiedMohrCoulombDruckerPrager;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<TrescaPlasticPotential<3>>>> mSmallStrainIsotropicDamage2DModifiedMohrCoulombTresca;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainIsotropicDamage2DTrescaVonMises;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<ModifiedMohrCoulombPlasticPotential<3>>>> mSmallStrainIsotropicDamage2DTrescaModifiedMohrCoulomb;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<DruckerPragerPlasticPotential<3>>>> mSmallStrainIsotropicDamage2DTrescaDruckerPrager;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<TrescaPlasticPotential<3>>>> mSmallStrainIsotropicDamage2DTrescaTresca;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainIsotropicDamage2DDruckerPragerVonMises;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<ModifiedMohrCoulombPlasticPotential<3>>>> mSmallStrainIsotropicDamage2DDruckerPragerModifiedMohrCoulomb;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<3>>>> mSmallStrainIsotropicDamage2DDruckerPragerDruckerPrager;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<TrescaPlasticPotential<3>>>> mSmallStrainIsotropicDamage2DDruckerPragerTresca;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainIsotropicDamage2DRankineVonMises;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<ModifiedMohrCoulombPlasticPotential<3>>>> mSmallStrainIsotropicDamage2DRankineModifiedMohrCoulomb;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<DruckerPragerPlasticPotential<3>>>> mSmallStrainIsotropicDamage2DRankineDruckerPrager;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<TrescaPlasticPotential<3>>>> mSmallStrainIsotropicDamage2DRankineTresca;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainIsotropicDamage2DSimoJuVonMises;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<ModifiedMohrCoulombPlasticPotential<3>>>> mSmallStrainIsotropicDamage2DSimoJuModifiedMohrCoulomb;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<DruckerPragerPlasticPotential<3>>>> mSmallStrainIsotropicDamage2DSimoJuDruckerPrager;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<TrescaPlasticPotential<3>>>> mSmallStrainIsotropicDamage2DSimoJuTresca;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<MohrCoulombPlasticPotential<3>>>> mSmallStrainIsotropicDamage2DVonMisesMohrCoulomb;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>> mSmallStrainIsotropicDamage2DMohrCoulombVonMises;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<MohrCoulombYieldSurface<MohrCoulombPlasticPotential<3>>>> mSmallStrainIsotropicDamage2DMohrCoulombMohrCoulomb;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<MohrCoulombYieldSurface<DruckerPragerPlasticPotential<3>>>> mSmallStrainIsotropicDamage2DMohrCoulombDruckerPrager;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<MohrCoulombYieldSurface<TrescaPlasticPotential<3>>>> mSmallStrainIsotropicDamage2DMohrCoulombTresca;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<MohrCoulombPlasticPotential<3>>>> mSmallStrainIsotropicDamage2DTrescaMohrCoulomb;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<MohrCoulombPlasticPotential<3>>>> mSmallStrainIsotropicDamage2DDruckerPragerMohrCoulomb;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<MohrCoulombPlasticPotential<3>>>> mSmallStrainIsotropicDamage2DRankineMohrCoulomb;
    const GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<MohrCoulombPlasticPotential<3>>>> mSmallStrainIsotropicDamage2DSimoJuMohrCoulomb;

    // HCF (High Cycle Fatigue)
    const GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainHighCycleFatigue3DLawVonMisesVonMises;
    const GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>> mSmallStrainHighCycleFatigue3DLawVonMisesModifiedMohrCoulomb;
    const GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<DruckerPragerPlasticPotential<6>>>> mSmallStrainHighCycleFatigue3DLawVonMisesDruckerPrager;
    const GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<TrescaPlasticPotential<6>>>> mSmallStrainHighCycleFatigue3DLawVonMisesTresca;
    const GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainHighCycleFatigue3DLawModifiedMohrCoulombVonMises;
    const GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>> mSmallStrainHighCycleFatigue3DLawModifiedMohrCoulombModifiedMohrCoulomb;
    const GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<DruckerPragerPlasticPotential<6>>>> mSmallStrainHighCycleFatigue3DLawModifiedMohrCoulombDruckerPrager;
    const GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<TrescaPlasticPotential<6>>>> mSmallStrainHighCycleFatigue3DLawModifiedMohrCoulombTresca;
    const GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainHighCycleFatigue3DLawTrescaVonMises;
    const GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>> mSmallStrainHighCycleFatigue3DLawTrescaModifiedMohrCoulomb;
    const GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<DruckerPragerPlasticPotential<6>>>> mSmallStrainHighCycleFatigue3DLawTrescaDruckerPrager;
    const GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<TrescaPlasticPotential<6>>>> mSmallStrainHighCycleFatigue3DLawTrescaTresca;
    const GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainHighCycleFatigue3DLawDruckerPragerVonMises;
    const GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>> mSmallStrainHighCycleFatigue3DLawDruckerPragerModifiedMohrCoulomb;
    const GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<6>>>> mSmallStrainHighCycleFatigue3DLawDruckerPragerDruckerPrager;
    const GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<TrescaPlasticPotential<6>>>> mSmallStrainHighCycleFatigue3DLawDruckerPragerTresca;
    const GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainHighCycleFatigue3DLawRankineVonMises;
    const GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>> mSmallStrainHighCycleFatigue3DLawRankineModifiedMohrCoulomb;
    const GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<DruckerPragerPlasticPotential<6>>>> mSmallStrainHighCycleFatigue3DLawRankineDruckerPrager;
    const GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<TrescaPlasticPotential<6>>>> mSmallStrainHighCycleFatigue3DLawRankineTresca;
    const GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>> mSmallStrainHighCycleFatigue3DLawSimoJuVonMises;
    const GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>> mSmallStrainHighCycleFatigue3DLawSimoJuModifiedMohrCoulomb;
    const GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<DruckerPragerPlasticPotential<6>>>> mSmallStrainHighCycleFatigue3DLawSimoJuDruckerPrager;
    const GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<TrescaPlasticPotential<6>>>> mSmallStrainHighCycleFatigue3DLawSimoJuTresca;

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
	const ParallelRuleOfMixturesLaw<2> mParallelRuleOfMixturesLaw2D;

    // Anisotropic law
    const GenericAnisotropic3DLaw mGenericAnisotropic3DLaw;

    const AssociativePlasticDamageModel <VonMisesYieldSurface<VonMisesPlasticPotential<6>>> mAssociativePlasticDamageModel3DVonMisesVonMises;
    const AssociativePlasticDamageModel <DruckerPragerYieldSurface<DruckerPragerPlasticPotential<6>>> mAssociativePlasticDamageModel3DDruckerPragerDruckerPrager;
    const AssociativePlasticDamageModel <ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>> mAssociativePlasticDamageModel3DModifiedMohrCoulombModifiedMohrCoulomb;
    const AssociativePlasticDamageModel <RankineYieldSurface<RankinePlasticPotential<6>>> mAssociativePlasticDamageModel3DRankineRankine;

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

#endif // KRATOS_CONSTITUTIVE_LAWS_APPLICATION_H_INCLUDED  defined
