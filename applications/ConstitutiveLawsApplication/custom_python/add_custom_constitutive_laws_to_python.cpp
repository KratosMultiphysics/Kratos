// KRATOS ___                _   _ _         _   _             __                       _
//       / __\___  _ __  ___| |_(_) |_ _   _| |_(_)_   _____  / /  __ ___      _____   /_\  _ __  _ __
//      / /  / _ \| '_ \/ __| __| | __| | | | __| \ \ / / _ \/ /  / _` \ \ /\ / / __| //_\\| '_ \| '_  |
//     / /__| (_) | | | \__ \ |_| | |_| |_| | |_| |\ V /  __/ /__| (_| |\ V  V /\__ \/  _  \ |_) | |_) |
//     \____/\___/|_| |_|___/\__|_|\__|\__,_|\__|_| \_/ \___\____/\__,_| \_/\_/ |___/\_/ \_/ .__/| .__/
//                                                                                         |_|   |_|
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Riccardo Rossi
//


// System includes

// External includes

// Project includes
#include "custom_python/add_custom_constitutive_laws_to_python.h"

// Elastic laws
#include "custom_constitutive/finite_strains/hyperelasticity/hyper_elastic_isotropic_ogden_1d.h"
#include "custom_constitutive/finite_strains/hyperelasticity/hyper_elastic_isotropic_henky_1d.h"
#include "custom_constitutive/small_strains/plasticity/truss_plasticity_constitutive_law.h"
#include "custom_constitutive/small_strains/linear/elastic_isotropic_plane_stress_uncoupled_shear.h"
#include "custom_constitutive/finite_strains/hyperelasticity/hyper_elastic_isotropic_kirchhoff_3d.h"
#include "custom_constitutive/finite_strains/hyperelasticity/hyper_elastic_isotropic_kirchhoff_plane_stress_2d.h"
#include "custom_constitutive/finite_strains/hyperelasticity/hyper_elastic_isotropic_kirchhoff_plane_strain_2d.h"
#include "custom_constitutive/finite_strains/hyperelasticity/hyper_elastic_isotropic_neo_hookean_3d.h"
#include "custom_constitutive/finite_strains/hyperelasticity/hyper_elastic_isotropic_q_incomp_isoch_neo_hook_3d.h"
#include "custom_constitutive/finite_strains/hyperelasticity/hyper_elastic_isotropic_neo_hookean_plane_strain_2d.h"
#include "custom_constitutive/small_strains/linear/linear_elastic_orthotropic_2D_law.h"
#include "custom_constitutive/small_strains/plasticity/small_strain_j2_plasticity_3d.h"
#include "custom_constitutive/small_strains/plasticity/small_strain_j2_plasticity_plane_strain_2d.h"
#include "custom_constitutive/small_strains/damage/small_strain_isotropic_damage_3d.h"
#include "custom_constitutive/small_strains/damage/small_strain_isotropic_damage_implex_3d.h"
#include "custom_constitutive/small_strains/damage/small_strain_isotropic_damage_traction_only_3d.h"
#include "custom_constitutive/small_strains/damage/small_strain_isotropic_damage_traction_only_implex_3d.h"
#include "custom_constitutive/small_strains/linear/wrinkling_linear_2d_law.h"
#include "custom_constitutive/small_strains/linear/multi_linear_elastic_1d_law.h"
#include "custom_constitutive/small_strains/linear/multi_linear_isotropic_plane_stress_2d.h"

// Plastic, damage laws and viscosities
#include "custom_constitutive/small_strains/damage/generic_small_strain_isotropic_damage_plane_stress.h"
#include "custom_constitutive/small_strains/damage/small_strain_isotropic_damage_plane_strain_2d.h"
#include "custom_constitutive/small_strains/plasticity/small_strain_isotropic_plasticity_factory.h"
#include "custom_constitutive/small_strains/plasticity/small_strain_kinematic_plasticity_factory.h"
#include "custom_constitutive/finite_strains/plasticity/finite_strain_isotropic_plasticity_factory.h"
#include "custom_constitutive/small_strains/damage/small_strain_isotropic_damage_factory.h"
#include "custom_constitutive/small_strains/viscous/viscous_generalized_maxwell.h"
#include "custom_constitutive/small_strains/viscous/viscous_generalized_kelvin.h"
#include "custom_constitutive/small_strains/viscoplasticity/generic_small_strain_viscoplasticity_3d.h"
#include "custom_constitutive/small_strains/plasticity/generic_small_strain_isotropic_plasticity.h"
#include "custom_constitutive/finite_strains/plasticity/generic_finite_strain_isotropic_plasticity.h"
#include "custom_constitutive/finite_strains/plasticity/generic_finite_strain_kinematic_plasticity.h"
#include "custom_constitutive/small_strains/damage/generic_small_strain_isotropic_damage.h"
#include "custom_constitutive/small_strains/damage/generic_small_strain_d_plus_d_minus_damage.h"
#include "custom_constitutive/small_strains/plasticity/generic_small_strain_kinematic_plasticity.h"
#include "custom_constitutive/small_strains/fatigue/generic_small_strain_high_cycle_fatigue_law.h"
#include "custom_constitutive/small_strains/plasticity/plasticity_isotropic_kinematic_j2.h"
#include "custom_constitutive/small_strains/damage/plane_stress_d_plus_d_minus_damage_masonry_2d.h"
#include "custom_constitutive/small_strains/damage/d_plus_d_minus_damage_masonry_3d.h"
#include "custom_constitutive/small_strains/plastic_damage/generic_small_strain_plastic_damage_model.h"
#include "custom_constitutive/small_strains/damage/generic_small_strain_orthotropic_damage.h"
#include "custom_constitutive/composites/serial_parallel_rule_of_mixtures_law.h"

// Integrators
#include "custom_constitutive/auxiliary_files/cl_integrators/generic_cl_integrator_damage.h"
#include "custom_constitutive/auxiliary_files/cl_integrators/generic_cl_integrator_plasticity.h"
#include "custom_constitutive/auxiliary_files/cl_integrators/generic_cl_integrator_kinematic_plasticity.h"
#include "custom_constitutive/auxiliary_files/cl_integrators/d+d-cl_integrators/generic_compression_cl_integrator.h"
#include "custom_constitutive/auxiliary_files/cl_integrators/d+d-cl_integrators/generic_tension_cl_integrator.h"

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

// Thermal CL's
#include "custom_constitutive/thermal/small_strains/elastic/thermal_elastic_isotropic_3d.h"
#include "custom_constitutive/thermal/small_strains/elastic/thermal_linear_plane_strain.h"
#include "custom_constitutive/thermal/small_strains/elastic/thermal_linear_plane_stress.h"
#include "custom_constitutive/thermal/small_strains/damage/generic_small_strain_thermal_isotropic_damage.h"
#include "custom_constitutive/thermal/small_strains/damage/generic_small_strain_thermal_isotropic_damage_plane_stress.h"

#include "custom_constitutive/thermal/auxiliary_files/thermal_yield_surfaces/thermal_von_mises_yield_surface.h"

namespace Kratos::Python {

void AddCustomConstitutiveLawsToPython(pybind11::module& m)
{
    namespace py = pybind11;

    py::class_< MultiLinearElastic1DLaw, typename MultiLinearElastic1DLaw::Pointer, TrussConstitutiveLaw >
    (m, "MultiLinearElastic1DLaw").def(py::init<>() )
    ;

    py::class_< MultiLinearIsotropicPlaneStress2D, typename MultiLinearIsotropicPlaneStress2D::Pointer, LinearPlaneStress >
    (m, "MultiLinearIsotropicPlaneStress2D").def(py::init<>() )
    ;

    py::class_< TractionSeparationLaw3D<3>, typename TractionSeparationLaw3D<3>::Pointer,  ConstitutiveLaw  >
    (m,"TractionSeparationLaw3D").def(py::init<>())
    ;

    py::class_< WrinklingLinear2DLaw, typename WrinklingLinear2DLaw::Pointer, ConstitutiveLaw >
    (m, "WrinklingLinear2DLaw").def(py::init<>() )
    ;

    py::class_< TrussPlasticityConstitutiveLaw, typename TrussPlasticityConstitutiveLaw::Pointer, ConstitutiveLaw >
    (m, "TrussPlasticityConstitutiveLaw").def(py::init<>() )
    ;

    py::class_< HyperElasticIsotropicOgden1D, typename HyperElasticIsotropicOgden1D::Pointer, ConstitutiveLaw >
    (m, "HyperElasticIsotropicOgden1D").def(py::init<>() )
    ;

    py::class_< HyperElasticIsotropicHenky1D, typename HyperElasticIsotropicHenky1D::Pointer, ConstitutiveLaw >
    (m, "HyperElasticIsotropicHenky1D").def(py::init<>() )
    ;

    py::class_< ElasticIsotropicPlaneStressUncoupledShear, typename ElasticIsotropicPlaneStressUncoupledShear::Pointer, ConstitutiveLaw >
    (m, "ElasticPlaneStressUncoupledShear2DLaw").def(py::init<>() )
    ;

    py::class_< HyperElasticIsotropicKirchhoff3D, typename HyperElasticIsotropicKirchhoff3D::Pointer, ConstitutiveLaw >
    (m, "KirchhoffSaintVenant3DLaw").def(py::init<>() )
    ;

    py::class_< HyperElasticIsotropicKirchhoffPlaneStress2D, typename HyperElasticIsotropicKirchhoffPlaneStress2D::Pointer, ConstitutiveLaw >
    (m, "KirchhoffSaintVenantPlaneStress2DLaw").def(py::init<>() )
    ;

    py::class_< HyperElasticIsotropicKirchhoffPlaneStrain2D, typename HyperElasticIsotropicKirchhoffPlaneStrain2D::Pointer, ConstitutiveLaw >
    (m, "KirchhoffSaintVenantPlaneStrain2DLaw").def(py::init<>() )
    ;

    py::class_< HyperElasticIsotropicNeoHookean3D, typename HyperElasticIsotropicNeoHookean3D::Pointer, ConstitutiveLaw >
    (m, "HyperElastic3DLaw").def(py::init<>() )
    ;

    py::class_< HyperElasticIsotropicQuasiIncompressibleIshochoricNeoHookean3D, typename HyperElasticIsotropicQuasiIncompressibleIshochoricNeoHookean3D::Pointer, ConstitutiveLaw >
    (m, "HyperElasticQuasiIncompressibleNeoHookean3DLaw").def(py::init<>() )
    ;

    py::class_< HyperElasticIsotropicNeoHookeanPlaneStrain2D, typename HyperElasticIsotropicNeoHookeanPlaneStrain2D::Pointer, ConstitutiveLaw >
    (m, "HyperElasticPlaneStrain2DLaw").def(py::init<>() )
    ;

    py::class_< LinearElasticOrthotropic2DLaw, typename LinearElasticOrthotropic2DLaw::Pointer, ConstitutiveLaw >
    (m,"LinearElasticOrthotropic2DLaw").def(py::init<>())
    ;

    py::class_< SmallStrainJ2PlasticityPlaneStrain2D, typename SmallStrainJ2PlasticityPlaneStrain2D::Pointer,  ConstitutiveLaw >
    (m,"SmallStrainJ2PlasticityPlaneStrain2DLaw").def(py::init<>())
    ;

    py::class_< SmallStrainJ2Plasticity3D, typename SmallStrainJ2Plasticity3D::Pointer,  ConstitutiveLaw >
    (m,"SmallStrainJ2Plasticity3DLaw").def(py::init<>())
    ;

    py::class_< SmallStrainIsotropicDamagePlaneStrain2D, typename SmallStrainIsotropicDamagePlaneStrain2D::Pointer,  ConstitutiveLaw  >
    (m,"SmallStrainIsotropicDamagePlaneStrain2DLaw").def(py::init<>())
    ;

    py::class_< SmallStrainIsotropicDamage3D, typename SmallStrainIsotropicDamage3D::Pointer,  ConstitutiveLaw  >
    (m,"SmallStrainIsotropicDamage3DLaw").def(py::init<>())
    ;

    py::class_< SmallStrainIsotropicDamageImplex3D, typename SmallStrainIsotropicDamageImplex3D::Pointer,  ConstitutiveLaw  >
    (m,"SmallStrainIsotropicDamageImplex3DLaw").def(py::init<>())
    ;

    py::class_< SmallStrainIsotropicDamageTractionOnly3D, typename SmallStrainIsotropicDamageTractionOnly3D::Pointer,  ConstitutiveLaw  >
    (m,"SmallStrainIsotropicDamageTractionOnly3DLaw").def(py::init<>())
    ;

    py::class_< SmallStrainIsotropicDamageTractionOnlyImplex3D, typename SmallStrainIsotropicDamageTractionOnlyImplex3D::Pointer,  ConstitutiveLaw  >
    (m,"SmallStrainIsotropicDamageTractionOnlyImplex3DLaw").def(py::init<>())
    ;

    py::class_< PlasticityIsotropicKinematicJ2, typename PlasticityIsotropicKinematicJ2::Pointer,  ConstitutiveLaw >
    (m,"PlasticityIsotropicKinematicJ2Law").def(py::init<>())
    ;

    py::class_< SmallStrainIsotropicDamageFactory, typename SmallStrainIsotropicDamageFactory::Pointer,  ConstitutiveLaw  >
    (m,"SmallStrainIsotropicDamageFactory").def(py::init<>())
    ;

    py::class_< SmallStrainIsotropicPlasticityFactory, typename SmallStrainIsotropicPlasticityFactory::Pointer,  ConstitutiveLaw  >
    (m,"SmallStrainIsotropicPlasticityFactory").def(py::init<>())
    ;

    py::class_< SmallStrainKinematicPlasticityFactory, typename SmallStrainKinematicPlasticityFactory::Pointer,  ConstitutiveLaw  >
    (m,"SmallStrainKinematicPlasticityFactory").def(py::init<>())
    ;

    py::class_< FiniteStrainIsotropicPlasticityFactory, typename FiniteStrainIsotropicPlasticityFactory::Pointer,  ConstitutiveLaw  >
    (m,"FiniteStrainIsotropicPlasticityFactory").def(py::init<>())
    ;

    py::class_< ViscousGeneralizedKelvin<ElasticIsotropic3D>, typename ViscousGeneralizedKelvin<ElasticIsotropic3D>::Pointer,  ConstitutiveLaw  >
    (m,"ViscousGeneralizedKelvin3D").def(py::init<>())
    ;

    py::class_< ViscousGeneralizedMaxwell<ElasticIsotropic3D>, typename ViscousGeneralizedMaxwell<ElasticIsotropic3D>::Pointer,  ConstitutiveLaw  >
    (m,"ViscousGeneralizedMaxwell3D").def(py::init<>())
    ;


    // Custom Constitutive Laws Registration
    // Isotropic Plasticity
    /* Small strain */

    py::class_< GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<VonMisesPlasticPotential<3>>>>,
    typename GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<VonMisesPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicPlasticityPlaneStrainVonMisesVonMises").def(py::init<>());

    py::class_< GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<ModifiedMohrCoulombPlasticPotential<3>>>>,
    typename GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<ModifiedMohrCoulombPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicPlasticityPlaneStrainVonMisesModifiedMohrCoulomb").def(py::init<>());

    py::class_< GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<DruckerPragerPlasticPotential<3>>>>,
    typename GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<DruckerPragerPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicPlasticityPlaneStrainVonMisesDruckerPrager").def(py::init<>());

    py::class_< GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<TrescaPlasticPotential<3>>>>,
    typename GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<TrescaPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicPlasticityPlaneStrainVonMisesTresca").def(py::init<>());

    py::class_< GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>>,
    typename GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicPlasticityPlaneStrainModifiedMohrCoulombVonMises").def(py::init<>());

    py::class_< GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<3>>>>,
    typename GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicPlasticityPlaneStrainModifiedMohrCoulombModifiedMohrCoulomb").def(py::init<>());

    py::class_< GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<DruckerPragerPlasticPotential<3>>>>,
    typename GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<DruckerPragerPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicPlasticityPlaneStrainModifiedMohrCoulombDruckerPrager").def(py::init<>());

    py::class_< GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<TrescaPlasticPotential<3>>>>,
    typename GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<TrescaPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicPlasticityPlaneStrainModifiedMohrCoulombTresca").def(py::init<>());

    py::class_< GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<VonMisesPlasticPotential<3>>>>,
    typename GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<VonMisesPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicPlasticityPlaneStrainTrescaVonMises").def(py::init<>());

    py::class_< GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<ModifiedMohrCoulombPlasticPotential<3>>>>,
    typename GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<ModifiedMohrCoulombPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicPlasticityPlaneStrainTrescaModifiedMohrCoulomb").def(py::init<>());

    py::class_< GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<DruckerPragerPlasticPotential<3>>>>,
    typename GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<DruckerPragerPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicPlasticityPlaneStrainTrescaDruckerPrager").def(py::init<>());

    py::class_< GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<TrescaPlasticPotential<3>>>>,
    typename GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<TrescaPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicPlasticityPlaneStrainTrescaTresca").def(py::init<>());

    py::class_< GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<VonMisesPlasticPotential<3>>>>,
    typename GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<VonMisesPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicPlasticityPlaneStrainDruckerPragerVonMises").def(py::init<>());

    py::class_< GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<ModifiedMohrCoulombPlasticPotential<3>>>>,
    typename GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<ModifiedMohrCoulombPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicPlasticityPlaneStrainDruckerPragerModifiedMohrCoulomb").def(py::init<>());

    py::class_< GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<3>>>>,
    typename GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicPlasticityPlaneStrainDruckerPragerDruckerPrager").def(py::init<>());

    py::class_< GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<TrescaPlasticPotential<3>>>>,
    typename GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<TrescaPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicPlasticityPlaneStrainDruckerPragerTresca").def(py::init<>());

    py::class_< GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<MohrCoulombPlasticPotential<3>>>>,
    typename GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<MohrCoulombPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicPlasticityPlaneStrainVonMisesMohrCoulomb").def(py::init<>());

	py::class_< GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>>,
    typename GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicPlasticityPlaneStrainMohrCoulombVonMises").def(py::init<>());

    py::class_< GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<MohrCoulombPlasticPotential<3>>>>,
    typename GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<MohrCoulombPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicPlasticityPlaneStrainMohrCoulombMohrCoulomb").def(py::init<>());

    py::class_< GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<DruckerPragerPlasticPotential<3>>>>,
    typename GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<DruckerPragerPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicPlasticityPlaneStrainMohrCoulombDruckerPrager").def(py::init<>());

    py::class_< GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<TrescaPlasticPotential<3>>>>,
    typename GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<TrescaPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicPlasticityPlaneStrainMohrCoulombTresca").def(py::init<>());

    py::class_< GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<MohrCoulombPlasticPotential<3>>>>,
    typename GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<MohrCoulombPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicPlasticityPlaneStrainTrescaMohrCoulomb").def(py::init<>());

    py::class_< GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<MohrCoulombPlasticPotential<3>>>>,
    typename GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<MohrCoulombPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicPlasticityPlaneStrainDruckerPragerMohrCoulomb").def(py::init<>());

    // Isotropic Plasticity 3D
    py::class_< GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicPlasticity3DVonMisesVonMises").def(py::init<>());

    py::class_< GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>,
    typename GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicPlasticity3DVonMisesModifiedMohrCoulomb").def(py::init<>());

    py::class_< GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<DruckerPragerPlasticPotential<6>>>>,
    typename GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<DruckerPragerPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicPlasticity3DVonMisesDruckerPrager").def(py::init<>());

    py::class_< GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<TrescaPlasticPotential<6>>>>,
    typename GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<TrescaPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicPlasticity3DVonMisesTresca").def(py::init<>());

    py::class_< GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicPlasticity3DModifiedMohrCoulombVonMises").def(py::init<>());

    py::class_< GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>,
    typename GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicPlasticity3DModifiedMohrCoulombModifiedMohrCoulomb").def(py::init<>());

    py::class_< GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<DruckerPragerPlasticPotential<6>>>>,
    typename GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<DruckerPragerPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicPlasticity3DModifiedMohrCoulombDruckerPrager").def(py::init<>());

    py::class_< GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<TrescaPlasticPotential<6>>>>,
    typename GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<TrescaPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicPlasticity3DModifiedMohrCoulombTresca").def(py::init<>());

    py::class_< GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicPlasticity3DTrescaVonMises").def(py::init<>());

    py::class_< GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>,
    typename GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicPlasticity3DTrescaModifiedMohrCoulomb").def(py::init<>());

    py::class_< GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<DruckerPragerPlasticPotential<6>>>>,
    typename GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<DruckerPragerPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicPlasticity3DTrescaDruckerPrager").def(py::init<>());

    py::class_< GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<TrescaPlasticPotential<6>>>>,
    typename GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<TrescaPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicPlasticity3DTrescaTresca").def(py::init<>());

    py::class_< GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicPlasticity3DDruckerPragerVonMises").def(py::init<>());

    py::class_< GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>,
    typename GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicPlasticity3DDruckerPragerModifiedMohrCoulomb").def(py::init<>());

    py::class_< GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<6>>>>,
    typename GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicPlasticity3DDruckerPragerDruckerPrager").def(py::init<>());

    py::class_< GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<TrescaPlasticPotential<6>>>>,
    typename GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<TrescaPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicPlasticity3DDruckerPragerTresca").def(py::init<>());

    py::class_< GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<MohrCoulombPlasticPotential<6>>>>,
    typename GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<MohrCoulombPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicPlasticity3DVonMisesMohrCoulomb").def(py::init<>());

	py::class_< GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicPlasticity3DMohrCoulombVonMises").def(py::init<>());

    py::class_< GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<MohrCoulombPlasticPotential<6>>>>,
    typename GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<MohrCoulombPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicPlasticity3DMohrCoulombMohrCoulomb").def(py::init<>());

    py::class_< GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<DruckerPragerPlasticPotential<6>>>>,
    typename GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<DruckerPragerPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicPlasticity3DMohrCoulombDruckerPrager").def(py::init<>());

    py::class_< GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<TrescaPlasticPotential<6>>>>,
    typename GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<TrescaPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicPlasticity3DMohrCoulombTresca").def(py::init<>());

    py::class_< GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<MohrCoulombPlasticPotential<6>>>>,
    typename GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<MohrCoulombPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicPlasticity3DTrescaMohrCoulomb").def(py::init<>());

    py::class_< GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<MohrCoulombPlasticPotential<6>>>>,
    typename GenericSmallStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<MohrCoulombPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicPlasticity3DDruckerPragerMohrCoulomb").def(py::init<>());

    // Plastic Damage Model
    py::class_< GenericSmallStrainPlasticDamageModel <GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>, GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainPlasticDamageModel <GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>, GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainPlasticDamageModel3DVonMisesVonMisesVonMises").def(py::init<>());

    py::class_< GenericSmallStrainPlasticDamageModel <GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>, GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainPlasticDamageModel <GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>, GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainPlasticDamageModel3DVonMisesVonMisesDruckerPrager").def(py::init<>());


    // Kinematic Plasticity
    /* Small strain */

	py::class_< GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<VonMisesPlasticPotential<3>>>>,
    typename GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<VonMisesPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainKinematicPlasticityPlaneStrainVonMisesVonMises").def(py::init<>());

	py::class_< GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<ModifiedMohrCoulombPlasticPotential<3>>>>,
    typename GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<ModifiedMohrCoulombPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainKinematicPlasticityPlaneStrainVonMisesModifiedMohrCoulomb").def(py::init<>());

	py::class_< GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<DruckerPragerPlasticPotential<3>>>>,
    typename GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<DruckerPragerPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainKinematicPlasticityPlaneStrainVonMisesDruckerPrager").def(py::init<>());

	py::class_< GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<TrescaPlasticPotential<3>>>>,
    typename GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<TrescaPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainKinematicPlasticityPlaneStrainVonMisesTresca").def(py::init<>());

	py::class_< GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>>,
    typename GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainKinematicPlasticityPlaneStrainModifiedMohrCoulombVonMises").def(py::init<>());

	py::class_< GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<3>>>>,
    typename GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainKinematicPlasticityPlaneStrainModifiedMohrCoulombModifiedMohrCoulomb").def(py::init<>());

	py::class_< GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<ModifiedMohrCoulombYieldSurface<DruckerPragerPlasticPotential<3>>>>,
    typename GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<ModifiedMohrCoulombYieldSurface<DruckerPragerPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainKinematicPlasticityPlaneStrainModifiedMohrCoulombDruckerPrager").def(py::init<>());

	py::class_< GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<ModifiedMohrCoulombYieldSurface<TrescaPlasticPotential<3>>>>,
    typename GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<ModifiedMohrCoulombYieldSurface<TrescaPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainKinematicPlasticityPlaneStrainModifiedMohrCoulombTresca").def(py::init<>());

	py::class_< GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<VonMisesPlasticPotential<3>>>>,
    typename GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<VonMisesPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainKinematicPlasticityPlaneStrainTrescaVonMises").def(py::init<>());

	py::class_< GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<ModifiedMohrCoulombPlasticPotential<3>>>>,
    typename GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<ModifiedMohrCoulombPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainKinematicPlasticityPlaneStrainTrescaModifiedMohrCoulomb").def(py::init<>());

	py::class_< GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<DruckerPragerPlasticPotential<3>>>>,
    typename GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<DruckerPragerPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainKinematicPlasticityPlaneStrainTrescaDruckerPrager").def(py::init<>());

	py::class_< GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<TrescaPlasticPotential<3>>>>,
    typename GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<TrescaPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainKinematicPlasticityPlaneStrainTrescaTresca").def(py::init<>());

	py::class_< GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<VonMisesPlasticPotential<3>>>>,
    typename GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<VonMisesPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainKinematicPlasticityPlaneStrainDruckerPragerVonMises").def(py::init<>());

	py::class_< GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<ModifiedMohrCoulombPlasticPotential<3>>>>,
    typename GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<ModifiedMohrCoulombPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainKinematicPlasticityPlaneStrainDruckerPragerModifiedMohrCoulomb").def(py::init<>());

	py::class_< GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<3>>>>,
    typename GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainKinematicPlasticityPlaneStrainDruckerPragerDruckerPrager").def(py::init<>());

	py::class_< GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<TrescaPlasticPotential<3>>>>,
    typename GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<TrescaPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainKinematicPlasticityPlaneStrainDruckerPragerTresca").def(py::init<>());


    // Kinematic isotropic 3D plasticity
	py::class_< GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainKinematicPlasticity3DVonMisesVonMises").def(py::init<>());

	py::class_< GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>,
    typename GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainKinematicPlasticity3DVonMisesModifiedMohrCoulomb").def(py::init<>());

	py::class_< GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<DruckerPragerPlasticPotential<6>>>>,
    typename GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<DruckerPragerPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainKinematicPlasticity3DVonMisesDruckerPrager").def(py::init<>());

	py::class_< GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<TrescaPlasticPotential<6>>>>,
    typename GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<TrescaPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainKinematicPlasticity3DVonMisesTresca").def(py::init<>());

	py::class_< GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainKinematicPlasticity3DModifiedMohrCoulombVonMises").def(py::init<>());

	py::class_< GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>,
    typename GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainKinematicPlasticity3DModifiedMohrCoulombModifiedMohrCoulomb").def(py::init<>());

	py::class_< GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<ModifiedMohrCoulombYieldSurface<DruckerPragerPlasticPotential<6>>>>,
    typename GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<ModifiedMohrCoulombYieldSurface<DruckerPragerPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainKinematicPlasticity3DModifiedMohrCoulombDruckerPrager").def(py::init<>());

	py::class_< GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<ModifiedMohrCoulombYieldSurface<TrescaPlasticPotential<6>>>>,
    typename GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<ModifiedMohrCoulombYieldSurface<TrescaPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainKinematicPlasticity3DModifiedMohrCoulombTresca").def(py::init<>());

	py::class_< GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainKinematicPlasticity3DTrescaVonMises").def(py::init<>());

	py::class_< GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>,
    typename GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainKinematicPlasticity3DTrescaModifiedMohrCoulomb").def(py::init<>());

	py::class_< GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<DruckerPragerPlasticPotential<6>>>>,
    typename GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<DruckerPragerPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainKinematicPlasticity3DTrescaDruckerPrager").def(py::init<>());

	py::class_< GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<TrescaPlasticPotential<6>>>>,
    typename GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<TrescaPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainKinematicPlasticity3DTrescaTresca").def(py::init<>());

	py::class_< GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainKinematicPlasticity3DDruckerPragerVonMises").def(py::init<>());

	py::class_< GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>,
    typename GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainKinematicPlasticity3DDruckerPragerModifiedMohrCoulomb").def(py::init<>());

	py::class_< GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<6>>>>,
    typename GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainKinematicPlasticity3DDruckerPragerDruckerPrager").def(py::init<>());

	py::class_< GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<TrescaPlasticPotential<6>>>>,
    typename GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<TrescaPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainKinematicPlasticity3DDruckerPragerTresca").def(py::init<>());

    // HCF
    py::class_<  GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw > (m,"SmallStrainHighCycleFatigue3DLawVonMises").def(py::init<>());

    py::class_<  GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>,
    typename GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw > (m,"SmallStrainHighCycleFatigue3DLawModifiedMohrCoulomb").def(py::init<>());

    py::class_<  GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<TrescaPlasticPotential<6>>>>,
    typename GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<TrescaPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw > (m,"SmallStrainHighCycleFatigue3DLawTresca").def(py::init<>());

    py::class_<  GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<6>>>>,
    typename GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw > (m,"SmallStrainHighCycleFatigue3DLawDruckerPrager").def(py::init<>());

    py::class_<  GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw > (m,"SmallStrainHighCycleFatigue3DLawSimoJu").def(py::init<>());

    py::class_<  GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<MohrCoulombYieldSurface<MohrCoulombPlasticPotential<6>>>>,
    typename GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<MohrCoulombYieldSurface<MohrCoulombPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw > (m,"SmallStrainHighCycleFatigue3DLawMohrCoulomb").def(py::init<>());

    py::class_<  GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<RankinePlasticPotential<6>>>>,
    typename GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<RankinePlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw > (m,"SmallStrainHighCycleFatigue3DLawRankine").def(py::init<>());

    // kinematic plasticity
    py::class_< GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<MohrCoulombPlasticPotential<6>>>>,
    typename GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<MohrCoulombPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainKinematicPlasticity3DVonMisesMohrCoulomb").def(py::init<>());

	py::class_< GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<MohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<MohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainKinematicPlasticity3DMohrCoulombVonMises").def(py::init<>());

    py::class_< GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<MohrCoulombYieldSurface<MohrCoulombPlasticPotential<6>>>>,
    typename GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<MohrCoulombYieldSurface<MohrCoulombPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainKinematicPlasticity3DMohrCoulombMohrCoulomb").def(py::init<>());

    py::class_< GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<MohrCoulombYieldSurface<DruckerPragerPlasticPotential<6>>>>,
    typename GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<MohrCoulombYieldSurface<DruckerPragerPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainKinematicPlasticity3DMohrCoulombDruckerPrager").def(py::init<>());

    py::class_< GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<MohrCoulombYieldSurface<TrescaPlasticPotential<6>>>>,
    typename GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<MohrCoulombYieldSurface<TrescaPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainKinematicPlasticity3DMohrCoulombTresca").def(py::init<>());

    py::class_< GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<MohrCoulombPlasticPotential<6>>>>,
    typename GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<MohrCoulombPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainKinematicPlasticity3DTrescaMohrCoulomb").def(py::init<>());

    py::class_< GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<MohrCoulombPlasticPotential<6>>>>,
    typename GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<MohrCoulombPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainKinematicPlasticity3DDruckerPragerMohrCoulomb").def(py::init<>());

    /* Finite strain */
    // Isotropic
    py::class_< GenericFiniteStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericFiniteStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"FiniteStrainIsotropicPlasticity3DVonMisesVonMises").def(py::init<>());

    py::class_< GenericFiniteStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>,
    typename GenericFiniteStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"FiniteStrainIsotropicPlasticity3DVonMisesModifiedMohrCoulomb").def(py::init<>());

    py::class_< GenericFiniteStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<DruckerPragerPlasticPotential<6>>>>,
    typename GenericFiniteStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<DruckerPragerPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"FiniteStrainIsotropicPlasticity3DVonMisesDruckerPrager").def(py::init<>());

    py::class_< GenericFiniteStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<TrescaPlasticPotential<6>>>>,
    typename GenericFiniteStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<TrescaPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"FiniteStrainIsotropicPlasticity3DVonMisesTresca").def(py::init<>());

    py::class_< GenericFiniteStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericFiniteStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"FiniteStrainIsotropicPlasticity3DModifiedMohrCoulombVonMises").def(py::init<>());

    py::class_< GenericFiniteStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>,
    typename GenericFiniteStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"FiniteStrainIsotropicPlasticity3DModifiedMohrCoulombModifiedMohrCoulomb").def(py::init<>());

    py::class_< GenericFiniteStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<DruckerPragerPlasticPotential<6>>>>,
    typename GenericFiniteStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<DruckerPragerPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"FiniteStrainIsotropicPlasticity3DModifiedMohrCoulombDruckerPrager").def(py::init<>());

    py::class_< GenericFiniteStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<TrescaPlasticPotential<6>>>>,
    typename GenericFiniteStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<TrescaPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"FiniteStrainIsotropicPlasticity3DModifiedMohrCoulombTresca").def(py::init<>());

    py::class_< GenericFiniteStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericFiniteStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"FiniteStrainIsotropicPlasticity3DTrescaVonMises").def(py::init<>());

    py::class_< GenericFiniteStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>,
    typename GenericFiniteStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"FiniteStrainIsotropicPlasticity3DTrescaModifiedMohrCoulomb").def(py::init<>());

    py::class_< GenericFiniteStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<DruckerPragerPlasticPotential<6>>>>,
    typename GenericFiniteStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<DruckerPragerPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"FiniteStrainIsotropicPlasticity3DTrescaDruckerPrager").def(py::init<>());

    py::class_< GenericFiniteStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<TrescaPlasticPotential<6>>>>,
    typename GenericFiniteStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<TrescaPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"FiniteStrainIsotropicPlasticity3DTrescaTresca").def(py::init<>());

    py::class_< GenericFiniteStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericFiniteStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"FiniteStrainIsotropicPlasticity3DDruckerPragerVonMises").def(py::init<>());

    py::class_< GenericFiniteStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>,
    typename GenericFiniteStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"FiniteStrainIsotropicPlasticity3DDruckerPragerModifiedMohrCoulomb").def(py::init<>());

    py::class_< GenericFiniteStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<6>>>>,
    typename GenericFiniteStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"FiniteStrainIsotropicPlasticity3DDruckerPragerDruckerPrager").def(py::init<>());

    py::class_< GenericFiniteStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<TrescaPlasticPotential<6>>>>,
    typename GenericFiniteStrainIsotropicPlasticity <GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<TrescaPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"FiniteStrainIsotropicPlasticity3DDruckerPragerTresca").def(py::init<>());

    /* Finite strain */
    // Kinematic
    py::class_< GenericFiniteStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericFiniteStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"FiniteStrainKinematicPlasticity3DVonMisesVonMises").def(py::init<>());

    py::class_< GenericFiniteStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>,
    typename GenericFiniteStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"FiniteStrainKinematicPlasticity3DVonMisesModifiedMohrCoulomb").def(py::init<>());

    py::class_< GenericFiniteStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<DruckerPragerPlasticPotential<6>>>>,
    typename GenericFiniteStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<DruckerPragerPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"FiniteStrainKinematicPlasticity3DVonMisesDruckerPrager").def(py::init<>());

    py::class_< GenericFiniteStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<TrescaPlasticPotential<6>>>>,
    typename GenericFiniteStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<TrescaPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"FiniteStrainKinematicPlasticity3DVonMisesTresca").def(py::init<>());

    py::class_< GenericFiniteStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericFiniteStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"FiniteStrainKinematicPlasticity3DModifiedMohrCoulombVonMises").def(py::init<>());

    py::class_< GenericFiniteStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>,
    typename GenericFiniteStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"FiniteStrainKinematicPlasticity3DModifiedMohrCoulombModifiedMohrCoulomb").def(py::init<>());

    py::class_< GenericFiniteStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<ModifiedMohrCoulombYieldSurface<DruckerPragerPlasticPotential<6>>>>,
    typename GenericFiniteStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<ModifiedMohrCoulombYieldSurface<DruckerPragerPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"FiniteStrainKinematicPlasticity3DModifiedMohrCoulombDruckerPrager").def(py::init<>());

    py::class_< GenericFiniteStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<ModifiedMohrCoulombYieldSurface<TrescaPlasticPotential<6>>>>,
    typename GenericFiniteStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<ModifiedMohrCoulombYieldSurface<TrescaPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"FiniteStrainKinematicPlasticity3DModifiedMohrCoulombTresca").def(py::init<>());

    py::class_< GenericFiniteStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericFiniteStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"FiniteStrainKinematicPlasticity3DTrescaVonMises").def(py::init<>());

    py::class_< GenericFiniteStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>,
    typename GenericFiniteStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"FiniteStrainKinematicPlasticity3DTrescaModifiedMohrCoulomb").def(py::init<>());

    py::class_< GenericFiniteStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<DruckerPragerPlasticPotential<6>>>>,
    typename GenericFiniteStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<DruckerPragerPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"FiniteStrainKinematicPlasticity3DTrescaDruckerPrager").def(py::init<>());

    py::class_< GenericFiniteStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<TrescaPlasticPotential<6>>>>,
    typename GenericFiniteStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<TrescaPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"FiniteStrainKinematicPlasticity3DTrescaTresca").def(py::init<>());

    py::class_< GenericFiniteStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericFiniteStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"FiniteStrainKinematicPlasticity3DDruckerPragerVonMises").def(py::init<>());

    py::class_< GenericFiniteStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>,
    typename GenericFiniteStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"FiniteStrainKinematicPlasticity3DDruckerPragerModifiedMohrCoulomb").def(py::init<>());

    py::class_< GenericFiniteStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<6>>>>,
    typename GenericFiniteStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"FiniteStrainKinematicPlasticity3DDruckerPragerDruckerPrager").def(py::init<>());

    py::class_< GenericFiniteStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<TrescaPlasticPotential<6>>>>,
    typename GenericFiniteStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<TrescaPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"FiniteStrainKinematicPlasticity3DDruckerPragerTresca").def(py::init<>());


    // Damage
    /* Small strain */

    // damage 3D
    py::class_<  GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw > (m,"SmallStrainIsotropicDamage3DVonMises").def(py::init<>());

    py::class_<  GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>,
    typename GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw > (m,"SmallStrainIsotropicDamage3DModifiedMohrCoulomb").def(py::init<>());

    py::class_<  GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<TrescaPlasticPotential<6>>>>,
    typename GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<TrescaPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw > (m,"SmallStrainIsotropicDamage3DTresca").def(py::init<>());

    py::class_<  GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<6>>>>,
    typename GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw > (m,"SmallStrainIsotropicDamage3DDruckerPrager").def(py::init<>());

    py::class_<  GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw > (m,"SmallStrainIsotropicDamage3DSimoJu").def(py::init<>());

    py::class_<  GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<MohrCoulombYieldSurface<MohrCoulombPlasticPotential<6>>>>,
    typename GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<MohrCoulombYieldSurface<MohrCoulombPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw > (m,"SmallStrainIsotropicDamage3DMohrCoulomb").def(py::init<>());

    py::class_<  GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<RankinePlasticPotential<6>>>>,
    typename GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<RankinePlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw > (m,"SmallStrainIsotropicDamage3DRankine").def(py::init<>());

    // damage plane strain
    py::class_<  GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<VonMisesPlasticPotential<3>>>>,
    typename GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<VonMisesPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw > (m,"SmallStrainIsotropicDamagePlaneStrainVonMises").def(py::init<>());

    py::class_<  GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<3>>>>,
    typename GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw > (m,"SmallStrainIsotropicDamagePlaneStrainModifiedMohrCoulomb").def(py::init<>());

    py::class_<  GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<TrescaPlasticPotential<3>>>>,
    typename GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<TrescaPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw > (m,"SmallStrainIsotropicDamagePlaneStrainTresca").def(py::init<>());

    py::class_<  GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<3>>>>,
    typename GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw > (m,"SmallStrainIsotropicDamagePlaneStrainDruckerPrager").def(py::init<>());

    py::class_<  GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<VonMisesPlasticPotential<3>>>>,
    typename GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<VonMisesPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw > (m,"SmallStrainIsotropicDamagePlaneStrainSimoJu").def(py::init<>());

    py::class_<  GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<MohrCoulombYieldSurface<MohrCoulombPlasticPotential<3>>>>,
    typename GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<MohrCoulombYieldSurface<MohrCoulombPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw > (m,"SmallStrainIsotropicDamagePlaneStrainMohrCoulomb").def(py::init<>());

    py::class_<  GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<RankinePlasticPotential<3>>>>,
    typename GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<RankinePlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw > (m,"SmallStrainIsotropicDamagePlaneStrainRankine").def(py::init<>());

    // damage plane stress
    py::class_<  GenericSmallStrainIsotropicDamagePlaneStress <GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<VonMisesPlasticPotential<3>>>>,
    typename GenericSmallStrainIsotropicDamagePlaneStress <GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<VonMisesPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw > (m,"SmallStrainIsotropicDamagePlaneStressVonMises").def(py::init<>());

    py::class_<  GenericSmallStrainIsotropicDamagePlaneStress <GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<3>>>>,
    typename GenericSmallStrainIsotropicDamagePlaneStress <GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw > (m,"SmallStrainIsotropicDamagePlaneStressModifiedMohrCoulomb").def(py::init<>());

    py::class_<  GenericSmallStrainIsotropicDamagePlaneStress <GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<TrescaPlasticPotential<3>>>>,
    typename GenericSmallStrainIsotropicDamagePlaneStress <GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<TrescaPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw > (m,"SmallStrainIsotropicDamagePlaneStressTresca").def(py::init<>());

    py::class_<  GenericSmallStrainIsotropicDamagePlaneStress <GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<3>>>>,
    typename GenericSmallStrainIsotropicDamagePlaneStress <GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw > (m,"SmallStrainIsotropicDamagePlaneStressDruckerPrager").def(py::init<>());

    py::class_<  GenericSmallStrainIsotropicDamagePlaneStress <GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<VonMisesPlasticPotential<3>>>>,
    typename GenericSmallStrainIsotropicDamagePlaneStress <GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<VonMisesPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw > (m,"SmallStrainIsotropicDamagePlaneStressSimoJu").def(py::init<>());

    py::class_<  GenericSmallStrainIsotropicDamagePlaneStress <GenericConstitutiveLawIntegratorDamage<MohrCoulombYieldSurface<MohrCoulombPlasticPotential<3>>>>,
    typename GenericSmallStrainIsotropicDamagePlaneStress <GenericConstitutiveLawIntegratorDamage<MohrCoulombYieldSurface<MohrCoulombPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw > (m,"SmallStrainIsotropicDamagePlaneStressMohrCoulomb").def(py::init<>());

    py::class_<  GenericSmallStrainIsotropicDamagePlaneStress <GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<RankinePlasticPotential<3>>>>,
    typename GenericSmallStrainIsotropicDamagePlaneStress <GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<RankinePlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw > (m,"SmallStrainIsotropicDamagePlaneStressRankine").def(py::init<>());


    py::class_< GenericSmallStrainViscoplasticity3D, typename GenericSmallStrainViscoplasticity3D::Pointer,  ConstitutiveLaw  >
    (m,"GenericSmallStrainViscoplasticity3D").def(py::init<>())
    ;

    py::class_< GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw>(m, "SmallStrainDplusDminusDamageModifiedMohrCoulombModifiedMohrCoulomb3D").def(py::init<>())
    ;

    py::class_< GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<RankineYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<RankineYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw>(m, "SmallStrainDplusDminusDamageModifiedMohrCoulombRankine3D").def(py::init<>())
    ;

    py::class_< GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw>(m, "SmallStrainDplusDminusDamageModifiedMohrCoulombSimoJu3D").def(py::init<>())
    ;

    py::class_< GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw>(m, "SmallStrainDplusDminusDamageModifiedMohrCoulombVonMises3D").def(py::init<>())
    ;

    py::class_< GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<TrescaYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<TrescaYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw>(m, "SmallStrainDplusDminusDamageModifiedMohrCoulombTresca3D").def(py::init<>())
    ;

    py::class_< GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw>(m, "SmallStrainDplusDminusDamageModifiedMohrCoulombDruckerPrager3D").def(py::init<>())
    ;

    py::class_< GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<RankineYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<RankineYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw>(m, "SmallStrainDplusDminusDamageRankineModifiedMohrCoulomb3D").def(py::init<>())
    ;

    py::class_< GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<RankineYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<RankineYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<RankineYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<RankineYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw>(m, "SmallStrainDplusDminusDamageRankineRankine3D").def(py::init<>())
    ;

    py::class_< GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<RankineYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<RankineYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw>(m, "SmallStrainDplusDminusDamageRankineSimoJu3D").def(py::init<>())
    ;

    py::class_< GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<RankineYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<RankineYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw>(m, "SmallStrainDplusDminusDamageRankineVonMises3D").def(py::init<>())
    ;

    py::class_< GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<RankineYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<TrescaYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<RankineYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<TrescaYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw>(m, "SmallStrainDplusDminusDamageRankineTresca3D").def(py::init<>())
    ;

    py::class_< GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<RankineYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<RankineYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw>(m, "SmallStrainDplusDminusDamageRankineDruckerPrager3D").def(py::init<>())
    ;

    py::class_< GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw>(m, "SmallStrainDplusDminusDamageSimoJuModifiedMohrCoulomb3D").def(py::init<>())
    ;

    py::class_< GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<RankineYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<RankineYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw>(m, "SmallStrainDplusDminusDamageSimoJuRankine3D").def(py::init<>())
    ;

    py::class_< GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw>(m, "SmallStrainDplusDminusDamageSimoJuSimoJu3D").def(py::init<>())
    ;

    py::class_< GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw>(m, "SmallStrainDplusDminusDamageSimoJuVonMises3D").def(py::init<>())
    ;

    py::class_< GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<TrescaYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<TrescaYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw>(m, "SmallStrainDplusDminusDamageSimoJuTresca3D").def(py::init<>())
    ;

    py::class_< GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw>(m, "SmallStrainDplusDminusDamageSimoJuDruckerPrager3D").def(py::init<>())
    ;

    py::class_< GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw>(m, "SmallStrainDplusDminusDamageVonMisesModifiedMohrCoulomb3D").def(py::init<>())
    ;

    py::class_< GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<RankineYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<RankineYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw>(m, "SmallStrainDplusDminusDamageVonMisesRankine3D").def(py::init<>())
    ;

    py::class_< GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw>(m, "SmallStrainDplusDminusDamageVonMisesSimoJu3D").def(py::init<>())
    ;

    py::class_< GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw>(m, "SmallStrainDplusDminusDamageVonMisesVonMises3D").def(py::init<>())
    ;

    py::class_< GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<TrescaYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<TrescaYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw>(m, "SmallStrainDplusDminusDamageVonMisesTresca3D").def(py::init<>())
    ;

    py::class_< GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw>(m, "SmallStrainDplusDminusDamageVonMisesDruckerPrager3D").def(py::init<>())
    ;

    py::class_< GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<TrescaYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<TrescaYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw>(m, "SmallStrainDplusDminusDamageTrescaModifiedMohrCoulomb3D").def(py::init<>())
    ;

    py::class_< GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<TrescaYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<RankineYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<TrescaYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<RankineYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw>(m, "SmallStrainDplusDminusDamageTrescaRankine3D").def(py::init<>())
    ;

    py::class_< GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<TrescaYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<TrescaYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw>(m, "SmallStrainDplusDminusDamageTrescaSimoJu3D").def(py::init<>())
    ;

    py::class_< GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<TrescaYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<TrescaYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw>(m, "SmallStrainDplusDminusDamageTrescaVonMises3D").def(py::init<>())
    ;

    py::class_< GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<TrescaYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<TrescaYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<TrescaYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<TrescaYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw>(m, "SmallStrainDplusDminusDamageTrescaTresca3D").def(py::init<>())
    ;

    py::class_< GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<TrescaYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<TrescaYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw>(m, "SmallStrainDplusDminusDamageTrescaDruckerPrager3D").def(py::init<>())
    ;

    py::class_< GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw>(m, "SmallStrainDplusDminusDamageDruckerPragerModifiedMohrCoulomb3D").def(py::init<>())
    ;

    py::class_< GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<RankineYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<RankineYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw>(m, "SmallStrainDplusDminusDamageDruckerPragerRankine3D").def(py::init<>())
    ;

    py::class_< GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw>(m, "SmallStrainDplusDminusDamageDruckerPragerSimoJu3D").def(py::init<>())
    ;

    py::class_< GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw>(m, "SmallStrainDplusDminusDamageDruckerPragerVonMises3D").def(py::init<>())
    ;

    py::class_< GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<TrescaYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<TrescaYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw>(m, "SmallStrainDplusDminusDamageDruckerPragerTresca3D").def(py::init<>())
    ;

    py::class_< GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw>(m, "SmallStrainDplusDminusDamageDruckerPragerDruckerPrager3D").def(py::init<>())
    ;

    py::class_< GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw>(m, "SmallStrainDplusDminusDamageMohrCoulombMohrCoulomb3D").def(py::init<>())
    ;

    py::class_< GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<RankineYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<RankineYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw>(m, "SmallStrainDplusDminusDamageMohrCoulombRankine3D").def(py::init<>())
    ;

    py::class_< GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw>(m, "SmallStrainDplusDminusDamageMohrCoulombSimoJu3D").def(py::init<>())
    ;

    py::class_< GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw>(m, "SmallStrainDplusDminusDamageMohrCoulombVonMises3D").def(py::init<>())
    ;

    py::class_< GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<TrescaYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<TrescaYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw>(m, "SmallStrainDplusDminusDamageMohrCoulombTresca3D").def(py::init<>())
    ;

    py::class_< GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw>(m, "SmallStrainDplusDminusDamageMohrCoulombDruckerPrager3D").def(py::init<>())
    ;

	py::class_< GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<RankineYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<RankineYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw>(m, "SmallStrainDplusDminusDamageRankineMohrCoulomb3D").def(py::init<>())
    ;

    py::class_< GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw>(m, "SmallStrainDplusDminusDamageSimoJuMohrCoulomb3D").def(py::init<>())
    ;

	py::class_< GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw>(m, "SmallStrainDplusDminusDamageVonMisesMohrCoulomb3D").def(py::init<>())
    ;

    py::class_< GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<TrescaYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<TrescaYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw>(m, "SmallStrainDplusDminusDamageTrescaMohrCoulomb3D").def(py::init<>())
    ;

    py::class_< GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw>(m, "SmallStrainDplusDminusDamageDruckerPragerMohrCoulomb3D").def(py::init<>())
    ;

    py::class_< GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>>,
    typename GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw>(m, "SmallStrainDplusDminusDamageMohrCoulombMohrCoulomb2D").def(py::init<>())
    ;

    py::class_< GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<RankineYieldSurface<VonMisesPlasticPotential<3>>>>,
    typename GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<RankineYieldSurface<VonMisesPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw>(m, "SmallStrainDplusDminusDamageMohrCoulombRankine2D").def(py::init<>())
    ;

    py::class_< GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<SimoJuYieldSurface<VonMisesPlasticPotential<3>>>>,
    typename GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<SimoJuYieldSurface<VonMisesPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw>(m, "SmallStrainDplusDminusDamageMohrCoulombSimoJu2D").def(py::init<>())
    ;

    py::class_< GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<VonMisesYieldSurface<VonMisesPlasticPotential<3>>>>,
    typename GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<VonMisesYieldSurface<VonMisesPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw>(m, "SmallStrainDplusDminusDamageMohrCoulombVonMises2D").def(py::init<>())
    ;

    py::class_< GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<TrescaYieldSurface<VonMisesPlasticPotential<3>>>>,
    typename GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<TrescaYieldSurface<VonMisesPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw>(m, "SmallStrainDplusDminusDamageMohrCoulombTresca2D").def(py::init<>())
    ;

    py::class_< GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<3>>>>,
    typename GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw>(m, "SmallStrainDplusDminusDamageMohrCoulombDruckerPrager2D").def(py::init<>())
    ;

	py::class_< GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<RankineYieldSurface<VonMisesPlasticPotential<3>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>>,
    typename GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<RankineYieldSurface<VonMisesPlasticPotential<3>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw>(m, "SmallStrainDplusDminusDamageRankineMohrCoulomb2D").def(py::init<>())
    ;

    py::class_< GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<SimoJuYieldSurface<VonMisesPlasticPotential<3>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>>,
    typename GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<SimoJuYieldSurface<VonMisesPlasticPotential<3>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw>(m, "SmallStrainDplusDminusDamageSimoJuMohrCoulomb2D").def(py::init<>())
    ;

	py::class_< GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<VonMisesYieldSurface<VonMisesPlasticPotential<3>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>>,
    typename GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<VonMisesYieldSurface<VonMisesPlasticPotential<3>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw>(m, "SmallStrainDplusDminusDamageVonMisesMohrCoulomb2D").def(py::init<>())
    ;

    py::class_< GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<TrescaYieldSurface<VonMisesPlasticPotential<3>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>>,
    typename GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<TrescaYieldSurface<VonMisesPlasticPotential<3>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw>(m, "SmallStrainDplusDminusDamageTrescaMohrCoulomb2D").def(py::init<>())
    ;

    py::class_< GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<3>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>>,
    typename GenericSmallStrainDplusDminusDamage<GenericTensionConstitutiveLawIntegratorDplusDminusDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<3>>>, GenericCompressionConstitutiveLawIntegratorDplusDminusDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw>(m, "SmallStrainDplusDminusDamageDruckerPragerMohrCoulomb2D").def(py::init<>())
    ;

    py::class_< DamageDPlusDMinusMasonry2DLaw, typename DamageDPlusDMinusMasonry2DLaw::Pointer, ConstitutiveLaw >
    (m, "DamageDPlusDMinusPlaneStressMasonry2DLaw").def(py::init<>())
    ;

	py::class_< DamageDPlusDMinusMasonry3DLaw, typename DamageDPlusDMinusMasonry3DLaw::Pointer, ConstitutiveLaw >
	(m, "DamageDPlusDMinusMasonry3DLaw").def(py::init<>())
	;

    py::class_<GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw>
    (m,"SmallStrainOrthotropicDamageRankine3D").def(py::init<>());

    py::class_<GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw>
    (m,"SmallStrainOrthotropicDamageVonMises3D").def(py::init<>());

    py::class_<GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw>
    (m,"SmallStrainOrthotropicDamageDruckerPrager3D").def(py::init<>());

    py::class_<GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw>
    (m,"SmallStrainOrthotropicDamageTresca3D").def(py::init<>());

    py::class_<GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw>
    (m,"SmallStrainOrthotropicDamageMohrCoulomb3D").def(py::init<>());

    py::class_<GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw>
    (m,"SmallStrainOrthotropicDamageModifiedMohrCoulomb3D").def(py::init<>());

    py::class_<GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw>
    (m,"SmallStrainOrthotropicDamageSimoJu3D").def(py::init<>());

    py::class_<GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<VonMisesPlasticPotential<3>>>>,
    typename GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<VonMisesPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw>
    (m,"SmallStrainOrthotropicDamageRankine2D").def(py::init<>());

    py::class_<GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<VonMisesPlasticPotential<3>>>>,
    typename GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<VonMisesPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw>
    (m,"SmallStrainOrthotropicDamageVonMises2D").def(py::init<>());

    py::class_<GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<3>>>>,
    typename GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw>
    (m,"SmallStrainOrthotropicDamageDruckerPrager2D").def(py::init<>());

    py::class_<GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<VonMisesPlasticPotential<3>>>>,
    typename GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<VonMisesPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw>
    (m,"SmallStrainOrthotropicDamageTresca2D").def(py::init<>());

    py::class_<GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>>,
    typename GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw>
    (m,"SmallStrainOrthotropicDamageMohrCoulomb2D").def(py::init<>());

    py::class_<GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>>,
    typename GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw>
    (m,"SmallStrainOrthotropicDamageModifiedMohrCoulomb2D").def(py::init<>());

    py::class_<GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<VonMisesPlasticPotential<3>>>>,
    typename GenericSmallStrainOrthotropicDamage<GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<VonMisesPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw>
    (m,"SmallStrainOrthotropicDamageSimoJu2D").def(py::init<>());

    py::class_< ParallelRuleOfMixturesLaw<3>, typename ParallelRuleOfMixturesLaw<3>::Pointer,  ConstitutiveLaw  >
    (m,"ParallelRuleOfMixturesLaw3D").def(py::init<>())
    ;

    py::class_< ParallelRuleOfMixturesLaw<2>, typename ParallelRuleOfMixturesLaw<2>::Pointer,  ConstitutiveLaw  >
    (m,"ParallelRuleOfMixturesLaw2D").def(py::init<>())
    ;

    py::class_< SerialParallelRuleOfMixturesLaw, typename SerialParallelRuleOfMixturesLaw::Pointer,  ConstitutiveLaw  >
    (m,"SerialParallelRuleOfMixturesLaw").def(py::init<>())
    ;

    py::class_< AssociativePlasticDamageModel <VonMisesYieldSurface<VonMisesPlasticPotential<6>>>,
    typename AssociativePlasticDamageModel <VonMisesYieldSurface<VonMisesPlasticPotential<6>>>::Pointer,
    ConstitutiveLaw >
    (m,"AssociativePlasticDamageModel3DVonMises").def(py::init<>());

    py::class_< AssociativePlasticDamageModel <DruckerPragerYieldSurface<DruckerPragerPlasticPotential<6>>>,
    typename AssociativePlasticDamageModel <DruckerPragerYieldSurface<DruckerPragerPlasticPotential<6>>>::Pointer,
    ConstitutiveLaw >
    (m,"AssociativePlasticDamageModel3DDruckerPrager").def(py::init<>());

    py::class_< AssociativePlasticDamageModel <ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>,
    typename AssociativePlasticDamageModel <ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>::Pointer,
    ConstitutiveLaw >
    (m,"AssociativePlasticDamageModel3DModifiedMohrCoulomb").def(py::init<>());

    py::class_< AssociativePlasticDamageModel <RankineYieldSurface<RankinePlasticPotential<6>>>,
    typename AssociativePlasticDamageModel <RankineYieldSurface<RankinePlasticPotential<6>>>::Pointer,
    ConstitutiveLaw >
    (m,"AssociativePlasticDamageModel3DRankine").def(py::init<>());

    // Thermal CL's
    py::class_< ThermalElasticIsotropic3D, typename ThermalElasticIsotropic3D::Pointer, ConstitutiveLaw >
    (m,"ThermalElasticIsotropic3D").def(py::init<>());

    py::class_< ThermalLinearPlaneStrain, typename ThermalLinearPlaneStrain::Pointer, ConstitutiveLaw >
    (m,"ThermalLinearPlaneStrain").def(py::init<>());

    py::class_< ThermalLinearPlaneStress, typename ThermalLinearPlaneStress::Pointer, ConstitutiveLaw >
    (m,"ThermalLinearPlaneStress").def(py::init<>());


    // damage 3D
    py::class_<  GenericSmallStrainThermalIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ThermalVonMisesYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainThermalIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ThermalVonMisesYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw > (m,"SmallStrainThermalIsotropicDamage3DVonMises").def(py::init<>());

    py::class_<  GenericSmallStrainThermalIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ThermalModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainThermalIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ThermalModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw > (m,"SmallStrainThermalIsotropicDamage3DModifiedMohrCoulomb").def(py::init<>());

    py::class_<  GenericSmallStrainThermalIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ThermalTrescaYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainThermalIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ThermalTrescaYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw > (m,"SmallStrainThermalIsotropicDamage3DTresca").def(py::init<>());

    py::class_<  GenericSmallStrainThermalIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ThermalDruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainThermalIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ThermalDruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw > (m,"SmallStrainThermalIsotropicDamage3DDruckerPrager").def(py::init<>());

    py::class_<  GenericSmallStrainThermalIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ThermalSimoJuYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainThermalIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ThermalSimoJuYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw > (m,"SmallStrainThermalIsotropicDamage3DSimoJu").def(py::init<>());

    py::class_<  GenericSmallStrainThermalIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ThermalMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainThermalIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ThermalMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw > (m,"SmallStrainThermalIsotropicDamage3DMohrCoulomb").def(py::init<>());

    py::class_<  GenericSmallStrainThermalIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ThermalRankineYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainThermalIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ThermalRankineYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw > (m,"SmallStrainThermalIsotropicDamage3DRankine").def(py::init<>());

    // damage plane strain
    py::class_<  GenericSmallStrainThermalIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ThermalVonMisesYieldSurface<VonMisesPlasticPotential<3>>>>,
    typename GenericSmallStrainThermalIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ThermalVonMisesYieldSurface<VonMisesPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw > (m,"SmallStrainThermalIsotropicDamagePlaneStrainVonMises").def(py::init<>());

    py::class_<  GenericSmallStrainThermalIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ThermalModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>>,
    typename GenericSmallStrainThermalIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ThermalModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw > (m,"SmallStrainThermalIsotropicDamagePlaneStrainModifiedMohrCoulomb").def(py::init<>());

    py::class_<  GenericSmallStrainThermalIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ThermalTrescaYieldSurface<VonMisesPlasticPotential<3>>>>,
    typename GenericSmallStrainThermalIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ThermalTrescaYieldSurface<VonMisesPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw > (m,"SmallStrainThermalIsotropicDamagePlaneStrainTresca").def(py::init<>());

    py::class_<  GenericSmallStrainThermalIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ThermalDruckerPragerYieldSurface<VonMisesPlasticPotential<3>>>>,
    typename GenericSmallStrainThermalIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ThermalDruckerPragerYieldSurface<VonMisesPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw > (m,"SmallStrainThermalIsotropicDamagePlaneStrainDruckerPrager").def(py::init<>());

    py::class_<  GenericSmallStrainThermalIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ThermalSimoJuYieldSurface<VonMisesPlasticPotential<3>>>>,
    typename GenericSmallStrainThermalIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ThermalSimoJuYieldSurface<VonMisesPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw > (m,"SmallStrainThermalIsotropicDamagePlaneStrainSimoJu").def(py::init<>());

    py::class_<  GenericSmallStrainThermalIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ThermalMohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>>,
    typename GenericSmallStrainThermalIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ThermalMohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw > (m,"SmallStrainThermalIsotropicDamagePlaneStrainMohrCoulomb").def(py::init<>());

    py::class_<  GenericSmallStrainThermalIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ThermalRankineYieldSurface<VonMisesPlasticPotential<3>>>>,
    typename GenericSmallStrainThermalIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ThermalRankineYieldSurface<VonMisesPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw > (m,"SmallStrainThermalIsotropicDamagePlaneStrainRankine").def(py::init<>());

    py::class_<  GenericSmallStrainThermalIsotropicDamagePlaneStress <GenericConstitutiveLawIntegratorDamage<ThermalModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>>,
    typename GenericSmallStrainThermalIsotropicDamagePlaneStress <GenericConstitutiveLawIntegratorDamage<ThermalModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw > (m,"SmallStrainThermalIsotropicDamagePlaneStressModifiedMohrCoulomb").def(py::init<>());

    // Plane stress thermal damage
    py::class_<  GenericSmallStrainThermalIsotropicDamagePlaneStress <GenericConstitutiveLawIntegratorDamage<ThermalTrescaYieldSurface<VonMisesPlasticPotential<3>>>>,
    typename GenericSmallStrainThermalIsotropicDamagePlaneStress <GenericConstitutiveLawIntegratorDamage<ThermalTrescaYieldSurface<VonMisesPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw > (m,"SmallStrainThermalIsotropicDamagePlaneStressTresca").def(py::init<>());

    py::class_<  GenericSmallStrainThermalIsotropicDamagePlaneStress <GenericConstitutiveLawIntegratorDamage<ThermalDruckerPragerYieldSurface<VonMisesPlasticPotential<3>>>>,
    typename GenericSmallStrainThermalIsotropicDamagePlaneStress <GenericConstitutiveLawIntegratorDamage<ThermalDruckerPragerYieldSurface<VonMisesPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw > (m,"SmallStrainThermalIsotropicDamagePlaneStressDruckerPrager").def(py::init<>());

    py::class_<  GenericSmallStrainThermalIsotropicDamagePlaneStress <GenericConstitutiveLawIntegratorDamage<ThermalSimoJuYieldSurface<VonMisesPlasticPotential<3>>>>,
    typename GenericSmallStrainThermalIsotropicDamagePlaneStress <GenericConstitutiveLawIntegratorDamage<ThermalSimoJuYieldSurface<VonMisesPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw > (m,"SmallStrainThermalIsotropicDamagePlaneStressSimoJu").def(py::init<>());

    py::class_<  GenericSmallStrainThermalIsotropicDamagePlaneStress <GenericConstitutiveLawIntegratorDamage<ThermalMohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>>,
    typename GenericSmallStrainThermalIsotropicDamagePlaneStress <GenericConstitutiveLawIntegratorDamage<ThermalMohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw > (m,"SmallStrainThermalIsotropicDamagePlaneStressMohrCoulomb").def(py::init<>());

    py::class_<  GenericSmallStrainThermalIsotropicDamagePlaneStress <GenericConstitutiveLawIntegratorDamage<ThermalRankineYieldSurface<VonMisesPlasticPotential<3>>>>,
    typename GenericSmallStrainThermalIsotropicDamagePlaneStress <GenericConstitutiveLawIntegratorDamage<ThermalRankineYieldSurface<VonMisesPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw > (m,"SmallStrainThermalIsotropicDamagePlaneStressRankine").def(py::init<>());

}

}  // namespace Kratos::Python.
