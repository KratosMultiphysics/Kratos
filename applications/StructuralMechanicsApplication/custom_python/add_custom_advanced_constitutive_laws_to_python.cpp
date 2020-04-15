// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Riccardo Rossi
//


// System includes

// External includes

// Project includes
#include "custom_python/add_custom_advanced_constitutive_laws_to_python.h"

// Elastic laws
#include "custom_advanced_constitutive/hyper_elastic_isotropic_ogden_1d.h"
#include "custom_advanced_constitutive/hyper_elastic_isotropic_henky_1d.h"
#include "custom_advanced_constitutive/truss_plasticity_constitutive_law.h"
#include "custom_advanced_constitutive/elastic_isotropic_plane_stress_uncoupled_shear.h"
#include "custom_advanced_constitutive/hyper_elastic_isotropic_kirchhoff_3d.h"
#include "custom_advanced_constitutive/hyper_elastic_isotropic_kirchhoff_plane_stress_2d.h"
#include "custom_advanced_constitutive/hyper_elastic_isotropic_kirchhoff_plane_strain_2d.h"
#include "custom_advanced_constitutive/hyper_elastic_isotropic_neo_hookean_3d.h"
#include "custom_advanced_constitutive/hyper_elastic_isotropic_neo_hookean_plane_strain_2d.h"
#include "custom_advanced_constitutive/linear_elastic_orthotropic_2D_law.h"
#include "custom_advanced_constitutive/small_strain_j2_plasticity_3d.h"
#include "custom_advanced_constitutive/small_strain_j2_plasticity_plane_strain_2d.h"
#include "custom_advanced_constitutive/small_strain_isotropic_damage_3d.h"
#include "custom_advanced_constitutive/small_strain_isotropic_damage_traction_only_3d.h"
#include "custom_advanced_constitutive/wrinkling_linear_2d_law.h"

// Plastic, damage laws and viscosities
#include "custom_advanced_constitutive/small_strain_isotropic_damage_plane_strain_2d.h"
#include "custom_advanced_constitutive/small_strain_isotropic_plasticity_factory.h"
#include "custom_advanced_constitutive/small_strain_kinematic_plasticity_factory.h"
#include "custom_advanced_constitutive/finite_strain_isotropic_plasticity_factory.h"
#include "custom_advanced_constitutive/small_strain_isotropic_damage_factory.h"
#include "custom_advanced_constitutive/viscous_generalized_maxwell.h"
#include "custom_advanced_constitutive/viscous_generalized_kelvin.h"
#include "custom_advanced_constitutive/generic_small_strain_viscoplasticity_3d.h"
#include "custom_advanced_constitutive/generic_small_strain_isotropic_plasticity.h"
#include "custom_advanced_constitutive/generic_finite_strain_isotropic_plasticity.h"
#include "custom_advanced_constitutive/generic_finite_strain_kinematic_plasticity.h"
#include "custom_advanced_constitutive/generic_small_strain_isotropic_damage.h"
#include "custom_advanced_constitutive/generic_small_strain_d_plus_d_minus_damage.h"
#include "custom_advanced_constitutive/generic_small_strain_kinematic_plasticity.h"
#include "custom_advanced_constitutive/generic_small_strain_high_cycle_fatigue_law.h"
#include "custom_advanced_constitutive/plasticity_isotropic_kinematic_j2.h"
#include "custom_advanced_constitutive/plane_stress_d_plus_d_minus_damage_masonry_2d.h"
#include "custom_advanced_constitutive/d_plus_d_minus_damage_masonry_3d.h"
#include "custom_advanced_constitutive/generic_small_strain_plastic_damage_model.h"
#include "custom_advanced_constitutive/generic_small_strain_orthotropic_damage.h"
#include "custom_advanced_constitutive/serial_parallel_rule_of_mixtures_law.h"

// Integrators
#include "custom_advanced_constitutive/constitutive_laws_integrators/generic_constitutive_law_integrator_damage.h"
#include "custom_advanced_constitutive/constitutive_laws_integrators/generic_constitutive_law_integrator_plasticity.h"
#include "custom_advanced_constitutive/constitutive_laws_integrators/generic_finite_strain_constitutive_law_integrator_plasticity.h"
#include "custom_advanced_constitutive/constitutive_laws_integrators/generic_finite_strain_constitutive_law_integrator_kinematic_plasticity.h"
#include "custom_advanced_constitutive/constitutive_laws_integrators/generic_constitutive_law_integrator_kinematic_plasticity.h"
#include "custom_advanced_constitutive/constitutive_laws_integrators/d+d-constitutive_law_integrators/generic_compression_constitutive_law_integrator.h"
#include "custom_advanced_constitutive/constitutive_laws_integrators/d+d-constitutive_law_integrators/generic_tension_constitutive_law_integrator.h"

// Yield surfaces
#include "custom_advanced_constitutive/yield_surfaces/generic_yield_surface.h"
#include "custom_advanced_constitutive/yield_surfaces/von_mises_yield_surface.h"
#include "custom_advanced_constitutive/yield_surfaces/modified_mohr_coulomb_yield_surface.h"
#include "custom_advanced_constitutive/yield_surfaces/mohr_coulomb_yield_surface.h"
#include "custom_advanced_constitutive/yield_surfaces/rankine_yield_surface.h"
#include "custom_advanced_constitutive/yield_surfaces/simo_ju_yield_surface.h"
#include "custom_advanced_constitutive/yield_surfaces/drucker_prager_yield_surface.h"
#include "custom_advanced_constitutive/yield_surfaces/tresca_yield_surface.h"

// Plastic potentials
#include "custom_advanced_constitutive/plastic_potentials/generic_plastic_potential.h"
#include "custom_advanced_constitutive/plastic_potentials/von_mises_plastic_potential.h"
#include "custom_advanced_constitutive/plastic_potentials/tresca_plastic_potential.h"
#include "custom_advanced_constitutive/plastic_potentials/modified_mohr_coulomb_plastic_potential.h"
#include "custom_advanced_constitutive/plastic_potentials/mohr_coulomb_plastic_potential.h"
#include "custom_advanced_constitutive/plastic_potentials/drucker_prager_plastic_potential.h"

// Rules of mixtures
#include "custom_advanced_constitutive/rule_of_mixtures_law.h"

namespace Kratos {
namespace Python {

void  AddCustomAdvancedConstitutiveLawsToPython(pybind11::module& m)
{
    namespace py = pybind11;

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

    py::class_< SmallStrainIsotropicDamageTractionOnly3D, typename SmallStrainIsotropicDamageTractionOnly3D::Pointer,  ConstitutiveLaw  >
    (m,"SmallStrainIsotropicDamageTractionOnly3DLaw").def(py::init<>())
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
    ConstitutiveLaw >
    (m,"SmallStrainHighCycleFatigue3DLawVonMisesVonMises").def(py::init<>());

    py::class_<  GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>,
    typename GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainHighCycleFatigue3DLawVonMisesModifiedMohrCoulomb").def(py::init<>());

    py::class_<  GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<DruckerPragerPlasticPotential<6>>>>,
    typename GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<DruckerPragerPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainHighCycleFatigue3DLawVonMisesDruckerPrager").def(py::init<>());

    py::class_<  GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<TrescaPlasticPotential<6>>>>,
    typename GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<TrescaPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainHighCycleFatigue3DLawVonMisesTresca").def(py::init<>());

    py::class_<  GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainHighCycleFatigue3DLawModifiedMohrCoulombVonMises").def(py::init<>());

    py::class_<  GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>,
    typename GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainHighCycleFatigue3DLawModifiedMohrCoulombModifiedMohrCoulomb").def(py::init<>());

    py::class_<  GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<DruckerPragerPlasticPotential<6>>>>,
    typename GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<DruckerPragerPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainHighCycleFatigue3DLawModifiedMohrCoulombDruckerPrager").def(py::init<>());

    py::class_<  GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<TrescaPlasticPotential<6>>>>,
    typename GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<TrescaPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainHighCycleFatigue3DLawModifiedMohrCoulombTresca").def(py::init<>());

    py::class_<  GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainHighCycleFatigue3DLawTrescaVonMises").def(py::init<>());

    py::class_<  GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>,
    typename GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainHighCycleFatigue3DLawTrescaModifiedMohrCoulomb").def(py::init<>());

    py::class_<  GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<DruckerPragerPlasticPotential<6>>>>,
    typename GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<DruckerPragerPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainHighCycleFatigue3DLawTrescaDruckerPrager").def(py::init<>());

    py::class_<  GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<TrescaPlasticPotential<6>>>>,
    typename GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<TrescaPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainHighCycleFatigue3DLawTrescaTresca").def(py::init<>());

    py::class_<  GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainHighCycleFatigue3DLawDruckerPragerVonMises").def(py::init<>());

    py::class_<  GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>,
    typename GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainHighCycleFatigue3DLawDruckerPragerModifiedMohrCoulomb").def(py::init<>());

    py::class_<  GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<6>>>>,
    typename GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainHighCycleFatigue3DLawDruckerPragerDruckerPrager").def(py::init<>());

    py::class_<  GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<TrescaPlasticPotential<6>>>>,
    typename GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<TrescaPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainHighCycleFatigue3DLawDruckerPragerTresca").def(py::init<>());

    py::class_<  GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainHighCycleFatigue3DLawRankineVonMises").def(py::init<>());

    py::class_<  GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>,
    typename GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainHighCycleFatigue3DLawRankineModifiedMohrCoulomb").def(py::init<>());

    py::class_<  GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<DruckerPragerPlasticPotential<6>>>>,
    typename GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<DruckerPragerPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainHighCycleFatigue3DLawRankineDruckerPrager").def(py::init<>());

    py::class_<  GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<TrescaPlasticPotential<6>>>>,
    typename GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<TrescaPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainHighCycleFatigue3DLawRankineTresca").def(py::init<>());

    py::class_<  GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainHighCycleFatigue3DLawSimoJuVonMises").def(py::init<>());

    py::class_<  GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>,
    typename GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainHighCycleFatigue3DLawSimoJuModifiedMohrCoulomb").def(py::init<>());

    py::class_<  GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<DruckerPragerPlasticPotential<6>>>>,
    typename GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<DruckerPragerPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainHighCycleFatigue3DLawSimoJuDruckerPrager").def(py::init<>());

    py::class_<  GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<TrescaPlasticPotential<6>>>>,
    typename GenericSmallStrainHighCycleFatigueLaw <GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<TrescaPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainHighCycleFatigue3DLawSimoJuTresca").def(py::init<>());

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
    // Kirchhoff
    py::class_< GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticIsotropicKirchhoffPlasticity3DVonMisesVonMises").def(py::init<>());

    py::class_< GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>,
    typename GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticIsotropicKirchhoffPlasticity3DVonMisesModifiedMohrCoulomb").def(py::init<>());

    py::class_< GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<DruckerPragerPlasticPotential<6>>>>,
    typename GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<DruckerPragerPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticIsotropicKirchhoffPlasticity3DVonMisesDruckerPrager").def(py::init<>());

    py::class_< GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<TrescaPlasticPotential<6>>>>,
    typename GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<TrescaPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticIsotropicKirchhoffPlasticity3DVonMisesTresca").def(py::init<>());

    py::class_< GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticIsotropicKirchhoffPlasticity3DModifiedMohrCoulombVonMises").def(py::init<>());

    py::class_< GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>,
    typename GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticIsotropicKirchhoffPlasticity3DModifiedMohrCoulombModifiedMohrCoulomb").def(py::init<>());

    py::class_< GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<DruckerPragerPlasticPotential<6>>>>,
    typename GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<DruckerPragerPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticIsotropicKirchhoffPlasticity3DModifiedMohrCoulombDruckerPrager").def(py::init<>());

    py::class_< GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<TrescaPlasticPotential<6>>>>,
    typename GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<TrescaPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticIsotropicKirchhoffPlasticity3DModifiedMohrCoulombTresca").def(py::init<>());

    py::class_< GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticIsotropicKirchhoffPlasticity3DTrescaVonMises").def(py::init<>());

    py::class_< GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>,
    typename GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticIsotropicKirchhoffPlasticity3DTrescaModifiedMohrCoulomb").def(py::init<>());

    py::class_< GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<DruckerPragerPlasticPotential<6>>>>,
    typename GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<DruckerPragerPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticIsotropicKirchhoffPlasticity3DTrescaDruckerPrager").def(py::init<>());

    py::class_< GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<TrescaPlasticPotential<6>>>>,
    typename GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<TrescaPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticIsotropicKirchhoffPlasticity3DTrescaTresca").def(py::init<>());

    py::class_< GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticIsotropicKirchhoffPlasticity3DDruckerPragerVonMises").def(py::init<>());

    py::class_< GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>,
    typename GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticIsotropicKirchhoffPlasticity3DDruckerPragerModifiedMohrCoulomb").def(py::init<>());

    py::class_< GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<6>>>>,
    typename GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticIsotropicKirchhoffPlasticity3DDruckerPragerDruckerPrager").def(py::init<>());

    py::class_< GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<TrescaPlasticPotential<6>>>>,
    typename GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<TrescaPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticIsotropicKirchhoffPlasticity3DDruckerPragerTresca").def(py::init<>());

    // Neo-Hookean
    py::class_< GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticIsotropicNeoHookeanPlasticity3DVonMisesVonMises").def(py::init<>());

    py::class_< GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>,
    typename GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticIsotropicNeoHookeanPlasticity3DVonMisesModifiedMohrCoulomb").def(py::init<>());

    py::class_< GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<DruckerPragerPlasticPotential<6>>>>,
    typename GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<DruckerPragerPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticIsotropicNeoHookeanPlasticity3DVonMisesDruckerPrager").def(py::init<>());

    py::class_< GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<TrescaPlasticPotential<6>>>>,
    typename GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<TrescaPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticIsotropicNeoHookeanPlasticity3DVonMisesTresca").def(py::init<>());

    py::class_< GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticIsotropicNeoHookeanPlasticity3DModifiedMohrCoulombVonMises").def(py::init<>());

    py::class_< GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>,
    typename GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticIsotropicNeoHookeanPlasticity3DModifiedMohrCoulombModifiedMohrCoulomb").def(py::init<>());

    py::class_< GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<DruckerPragerPlasticPotential<6>>>>,
    typename GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<DruckerPragerPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticIsotropicNeoHookeanPlasticity3DModifiedMohrCoulombDruckerPrager").def(py::init<>());

    py::class_< GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<TrescaPlasticPotential<6>>>>,
    typename GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<TrescaPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticIsotropicNeoHookeanPlasticity3DModifiedMohrCoulombTresca").def(py::init<>());

    py::class_< GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticIsotropicNeoHookeanPlasticity3DTrescaVonMises").def(py::init<>());

    py::class_< GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>,
    typename GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticIsotropicNeoHookeanPlasticity3DTrescaModifiedMohrCoulomb").def(py::init<>());

    py::class_< GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<DruckerPragerPlasticPotential<6>>>>,
    typename GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<DruckerPragerPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticIsotropicNeoHookeanPlasticity3DTrescaDruckerPrager").def(py::init<>());

    py::class_< GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<TrescaPlasticPotential<6>>>>,
    typename GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<TrescaPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticIsotropicNeoHookeanPlasticity3DTrescaTresca").def(py::init<>());

    py::class_< GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticIsotropicNeoHookeanPlasticity3DDruckerPragerVonMises").def(py::init<>());

    py::class_< GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>,
    typename GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticIsotropicNeoHookeanPlasticity3DDruckerPragerModifiedMohrCoulomb").def(py::init<>());

    py::class_< GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<6>>>>,
    typename GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticIsotropicNeoHookeanPlasticity3DDruckerPragerDruckerPrager").def(py::init<>());

    py::class_< GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<TrescaPlasticPotential<6>>>>,
    typename GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<TrescaPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticIsotropicNeoHookeanPlasticity3DDruckerPragerTresca").def(py::init<>());

    py::class_< GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<MohrCoulombPlasticPotential<6>>>>,
    typename GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<MohrCoulombPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticIsotropicKirchhoffPlasticity3DVonMisesMohrCoulomb").def(py::init<>());

	py::class_< GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticIsotropicKirchhoffPlasticity3DMohrCoulombVonMises").def(py::init<>());

    py::class_< GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<MohrCoulombPlasticPotential<6>>>>,
    typename GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<MohrCoulombPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticIsotropicKirchhoffPlasticity3DMohrCoulombMohrCoulomb").def(py::init<>());

    py::class_< GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<DruckerPragerPlasticPotential<6>>>>,
    typename GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<DruckerPragerPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticIsotropicKirchhoffPlasticity3DMohrCoulombDruckerPrager").def(py::init<>());

    py::class_< GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<TrescaPlasticPotential<6>>>>,
    typename GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<TrescaPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticIsotropicKirchhoffPlasticity3DMohrCoulombTresca").def(py::init<>());

    py::class_< GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<MohrCoulombPlasticPotential<6>>>>,
    typename GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<MohrCoulombPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticIsotropicKirchhoffPlasticity3DTrescaMohrCoulomb").def(py::init<>());

    py::class_< GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<MohrCoulombPlasticPotential<6>>>>,
    typename GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<MohrCoulombPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticIsotropicKirchhoffPlasticity3DDruckerPragerMohrCoulomb").def(py::init<>());

    py::class_< GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<MohrCoulombPlasticPotential<6>>>>,
    typename GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<MohrCoulombPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticIsotropicNeoHookeanPlasticity3DVonMisesMohrCoulomb").def(py::init<>());

    py::class_< GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticIsotropicNeoHookeanPlasticity3DMohrCoulombVonMises").def(py::init<>());

    py::class_< GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<MohrCoulombPlasticPotential<6>>>>,
    typename GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<MohrCoulombPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticIsotropicNeoHookeanPlasticity3DMohrCoulombMohrCoulomb").def(py::init<>());

    py::class_< GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<DruckerPragerPlasticPotential<6>>>>,
    typename GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<DruckerPragerPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticIsotropicNeoHookeanPlasticity3DMohrCoulombDruckerPrager").def(py::init<>());

    py::class_< GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<TrescaPlasticPotential<6>>>>,
    typename GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<TrescaPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticIsotropicNeoHookeanPlasticity3DMohrCoulombTresca").def(py::init<>());

    py::class_< GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<MohrCoulombPlasticPotential<6>>>>,
    typename GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<MohrCoulombPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticIsotropicNeoHookeanPlasticity3DTrescaMohrCoulomb").def(py::init<>());

    py::class_< GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<MohrCoulombPlasticPotential<6>>>>,
    typename GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<MohrCoulombPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticIsotropicNeoHookeanPlasticity3DDruckerPragerMohrCoulomb").def(py::init<>());

    /* Finite strain */
    // Kirchhoff
    py::class_< GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticKirchhoffKinematicPlasticity3DVonMisesVonMises").def(py::init<>());

    py::class_< GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>,
    typename GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticKirchhoffKinematicPlasticity3DVonMisesModifiedMohrCoulomb").def(py::init<>());

    py::class_< GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<DruckerPragerPlasticPotential<6>>>>,
    typename GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<DruckerPragerPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticKirchhoffKinematicPlasticity3DVonMisesDruckerPrager").def(py::init<>());

    py::class_< GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<TrescaPlasticPotential<6>>>>,
    typename GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<TrescaPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticKirchhoffKinematicPlasticity3DVonMisesTresca").def(py::init<>());

    py::class_< GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticKirchhoffKinematicPlasticity3DModifiedMohrCoulombVonMises").def(py::init<>());

    py::class_< GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>,
    typename GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticKirchhoffKinematicPlasticity3DModifiedMohrCoulombModifiedMohrCoulomb").def(py::init<>());

    py::class_< GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<ModifiedMohrCoulombYieldSurface<DruckerPragerPlasticPotential<6>>>>,
    typename GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<ModifiedMohrCoulombYieldSurface<DruckerPragerPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticKirchhoffKinematicPlasticity3DModifiedMohrCoulombDruckerPrager").def(py::init<>());

    py::class_< GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<ModifiedMohrCoulombYieldSurface<TrescaPlasticPotential<6>>>>,
    typename GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<ModifiedMohrCoulombYieldSurface<TrescaPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticKirchhoffKinematicPlasticity3DModifiedMohrCoulombTresca").def(py::init<>());

    py::class_< GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticKirchhoffKinematicPlasticity3DTrescaVonMises").def(py::init<>());

    py::class_< GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>,
    typename GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticKirchhoffKinematicPlasticity3DTrescaModifiedMohrCoulomb").def(py::init<>());

    py::class_< GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<DruckerPragerPlasticPotential<6>>>>,
    typename GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<DruckerPragerPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticKirchhoffKinematicPlasticity3DTrescaDruckerPrager").def(py::init<>());

    py::class_< GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<TrescaPlasticPotential<6>>>>,
    typename GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<TrescaPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticKirchhoffKinematicPlasticity3DTrescaTresca").def(py::init<>());

    py::class_< GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticKirchhoffKinematicPlasticity3DDruckerPragerVonMises").def(py::init<>());

    py::class_< GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>,
    typename GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticKirchhoffKinematicPlasticity3DDruckerPragerModifiedMohrCoulomb").def(py::init<>());

    py::class_< GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<6>>>>,
    typename GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticKirchhoffKinematicPlasticity3DDruckerPragerDruckerPrager").def(py::init<>());

    py::class_< GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<TrescaPlasticPotential<6>>>>,
    typename GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<TrescaPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticKirchhoffKinematicPlasticity3DDruckerPragerTresca").def(py::init<>());

    // Neo-Hookean
    py::class_< GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticNeoHookeanKinematicPlasticity3DVonMisesVonMises").def(py::init<>());

    py::class_< GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>,
    typename GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticNeoHookeanKinematicPlasticity3DVonMisesModifiedMohrCoulomb").def(py::init<>());

    py::class_< GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<DruckerPragerPlasticPotential<6>>>>,
    typename GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<DruckerPragerPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticNeoHookeanKinematicPlasticity3DVonMisesDruckerPrager").def(py::init<>());

    py::class_< GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<TrescaPlasticPotential<6>>>>,
    typename GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<TrescaPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticNeoHookeanKinematicPlasticity3DVonMisesTresca").def(py::init<>());

    py::class_< GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticNeoHookeanKinematicPlasticity3DModifiedMohrCoulombVonMises").def(py::init<>());

    py::class_< GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>,
    typename GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticNeoHookeanKinematicPlasticity3DModifiedMohrCoulombModifiedMohrCoulomb").def(py::init<>());

    py::class_< GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<ModifiedMohrCoulombYieldSurface<DruckerPragerPlasticPotential<6>>>>,
    typename GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<ModifiedMohrCoulombYieldSurface<DruckerPragerPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticNeoHookeanKinematicPlasticity3DModifiedMohrCoulombDruckerPrager").def(py::init<>());

    py::class_< GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<ModifiedMohrCoulombYieldSurface<TrescaPlasticPotential<6>>>>,
    typename GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<ModifiedMohrCoulombYieldSurface<TrescaPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticNeoHookeanKinematicPlasticity3DModifiedMohrCoulombTresca").def(py::init<>());

    py::class_< GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticNeoHookeanKinematicPlasticity3DTrescaVonMises").def(py::init<>());

    py::class_< GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>,
    typename GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticNeoHookeanKinematicPlasticity3DTrescaModifiedMohrCoulomb").def(py::init<>());

    py::class_< GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<DruckerPragerPlasticPotential<6>>>>,
    typename GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<DruckerPragerPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticNeoHookeanKinematicPlasticity3DTrescaDruckerPrager").def(py::init<>());

    py::class_< GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<TrescaPlasticPotential<6>>>>,
    typename GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<TrescaPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticNeoHookeanKinematicPlasticity3DTrescaTresca").def(py::init<>());

    py::class_< GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticNeoHookeanKinematicPlasticity3DDruckerPragerVonMises").def(py::init<>());

    py::class_< GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>,
    typename GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticNeoHookeanKinematicPlasticity3DDruckerPragerModifiedMohrCoulomb").def(py::init<>());

    py::class_< GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<6>>>>,
    typename GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticNeoHookeanKinematicPlasticity3DDruckerPragerDruckerPrager").def(py::init<>());

    py::class_< GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<TrescaPlasticPotential<6>>>>,
    typename GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<TrescaPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticNeoHookeanKinematicPlasticity3DDruckerPragerTresca").def(py::init<>());

    py::class_< GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<MohrCoulombPlasticPotential<6>>>>,
    typename GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<MohrCoulombPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticKirchhoffKinematicPlasticity3DVonMisesMohrCoulomb").def(py::init<>());

	py::class_< GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<MohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<MohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticKirchhoffKinematicPlasticity3DMohrCoulombVonMises").def(py::init<>());

    py::class_< GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<MohrCoulombYieldSurface<MohrCoulombPlasticPotential<6>>>>,
    typename GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<MohrCoulombYieldSurface<MohrCoulombPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticKirchhoffKinematicPlasticity3DMohrCoulombMohrCoulomb").def(py::init<>());

    py::class_< GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<MohrCoulombYieldSurface<DruckerPragerPlasticPotential<6>>>>,
    typename GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<MohrCoulombYieldSurface<DruckerPragerPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticKirchhoffKinematicPlasticity3DMohrCoulombDruckerPrager").def(py::init<>());

    py::class_< GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<MohrCoulombYieldSurface<TrescaPlasticPotential<6>>>>,
    typename GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<MohrCoulombYieldSurface<TrescaPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticKirchhoffKinematicPlasticity3DMohrCoulombTresca").def(py::init<>());

    py::class_< GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<MohrCoulombPlasticPotential<6>>>>,
    typename GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<MohrCoulombPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticKirchhoffKinematicPlasticity3DTrescaMohrCoulomb").def(py::init<>());

    py::class_< GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<MohrCoulombPlasticPotential<6>>>>,
    typename GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<MohrCoulombPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticKirchhoffKinematicPlasticity3DDruckerPragerMohrCoulomb").def(py::init<>());

    py::class_< GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<MohrCoulombPlasticPotential<6>>>>,
    typename GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<MohrCoulombPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticNeoHookeanKinematicPlasticity3DVonMisesMohrCoulomb").def(py::init<>());

    py::class_< GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<MohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<MohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticNeoHookeanKinematicPlasticity3DMohrCoulombVonMises").def(py::init<>());

    py::class_< GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<MohrCoulombYieldSurface<MohrCoulombPlasticPotential<6>>>>,
    typename GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<MohrCoulombYieldSurface<MohrCoulombPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticNeoHookeanKinematicPlasticity3DMohrCoulombMohrCoulomb").def(py::init<>());

    py::class_< GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<MohrCoulombYieldSurface<DruckerPragerPlasticPotential<6>>>>,
    typename GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<MohrCoulombYieldSurface<DruckerPragerPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticNeoHookeanKinematicPlasticity3DMohrCoulombDruckerPrager").def(py::init<>());

    py::class_< GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<MohrCoulombYieldSurface<TrescaPlasticPotential<6>>>>,
    typename GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<MohrCoulombYieldSurface<TrescaPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticNeoHookeanKinematicPlasticity3DMohrCoulombTresca").def(py::init<>());

    py::class_< GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<MohrCoulombPlasticPotential<6>>>>,
    typename GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<MohrCoulombPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticNeoHookeanKinematicPlasticity3DTrescaMohrCoulomb").def(py::init<>());

    py::class_< GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<MohrCoulombPlasticPotential<6>>>>,
    typename GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<MohrCoulombPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"HyperElasticNeoHookeanKinematicPlasticity3DDruckerPragerMohrCoulomb").def(py::init<>());

    // Damage
    /* Small strain */
    py::class_<  GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage3DVonMisesVonMises").def(py::init<>());

    py::class_<  GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>,
    typename GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage3DVonMisesModifiedMohrCoulomb").def(py::init<>());

    py::class_<  GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<DruckerPragerPlasticPotential<6>>>>,
    typename GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<DruckerPragerPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage3DVonMisesDruckerPrager").def(py::init<>());

    py::class_<  GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<TrescaPlasticPotential<6>>>>,
    typename GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<TrescaPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage3DVonMisesTresca").def(py::init<>());

    py::class_<  GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage3DModifiedMohrCoulombVonMises").def(py::init<>());

    py::class_<  GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>,
    typename GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage3DModifiedMohrCoulombModifiedMohrCoulomb").def(py::init<>());

    py::class_<  GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<DruckerPragerPlasticPotential<6>>>>,
    typename GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<DruckerPragerPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage3DModifiedMohrCoulombDruckerPrager").def(py::init<>());

    py::class_<  GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<TrescaPlasticPotential<6>>>>,
    typename GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<TrescaPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage3DModifiedMohrCoulombTresca").def(py::init<>());

    py::class_<  GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage3DTrescaVonMises").def(py::init<>());

    py::class_<  GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>,
    typename GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage3DTrescaModifiedMohrCoulomb").def(py::init<>());

    py::class_<  GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<DruckerPragerPlasticPotential<6>>>>,
    typename GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<DruckerPragerPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage3DTrescaDruckerPrager").def(py::init<>());

    py::class_<  GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<TrescaPlasticPotential<6>>>>,
    typename GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<TrescaPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage3DTrescaTresca").def(py::init<>());

    py::class_<  GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage3DDruckerPragerVonMises").def(py::init<>());

    py::class_<  GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>,
    typename GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage3DDruckerPragerModifiedMohrCoulomb").def(py::init<>());

    py::class_<  GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<6>>>>,
    typename GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage3DDruckerPragerDruckerPrager").def(py::init<>());

    py::class_<  GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<TrescaPlasticPotential<6>>>>,
    typename GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<TrescaPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage3DDruckerPragerTresca").def(py::init<>());

    py::class_<  GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage3DRankineVonMises").def(py::init<>());

    py::class_<  GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>,
    typename GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage3DRankineModifiedMohrCoulomb").def(py::init<>());

    py::class_<  GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<DruckerPragerPlasticPotential<6>>>>,
    typename GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<DruckerPragerPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage3DRankineDruckerPrager").def(py::init<>());

    py::class_<  GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<TrescaPlasticPotential<6>>>>,
    typename GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<TrescaPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage3DRankineTresca").def(py::init<>());

    py::class_<  GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage3DSimoJuVonMises").def(py::init<>());

    py::class_<  GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>,
    typename GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage3DSimoJuModifiedMohrCoulomb").def(py::init<>());

    py::class_<  GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<DruckerPragerPlasticPotential<6>>>>,
    typename GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<DruckerPragerPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage3DSimoJuDruckerPrager").def(py::init<>());

    py::class_<  GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<TrescaPlasticPotential<6>>>>,
    typename GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<TrescaPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage3DSimoJuTresca").def(py::init<>());

    py::class_<  GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<VonMisesPlasticPotential<3>>>>,
    typename GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<VonMisesPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage2DVonMisesVonMises").def(py::init<>());

    py::class_<  GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<ModifiedMohrCoulombPlasticPotential<3>>>>,
    typename GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<ModifiedMohrCoulombPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage2DVonMisesModifiedMohrCoulomb").def(py::init<>());

    py::class_<  GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<DruckerPragerPlasticPotential<3>>>>,
    typename GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<DruckerPragerPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage2DVonMisesDruckerPrager").def(py::init<>());

    py::class_<  GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<TrescaPlasticPotential<3>>>>,
    typename GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<TrescaPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage2DVonMisesTresca").def(py::init<>());

    py::class_<  GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>>,
    typename GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage2DModifiedMohrCoulombVonMises").def(py::init<>());

    py::class_<  GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<3>>>>,
    typename GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage2DModifiedMohrCoulombModifiedMohrCoulomb").def(py::init<>());

    py::class_<  GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<DruckerPragerPlasticPotential<3>>>>,
    typename GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<DruckerPragerPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage2DModifiedMohrCoulombDruckerPrager").def(py::init<>());

    py::class_<  GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<TrescaPlasticPotential<3>>>>,
    typename GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<TrescaPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage2DModifiedMohrCoulombTresca").def(py::init<>());

    py::class_<  GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<VonMisesPlasticPotential<3>>>>,
    typename GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<VonMisesPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage2DTrescaVonMises").def(py::init<>());

    py::class_<  GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<ModifiedMohrCoulombPlasticPotential<3>>>>,
    typename GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<ModifiedMohrCoulombPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage2DTrescaModifiedMohrCoulomb").def(py::init<>());

    py::class_<  GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<DruckerPragerPlasticPotential<3>>>>,
    typename GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<DruckerPragerPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage2DTrescaDruckerPrager").def(py::init<>());

    py::class_<  GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<TrescaPlasticPotential<3>>>>,
    typename GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<TrescaPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage2DTrescaTresca").def(py::init<>());

    py::class_<  GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<3>>>>,
    typename GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage2DDruckerPragerVonMises").def(py::init<>());

    py::class_<  GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<ModifiedMohrCoulombPlasticPotential<3>>>>,
    typename GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<ModifiedMohrCoulombPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage2DDruckerPragerModifiedMohrCoulomb").def(py::init<>());

    py::class_<  GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<3>>>>,
    typename GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage2DDruckerPragerDruckerPrager").def(py::init<>());

    py::class_<  GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<TrescaPlasticPotential<3>>>>,
    typename GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<TrescaPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage2DDruckerPragerTresca").def(py::init<>());

    py::class_<  GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<VonMisesPlasticPotential<3>>>>,
    typename GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<VonMisesPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage2DRankineVonMises").def(py::init<>());

    py::class_<  GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<ModifiedMohrCoulombPlasticPotential<3>>>>,
    typename GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<ModifiedMohrCoulombPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage2DRankineModifiedMohrCoulomb").def(py::init<>());

    py::class_<  GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<DruckerPragerPlasticPotential<3>>>>,
    typename GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<DruckerPragerPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage2DRankineDruckerPrager").def(py::init<>());

    py::class_<  GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<TrescaPlasticPotential<3>>>>,
    typename GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<TrescaPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage2DRankineTresca").def(py::init<>());

    py::class_<  GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<VonMisesPlasticPotential<3>>>>,
    typename GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<VonMisesPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage2DSimoJuVonMises").def(py::init<>());

    py::class_<  GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<ModifiedMohrCoulombPlasticPotential<3>>>>,
    typename GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<ModifiedMohrCoulombPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage2DSimoJuModifiedMohrCoulomb").def(py::init<>());

    py::class_<  GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<DruckerPragerPlasticPotential<3>>>>,
    typename GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<DruckerPragerPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage2DSimoJuDruckerPrager").def(py::init<>());

    py::class_<  GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<TrescaPlasticPotential<3>>>>,
    typename GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<TrescaPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage2DSimoJuTresca").def(py::init<>());

    py::class_<  GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<MohrCoulombPlasticPotential<3>>>>,
    typename GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<MohrCoulombPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage2DVonMisesMohrCoulomb").def(py::init<>());

    py::class_<  GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>>,
    typename GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage2DMohrCoulombVonMises").def(py::init<>());

    py::class_<  GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<MohrCoulombYieldSurface<MohrCoulombPlasticPotential<3>>>>,
    typename GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<MohrCoulombYieldSurface<MohrCoulombPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage2DMohrCoulombMohrCoulomb").def(py::init<>());

    py::class_<  GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<MohrCoulombYieldSurface<DruckerPragerPlasticPotential<3>>>>,
    typename GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<MohrCoulombYieldSurface<DruckerPragerPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage2DMohrCoulombDruckerPrager").def(py::init<>());

    py::class_<  GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<MohrCoulombYieldSurface<TrescaPlasticPotential<3>>>>,
    typename GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<MohrCoulombYieldSurface<TrescaPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage2DMohrCoulombTresca").def(py::init<>());

    py::class_<  GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<MohrCoulombPlasticPotential<3>>>>,
    typename GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<MohrCoulombPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage2DTrescaMohrCoulomb").def(py::init<>());

    py::class_<  GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<MohrCoulombPlasticPotential<3>>>>,
    typename GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<MohrCoulombPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage2DDruckerPragerMohrCoulomb").def(py::init<>());

    py::class_<  GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<MohrCoulombPlasticPotential<3>>>>,
    typename GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<MohrCoulombPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage2DRankineMohrCoulomb").def(py::init<>());

    py::class_<  GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<MohrCoulombPlasticPotential<3>>>>,
    typename GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<MohrCoulombPlasticPotential<3>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage2DSimoJuMohrCoulomb").def(py::init<>());

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

    py::class_<  GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<MohrCoulombPlasticPotential<6>>>>,
    typename GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<MohrCoulombPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage3DVonMisesMohrCoulomb").def(py::init<>());

    py::class_<  GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>,
    typename GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage3DMohrCoulombVonMises").def(py::init<>());

    py::class_<  GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<MohrCoulombYieldSurface<MohrCoulombPlasticPotential<6>>>>,
    typename GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<MohrCoulombYieldSurface<MohrCoulombPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage3DMohrCoulombMohrCoulomb").def(py::init<>());

    py::class_<  GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<MohrCoulombYieldSurface<DruckerPragerPlasticPotential<6>>>>,
    typename GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<MohrCoulombYieldSurface<DruckerPragerPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage3DMohrCoulombDruckerPrager").def(py::init<>());

    py::class_<  GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<MohrCoulombYieldSurface<TrescaPlasticPotential<6>>>>,
    typename GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<MohrCoulombYieldSurface<TrescaPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage3DMohrCoulombTresca").def(py::init<>());

    py::class_<  GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<MohrCoulombPlasticPotential<6>>>>,
    typename GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<MohrCoulombPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage3DTrescaMohrCoulomb").def(py::init<>());

    py::class_<  GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<MohrCoulombPlasticPotential<6>>>>,
    typename GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<MohrCoulombPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage3DDruckerPragerMohrCoulomb").def(py::init<>());

    py::class_<  GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<MohrCoulombPlasticPotential<6>>>>,
    typename GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<MohrCoulombPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage3DRankineMohrCoulomb").def(py::init<>());

    py::class_<  GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<MohrCoulombPlasticPotential<6>>>>,
    typename GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<MohrCoulombPlasticPotential<6>>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage3DSimoJuMohrCoulomb").def(py::init<>());

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

    py::class_< RuleOfMixturesLaw, typename RuleOfMixturesLaw::Pointer,  ConstitutiveLaw  >
    (m,"RuleOfMixturesLaw").def(py::init<>())
    ;

    py::class_< SerialParallelRuleOfMixturesLaw, typename SerialParallelRuleOfMixturesLaw::Pointer,  ConstitutiveLaw  >
    (m,"SerialParallelRuleOfMixturesLaw").def(py::init<>())
    ;
}

}  // namespace Python.
}  // namespace Kratos.
