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
#include "custom_python/add_custom_constitutive_laws_to_python.h"

#include "custom_constitutive/truss_constitutive_law.h"
#include "custom_constitutive/truss_plasticity_constitutive_law.h"
#include "custom_constitutive/beam_constitutive_law.h"
#include "custom_constitutive/elastic_isotropic_3d.h"
#include "custom_constitutive/axisym_elastic_isotropic.h"
#include "custom_constitutive/linear_plane_stress.h"
#include "custom_constitutive/linear_plane_strain.h"
#include "custom_constitutive/elastic_isotropic_plane_stress_uncoupled_shear.h"
#include "custom_constitutive/hyper_elastic_isotropic_kirchhoff_3d.h"
#include "custom_constitutive/hyper_elastic_isotropic_kirchhoff_plane_stress_2d.h"
#include "custom_constitutive/hyper_elastic_isotropic_kirchhoff_plane_strain_2d.h"
#include "custom_constitutive/hyper_elastic_isotropic_neo_hookean_3d.h"
#include "custom_constitutive/hyper_elastic_isotropic_neo_hookean_plane_strain_2d.h"
#include "custom_constitutive/linear_elastic_orthotropic_2D_law.h"
#include "custom_constitutive/linear_j2_plasticity_3d.h"
#include "custom_constitutive/linear_j2_plasticity_plane_strain_2d.h"
#include "custom_constitutive/linear_isotropic_damage_3D_law.h"

// CL
#include "custom_constitutive/small_strain_isotropic_plasticity_factory_3d.h"
#include "custom_constitutive/small_strain_isotropic_damage_factory_3d.h"
#include "custom_constitutive/viscous_generalized_maxwell_3d.h"
#include "custom_constitutive/viscous_generalized_kelvin_3d.h"
#include "custom_constitutive/generic_small_strain_viscoplasticity_3d.h"
#include "custom_constitutive/generic_small_strain_isotropic_plasticity_3d.h"
#include "custom_constitutive/generic_small_strain_isotropic_damage_3d.h"

namespace Kratos
{
namespace Python
{

using namespace pybind11;


void  AddCustomConstitutiveLawsToPython(pybind11::module& m)
{

    class_< TrussConstitutiveLaw, typename TrussConstitutiveLaw::Pointer, ConstitutiveLaw >
    (m, "TrussConstitutiveLaw").def(init<>() )
    ;

    class_< TrussPlasticityConstitutiveLaw, typename TrussPlasticityConstitutiveLaw::Pointer, ConstitutiveLaw >
    (m, "TrussPlasticityConstitutiveLaw").def(init<>() )
    ;

    class_< BeamConstitutiveLaw, typename BeamConstitutiveLaw::Pointer, ConstitutiveLaw >
    (m, "BeamConstitutiveLaw").def(init<>() )
    ;

    class_< LinearPlaneStress, typename LinearPlaneStress::Pointer, ConstitutiveLaw >
    (m, "LinearElasticPlaneStress2DLaw").def(init<>() )
    ;

    class_< LinearPlaneStrain, typename LinearPlaneStrain::Pointer, ConstitutiveLaw >
    (m, "LinearElasticPlaneStrain2DLaw").def(init<>() )
    ;

    class_< ElasticIsotropic3D, typename ElasticIsotropic3D::Pointer, ConstitutiveLaw >
    (m, "LinearElastic3DLaw").def(init<>() )
    ;

    class_< AxisymElasticIsotropic, typename AxisymElasticIsotropic::Pointer, ConstitutiveLaw >
    (m, "LinearElasticAxisym2DLaw").def(init<>() )
    ;

    class_< ElasticIsotropicPlaneStressUncoupledShear, typename ElasticIsotropicPlaneStressUncoupledShear::Pointer, ConstitutiveLaw >
    (m, "ElasticPlaneStressUncoupledShear2DLaw").def(init<>() )
    ;

    class_< HyperElasticIsotropicKirchhoff3D, typename HyperElasticIsotropicKirchhoff3D::Pointer, ConstitutiveLaw >
    (m, "KirchhoffSaintVenant3DLaw").def(init<>() )
    ;

    class_< HyperElasticIsotropicKirchhoffPlaneStress2D, typename HyperElasticIsotropicKirchhoffPlaneStress2D::Pointer, ConstitutiveLaw >
    (m, "KirchhoffSaintVenantPlaneStress2DLaw").def(init<>() )
    ;

    class_< HyperElasticIsotropicKirchhoffPlaneStrain2D, typename HyperElasticIsotropicKirchhoffPlaneStrain2D::Pointer, ConstitutiveLaw >
    (m, "KirchhoffSaintVenantPlaneStrain2DLaw").def(init<>() )
    ;

    class_< HyperElasticIsotropicNeoHookean3D, typename HyperElasticIsotropicNeoHookean3D::Pointer, ConstitutiveLaw >
    (m, "HyperElastic3DLaw").def(init<>() )
    ;

    class_< HyperElasticIsotropicNeoHookeanPlaneStrain2D, typename HyperElasticIsotropicNeoHookeanPlaneStrain2D::Pointer, ConstitutiveLaw >
    (m, "HyperElasticPlaneStrain2DLaw").def(init<>() )
    ;

    class_< LinearElasticOrthotropic2DLaw, typename LinearElasticOrthotropic2DLaw::Pointer, ConstitutiveLaw >
    (m,"LinearElasticOrthotropic2DLaw").def( init<>())
    ;

    class_< LinearJ2PlasticityPlaneStrain2D, typename LinearJ2PlasticityPlaneStrain2D::Pointer,  ConstitutiveLaw >
    (m,"LinearJ2PlasticityPlaneStrain2DLaw").def( init<>())
    ;

    class_< LinearJ2Plasticity3D, typename LinearJ2Plasticity3D::Pointer,  ConstitutiveLaw >
    (m,"LinearJ2Plasticity3DLaw").def( init<>())
    ;

    class_< LinearIsotropicDamage3D, typename LinearIsotropicDamage3D::Pointer,  ConstitutiveLaw  >
    (m,"LinearIsotropicDamage3DLaw").def( init<>())
    ;
	
    
    class_< SmallStrainIsotropicDamageFactory3D, typename SmallStrainIsotropicDamageFactory3D::Pointer,  ConstitutiveLaw  >
    (m,"SmallStrainIsotropicDamageFactory3D").def( init<>())
    ;

    class_< SmallStrainIsotropicPlasticityFactory3D, typename SmallStrainIsotropicPlasticityFactory3D::Pointer,  ConstitutiveLaw  >
    (m,"SmallStrainIsotropicPlasticityFactory3D").def( init<>())
    ;

    class_< ViscousGeneralizedKelvin3D, typename ViscousGeneralizedKelvin3D::Pointer,  ConstitutiveLaw  >
    (m,"ViscousGeneralizedKelvin3D").def( init<>())
    ;

    class_< ViscousGeneralizedMaxwell3D, typename ViscousGeneralizedMaxwell3D::Pointer,  ConstitutiveLaw  >
    (m,"ViscousGeneralizedMaxwell3D").def( init<>())
    ;

	// Custom Constitutive Laws Registration
	// Plasticity
    class_< GenericSmallStrainIsotropicPlasticity3D <GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<VonMisesPlasticPotential>>>,
    typename GenericSmallStrainIsotropicPlasticity3D <GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<VonMisesPlasticPotential>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicPlasticity3DVonMisesVonMises").def( init<>());

    class_< GenericSmallStrainIsotropicPlasticity3D <GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<ModifiedMohrCoulombPlasticPotential>>>,
    typename GenericSmallStrainIsotropicPlasticity3D <GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<ModifiedMohrCoulombPlasticPotential>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicPlasticity3DVonMisesModifiedMohrCoulomb").def( init<>());

    class_< GenericSmallStrainIsotropicPlasticity3D <GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<DruckerPragerPlasticPotential>>>,
    typename GenericSmallStrainIsotropicPlasticity3D <GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<DruckerPragerPlasticPotential>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicPlasticity3DVonMisesDruckerPrager").def( init<>());

    class_< GenericSmallStrainIsotropicPlasticity3D <GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<TrescaPlasticPotential>>>,
    typename GenericSmallStrainIsotropicPlasticity3D <GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<TrescaPlasticPotential>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicPlasticity3DVonMisesTresca").def( init<>());

    class_< GenericSmallStrainIsotropicPlasticity3D <GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential>>>,
    typename GenericSmallStrainIsotropicPlasticity3D <GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicPlasticity3DModifiedMohrCoulombVonMises").def( init<>());

    class_< GenericSmallStrainIsotropicPlasticity3D <GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential>>>,
    typename GenericSmallStrainIsotropicPlasticity3D <GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicPlasticity3DModifiedMohrCoulombModifiedMohrCoulomb").def( init<>());

    class_< GenericSmallStrainIsotropicPlasticity3D <GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<DruckerPragerPlasticPotential>>>,
    typename GenericSmallStrainIsotropicPlasticity3D <GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<DruckerPragerPlasticPotential>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicPlasticity3DModifiedMohrCoulombDruckerPrager").def( init<>());

    class_< GenericSmallStrainIsotropicPlasticity3D <GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<TrescaPlasticPotential>>>,
    typename GenericSmallStrainIsotropicPlasticity3D <GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<TrescaPlasticPotential>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicPlasticity3DModifiedMohrCoulombTresca").def( init<>());

    class_< GenericSmallStrainIsotropicPlasticity3D <GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<VonMisesPlasticPotential>>>,
    typename GenericSmallStrainIsotropicPlasticity3D <GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<VonMisesPlasticPotential>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicPlasticity3DTrescaVonMises").def( init<>());

    class_< GenericSmallStrainIsotropicPlasticity3D <GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<ModifiedMohrCoulombPlasticPotential>>>,
    typename GenericSmallStrainIsotropicPlasticity3D <GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<ModifiedMohrCoulombPlasticPotential>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicPlasticity3DTrescaModifiedMohrCoulomb").def( init<>());

    class_< GenericSmallStrainIsotropicPlasticity3D <GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<DruckerPragerPlasticPotential>>>,
    typename GenericSmallStrainIsotropicPlasticity3D <GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<DruckerPragerPlasticPotential>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicPlasticity3DTrescaDruckerPrager").def( init<>());

    class_< GenericSmallStrainIsotropicPlasticity3D <GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<TrescaPlasticPotential>>>,
    typename GenericSmallStrainIsotropicPlasticity3D <GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<TrescaPlasticPotential>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicPlasticity3DTrescaTresca").def( init<>());

    class_< GenericSmallStrainIsotropicPlasticity3D <GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<VonMisesPlasticPotential>>>,
    typename GenericSmallStrainIsotropicPlasticity3D <GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<VonMisesPlasticPotential>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicPlasticity3DDruckerPragerVonMises").def( init<>());

    class_< GenericSmallStrainIsotropicPlasticity3D <GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<ModifiedMohrCoulombPlasticPotential>>>,
    typename GenericSmallStrainIsotropicPlasticity3D <GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<ModifiedMohrCoulombPlasticPotential>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicPlasticity3DDruckerPragerModifiedMohrCoulomb").def( init<>());

    class_< GenericSmallStrainIsotropicPlasticity3D <GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<DruckerPragerPlasticPotential>>>,
    typename GenericSmallStrainIsotropicPlasticity3D <GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<DruckerPragerPlasticPotential>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicPlasticity3DDruckerPragerDruckerPrager").def( init<>());

    class_< GenericSmallStrainIsotropicPlasticity3D <GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<TrescaPlasticPotential>>>,
    typename GenericSmallStrainIsotropicPlasticity3D <GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<TrescaPlasticPotential>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicPlasticity3DDruckerPragerTresca").def( init<>());

	// Damage
    class_<  GenericSmallStrainIsotropicDamage3D <GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<VonMisesPlasticPotential>>>,
    typename GenericSmallStrainIsotropicDamage3D <GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<VonMisesPlasticPotential>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage3DVonMisesVonMises").def( init<>());

    class_<  GenericSmallStrainIsotropicDamage3D <GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<ModifiedMohrCoulombPlasticPotential>>>,
    typename GenericSmallStrainIsotropicDamage3D <GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<ModifiedMohrCoulombPlasticPotential>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage3DVonMisesModifiedMohrCoulomb").def( init<>());

    class_<  GenericSmallStrainIsotropicDamage3D <GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<DruckerPragerPlasticPotential>>>,
    typename GenericSmallStrainIsotropicDamage3D <GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<DruckerPragerPlasticPotential>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage3DVonMisesDruckerPrager").def( init<>());

    class_<  GenericSmallStrainIsotropicDamage3D <GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<TrescaPlasticPotential>>>,
    typename GenericSmallStrainIsotropicDamage3D <GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<TrescaPlasticPotential>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage3DVonMisesTresca").def( init<>());

    class_<  GenericSmallStrainIsotropicDamage3D <GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential>>>,
    typename GenericSmallStrainIsotropicDamage3D <GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage3DModifiedMohrCoulombVonMises").def( init<>());

    class_<  GenericSmallStrainIsotropicDamage3D <GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential>>>,
    typename GenericSmallStrainIsotropicDamage3D <GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage3DModifiedMohrCoulombModifiedMohrCoulomb").def( init<>());

    class_<  GenericSmallStrainIsotropicDamage3D <GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<DruckerPragerPlasticPotential>>>,
    typename GenericSmallStrainIsotropicDamage3D <GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<DruckerPragerPlasticPotential>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage3DModifiedMohrCoulombDruckerPrager").def( init<>());

    class_<  GenericSmallStrainIsotropicDamage3D <GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<TrescaPlasticPotential>>>,
    typename GenericSmallStrainIsotropicDamage3D <GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<TrescaPlasticPotential>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage3DModifiedMohrCoulombTresca").def( init<>());

    class_<  GenericSmallStrainIsotropicDamage3D <GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<VonMisesPlasticPotential>>>,
    typename GenericSmallStrainIsotropicDamage3D <GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<VonMisesPlasticPotential>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage3DTrescaVonMises").def( init<>());

    class_<  GenericSmallStrainIsotropicDamage3D <GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<ModifiedMohrCoulombPlasticPotential>>>,
    typename GenericSmallStrainIsotropicDamage3D <GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<ModifiedMohrCoulombPlasticPotential>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage3DTrescaModifiedMohrCoulomb").def( init<>());

    class_<  GenericSmallStrainIsotropicDamage3D <GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<DruckerPragerPlasticPotential>>>,
    typename GenericSmallStrainIsotropicDamage3D <GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<DruckerPragerPlasticPotential>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage3DTrescaDruckerPrager").def( init<>());

    class_<  GenericSmallStrainIsotropicDamage3D <GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<TrescaPlasticPotential>>>,
    typename GenericSmallStrainIsotropicDamage3D <GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<TrescaPlasticPotential>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage3DTrescaTresca").def( init<>());

    class_<  GenericSmallStrainIsotropicDamage3D <GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential>>>,
    typename GenericSmallStrainIsotropicDamage3D <GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage3DDruckerPragerVonMises").def( init<>());

    class_<  GenericSmallStrainIsotropicDamage3D <GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<ModifiedMohrCoulombPlasticPotential>>>,
    typename GenericSmallStrainIsotropicDamage3D <GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<ModifiedMohrCoulombPlasticPotential>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage3DDruckerPragerModifiedMohrCoulomb").def( init<>());

    class_<  GenericSmallStrainIsotropicDamage3D <GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<DruckerPragerPlasticPotential>>>,
    typename GenericSmallStrainIsotropicDamage3D <GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<DruckerPragerPlasticPotential>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage3DDruckerPragerDruckerPrager").def( init<>());

    class_<  GenericSmallStrainIsotropicDamage3D <GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<TrescaPlasticPotential>>>,
    typename GenericSmallStrainIsotropicDamage3D <GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<TrescaPlasticPotential>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage3DDruckerPragerTresca").def( init<>());

    class_<  GenericSmallStrainIsotropicDamage3D <GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<VonMisesPlasticPotential>>>,
    typename GenericSmallStrainIsotropicDamage3D <GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<VonMisesPlasticPotential>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage3DRankineVonMises").def( init<>());

    class_<  GenericSmallStrainIsotropicDamage3D <GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<ModifiedMohrCoulombPlasticPotential>>>,
    typename GenericSmallStrainIsotropicDamage3D <GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<ModifiedMohrCoulombPlasticPotential>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage3DRankineModifiedMohrCoulomb").def( init<>());

    class_<  GenericSmallStrainIsotropicDamage3D <GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<DruckerPragerPlasticPotential>>>,
    typename GenericSmallStrainIsotropicDamage3D <GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<DruckerPragerPlasticPotential>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage3DRankineDruckerPrager").def( init<>());

    class_<  GenericSmallStrainIsotropicDamage3D <GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<TrescaPlasticPotential>>>,
    typename GenericSmallStrainIsotropicDamage3D <GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<TrescaPlasticPotential>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage3DRankineTresca").def( init<>());

    class_<  GenericSmallStrainIsotropicDamage3D <GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<VonMisesPlasticPotential>>>,
    typename GenericSmallStrainIsotropicDamage3D <GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<VonMisesPlasticPotential>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage3DSimoJuVonMises").def( init<>());

    class_<  GenericSmallStrainIsotropicDamage3D <GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<ModifiedMohrCoulombPlasticPotential>>>,
    typename GenericSmallStrainIsotropicDamage3D <GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<ModifiedMohrCoulombPlasticPotential>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage3DSimoJuModifiedMohrCoulomb").def( init<>());

    class_<  GenericSmallStrainIsotropicDamage3D <GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<DruckerPragerPlasticPotential>>>,
    typename GenericSmallStrainIsotropicDamage3D <GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<DruckerPragerPlasticPotential>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage3DSimoJuDruckerPrager").def( init<>());

    class_<  GenericSmallStrainIsotropicDamage3D <GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<TrescaPlasticPotential>>>,
    typename GenericSmallStrainIsotropicDamage3D <GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<TrescaPlasticPotential>>>::Pointer,
    ConstitutiveLaw >
    (m,"SmallStrainIsotropicDamage3DSimoJuTresca").def( init<>());

    class_< GenericSmallStrainViscoplasticity3D, typename GenericSmallStrainViscoplasticity3D::Pointer,  ConstitutiveLaw  >
    (m,"GenericSmallStrainViscoplasticity3D").def( init<>())
    ;
}

}  // namespace Python.
}  // namespace Kratos.
