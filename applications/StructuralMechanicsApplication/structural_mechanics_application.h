// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                     license: structural_mechanics_application/license.txt
//
//  Main authors:    Riccardo Rossi
//    Co-authors:    Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_STRUCTURAL_MECHANICS_APPLICATION_H_INCLUDED )
#define  KRATOS_STRUCTURAL_MECHANICS_APPLICATION_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"

/* ELEMENTS */

/* Adding truss element */
#include "custom_elements/truss_element_3D2N.hpp"
#include "custom_elements/truss_element_linear_3D2N.hpp"
#include "custom_elements/cable_element_3D2N.hpp"

/* Adding beam element */
#include "custom_elements/cr_beam_element_3D2N.hpp"
#include "custom_elements/cr_beam_element_linear_3D2N.hpp"
#include "custom_elements/cr_beam_element_2D2N.hpp"
#include "custom_elements/cr_beam_element_linear_2D2N.hpp"

/* Adding the adjoint elements */
#include "custom_response_functions/adjoint_elements/adjoint_finite_difference_shell_element.h"
#include "custom_response_functions/adjoint_elements/adjoint_finite_difference_cr_beam_element_3D2N.h"
#include "custom_response_functions/adjoint_elements/adjoint_finite_difference_truss_element_3D2N.h"
#include "custom_response_functions/adjoint_elements/adjoint_finite_difference_truss_element_linear_3D2N.h"
#include "custom_response_functions/adjoint_elements/adjoint_solid_element.h"
#include "custom_response_functions/adjoint_elements/adjoint_finite_difference_small_displacement_element.h"

/* Adding shells and membranes elements */
#include "custom_elements/isotropic_shell_element.hpp"
#include "custom_elements/prestress_membrane_element.hpp"
#include "custom_elements/shell_thick_element_3D4N.hpp"
#include "custom_elements/shell_thin_element_3D4N.hpp"
#include "custom_elements/shell_thin_element_3D3N.hpp"
#include "custom_elements/shell_thick_element_3D3N.hpp"
#include "custom_elements/nodal_concentrated_element.hpp"

/* Adding the spring damper element */
#include "custom_elements/spring_damper_element_3D2N.hpp"

/* Adding the SPRISM element */
#include "custom_elements/solid_shell_element_sprism_3D6N.h"

/* Adding solid elements */
#include "custom_elements/small_displacement.h"
#include "custom_elements/axisym_small_displacement.h"
#include "custom_elements/total_lagrangian.h"
#include "custom_elements/axisym_total_lagrangian.h"
#include "custom_elements/updated_lagrangian.h"
#include "custom_elements/axisym_updated_lagrangian.h"
#include "custom_elements/small_displacement_bbar.h"

/* CONDITIONS */
#include "custom_conditions/base_load_condition.h"
#include "custom_conditions/point_load_condition.h"
#include "custom_conditions/point_contact_condition.h"
#include "custom_conditions/axisym_point_load_condition.h"
#include "custom_conditions/line_load_condition.h"
#include "custom_conditions/line_load_condition_2d.h"
#include "custom_conditions/axisym_line_load_condition_2d.h"
#include "custom_conditions/surface_load_condition_3d.h"
#include "custom_conditions/point_moment_condition_3d.h"

/* Adding the adjoint conditions */
#include "custom_response_functions/adjoint_conditions/adjoint_semi_analytic_point_load_condition.h"

/* CONSTITUTIVE LAWS */
#include "custom_constitutive/truss_plasticity_constitutive_law.h"
#include "custom_constitutive/truss_constitutive_law.h"
#include "custom_constitutive/beam_constitutive_law.h"
#include "custom_constitutive/elastic_isotropic_3d.h"
#include "custom_constitutive/axisym_elastic_isotropic.h"
#include "custom_constitutive/linear_plane_strain.h"
#include "custom_constitutive/linear_plane_stress.h"
#include "custom_constitutive/elastic_isotropic_plane_stress_uncoupled_shear.h"
#include "custom_constitutive/hyper_elastic_isotropic_kirchhoff_3d.h"
#include "custom_constitutive/hyper_elastic_isotropic_kirchhoff_plane_stress_2d.h"
#include "custom_constitutive/hyper_elastic_isotropic_kirchhoff_plane_strain_2d.h"
#include "custom_constitutive/hyper_elastic_isotropic_neo_hookean_3d.h"
#include "custom_constitutive/hyper_elastic_isotropic_neo_hookean_plane_strain_2d.h"
#include "custom_constitutive/linear_elastic_orthotropic_2D_law.h"
#include "custom_constitutive/small_strain_j2_plasticity_plane_strain_2d.h"
#include "custom_constitutive/small_strain_j2_plasticity_3d.h"
#include "custom_constitutive/small_strain_isotropic_damage_3d.h"
#include "custom_constitutive/small_strain_isotropic_damage_plane_strain_2d.h"
#include "custom_constitutive/small_strain_isotropic_damage_traction_only_3d.h"
#include "custom_constitutive/plane_stress_d_plus_d_minus_damage_masonry_2d.h"
#include "custom_constitutive/d_plus_d_minus_damage_masonry_3d.h"

// Advanced Constitutive laws
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

// Integrators
#include "custom_constitutive/constitutive_laws_integrators/generic_constitutive_law_integrator_damage.h"
#include "custom_constitutive/constitutive_laws_integrators/generic_constitutive_law_integrator_plasticity.h"
#include "custom_constitutive/constitutive_laws_integrators/generic_constitutive_law_integrator_kinematic_plasticity.h"
#include "custom_constitutive/constitutive_laws_integrators/generic_finite_strain_constitutive_law_integrator_plasticity.h"
#include "custom_constitutive/constitutive_laws_integrators/generic_finite_strain_constitutive_law_integrator_kinematic_plasticity.h"
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

namespace Kratos
{

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

/**
 * @class KratosStructuralMechanicsApplication
 * @ingroup StructuralMechanicsApplication
 * @brief This application features Elements, Conditions, Constitutive laws and Utilities for structural analysis problems
 * @author Riccardo Rossi
 */
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) KratosStructuralMechanicsApplication : public KratosApplication
{
public:
    ///@name Type Definitions
    ///@{


    /// Pointer definition of KratosStructuralMechanicsApplication
    KRATOS_CLASS_POINTER_DEFINITION(KratosStructuralMechanicsApplication);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    KratosStructuralMechanicsApplication();

    /// Destructor.
    ~KratosStructuralMechanicsApplication() override {}


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
        return "KratosStructuralMechanicsApplication";
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
        KRATOS_INFO("KratosStructuralMechanicsApplication") << "Has the following number of components: " << KratosComponents<VariableData>::GetComponents().size() << std::endl;
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

//     static const ApplicationCondition  msApplicationCondition;

    ///@}
    ///@name Member Variables
    ///@{


    /* ELEMENTS */

    // Adding the truss element
    const TrussElement3D2N mTrussElement3D2N;
    const TrussElementLinear3D2N mTrussLinearElement3D2N;
    const CableElement3D2N mCableElement3D2N;

    // Adding the beam element
    const CrBeamElement3D2N mCrBeamElement3D2N;
    const CrBeamElementLinear3D2N mCrLinearBeamElement3D2N;
    const CrBeamElement2D2N mCrBeamElement2D2N;
    const CrBeamElementLinear2D2N mCrLinearBeamElement2D2N;


    // Adding the shells elements
    const IsotropicShellElement mIsotropicShellElement3D3N;
    const ShellThickElement3D4N mShellThickElement3D4N;
    const ShellThickElement3D4N mShellThickCorotationalElement3D4N;
    const ShellThinElement3D4N   mShellThinCorotationalElement3D4N;
    const ShellThinElement3D3N mShellThinElement3D3N;
    const ShellThinElement3D3N mShellThinCorotationalElement3D3N;
    const ShellThickElement3D3N  mShellThickCorotationalElement3D3N;

    // Adding the membrane element
    const PrestressMembraneElement mPreStressMembraneElement3D3N;
    const PrestressMembraneElement mPreStressMembraneElement3D4N;

    // Adding the SPRISM element
    const SolidShellElementSprism3D6N mSolidShellElementSprism3D6N;

    // Adding the nodal concentrated element
    const NodalConcentratedElement mNodalConcentratedElement2D1N;
    const NodalConcentratedElement mNodalConcentratedDampedElement2D1N;
    const NodalConcentratedElement mNodalConcentratedElement3D1N;
    const NodalConcentratedElement mNodalConcentratedDampedElement3D1N;

    // Linear kinematic elements
    const SmallDisplacement mSmallDisplacement2D3N;
    const SmallDisplacement mSmallDisplacement2D4N;
    const SmallDisplacement mSmallDisplacement2D6N;
    const SmallDisplacement mSmallDisplacement2D8N;
    const SmallDisplacement mSmallDisplacement2D9N;
    const SmallDisplacement mSmallDisplacement3D4N;
    const SmallDisplacement mSmallDisplacement3D6N;
    const SmallDisplacement mSmallDisplacement3D8N;
    const SmallDisplacement mSmallDisplacement3D10N;
    const SmallDisplacement mSmallDisplacement3D15N;
    const SmallDisplacement mSmallDisplacement3D20N;
    const SmallDisplacement mSmallDisplacement3D27N;

    const SmallDisplacementBbar mSmallDisplacementBbar2D4N;
    const SmallDisplacementBbar mSmallDisplacementBbar3D8N;

    const AxisymSmallDisplacement mAxisymSmallDisplacement2D3N;
    const AxisymSmallDisplacement mAxisymSmallDisplacement2D4N;
    const AxisymSmallDisplacement mAxisymSmallDisplacement2D6N;
    const AxisymSmallDisplacement mAxisymSmallDisplacement2D8N;
    const AxisymSmallDisplacement mAxisymSmallDisplacement2D9N;

    // Total lagrangian
    const TotalLagrangian mTotalLagrangian2D3N;
    const TotalLagrangian mTotalLagrangian2D4N;
    const TotalLagrangian mTotalLagrangian2D6N;
    const TotalLagrangian mTotalLagrangian2D8N;
    const TotalLagrangian mTotalLagrangian2D9N;
    const TotalLagrangian mTotalLagrangian3D4N;
    const TotalLagrangian mTotalLagrangian3D6N;
    const TotalLagrangian mTotalLagrangian3D8N;
    const TotalLagrangian mTotalLagrangian3D10N;
    const TotalLagrangian mTotalLagrangian3D15N;
    const TotalLagrangian mTotalLagrangian3D20N;
    const TotalLagrangian mTotalLagrangian3D27N;

    const AxisymTotalLagrangian mAxisymTotalLagrangian2D3N;
    const AxisymTotalLagrangian mAxisymTotalLagrangian2D4N;
    const AxisymTotalLagrangian mAxisymTotalLagrangian2D6N;
    const AxisymTotalLagrangian mAxisymTotalLagrangian2D8N;
    const AxisymTotalLagrangian mAxisymTotalLagrangian2D9N;

    // Updated lagrangian
    const UpdatedLagrangian mUpdatedLagrangian2D3N;
    const UpdatedLagrangian mUpdatedLagrangian2D4N;
    const UpdatedLagrangian mUpdatedLagrangian2D6N;
    const UpdatedLagrangian mUpdatedLagrangian2D8N;
    const UpdatedLagrangian mUpdatedLagrangian2D9N;
    const UpdatedLagrangian mUpdatedLagrangian3D4N;
    const UpdatedLagrangian mUpdatedLagrangian3D6N;
    const UpdatedLagrangian mUpdatedLagrangian3D8N;
    const UpdatedLagrangian mUpdatedLagrangian3D10N;
    const UpdatedLagrangian mUpdatedLagrangian3D15N;
    const UpdatedLagrangian mUpdatedLagrangian3D20N;
    const UpdatedLagrangian mUpdatedLagrangian3D27N;

    const AxisymUpdatedLagrangian mAxisymUpdatedLagrangian2D3N;
    const AxisymUpdatedLagrangian mAxisymUpdatedLagrangian2D4N;
    const AxisymUpdatedLagrangian mAxisymUpdatedLagrangian2D6N;
    const AxisymUpdatedLagrangian mAxisymUpdatedLagrangian2D8N;
    const AxisymUpdatedLagrangian mAxisymUpdatedLagrangian2D9N;

    // Adding the spring damper element
    const SpringDamperElement3D2N mSpringDamperElement3D2N;

    // Adding adjoint elements
    const AdjointFiniteDifferencingShellElement<ShellThinElement3D3N> mAdjointFiniteDifferencingShellThinElement3D3N;
    const AdjointFiniteDifferenceCrBeamElement<CrBeamElementLinear3D2N> mAdjointFiniteDifferenceCrBeamElementLinear3D2N;
    const AdjointFiniteDifferenceTrussElement<TrussElement3D2N> mAdjointFiniteDifferenceTrussElement3D2N;
    const AdjointFiniteDifferenceTrussElementLinear<TrussElementLinear3D2N> mAdjointFiniteDifferenceTrussLinearElement3D2N;
    const AdjointSolidElement<TotalLagrangian> mTotalLagrangianAdjoint2D3N;
    const AdjointSolidElement<TotalLagrangian> mTotalLagrangianAdjoint2D4N;
    const AdjointSolidElement<TotalLagrangian> mTotalLagrangianAdjoint2D6N;
    const AdjointSolidElement<TotalLagrangian> mTotalLagrangianAdjoint3D4N;
    const AdjointSolidElement<TotalLagrangian> mTotalLagrangianAdjoint3D8N;
    const AdjointFiniteDifferencingSmallDisplacementElement<SmallDisplacement> mAdjointFiniteDifferencingSmallDisplacementElement3D4N;
    const AdjointFiniteDifferencingSmallDisplacementElement<SmallDisplacement> mAdjointFiniteDifferencingSmallDisplacementElement3D6N;
    const AdjointFiniteDifferencingSmallDisplacementElement<SmallDisplacement> mAdjointFiniteDifferencingSmallDisplacementElement3D8N;

    /* CONDITIONS*/
    // Point load
    const PointLoadCondition mPointLoadCondition2D1N;
    const PointLoadCondition mPointLoadCondition3D1N;
    const PointContactCondition mPointContactCondition2D1N;
    const PointContactCondition mPointContactCondition3D1N;

    const AxisymPointLoadCondition mAxisymPointLoadCondition2D1N;

    // Line load
    const LineLoadCondition<2> mLineLoadCondition2D2N;
    const LineLoadCondition<2> mLineLoadCondition2D3N;
    const LineLoadCondition<3> mLineLoadCondition3D2N;
    const LineLoadCondition<3> mLineLoadCondition3D3N;

    const AxisymLineLoadCondition2D mAxisymLineLoadCondition2D2N;
    const AxisymLineLoadCondition2D mAxisymLineLoadCondition2D3N;

    // Surface load
    const SurfaceLoadCondition3D mSurfaceLoadCondition3D3N;
    const SurfaceLoadCondition3D mSurfaceLoadCondition3D4N;
    const SurfaceLoadCondition3D mSurfaceLoadCondition3D6N;
    const SurfaceLoadCondition3D mSurfaceLoadCondition3D8N;
    const SurfaceLoadCondition3D mSurfaceLoadCondition3D9N;

    // Point moment
    const PointMomentCondition3D mPointMomentCondition3D1N;

    // Adjoint Conditions
    const AdjointSemiAnalyticPointLoadCondition<PointLoadCondition> mAdjointSemiAnalyticPointLoadCondition2D1N;
    const AdjointSemiAnalyticPointLoadCondition<PointLoadCondition> mAdjointSemiAnalyticPointLoadCondition3D1N;

    /* CONSTITUTIVE LAWS */
    // Linear elastics laws
    const TrussConstitutiveLaw mTrussConstitutiveLaw;
    const TrussPlasticityConstitutiveLaw mTrussPlasticityConstitutiveLaw;
    const BeamConstitutiveLaw mBeamConstitutiveLaw;
    const ElasticIsotropic3D mElasticIsotropic3D;
    const AxisymElasticIsotropic mAxisymElasticIsotropic;
    const LinearPlaneStrain  mLinearPlaneStrain;
    const LinearPlaneStress  mLinearPlaneStress;
    const ElasticIsotropicPlaneStressUncoupledShear  mElasticIsotropicPlaneStressUncoupledShear;
    const HyperElasticIsotropicKirchhoff3D  mHyperElasticIsotropicKirchhoff3D;
    const HyperElasticIsotropicKirchhoffPlaneStress2D  mHyperElasticIsotropicKirchhoffPlaneStress2D;
    const HyperElasticIsotropicKirchhoffPlaneStrain2D  mHyperElasticIsotropicKirchhoffPlaneStrain2D;
    const HyperElasticIsotropicNeoHookean3D  mHyperElasticIsotropicNeoHookean3D;
    const HyperElasticIsotropicNeoHookeanPlaneStrain2D  mHyperElasticIsotropicNeoHookeanPlaneStrain2D;
    const LinearElasticOrthotropic2DLaw mLinearElasticOrthotropic2DLaw;

    const SmallStrainJ2Plasticity3D mSmallStrainJ2Plasticity3D;
    const SmallStrainJ2PlasticityPlaneStrain2D mSmallStrainJ2PlasticityPlaneStrain2D;
    const SmallStrainIsotropicDamage3D mSmallStrainIsotropicDamage3D;
    const SmallStrainIsotropicDamagePlaneStrain2D mSmallStrainIsotropicDamagePlaneStrain2D;
    const SmallStrainIsotropicDamageTractionOnly3D mSmallStrainIsotropicDamageTractionOnly3D;

    // Damage and plasticity laws
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


    /* Finite strain */
    // Kirchhoff
    const GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>> mHyperElasticIsotropicKirchhoffPlasticity3DVonMisesVonMises;
    const GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>> mHyperElasticIsotropicKirchhoffPlasticity3DVonMisesModifiedMohrCoulomb;
    const GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<DruckerPragerPlasticPotential<6>>>> mHyperElasticIsotropicKirchhoffPlasticity3DVonMisesDruckerPrager;
    const GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<TrescaPlasticPotential<6>>>> mHyperElasticIsotropicKirchhoffPlasticity3DVonMisesTresca;
    const GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>> mHyperElasticIsotropicKirchhoffPlasticity3DModifiedMohrCoulombVonMises;
    const GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>> mHyperElasticIsotropicKirchhoffPlasticity3DModifiedMohrCoulombModifiedMohrCoulomb;
    const GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<DruckerPragerPlasticPotential<6>>>> mHyperElasticIsotropicKirchhoffPlasticity3DModifiedMohrCoulombDruckerPrager;
    const GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<TrescaPlasticPotential<6>>>> mHyperElasticIsotropicKirchhoffPlasticity3DModifiedMohrCoulombTresca;
    const GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<VonMisesPlasticPotential<6>>>> mHyperElasticIsotropicKirchhoffPlasticity3DTrescaVonMises;
    const GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>> mHyperElasticIsotropicKirchhoffPlasticity3DTrescaModifiedMohrCoulomb;
    const GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<DruckerPragerPlasticPotential<6>>>> mHyperElasticIsotropicKirchhoffPlasticity3DTrescaDruckerPrager;
    const GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<TrescaPlasticPotential<6>>>> mHyperElasticIsotropicKirchhoffPlasticity3DTrescaTresca;
    const GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>> mHyperElasticIsotropicKirchhoffPlasticity3DDruckerPragerVonMises;
    const GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>> mHyperElasticIsotropicKirchhoffPlasticity3DDruckerPragerModifiedMohrCoulomb;
    const GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<6>>>> mHyperElasticIsotropicKirchhoffPlasticity3DDruckerPragerDruckerPrager;
    const GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<TrescaPlasticPotential<6>>>> mHyperElasticIsotropicKirchhoffPlasticity3DDruckerPragerTresca;
    const GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<MohrCoulombPlasticPotential<6>>>> mHyperElasticIsotropicKirchhoffPlasticity3DVonMisesMohrCoulomb;
    const GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>> mHyperElasticIsotropicKirchhoffPlasticity3DMohrCoulombVonMises;
    const GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<MohrCoulombPlasticPotential<6>>>> mHyperElasticIsotropicKirchhoffPlasticity3DMohrCoulombMohrCoulomb;
    const GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<DruckerPragerPlasticPotential<6>>>> mHyperElasticIsotropicKirchhoffPlasticity3DMohrCoulombDruckerPrager;
    const GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<TrescaPlasticPotential<6>>>> mHyperElasticIsotropicKirchhoffPlasticity3DMohrCoulombTresca;
    const GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<MohrCoulombPlasticPotential<6>>>> mHyperElasticIsotropicKirchhoffPlasticity3DTrescaMohrCoulomb;
    const GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<MohrCoulombPlasticPotential<6>>>> mHyperElasticIsotropicKirchhoffPlasticity3DDruckerPragerMohrCoulomb;
    // Neo-Hookean
    const GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>> mHyperElasticIsotropicNeoHookeanPlasticity3DVonMisesVonMises;
    const GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>> mHyperElasticIsotropicNeoHookeanPlasticity3DVonMisesModifiedMohrCoulomb;
    const GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<DruckerPragerPlasticPotential<6>>>> mHyperElasticIsotropicNeoHookeanPlasticity3DVonMisesDruckerPrager;
    const GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<TrescaPlasticPotential<6>>>> mHyperElasticIsotropicNeoHookeanPlasticity3DVonMisesTresca;
    const GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>> mHyperElasticIsotropicNeoHookeanPlasticity3DModifiedMohrCoulombVonMises;
    const GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>> mHyperElasticIsotropicNeoHookeanPlasticity3DModifiedMohrCoulombModifiedMohrCoulomb;
    const GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<DruckerPragerPlasticPotential<6>>>> mHyperElasticIsotropicNeoHookeanPlasticity3DModifiedMohrCoulombDruckerPrager;
    const GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<TrescaPlasticPotential<6>>>> mHyperElasticIsotropicNeoHookeanPlasticity3DModifiedMohrCoulombTresca;
    const GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<VonMisesPlasticPotential<6>>>> mHyperElasticIsotropicNeoHookeanPlasticity3DTrescaVonMises;
    const GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>> mHyperElasticIsotropicNeoHookeanPlasticity3DTrescaModifiedMohrCoulomb;
    const GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<DruckerPragerPlasticPotential<6>>>> mHyperElasticIsotropicNeoHookeanPlasticity3DTrescaDruckerPrager;
    const GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<TrescaPlasticPotential<6>>>> mHyperElasticIsotropicNeoHookeanPlasticity3DTrescaTresca;
    const GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>> mHyperElasticIsotropicNeoHookeanPlasticity3DDruckerPragerVonMises;
    const GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>> mHyperElasticIsotropicNeoHookeanPlasticity3DDruckerPragerModifiedMohrCoulomb;
    const GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<6>>>> mHyperElasticIsotropicNeoHookeanPlasticity3DDruckerPragerDruckerPrager;
    const GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<TrescaPlasticPotential<6>>>> mHyperElasticIsotropicNeoHookeanPlasticity3DDruckerPragerTresca;
    const GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<MohrCoulombPlasticPotential<6>>>> mHyperElasticIsotropicNeoHookeanPlasticity3DVonMisesMohrCoulomb;
    const GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>> mHyperElasticIsotropicNeoHookeanPlasticity3DMohrCoulombVonMises;
    const GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<MohrCoulombPlasticPotential<6>>>> mHyperElasticIsotropicNeoHookeanPlasticity3DMohrCoulombMohrCoulomb;
    const GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<DruckerPragerPlasticPotential<6>>>> mHyperElasticIsotropicNeoHookeanPlasticity3DMohrCoulombDruckerPrager;
    const GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<TrescaPlasticPotential<6>>>> mHyperElasticIsotropicNeoHookeanPlasticity3DMohrCoulombTresca;
    const GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<MohrCoulombPlasticPotential<6>>>> mHyperElasticIsotropicNeoHookeanPlasticity3DTrescaMohrCoulomb;
    const GenericFiniteStrainIsotropicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<MohrCoulombPlasticPotential<6>>>> mHyperElasticIsotropicNeoHookeanPlasticity3DDruckerPragerMohrCoulomb;

    /* Finite strain */
    // Kirchhoff
    const GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>> mHyperElasticKirchhoffKinematicPlasticity3DVonMisesVonMises;
    const GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>> mHyperElasticKirchhoffKinematicPlasticity3DVonMisesModifiedMohrCoulomb;
    const GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<DruckerPragerPlasticPotential<6>>>> mHyperElasticKirchhoffKinematicPlasticity3DVonMisesDruckerPrager;
    const GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<TrescaPlasticPotential<6>>>> mHyperElasticKirchhoffKinematicPlasticity3DVonMisesTresca;
    const GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>> mHyperElasticKirchhoffKinematicPlasticity3DModifiedMohrCoulombVonMises;
    const GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>> mHyperElasticKirchhoffKinematicPlasticity3DModifiedMohrCoulombModifiedMohrCoulomb;
    const GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<ModifiedMohrCoulombYieldSurface<DruckerPragerPlasticPotential<6>>>> mHyperElasticKirchhoffKinematicPlasticity3DModifiedMohrCoulombDruckerPrager;
    const GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<ModifiedMohrCoulombYieldSurface<TrescaPlasticPotential<6>>>> mHyperElasticKirchhoffKinematicPlasticity3DModifiedMohrCoulombTresca;
    const GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<VonMisesPlasticPotential<6>>>> mHyperElasticKirchhoffKinematicPlasticity3DTrescaVonMises;
    const GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>> mHyperElasticKirchhoffKinematicPlasticity3DTrescaModifiedMohrCoulomb;
    const GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<DruckerPragerPlasticPotential<6>>>> mHyperElasticKirchhoffKinematicPlasticity3DTrescaDruckerPrager;
    const GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<TrescaPlasticPotential<6>>>> mHyperElasticKirchhoffKinematicPlasticity3DTrescaTresca;
    const GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>> mHyperElasticKirchhoffKinematicPlasticity3DDruckerPragerVonMises;
    const GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>> mHyperElasticKirchhoffKinematicPlasticity3DDruckerPragerModifiedMohrCoulomb;
    const GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<6>>>> mHyperElasticKirchhoffKinematicPlasticity3DDruckerPragerDruckerPrager;
    const GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<TrescaPlasticPotential<6>>>> mHyperElasticKirchhoffKinematicPlasticity3DDruckerPragerTresca;
    const GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<MohrCoulombPlasticPotential<6>>>> mHyperElasticKirchhoffKinematicPlasticity3DVonMisesMohrCoulomb;
    const GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<MohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>> mHyperElasticKirchhoffKinematicPlasticity3DMohrCoulombVonMises;
    const GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<MohrCoulombYieldSurface<MohrCoulombPlasticPotential<6>>>> mHyperElasticKirchhoffKinematicPlasticity3DMohrCoulombMohrCoulomb;
    const GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<MohrCoulombYieldSurface<DruckerPragerPlasticPotential<6>>>> mHyperElasticKirchhoffKinematicPlasticity3DMohrCoulombDruckerPrager;
    const GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<MohrCoulombYieldSurface<TrescaPlasticPotential<6>>>> mHyperElasticKirchhoffKinematicPlasticity3DMohrCoulombTresca;
    const GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<MohrCoulombPlasticPotential<6>>>> mHyperElasticKirchhoffKinematicPlasticity3DTrescaMohrCoulomb;
    const GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<MohrCoulombPlasticPotential<6>>>> mHyperElasticKirchhoffKinematicPlasticity3DDruckerPragerMohrCoulomb;
    // Neo-Hookean
    const GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>> mHyperElasticNeoHookeanKinematicPlasticity3DVonMisesVonMises;
    const GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>> mHyperElasticNeoHookeanKinematicPlasticity3DVonMisesModifiedMohrCoulomb;
    const GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<DruckerPragerPlasticPotential<6>>>> mHyperElasticNeoHookeanKinematicPlasticity3DVonMisesDruckerPrager;
    const GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<TrescaPlasticPotential<6>>>> mHyperElasticNeoHookeanKinematicPlasticity3DVonMisesTresca;
    const GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>> mHyperElasticNeoHookeanKinematicPlasticity3DModifiedMohrCoulombVonMises;
    const GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>> mHyperElasticNeoHookeanKinematicPlasticity3DModifiedMohrCoulombModifiedMohrCoulomb;
    const GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<ModifiedMohrCoulombYieldSurface<DruckerPragerPlasticPotential<6>>>> mHyperElasticNeoHookeanKinematicPlasticity3DModifiedMohrCoulombDruckerPrager;
    const GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<ModifiedMohrCoulombYieldSurface<TrescaPlasticPotential<6>>>> mHyperElasticNeoHookeanKinematicPlasticity3DModifiedMohrCoulombTresca;
    const GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<VonMisesPlasticPotential<6>>>> mHyperElasticNeoHookeanKinematicPlasticity3DTrescaVonMises;
    const GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>> mHyperElasticNeoHookeanKinematicPlasticity3DTrescaModifiedMohrCoulomb;
    const GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<DruckerPragerPlasticPotential<6>>>> mHyperElasticNeoHookeanKinematicPlasticity3DTrescaDruckerPrager;
    const GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<TrescaPlasticPotential<6>>>> mHyperElasticNeoHookeanKinematicPlasticity3DTrescaTresca;
    const GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>> mHyperElasticNeoHookeanKinematicPlasticity3DDruckerPragerVonMises;
    const GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>> mHyperElasticNeoHookeanKinematicPlasticity3DDruckerPragerModifiedMohrCoulomb;
    const GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<6>>>> mHyperElasticNeoHookeanKinematicPlasticity3DDruckerPragerDruckerPrager;
    const GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<TrescaPlasticPotential<6>>>> mHyperElasticNeoHookeanKinematicPlasticity3DDruckerPragerTresca;
    const GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<MohrCoulombPlasticPotential<6>>>> mHyperElasticNeoHookeanKinematicPlasticity3DVonMisesMohrCoulomb;
    const GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<MohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>> mHyperElasticNeoHookeanKinematicPlasticity3DMohrCoulombVonMises;
    const GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<MohrCoulombYieldSurface<MohrCoulombPlasticPotential<6>>>> mHyperElasticNeoHookeanKinematicPlasticity3DMohrCoulombMohrCoulomb;
    const GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<MohrCoulombYieldSurface<DruckerPragerPlasticPotential<6>>>> mHyperElasticNeoHookeanKinematicPlasticity3DMohrCoulombDruckerPrager;
    const GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<MohrCoulombYieldSurface<TrescaPlasticPotential<6>>>> mHyperElasticNeoHookeanKinematicPlasticity3DMohrCoulombTresca;
    const GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<MohrCoulombPlasticPotential<6>>>> mHyperElasticNeoHookeanKinematicPlasticity3DTrescaMohrCoulomb;
    const GenericFiniteStrainKinematicPlasticity <HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<MohrCoulombPlasticPotential<6>>>> mHyperElasticNeoHookeanKinematicPlasticity3DDruckerPragerMohrCoulomb;

    /// Damage
    /* Small strain */
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

    // d+d- laws
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
    KratosStructuralMechanicsApplication& operator=(KratosStructuralMechanicsApplication const& rOther);

    /// Copy constructor.
    KratosStructuralMechanicsApplication(KratosStructuralMechanicsApplication const& rOther);


    ///@}

}; // Class KratosStructuralMechanicsApplication

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // KRATOS_STRUCTURAL_MECHANICS_APPLICATION_H_INCLUDED  defined


