// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Riccardo Rossi
//    Co-authors:    Vicente Mataix Ferrándiz
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

#include "structural_mechanics_application_variables.h"

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
#include "custom_response_functions/adjoint_elements/adjoint_finite_difference_base_element.h"
#include "custom_response_functions/adjoint_elements/adjoint_finite_difference_shell_element.h"
#include "custom_response_functions/adjoint_elements/adjoint_finite_difference_cr_beam_element_3D2N.h"

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
#include "custom_constitutive/linear_j2_plasticity_plane_strain_2d.h"
#include "custom_constitutive/linear_j2_plasticity_3d.h"
#include "custom_constitutive/linear_isotropic_damage_3D_law.h"

// Custom Constitutive laws
#include "custom_constitutive/small_strain_isotropic_plasticity_factory_3d.h"
#include "custom_constitutive/generic_small_strain_isotropic_plasticity_3d.h"
#include "custom_constitutive/generic_small_strain_isotropic_damage_3d.h"
#include "custom_constitutive/small_strain_isotropic_damage_factory_3d.h"
#include "custom_constitutive/viscous_generalized_kelvin_3d.h"
#include "custom_constitutive/generic_small_strain_viscoplasticity_3d.h"
#include "custom_constitutive/viscous_generalized_maxwell_3d.h"


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

/// Short class definition.
/**
 * This application features Elements, Conditions, Constitutive laws and Utilities
 * for structural analysis problems
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
    const AdjointFiniteDifferencingBaseElement mAdjointFiniteDifferencingBaseElement;
    const AdjointFiniteDifferencingShellElement mAdjointFiniteDifferencingShellElement;
    const AdjointFiniteDifferenceCrBeamElement mAdjointFiniteDifferenceCrBeamElement;

    /* CONDITIONS*/
    // Point load
    const PointLoadCondition mPointLoadCondition2D1N;
    const PointLoadCondition mPointLoadCondition3D1N;
    const PointContactCondition mPointContactCondition2D1N;
    const PointContactCondition mPointContactCondition3D1N;

    const AxisymPointLoadCondition mAxisymPointLoadCondition2D1N;

    // Line load
    const LineLoadCondition2D mLineLoadCondition2D2N;
    const LineLoadCondition2D mLineLoadCondition2D3N;

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
    const AdjointSemiAnalyticPointLoadCondition mAdjointSemiAnalyticPointLoadCondition2D1N;
    const AdjointSemiAnalyticPointLoadCondition mAdjointSemiAnalyticPointLoadCondition3D1N;

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

    const LinearJ2Plasticity3D mLinearJ2Plasticity3D;
    const LinearJ2PlasticityPlaneStrain2D mLinearJ2PlasticityPlaneStrain2D;
    const LinearIsotropicDamage3D mLinearIsotropicDamage3D;

    // Damage and plasticity laws
    const SmallStrainIsotropicPlasticityFactory3D mSmallStrainIsotropicPlasticityFactory3D;
    const SmallStrainIsotropicDamageFactory3D mSmallStrainIsotropicDamageFactory3D;
    const ViscousGeneralizedKelvin3D mViscousGeneralizedKelvin3D;
    const ViscousGeneralizedMaxwell3D mViscousGeneralizedMaxwell3D;
    const GenericSmallStrainViscoplasticity3D mGenericSmallStrainViscoplasticity3D;

	const GenericSmallStrainIsotropicPlasticity3D <GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<VonMisesPlasticPotential>>> mSmallStrainIsotropicPlasticity3DVonMisesVonMises;
    const GenericSmallStrainIsotropicPlasticity3D <GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<ModifiedMohrCoulombPlasticPotential>>> mSmallStrainIsotropicPlasticity3DVonMisesModifiedMohrCoulomb;
    const GenericSmallStrainIsotropicPlasticity3D <GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<DruckerPragerPlasticPotential>>> mSmallStrainIsotropicPlasticity3DVonMisesDruckerPrager;
    const GenericSmallStrainIsotropicPlasticity3D <GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<TrescaPlasticPotential>>> mSmallStrainIsotropicPlasticity3DVonMisesTresca;
    const GenericSmallStrainIsotropicPlasticity3D <GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential>>> mSmallStrainIsotropicPlasticity3DModifiedMohrCoulombVonMises;
    const GenericSmallStrainIsotropicPlasticity3D <GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential>>> mSmallStrainIsotropicPlasticity3DModifiedMohrCoulombModifiedMohrCoulomb;
    const GenericSmallStrainIsotropicPlasticity3D <GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<DruckerPragerPlasticPotential>>> mSmallStrainIsotropicPlasticity3DModifiedMohrCoulombDruckerPrager;
    const GenericSmallStrainIsotropicPlasticity3D <GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<TrescaPlasticPotential>>> mSmallStrainIsotropicPlasticity3DModifiedMohrCoulombTresca;
    const GenericSmallStrainIsotropicPlasticity3D <GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<VonMisesPlasticPotential>>> mSmallStrainIsotropicPlasticity3DTrescaVonMises;
    const GenericSmallStrainIsotropicPlasticity3D <GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<ModifiedMohrCoulombPlasticPotential>>> mSmallStrainIsotropicPlasticity3DTrescaModifiedMohrCoulomb;
    const GenericSmallStrainIsotropicPlasticity3D <GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<DruckerPragerPlasticPotential>>> mSmallStrainIsotropicPlasticity3DTrescaDruckerPrager;
    const GenericSmallStrainIsotropicPlasticity3D <GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<TrescaPlasticPotential>>> mSmallStrainIsotropicPlasticity3DTrescaTresca;
    const GenericSmallStrainIsotropicPlasticity3D <GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<VonMisesPlasticPotential>>> mSmallStrainIsotropicPlasticity3DDruckerPragerVonMises;
    const GenericSmallStrainIsotropicPlasticity3D <GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<ModifiedMohrCoulombPlasticPotential>>> mSmallStrainIsotropicPlasticity3DDruckerPragerModifiedMohrCoulomb;
    const GenericSmallStrainIsotropicPlasticity3D <GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<DruckerPragerPlasticPotential>>> mSmallStrainIsotropicPlasticity3DDruckerPragerDruckerPrager;
    const GenericSmallStrainIsotropicPlasticity3D <GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<TrescaPlasticPotential>>> mSmallStrainIsotropicPlasticity3DDruckerPragerTresca;
    const GenericSmallStrainIsotropicDamage3D <GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<VonMisesPlasticPotential>>> mSmallStrainIsotropicDamage3DVonMisesVonMises;
    const GenericSmallStrainIsotropicDamage3D <GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<ModifiedMohrCoulombPlasticPotential>>> mSmallStrainIsotropicDamage3DVonMisesModifiedMohrCoulomb;
    const GenericSmallStrainIsotropicDamage3D <GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<DruckerPragerPlasticPotential>>> mSmallStrainIsotropicDamage3DVonMisesDruckerPrager;
    const GenericSmallStrainIsotropicDamage3D <GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<TrescaPlasticPotential>>> mSmallStrainIsotropicDamage3DVonMisesTresca;
    const GenericSmallStrainIsotropicDamage3D <GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential>>> mSmallStrainIsotropicDamage3DModifiedMohrCoulombVonMises;
    const GenericSmallStrainIsotropicDamage3D <GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential>>> mSmallStrainIsotropicDamage3DModifiedMohrCoulombModifiedMohrCoulomb;
    const GenericSmallStrainIsotropicDamage3D <GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<DruckerPragerPlasticPotential>>> mSmallStrainIsotropicDamage3DModifiedMohrCoulombDruckerPrager;
    const GenericSmallStrainIsotropicDamage3D <GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<TrescaPlasticPotential>>> mSmallStrainIsotropicDamage3DModifiedMohrCoulombTresca;
    const GenericSmallStrainIsotropicDamage3D <GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<VonMisesPlasticPotential>>> mSmallStrainIsotropicDamage3DTrescaVonMises;
    const GenericSmallStrainIsotropicDamage3D <GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<ModifiedMohrCoulombPlasticPotential>>> mSmallStrainIsotropicDamage3DTrescaModifiedMohrCoulomb;
    const GenericSmallStrainIsotropicDamage3D <GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<DruckerPragerPlasticPotential>>> mSmallStrainIsotropicDamage3DTrescaDruckerPrager;
    const GenericSmallStrainIsotropicDamage3D <GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<TrescaPlasticPotential>>> mSmallStrainIsotropicDamage3DTrescaTresca;
    const GenericSmallStrainIsotropicDamage3D <GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential>>> mSmallStrainIsotropicDamage3DDruckerPragerVonMises;
    const GenericSmallStrainIsotropicDamage3D <GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<ModifiedMohrCoulombPlasticPotential>>> mSmallStrainIsotropicDamage3DDruckerPragerModifiedMohrCoulomb;
    const GenericSmallStrainIsotropicDamage3D <GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<DruckerPragerPlasticPotential>>> mSmallStrainIsotropicDamage3DDruckerPragerDruckerPrager;
    const GenericSmallStrainIsotropicDamage3D <GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<TrescaPlasticPotential>>> mSmallStrainIsotropicDamage3DDruckerPragerTresca;
    const GenericSmallStrainIsotropicDamage3D <GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<VonMisesPlasticPotential>>> mSmallStrainIsotropicDamage3DRankineVonMises;
    const GenericSmallStrainIsotropicDamage3D <GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<ModifiedMohrCoulombPlasticPotential>>> mSmallStrainIsotropicDamage3DRankineModifiedMohrCoulomb;
    const GenericSmallStrainIsotropicDamage3D <GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<DruckerPragerPlasticPotential>>> mSmallStrainIsotropicDamage3DRankineDruckerPrager;
    const GenericSmallStrainIsotropicDamage3D <GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<TrescaPlasticPotential>>> mSmallStrainIsotropicDamage3DRankineTresca;
    const GenericSmallStrainIsotropicDamage3D <GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<VonMisesPlasticPotential>>> mSmallStrainIsotropicDamage3DSimoJuVonMises;
    const GenericSmallStrainIsotropicDamage3D <GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<ModifiedMohrCoulombPlasticPotential>>> mSmallStrainIsotropicDamage3DSimoJuModifiedMohrCoulomb;
    const GenericSmallStrainIsotropicDamage3D <GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<DruckerPragerPlasticPotential>>> mSmallStrainIsotropicDamage3DSimoJuDruckerPrager;
    const GenericSmallStrainIsotropicDamage3D <GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<TrescaPlasticPotential>>> mSmallStrainIsotropicDamage3DSimoJuTresca;

	
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


