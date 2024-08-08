// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Riccardo Rossi
//  Co-authors:      Vicente Mataix Ferrandiz
//

#pragma once

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"

/* ELEMENTS */

/* 0D elements */
#include "custom_elements/nodal_concentrated_element.hpp"

/* Mass elements */
#include "custom_elements/mass_element.h"

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
#include "custom_response_functions/adjoint_elements/adjoint_finite_difference_spring_damper_element_3D2N.h"

/* Adding shells and membranes elements */
#include "custom_elements/isotropic_shell_element.hpp"
#include "custom_elements/membrane_element.hpp"
#include "custom_elements/membrane_element_2D2N.h"
#include "custom_elements/shell_thick_element_3D4N.hpp"
#include "custom_elements/shell_thin_element_3D4N.hpp"
#include "custom_elements/shell_thin_element_3D3N.hpp"
#include "custom_elements/shell_thick_element_3D3N.hpp"


/* Adding the bushing element */
#include "custom_elements/bushing_element.h"

/* Adding the spring damper element */
#include "custom_elements/spring_damper_element.hpp"

/* Adding the SPRISM element */
#include "custom_elements/solid_shell_element_sprism_3D6N.h"

/* Adding solid elements */
#include "custom_elements/small_displacement.h"
#include "custom_elements/axisym_small_displacement.h"
#include "custom_elements/z_strain_driven_2p5_small_displacement.h"
#include "custom_elements/total_lagrangian.h"
#include "custom_elements/axisym_total_lagrangian.h"
#include "custom_elements/updated_lagrangian.h"
#include "custom_elements/axisym_updated_lagrangian.h"
#include "custom_elements/small_displacement_bbar.h"
#include "custom_elements/small_displacement_shifted_boundary_element.h"

/* Adding the mixed solid elements */
#include "custom_elements/small_displacement_mixed_volumetric_strain_element.h"
#include "custom_elements/small_displacement_mixed_volumetric_strain_oss_element.h"
#include "custom_elements/small_displacement_mixed_volumetric_strain_oss_non_linear_element.h"
#include "custom_elements/total_lagrangian_mixed_volumetric_strain_element.h"
#include "custom_elements/total_lagrangian_q1p0_mixed_element.h"
#include "custom_elements/timoshenko_beam_element_2D2N.h"
#include "custom_elements/timoshenko_beam_element_2D3N.h"
#include "custom_elements/timoshenko_curved_beam_element_2D3N.h"

/* Conditions */
#include "custom_conditions/base_load_condition.h"
#include "custom_conditions/point_load_condition.h"
#include "custom_conditions/point_contact_condition.h"
#include "custom_conditions/axisym_point_load_condition.h"
#include "custom_conditions/line_load_condition.h"
#include "custom_conditions/small_displacement_line_load_condition.h"
#include "custom_conditions/axisym_line_load_condition_2d.h"
#include "custom_conditions/surface_load_condition_3d.h"
#include "custom_conditions/small_displacement_surface_load_condition_3d.h"
#include "custom_conditions/point_moment_condition_3d.h"
#include "custom_conditions/displacement_control_condition.h"
#include "custom_conditions/moving_load_condition.h"

/* Adding the displacement-based SBM condition */
#include "custom_conditions/displacement_shifted_boundary_condition.h"

/* Adding the adjoint conditions */
#include "custom_response_functions/adjoint_conditions/adjoint_semi_analytic_point_load_condition.h"
#include "custom_response_functions/adjoint_conditions/adjoint_semi_analytic_base_condition.h"

/* Constitutive Laws */
#include "custom_constitutive/truss_constitutive_law.h"
#include "custom_constitutive/beam_constitutive_law.h"
#include "custom_constitutive/elastic_isotropic_3d.h"
#include "custom_constitutive/axisym_elastic_isotropic.h"
#include "custom_constitutive/linear_plane_strain.h"
#include "custom_constitutive/linear_plane_stress.h"
#include "custom_constitutive/user_provided_linear_elastic_law.h"
// Constitutive laws for the Timoshenko beams
#include "custom_constitutive/timoshenko_beam_elastic_constitutive_law.h"


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
    const LinearTimoshenkoBeamElement2D2N mLinearTimoshenkoBeamElement2D2N;
    const LinearTimoshenkoBeamElement2D3N mLinearTimoshenkoBeamElement2D3N;
    const LinearTimoshenkoCurvedBeamElement2D3N mLinearTimoshenkoCurvedBeamElement2D3N;


    // Adding the shells elements
    const IsotropicShellElement mIsotropicShellElement3D3N;
    const ShellThickElement3D4N<ShellKinematics::LINEAR>                 mShellThickElement3D4N;
    const ShellThickElement3D4N<ShellKinematics::NONLINEAR_COROTATIONAL> mShellThickCorotationalElement3D4N;
    const ShellThinElement3D4N<ShellKinematics::NONLINEAR_COROTATIONAL>  mShellThinCorotationalElement3D4N;
    const ShellThinElement3D3N<ShellKinematics::LINEAR>                  mShellThinElement3D3N;
    const ShellThinElement3D3N<ShellKinematics::NONLINEAR_COROTATIONAL>  mShellThinCorotationalElement3D3N;
    const ShellThickElement3D3N<ShellKinematics::NONLINEAR_COROTATIONAL> mShellThickCorotationalElement3D3N;

    // Adding the membrane elements
    const MembraneElement mMembraneElement3D4N;
    const MembraneElement mMembraneElement3D3N;
    const MembraneElement2D2N mMembraneElement2D2N;

    // Adding the SPRISM element
    const SolidShellElementSprism3D6N mSolidShellElementSprism3D6N;

    // Adding the nodal concentrated element
    const NodalConcentratedElement mNodalConcentratedElement2D1N;
    const NodalConcentratedElement mNodalConcentratedDampedElement2D1N;
    const NodalConcentratedElement mNodalConcentratedElement3D1N;
    const NodalConcentratedElement mNodalConcentratedDampedElement3D1N;

    // Adding the mass elements
    const MassElement mLineMassElement3D2N;
    const MassElement mSurfaceMassElement3D3N;
    const MassElement mSurfaceMassElement3D4N;

    // Linear kinematic elements
    const SmallDisplacement mSmallDisplacement2D3N;
    const SmallDisplacement mSmallDisplacement2D4N;
    const SmallDisplacement mSmallDisplacement2D6N;
    const SmallDisplacement mSmallDisplacement2D8N;
    const SmallDisplacement mSmallDisplacement2D9N;
    const SmallDisplacement mSmallDisplacement2D10N;
    const SmallDisplacement mSmallDisplacement2D15N;
    const SmallDisplacement mSmallDisplacement3D4N;
    const SmallDisplacement mSmallDisplacement3D5N;
    const SmallDisplacement mSmallDisplacement3D6N;
    const SmallDisplacement mSmallDisplacement3D8N;
    const SmallDisplacement mSmallDisplacement3D10N;
    const SmallDisplacement mSmallDisplacement3D13N;
    const SmallDisplacement mSmallDisplacement3D15N;
    const SmallDisplacement mSmallDisplacement3D20N;
    const SmallDisplacement mSmallDisplacement3D27N;

    const SmallDisplacementBbar mSmallDisplacementBbar2D4N;
    const SmallDisplacementBbar mSmallDisplacementBbar3D8N;

    const SmallDisplacementShiftedBoundaryElement<2> mSmallDisplacementShiftedBoundaryElement2D3N;
    const SmallDisplacementShiftedBoundaryElement<3> mSmallDisplacementShiftedBoundaryElement3D4N;

    const SmallDisplacementMixedVolumetricStrainElement mSmallDisplacementMixedVolumetricStrainElement2D3N;
    const SmallDisplacementMixedVolumetricStrainElement mSmallDisplacementMixedVolumetricStrainElement2D4N;
    const SmallDisplacementMixedVolumetricStrainElement mSmallDisplacementMixedVolumetricStrainElement3D4N;
    const SmallDisplacementMixedVolumetricStrainElement mSmallDisplacementMixedVolumetricStrainElement3D8N;

    const SmallDisplacementMixedVolumetricStrainOssElement mSmallDisplacementMixedVolumetricStrainOssElement2D3N;
    const SmallDisplacementMixedVolumetricStrainOssElement mSmallDisplacementMixedVolumetricStrainOssElement2D4N;
    const SmallDisplacementMixedVolumetricStrainOssElement mSmallDisplacementMixedVolumetricStrainOssElement3D4N;
    const SmallDisplacementMixedVolumetricStrainOssElement mSmallDisplacementMixedVolumetricStrainOssElement3D8N;

    const SmallDisplacementMixedVolumetricStrainOssNonLinearElement mSmallDisplacementMixedVolumetricStrainOssNonLinearElement2D3N;
    const SmallDisplacementMixedVolumetricStrainOssNonLinearElement mSmallDisplacementMixedVolumetricStrainOssNonLinearElement2D4N;
    const SmallDisplacementMixedVolumetricStrainOssNonLinearElement mSmallDisplacementMixedVolumetricStrainOssNonLinearElement3D4N;
    const SmallDisplacementMixedVolumetricStrainOssNonLinearElement mSmallDisplacementMixedVolumetricStrainOssNonLinearElement3D8N;

    const TotalLagrangianMixedVolumetricStrainElement<2> mTotalLagrangianMixedVolumetricStrainElement2D3N;
    const TotalLagrangianMixedVolumetricStrainElement<3> mTotalLagrangianMixedVolumetricStrainElement3D4N;

    const TotalLagrangianQ1P0MixedElement mTotalLagrangianQ1P0MixedElement2D3N;
    const TotalLagrangianQ1P0MixedElement mTotalLagrangianQ1P0MixedElement2D4N;
    const TotalLagrangianQ1P0MixedElement mTotalLagrangianQ1P0MixedElement3D4N;
    const TotalLagrangianQ1P0MixedElement mTotalLagrangianQ1P0MixedElement3D10N;
    const TotalLagrangianQ1P0MixedElement mTotalLagrangianQ1P0MixedElement3D8N;
    const TotalLagrangianQ1P0MixedElement mTotalLagrangianQ1P0MixedElement3D20N;

    const AxisymSmallDisplacement mAxisymSmallDisplacement2D3N;
    const AxisymSmallDisplacement mAxisymSmallDisplacement2D4N;
    const AxisymSmallDisplacement mAxisymSmallDisplacement2D6N;
    const AxisymSmallDisplacement mAxisymSmallDisplacement2D8N;
    const AxisymSmallDisplacement mAxisymSmallDisplacement2D9N;

    const ZStrainDriven2p5DSmallDisplacement mZStrainDriven2p5DSmallDisplacement2D3N;
    const ZStrainDriven2p5DSmallDisplacement mZStrainDriven2p5DSmallDisplacement2D4N;
    const ZStrainDriven2p5DSmallDisplacement mZStrainDriven2p5DSmallDisplacement2D6N;
    const ZStrainDriven2p5DSmallDisplacement mZStrainDriven2p5DSmallDisplacement2D8N;
    const ZStrainDriven2p5DSmallDisplacement mZStrainDriven2p5DSmallDisplacement2D9N;

    // Total lagrangian
    const TotalLagrangian mTotalLagrangian2D3N;
    const TotalLagrangian mTotalLagrangian2D4N;
    const TotalLagrangian mTotalLagrangian2D6N;
    const TotalLagrangian mTotalLagrangian2D8N;
    const TotalLagrangian mTotalLagrangian2D9N;
    const TotalLagrangian mTotalLagrangian2D10N;
    const TotalLagrangian mTotalLagrangian2D15N;
    const TotalLagrangian mTotalLagrangian3D4N;
    const TotalLagrangian mTotalLagrangian3D5N;
    const TotalLagrangian mTotalLagrangian3D6N;
    const TotalLagrangian mTotalLagrangian3D8N;
    const TotalLagrangian mTotalLagrangian3D10N;
    const TotalLagrangian mTotalLagrangian3D13N;
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
    const UpdatedLagrangian mUpdatedLagrangian2D10N;
    const UpdatedLagrangian mUpdatedLagrangian2D15N;
    const UpdatedLagrangian mUpdatedLagrangian3D4N;
    const UpdatedLagrangian mUpdatedLagrangian3D5N;
    const UpdatedLagrangian mUpdatedLagrangian3D6N;
    const UpdatedLagrangian mUpdatedLagrangian3D8N;
    const UpdatedLagrangian mUpdatedLagrangian3D10N;
    const UpdatedLagrangian mUpdatedLagrangian3D13N;
    const UpdatedLagrangian mUpdatedLagrangian3D15N;
    const UpdatedLagrangian mUpdatedLagrangian3D20N;
    const UpdatedLagrangian mUpdatedLagrangian3D27N;

    const AxisymUpdatedLagrangian mAxisymUpdatedLagrangian2D3N;
    const AxisymUpdatedLagrangian mAxisymUpdatedLagrangian2D4N;
    const AxisymUpdatedLagrangian mAxisymUpdatedLagrangian2D6N;
    const AxisymUpdatedLagrangian mAxisymUpdatedLagrangian2D8N;
    const AxisymUpdatedLagrangian mAxisymUpdatedLagrangian2D9N;

    // Adding the bushing element
    const BushingElement mBushingElement3D2N;    

    // Adding the spring damper element
    const SpringDamperElement<2> mSpringDamperElement2D;
    const SpringDamperElement<3> mSpringDamperElement3D;

    // Adding adjoint elements
    const AdjointFiniteDifferencingShellElement<ShellThinElement3D3N<ShellKinematics::LINEAR>> mAdjointFiniteDifferencingShellThinElement3D3N;
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
    const AdjointFiniteDifferenceSpringDamperElement<SpringDamperElement<3>>  mAdjointFiniteDifferenceSpringDamperElement3D2N;

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
    const LineLoadCondition<2> mLineLoadCondition2D4N;
    const LineLoadCondition<2> mLineLoadCondition2D5N;
    const LineLoadCondition<3> mLineLoadCondition3D2N;
    const LineLoadCondition<3> mLineLoadCondition3D3N;

    const SmallDisplacementLineLoadCondition<2> mSmallDisplacementLineLoadCondition2D2N;
    const SmallDisplacementLineLoadCondition<2> mSmallDisplacementLineLoadCondition2D3N;
    const SmallDisplacementLineLoadCondition<2> mSmallDisplacementLineLoadCondition2D4N;
    const SmallDisplacementLineLoadCondition<2> mSmallDisplacementLineLoadCondition2D5N;
    const SmallDisplacementLineLoadCondition<3> mSmallDisplacementLineLoadCondition3D2N;
    const SmallDisplacementLineLoadCondition<3> mSmallDisplacementLineLoadCondition3D3N;

    const AxisymLineLoadCondition2D mAxisymLineLoadCondition2D2N;
    const AxisymLineLoadCondition2D mAxisymLineLoadCondition2D3N;

    // Surface load
    const SurfaceLoadCondition3D mSurfaceLoadCondition3D3N;
    const SurfaceLoadCondition3D mSurfaceLoadCondition3D4N;
    const SurfaceLoadCondition3D mSurfaceLoadCondition3D6N;
    const SurfaceLoadCondition3D mSurfaceLoadCondition3D8N;
    const SurfaceLoadCondition3D mSurfaceLoadCondition3D9N;

    const SmallDisplacementSurfaceLoadCondition3D mSmallDisplacementSurfaceLoadCondition3D3N;
    const SmallDisplacementSurfaceLoadCondition3D mSmallDisplacementSurfaceLoadCondition3D4N;
    const SmallDisplacementSurfaceLoadCondition3D mSmallDisplacementSurfaceLoadCondition3D6N;
    const SmallDisplacementSurfaceLoadCondition3D mSmallDisplacementSurfaceLoadCondition3D8N;
    const SmallDisplacementSurfaceLoadCondition3D mSmallDisplacementSurfaceLoadCondition3D9N;

    // Point moment
    const PointMomentCondition3D mPointMomentCondition3D1N;

    // Adjoint Conditions
    const AdjointSemiAnalyticPointLoadCondition<PointLoadCondition> mAdjointSemiAnalyticPointLoadCondition2D1N;
    const AdjointSemiAnalyticPointLoadCondition<PointLoadCondition> mAdjointSemiAnalyticPointLoadCondition3D1N;
    const AdjointSemiAnalyticBaseCondition<SurfaceLoadCondition3D> mAdjointSemiAnalyticSurfaceLoadCondition3D3N;
    const AdjointSemiAnalyticBaseCondition<SurfaceLoadCondition3D> mAdjointSemiAnalyticSurfaceLoadCondition3D4N;
    const AdjointSemiAnalyticBaseCondition<SmallDisplacementSurfaceLoadCondition3D> mAdjointSemiAnalyticSmallDisplacementSurfaceLoadCondition3D3N;
    const AdjointSemiAnalyticBaseCondition<SmallDisplacementSurfaceLoadCondition3D> mAdjointSemiAnalyticSmallDisplacementSurfaceLoadCondition3D4N;
    const AdjointSemiAnalyticBaseCondition<LineLoadCondition<3>> mAdjointSemiAnalyticLineLoadCondition3D2N;
    const AdjointSemiAnalyticBaseCondition<SmallDisplacementLineLoadCondition<3>> mAdjointSemiAnalyticSmallDisplacementLineLoadCondition3D2N;

    // Displacement-Control Conditions
    const DisplacementControlCondition mDisplacementControlCondition3D1N;

    // SBM displacement conditions
    const DisplacementShiftedBoundaryCondition mDisplacementShiftedBoundaryCondition;

    // Moving load
    const MovingLoadCondition<2,2> mMovingLoadCondition2D2N;
    const MovingLoadCondition<2, 3> mMovingLoadCondition2D3N;
    const MovingLoadCondition<3, 2> mMovingLoadCondition3D2N;
    const MovingLoadCondition<3, 3> mMovingLoadCondition3D3N;

    /* CONSTITUTIVE LAWS */
    // Linear elastics laws
    const TrussConstitutiveLaw mTrussConstitutiveLaw;
    const BeamConstitutiveLaw mBeamConstitutiveLaw;
    const ElasticIsotropic3D mElasticIsotropic3D;
    const AxisymElasticIsotropic mAxisymElasticIsotropic;
    const LinearPlaneStrain  mLinearPlaneStrain;
    const LinearPlaneStress  mLinearPlaneStress;
    const UserProvidedLinearElasticLaw<2> mUserProvidedLinearElastic2DLaw;
    const UserProvidedLinearElasticLaw<3> mUserProvidedLinearElastic3DLaw;
    const TimoshenkoBeamElasticConstitutiveLaw mTimoshenkoBeamElasticConstitutiveLaw;

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
