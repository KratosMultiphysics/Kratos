//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

#if !defined(KRATOS_RANS_APPLICATION_H_INCLUDED)
#define KRATOS_RANS_APPLICATION_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/kratos_application.h"

// Application includes

// stabilization validaton elements
#include "custom_elements/data_containers/stabilization_validation/circular_convection_element_data.h"
#include "custom_elements/data_containers/stabilization_validation/body_force_governed_cdr_element_data.h"
#include "custom_elements/data_containers/stabilization_validation/diffusion_element_data.h"

// incompressible potential flow elements
#include "custom_elements/incompressible_potential_flow_velocity_element.h"

// stabilized generic convection diffusion reaction elements
#include "custom_elements/convection_diffusion_reaction_cross_wind_stabilized_element.h"
#include "custom_elements/convection_diffusion_reaction_element.h"
#include "custom_elements/convection_diffusion_reaction_residual_based_flux_corrected_element.h"

// k-epsilon turbulence model element data
#include "custom_elements/data_containers/k_epsilon/epsilon_element_data.h"
#include "custom_elements/data_containers/k_epsilon/k_element_data.h"

// k-omega turbulence model element data
#include "custom_elements/data_containers/k_omega/k_element_data.h"
#include "custom_elements/data_containers/k_omega/omega_element_data.h"

// k-omega-sst turbulence model element data
#include "custom_elements/data_containers/k_omega_sst/k_element_data.h"
#include "custom_elements/data_containers/k_omega_sst/omega_element_data.h"

// vms monolithic wall conditions
#include "custom_conditions/vms_monolithic_u_based_wall_condition.h"
#include "custom_conditions/vms_monolithic_k_based_wall_condition.h"

// fractional step wall conditions
#include "custom_conditions/fractional_step_k_based_wall_condition.h"

// incompressible potential flow conditions
#include "custom_conditions/incompressible_potential_flow_velocity_inlet_condition.h"

// generic scalar wall flux condition
#include "custom_conditions/scalar_wall_flux_condition.h"

// k-epsilon turbulence model condition data
#include "custom_conditions/data_containers/k_epsilon/k_vis_based_wall_condition_data.h"
#include "custom_conditions/data_containers/k_epsilon/epsilon_k_based_wall_condition_data.h"
#include "custom_conditions/data_containers/k_epsilon/epsilon_u_based_wall_condition_data.h"

// k-omega turbulence model condition data
#include "custom_conditions/data_containers/k_omega/omega_k_based_wall_condition_data.h"
#include "custom_conditions/data_containers/k_omega/omega_u_based_wall_condition_data.h"
#include "custom_conditions/data_containers/k_omega/omega_vis_log_law_based_wall_condition_data.h"

// k-omega-sst turbulence model condition data
#include "custom_conditions/data_containers/k_omega_sst/omega_k_based_wall_condition_data.h"
#include "custom_conditions/data_containers/k_omega_sst/omega_u_based_wall_condition_data.h"

// constitutive laws
#include "custom_constitutive/newtonian_2d_law.h"
#include "custom_constitutive/newtonian_3d_law.h"
#include "custom_constitutive/rans_k_epsilon_newtonian_law.h"
#include "custom_constitutive/rans_k_omega_newtonian_law.h"
#include "custom_constitutive/rans_k_omega_sst_newtonian_law.h"
#include "custom_constitutive/rans_frozen_turbulence_newtonian_law.h"

// adjoint elements
#include "custom_elements/scalar_equation_adjoint_element.h"
#include "custom_elements/two_equation_turbulence_model_adjoint_element.h"

// adjoint element data containers
#include "custom_elements/data_containers/stabilization_validation/circular_convection_rfc_adjoint_element_data.h"
#include "custom_elements/data_containers/stabilization_validation/diffusion_rfc_adjoint_element_data.h"
#include "custom_elements/data_containers/k_epsilon/qsvms_rfc_adjoint_element_data.h"
#include "custom_elements/data_containers/k_omega/qsvms_rfc_adjoint_element_data.h"
#include "custom_elements/data_containers/k_omega_sst/qsvms_rfc_adjoint_element_data.h"

// adjoint conditions
#include "custom_conditions/two_equation_turbulence_model_adjoint_condition.h"
#include "custom_conditions/scalar_equation_adjoint_condition.h"

// adjoint condition data containers
// k-epsilon condition data containers
#include "custom_conditions/data_containers/k_epsilon/vms_monolithic_kk_based_epsilon_wall_condition_data.h"
#include "custom_conditions/data_containers/k_epsilon/vms_monolithic_ku_based_epsilon_wall_condition_data.h"
// k-omega condition data containers
#include "custom_conditions/data_containers/k_omega/vms_monolithic_kk_based_omega_wall_condition_data.h"
#include "custom_conditions/data_containers/k_omega/vms_monolithic_ku_based_omega_wall_condition_data.h"
// k-omega-sst condition data containers
#include "custom_conditions/data_containers/k_omega_sst/vms_monolithic_kk_based_omega_wall_condition_data.h"
#include "custom_conditions/data_containers/k_omega_sst/vms_monolithic_ku_based_omega_wall_condition_data.h"

namespace Kratos
{
///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
 */
class KRATOS_API(RANS_APPLICATION) KratosRANSApplication : public KratosApplication
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of KratosRANSApplication
    KRATOS_CLASS_POINTER_DEFINITION(KratosRANSApplication);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    KratosRANSApplication();

    /// Destructor.
    ~KratosRANSApplication() override = default;

    ///@}
    ///@name Operations
    ///@{

    void Register() override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "KratosRANSApplication";
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
        KRATOS_WATCH(KratosComponents<VariableData>::GetComponents().size());

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

private:
    ///@name Member Variables
    ///@{

    /// stabilization validation elements
    const ConvectionDiffusionReactionElement<2, 3, StabilizationValidationElementData::CircularConvectionElementData> mRansCircularConvectionAFC2D;
    const ConvectionDiffusionReactionCrossWindStabilizedElement<2, 3, StabilizationValidationElementData::CircularConvectionElementData> mRansCircularConvectionCWD2D;
    const ConvectionDiffusionReactionResidualBasedFluxCorrectedElement<2, 3, StabilizationValidationElementData::CircularConvectionElementData> mRansCircularConvectionRFC2D;

    const ConvectionDiffusionReactionElement<2, 3, StabilizationValidationElementData::BodyForceGovernedCDRElementData> mRansBodyForceGovernedCDRAFC2D;
    const ConvectionDiffusionReactionCrossWindStabilizedElement<2, 3, StabilizationValidationElementData::BodyForceGovernedCDRElementData> mRansBodyForceGovernedCDRCWD2D;
    const ConvectionDiffusionReactionResidualBasedFluxCorrectedElement<2, 3, StabilizationValidationElementData::BodyForceGovernedCDRElementData> mRansBodyForceGovernedCDRRFC2D;

    const ConvectionDiffusionReactionElement<2, 3, StabilizationValidationElementData::DiffusionElementData> mRansDiffusionAFC2D;
    const ConvectionDiffusionReactionCrossWindStabilizedElement<2, 3, StabilizationValidationElementData::DiffusionElementData> mRansDiffusionCWD2D;
    const ConvectionDiffusionReactionResidualBasedFluxCorrectedElement<2, 3, StabilizationValidationElementData::DiffusionElementData> mRansDiffusionRFC2D;

    /// incompressible potential flow elements
    const IncompressiblePotentialFlowVelocityElement<2, 3> mIncompressiblePotentialFlowVelocity2D;
    const IncompressiblePotentialFlowVelocityElement<3, 4> mIncompressiblePotentialFlowVelocity3D;

    /// k-epsilon turbulence model elements
    /// Algebraic flux correction based elements
    const ConvectionDiffusionReactionElement<2, 3, KEpsilonElementData::KElementData<2>> mRansKEpsilonKAFC2D;
    const ConvectionDiffusionReactionElement<3, 4, KEpsilonElementData::KElementData<3>> mRansKEpsilonKAFC3D;

    const ConvectionDiffusionReactionElement<2, 3, KEpsilonElementData::EpsilonElementData<2>> mRansKEpsilonEpsilonAFC2D;
    const ConvectionDiffusionReactionElement<3, 4, KEpsilonElementData::EpsilonElementData<3>> mRansKEpsilonEpsilonAFC3D;

    /// Residual based flux corrected elements
    const ConvectionDiffusionReactionResidualBasedFluxCorrectedElement<2, 3, KEpsilonElementData::KElementData<2>> mRansKEpsilonKRFC2D;
    const ConvectionDiffusionReactionResidualBasedFluxCorrectedElement<3, 4, KEpsilonElementData::KElementData<3>> mRansKEpsilonKRFC3D;

    const ConvectionDiffusionReactionResidualBasedFluxCorrectedElement<2, 3, KEpsilonElementData::EpsilonElementData<2>> mRansKEpsilonEpsilonRFC2D;
    const ConvectionDiffusionReactionResidualBasedFluxCorrectedElement<3, 4, KEpsilonElementData::EpsilonElementData<3>> mRansKEpsilonEpsilonRFC3D;

    /// Cross wind stabilization based elements
    const ConvectionDiffusionReactionCrossWindStabilizedElement<2, 3, KEpsilonElementData::KElementData<2>> mRansKEpsilonKCWD2D;
    const ConvectionDiffusionReactionCrossWindStabilizedElement<3, 4, KEpsilonElementData::KElementData<3>> mRansKEpsilonKCWD3D;

    const ConvectionDiffusionReactionCrossWindStabilizedElement<2, 3, KEpsilonElementData::EpsilonElementData<2>> mRansKEpsilonEpsilonCWD2D;
    const ConvectionDiffusionReactionCrossWindStabilizedElement<3, 4, KEpsilonElementData::EpsilonElementData<3>> mRansKEpsilonEpsilonCWD3D;

    /// k-omega turbulence model elements
    /// Algebraic flux correction based elements
    const ConvectionDiffusionReactionElement<2, 3, KOmegaElementData::KElementData<2>> mRansKOmegaKAFC2D;
    const ConvectionDiffusionReactionElement<3, 4, KOmegaElementData::KElementData<3>> mRansKOmegaKAFC3D;

    const ConvectionDiffusionReactionElement<2, 3, KOmegaElementData::OmegaElementData<2>> mRansKOmegaOmegaAFC2D;
    const ConvectionDiffusionReactionElement<3, 4, KOmegaElementData::OmegaElementData<3>> mRansKOmegaOmegaAFC3D;

    /// Residual based flux corrected elements
    const ConvectionDiffusionReactionResidualBasedFluxCorrectedElement<2, 3, KOmegaElementData::KElementData<2>> mRansKOmegaKRFC2D;
    const ConvectionDiffusionReactionResidualBasedFluxCorrectedElement<3, 4, KOmegaElementData::KElementData<3>> mRansKOmegaKRFC3D;

    const ConvectionDiffusionReactionResidualBasedFluxCorrectedElement<2, 3, KOmegaElementData::OmegaElementData<2>> mRansKOmegaOmegaRFC2D;
    const ConvectionDiffusionReactionResidualBasedFluxCorrectedElement<3, 4, KOmegaElementData::OmegaElementData<3>> mRansKOmegaOmegaRFC3D;

    /// Cross wind stabilization based elements
    const ConvectionDiffusionReactionCrossWindStabilizedElement<2, 3, KOmegaElementData::KElementData<2>> mRansKOmegaKCWD2D;
    const ConvectionDiffusionReactionCrossWindStabilizedElement<3, 4, KOmegaElementData::KElementData<3>> mRansKOmegaKCWD3D;

    const ConvectionDiffusionReactionCrossWindStabilizedElement<2, 3, KOmegaElementData::OmegaElementData<2>> mRansKOmegaOmegaCWD2D;
    const ConvectionDiffusionReactionCrossWindStabilizedElement<3, 4, KOmegaElementData::OmegaElementData<3>> mRansKOmegaOmegaCWD3D;

    /// k-omega-sst turbulence model elements
    /// Algebraic flux correction based elements
    const ConvectionDiffusionReactionElement<2, 3, KOmegaSSTElementData::KElementData<2>> mRansKOmegaSSTKAFC2D;
    const ConvectionDiffusionReactionElement<3, 4, KOmegaSSTElementData::KElementData<3>> mRansKOmegaSSTKAFC3D;

    const ConvectionDiffusionReactionElement<2, 3, KOmegaSSTElementData::OmegaElementData<2>> mRansKOmegaSSTOmegaAFC2D;
    const ConvectionDiffusionReactionElement<3, 4, KOmegaSSTElementData::OmegaElementData<3>> mRansKOmegaSSTOmegaAFC3D;

    /// Residual based flux corrected elements
    const ConvectionDiffusionReactionResidualBasedFluxCorrectedElement<2, 3, KOmegaSSTElementData::KElementData<2>> mRansKOmegaSSTKRFC2D;
    const ConvectionDiffusionReactionResidualBasedFluxCorrectedElement<3, 4, KOmegaSSTElementData::KElementData<3>> mRansKOmegaSSTKRFC3D;

    const ConvectionDiffusionReactionResidualBasedFluxCorrectedElement<2, 3, KOmegaSSTElementData::OmegaElementData<2>> mRansKOmegaSSTOmegaRFC2D;
    const ConvectionDiffusionReactionResidualBasedFluxCorrectedElement<3, 4, KOmegaSSTElementData::OmegaElementData<3>> mRansKOmegaSSTOmegaRFC3D;

    /// Cross wind stabilization based elements
    const ConvectionDiffusionReactionCrossWindStabilizedElement<2, 3, KOmegaSSTElementData::KElementData<2>> mRansKOmegaSSTKCWD2D;
    const ConvectionDiffusionReactionCrossWindStabilizedElement<3, 4, KOmegaSSTElementData::KElementData<3>> mRansKOmegaSSTKCWD3D;

    const ConvectionDiffusionReactionCrossWindStabilizedElement<2, 3, KOmegaSSTElementData::OmegaElementData<2>> mRansKOmegaSSTOmegaCWD2D;
    const ConvectionDiffusionReactionCrossWindStabilizedElement<3, 4, KOmegaSSTElementData::OmegaElementData<3>> mRansKOmegaSSTOmegaCWD3D;

    // vms monolithic k based wall conditions
    const VMSMonolithicKBasedWallCondition<2> mRansVMSMonolithicKBasedWall2D2N;
    const VMSMonolithicKBasedWallCondition<3> mRansVMSMonolithicKBasedWall3D3N;

    const VMSMonolithicUBasedWallCondition<2> mRansVMSMonolithicUBasedWall2D2N;
    const VMSMonolithicUBasedWallCondition<3> mRansVMSMonolithicUBasedWall3D3N;

    // fractional step wall conditions
    const FractionalStepKBasedWallCondition<2, 2> mFractionalStepKBasedWall2D2N;
    const FractionalStepKBasedWallCondition<3, 3> mFractionalStepKBasedWall3D3N;

    // incompressible potential flow conditions
    const IncompressiblePotentialFlowVelocityInletCondition<2, 2> mIncompressiblePotentialFlowVelocityInlet2D2N;
    const IncompressiblePotentialFlowVelocityInletCondition<3, 3> mIncompressiblePotentialFlowVelocityInlet3D3N;

    // k-epsilon turbulence model conditions
    const ScalarWallFluxCondition<2, 2, KEpsilonWallConditionData::KVisBasedWallConditionData> mRansKEpsilonKVisBasedWall2D2N;
    const ScalarWallFluxCondition<3, 3, KEpsilonWallConditionData::KVisBasedWallConditionData> mRansKEpsilonKVisBasedWall3D3N;

    const ScalarWallFluxCondition<2, 2, KEpsilonWallConditionData::EpsilonKBasedWallConditionData> mRansKEpsilonEpsilonKBasedWall2D2N;
    const ScalarWallFluxCondition<3, 3, KEpsilonWallConditionData::EpsilonKBasedWallConditionData> mRansKEpsilonEpsilonKBasedWall3D3N;

    const ScalarWallFluxCondition<2, 2, KEpsilonWallConditionData::EpsilonUBasedWallConditionData> mRansKEpsilonEpsilonUBasedWall2D2N;
    const ScalarWallFluxCondition<3, 3, KEpsilonWallConditionData::EpsilonUBasedWallConditionData> mRansKEpsilonEpsilonUBasedWall3D3N;

    // k-omega turbulence model conditions
    const ScalarWallFluxCondition<2, 2, KOmegaWallConditionData::OmegaKBasedWallConditionData> mRansKOmegaOmegaKBasedWall2D2N;
    const ScalarWallFluxCondition<3, 3, KOmegaWallConditionData::OmegaKBasedWallConditionData> mRansKOmegaOmegaKBasedWall3D3N;

    const ScalarWallFluxCondition<2, 2, KOmegaWallConditionData::OmegaUBasedWallConditionData> mRansKOmegaOmegaUBasedWall2D2N;
    const ScalarWallFluxCondition<3, 3, KOmegaWallConditionData::OmegaUBasedWallConditionData> mRansKOmegaOmegaUBasedWall3D3N;

    const ScalarWallFluxCondition<2, 2, KOmegaWallConditionData::OmegaVisLogBasedWallConditionData> mRansKOmegaOmegaVisLogBasedWall2D2N;
    const ScalarWallFluxCondition<3, 3, KOmegaWallConditionData::OmegaVisLogBasedWallConditionData> mRansKOmegaOmegaVisLogBasedWall3D3N;

    // k-omega-sst turbulence model conditions
    const ScalarWallFluxCondition<2, 2, KOmegaSSTWallConditionData::OmegaKBasedWallConditionData> mRansKOmegaSSTOmegaKBasedWall2D2N;
    const ScalarWallFluxCondition<3, 3, KOmegaSSTWallConditionData::OmegaKBasedWallConditionData> mRansKOmegaSSTOmegaKBasedWall3D3N;

    const ScalarWallFluxCondition<2, 2, KOmegaSSTWallConditionData::OmegaUBasedWallConditionData> mRansKOmegaSSTOmegaUBasedWall2D2N;
    const ScalarWallFluxCondition<3, 3, KOmegaSSTWallConditionData::OmegaUBasedWallConditionData> mRansKOmegaSSTOmegaUBasedWall3D3N;

    // constitutive laws
    const RansKEpsilonNewtonianLaw<Newtonian2DLaw> mRansKEpsilonNewtonian2DLaw;
    const RansKEpsilonNewtonianLaw<Newtonian3DLaw> mRansKEpsilonNewtonian3DLaw;

    const RansKOmegaNewtonianLaw<Newtonian2DLaw> mRansKOmegaNewtonian2DLaw;
    const RansKOmegaNewtonianLaw<Newtonian3DLaw> mRansKOmegaNewtonian3DLaw;

    const RansKOmegaSSTNewtonianLaw<2, Newtonian2DLaw> mRansKOmegaSSTNewtonian2DLaw;
    const RansKOmegaSSTNewtonianLaw<3, Newtonian3DLaw> mRansKOmegaSSTNewtonian3DLaw;

    const RansFrozenTurbulenceNewtonianLaw<RansKEpsilonNewtonianLaw<Newtonian2DLaw>> mRansFrozenTurbulenceKEpsilonNewtonian2DLaw;
    const RansFrozenTurbulenceNewtonianLaw<RansKEpsilonNewtonianLaw<Newtonian3DLaw>> mRansFrozenTurbulenceKEpsilonNewtonian3DLaw;

    const RansFrozenTurbulenceNewtonianLaw<RansKOmegaNewtonianLaw<Newtonian2DLaw>> mRansFrozenTurbulenceKOmegaNewtonian2DLaw;
    const RansFrozenTurbulenceNewtonianLaw<RansKOmegaNewtonianLaw<Newtonian3DLaw>> mRansFrozenTurbulenceKOmegaNewtonian3DLaw;

    const RansFrozenTurbulenceNewtonianLaw<RansKOmegaSSTNewtonianLaw<2, Newtonian2DLaw>> mRansFrozenTurbulenceKOmegaSSTNewtonian2DLaw;
    const RansFrozenTurbulenceNewtonianLaw<RansKOmegaSSTNewtonianLaw<3, Newtonian3DLaw>> mRansFrozenTurbulenceKOmegaSSTNewtonian3DLaw;

    // stabilization validation adjoint elements
    const ScalarEquationAdjointElement<2, 3, StabilizationValidationElementData::CircularConvectionRFCAdjointElementData> mRansCircularConvectionRFCAdjoint2D3N;
    const ScalarEquationAdjointElement<2, 3, StabilizationValidationElementData::DiffusionRFCAdjointElementData> mRansDiffusionRFCAdjoint2D3N;

    // stabilization validation adjoint conditions
    const ScalarEquationAdjointCondition<2, 2> mRansScalarEquationAdjoint2D2N;

    // k-epsilon turbulence model adjoint elements
    const TwoEquationTurbulenceModelAdjointElement<2, 3, KEpsilonElementData::QSVMSRFCAdjointElementData<2, 3>> mRansKEpsilonQSVMSRFCAdjoint2D3N;
    const TwoEquationTurbulenceModelAdjointElement<3, 4, KEpsilonElementData::QSVMSRFCAdjointElementData<3, 4>> mRansKEpsilonQSVMSRFCAdjoint3D4N;

    // k-omega turbulence model adjoint elements
    const TwoEquationTurbulenceModelAdjointElement<2, 3, KOmegaElementData::QSVMSRFCAdjointElementData<2, 3>> mRansKOmegaQSVMSRFCAdjoint2D3N;
    const TwoEquationTurbulenceModelAdjointElement<3, 4, KOmegaElementData::QSVMSRFCAdjointElementData<3, 4>> mRansKOmegaQSVMSRFCAdjoint3D4N;

    // k-omega-sst turbulence model adjoint elements
    const TwoEquationTurbulenceModelAdjointElement<2, 3, KOmegaSSTElementData::QSVMSRFCAdjointElementData<2, 3>> mRansKOmegaSSTQSVMSRFCAdjoint2D3N;
    const TwoEquationTurbulenceModelAdjointElement<3, 4, KOmegaSSTElementData::QSVMSRFCAdjointElementData<3, 4>> mRansKOmegaSSTQSVMSRFCAdjoint3D4N;

    // k-epsilon turbulence model adjoint conditions
    const TwoEquationTurbulenceModelAdjointCondition<2, 2, KEpsilonWallConditionData::VMSMonolithicKBasedEpsilonKBasedWallConditionData<2, 2>> mRansKEpsilonVMSKBasedEpsilonKBasedWallAdjoint2D2N;
    const TwoEquationTurbulenceModelAdjointCondition<3, 3, KEpsilonWallConditionData::VMSMonolithicKBasedEpsilonKBasedWallConditionData<3, 3>> mRansKEpsilonVMSKBasedEpsilonKBasedWallAdjoint3D3N;

    const TwoEquationTurbulenceModelAdjointCondition<2, 2, KEpsilonWallConditionData::VMSMonolithicKBasedEpsilonUBasedWallConditionData<2, 2>> mRansKEpsilonVMSKBasedEpsilonUBasedWallAdjoint2D2N;
    const TwoEquationTurbulenceModelAdjointCondition<3, 3, KEpsilonWallConditionData::VMSMonolithicKBasedEpsilonUBasedWallConditionData<3, 3>> mRansKEpsilonVMSKBasedEpsilonUBasedWallAdjoint3D3N;

    // k-omega turbulence model adjoint conditions
    const TwoEquationTurbulenceModelAdjointCondition<2, 2, KOmegaWallConditionData::VMSMonolithicKBasedOmegaKBasedWallConditionData<2, 2>> mRansKOmegaVMSKBasedOmegaKBasedWallAdjoint2D2N;
    const TwoEquationTurbulenceModelAdjointCondition<3, 3, KOmegaWallConditionData::VMSMonolithicKBasedOmegaKBasedWallConditionData<3, 3>> mRansKOmegaVMSKBasedOmegaKBasedWallAdjoint3D3N;

    const TwoEquationTurbulenceModelAdjointCondition<2, 2, KOmegaWallConditionData::VMSMonolithicKBasedOmegaUBasedWallConditionData<2, 2>> mRansKOmegaVMSKBasedOmegaUBasedWallAdjoint2D2N;
    const TwoEquationTurbulenceModelAdjointCondition<3, 3, KOmegaWallConditionData::VMSMonolithicKBasedOmegaUBasedWallConditionData<3, 3>> mRansKOmegaVMSKBasedOmegaUBasedWallAdjoint3D3N;

    // k-omega-sst turbulence model adjoint conditions
    const TwoEquationTurbulenceModelAdjointCondition<2, 2, KOmegaSSTWallConditionData::VMSMonolithicKBasedOmegaKBasedWallConditionData<2, 2>> mRansKOmegaSSTVMSKBasedOmegaKBasedWallAdjoint2D2N;
    const TwoEquationTurbulenceModelAdjointCondition<3, 3, KOmegaSSTWallConditionData::VMSMonolithicKBasedOmegaKBasedWallConditionData<3, 3>> mRansKOmegaSSTVMSKBasedOmegaKBasedWallAdjoint3D3N;

    const TwoEquationTurbulenceModelAdjointCondition<2, 2, KOmegaSSTWallConditionData::VMSMonolithicKBasedOmegaUBasedWallConditionData<2, 2>> mRansKOmegaSSTVMSKBasedOmegaUBasedWallAdjoint2D2N;
    const TwoEquationTurbulenceModelAdjointCondition<3, 3, KOmegaSSTWallConditionData::VMSMonolithicKBasedOmegaUBasedWallConditionData<3, 3>> mRansKOmegaSSTVMSKBasedOmegaUBasedWallAdjoint3D3N;

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    KratosRANSApplication& operator=(KratosRANSApplication const& rOther);

    /// Copy constructor.
    KratosRANSApplication(KratosRANSApplication const& rOther);

    ///@}

}; // Class KratosRANSApplication

///@}

} // namespace Kratos.

#endif // KRATOS_RANS_APPLICATION_H_INCLUDED  defined
