//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya (https://github.com/sunethwarna)
//

#if !defined(KRATOS_RANS_APPLICATION_H_INCLUDED)
#define KRATOS_RANS_APPLICATION_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/kratos_application.h"

// Element includes
// incompressible potential flow
#include "custom_elements/incompressible_potential_flow/incompressible_potential_flow_velocity_element.h"
#include "custom_elements/incompressible_potential_flow/incompressible_potential_flow_pressure_element.h"

// fractional step
#include "custom_elements/fractional_step/rans_fractional_step.h"

// k-epsilon turbulence model
#include "custom_elements/evm_k_epsilon/rans_evm_k_epsilon_epsilon.h"
#include "custom_elements/evm_k_epsilon/rans_evm_k_epsilon_k.h"
#include "custom_elements/evm_k_epsilon/rans_evm_k_epsilon_low_re_epsilon.h"
#include "custom_elements/evm_k_epsilon/rans_evm_k_epsilon_low_re_k.h"
#include "custom_elements/evm_k_epsilon/rans_evm_k_epsilon_k_afc_element.h"
#include "custom_elements/evm_k_epsilon/rans_evm_k_epsilon_epsilon_afc_element.h"
#include "custom_elements/evm_k_epsilon/rans_evm_k_epsilon_k_residual_based_fc_element.h"
#include "custom_elements/evm_k_epsilon/rans_evm_k_epsilon_epsilon_residual_based_fc_element.h"
#include "custom_elements/evm_k_epsilon/rans_evm_k_epsilon_k_cross_wind_stabilized_element.h"
#include "custom_elements/evm_k_epsilon/rans_evm_k_epsilon_epsilon_cross_wind_stabilized_element.h"

// Condition includes
#include "custom_conditions/evm_k_epsilon/rans_evm_k_epsilon_epsilon_velocity_based_wall_condition.h"
#include "custom_conditions/monolithic/rans_vms_monolithic_k_based_wall_condition.h"
#include "custom_conditions/incompressible_potential_flow/incompressible_potential_flow_velocity_condition.h"
#include "custom_conditions/incompressible_potential_flow/incompressible_potential_flow_pressure_condition.h"
#include "custom_conditions/fractional_step/fs_high_re_k_wall_condition.h"

// Adjoint element includes
#include "custom_elements/evm_k_epsilon/rans_evm_epsilon_adjoint.h"
#include "custom_elements/evm_k_epsilon/rans_evm_k_adjoint.h"
#include "custom_elements/evm_k_epsilon/rans_evm_k_epsilon_vms_adjoint.h"
#include "custom_elements/evm_k_epsilon/rans_evm_monolithic_k_epsilon_vms_adjoint.h"

// Adjoint condition includes
#include "custom_conditions/evm_k_epsilon/rans_evm_epsilon_adjoint_wall_condition.h"
#include "custom_conditions/evm_k_epsilon/rans_evm_vms_monolithic_adjoint_wall_condition.h"
#include "custom_conditions/evm_k_epsilon/rans_evm_monolithic_k_epsilon_vms_adjoint_wall_condition.h"

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

    /// incompressible potential flow elements
    const IncompressiblePotentialFlowVelocityElement<2, 3> mIncompressiblePotentialFlowVelocity2D;
    const IncompressiblePotentialFlowVelocityElement<3, 4> mIncompressiblePotentialFlowVelocity3D;
    const IncompressiblePotentialFlowPressureElement<2, 3> mIncompressiblePotentialFlowPressure2D;
    const IncompressiblePotentialFlowPressureElement<3, 4> mIncompressiblePotentialFlowPressure3D;

    /// incompressible potential flow conditions
    const IncompressiblePotentialFlowVelocityCondition<2, 2> mIncompressiblePotentialFlowVelocityCondition2D2N;
    const IncompressiblePotentialFlowVelocityCondition<3, 3> mIncompressiblePotentialFlowVelocityCondition3D3N;
    const IncompressiblePotentialFlowPressureCondition<2, 2> mIncompressiblePotentialFlowPressureCondition2D2N;
    const IncompressiblePotentialFlowPressureCondition<3, 3> mIncompressiblePotentialFlowPressureCondition3D3N;

    // Fractionalstep elements
    const RansFractionalStep<2> mRansFractionalStep2D;
    const RansFractionalStep<3> mRansFractionalStep3D;

    // FractionalStep wall conditions
    const FSHighReKWallCondition<2, 2> mFSHighReKWallCondition2D2N;
    const FSHighReKWallCondition<3, 3> mFSHighReKWallCondition3D3N;

    /// k-epsilon turbulence model elements
    /// Algebraic flux correction based elements
    const RansEvmKEpsilonKAFCElement<2, 3> mRansEvmKEpsilonKAFC2D;
    const RansEvmKEpsilonKAFCElement<3, 4> mRansEvmKEpsilonKAFC3D;

    const RansEvmKEpsilonEpsilonAFCElement<2, 3> mRansEvmKEpsilonEpsilonAFC2D;
    const RansEvmKEpsilonEpsilonAFCElement<3, 4> mRansEvmKEpsilonEpsilonAFC3D;

    /// Residual based flux corrected elements
    const RansEvmKEpsilonKResidualBasedFCElement<2, 3> mRansEvmKEpsilonKResidualBasedFC2D;
    const RansEvmKEpsilonKResidualBasedFCElement<3, 4> mRansEvmKEpsilonKResidualBasedFC3D;

    const RansEvmKEpsilonEpsilonResidualBasedFCElement<2, 3> mRansEvmKEpsilonEpsilonResidualBasedFC2D;
    const RansEvmKEpsilonEpsilonResidualBasedFCElement<3, 4> mRansEvmKEpsilonEpsilonResidualBasedFC3D;

    /// Cross wind stabilization based elements
    const RansEvmKEpsilonKCrossWindStabilizedElement<2, 3> mRansEvmKEpsilonKCrossWindStabilized2D;
    const RansEvmKEpsilonKCrossWindStabilizedElement<3, 4> mRansEvmKEpsilonKCrossWindStabilized3D;

    const RansEvmKEpsilonEpsilonCrossWindStabilizedElement<2, 3> mRansEvmKEpsilonEpsilonCrossWindStabilized2D;
    const RansEvmKEpsilonEpsilonCrossWindStabilizedElement<3, 4> mRansEvmKEpsilonEpsilonCrossWindStabilized3D;


    const RansEvmKEpsilonLowReKElement<2, 3> mRansEvmKEpsilonLowReK2D;
    const RansEvmKEpsilonLowReKElement<3, 4> mRansEvmKEpsilonLowReK3D;

    const RansEvmKEpsilonLowReEpsilonElement<2, 3> mRansEvmKEpsilonLowReEpsilon2D;
    const RansEvmKEpsilonLowReEpsilonElement<3, 4> mRansEvmKEpsilonLowReEpsilon3D;

    const RansEvmKEpsilonKElement<2, 3> mRansEvmKEpsilonK2D;
    const RansEvmKEpsilonKElement<3, 4> mRansEvmKEpsilonK3D;

    const RansEvmKEpsilonEpsilon<2, 3> mRansEvmKEpsilonEpsilon2D;
    const RansEvmKEpsilonEpsilon<3, 4> mRansEvmKEpsilonEpsilon3D;

    /// k-epsilon turbulence model conditions
    const RansEvmKEpsilonEpsilonVelocityBasedWallCondition<2> mRansEvmKEpsilonEpsilonVelocityBasedWallCondition2D2N;
    const RansEvmKEpsilonEpsilonVelocityBasedWallCondition<3> mRansEvmKEpsilonEpsilonVelocityBasedWallCondition3D3N;

    const RansVMSMonolithicKBasedWallCondition<2> mRansVMSMonolithicKBasedWallCondition2D2N;
    const RansVMSMonolithicKBasedWallCondition<3> mRansVMSMonolithicKBasedWallCondition3D3N;

    // k-epsilon adjoint elements
    const RansEvmEpsilonAdjoint<2, 3> mRansEvmEpsilonAdjoint2D3N;
    const RansEvmEpsilonAdjoint<3, 4> mRansEvmEpsilonAdjoint3D4N;

    const RansEvmKAdjoint<2, 3> mRansEvmKAdjoint2D3N;
    const RansEvmKAdjoint<3, 4> mRansEvmKAdjoint3D4N;

    const RansEvmKEpsilonVMSAdjoint<2> mRansEvmKEpsilonVMSAdjoint2D3N;
    const RansEvmKEpsilonVMSAdjoint<3> mRansEvmKEpsilonVMSAdjoint3D4N;

    const RansEvmMonolithicKEpsilonVMSAdjoint<2> mRansEvmMonolithicKEpsilonVMSAdjoint2D;
    const RansEvmMonolithicKEpsilonVMSAdjoint<3> mRansEvmMonolithicKEpsilonVMSAdjoint3D;

    // k-epsilon adjoint conditions
    const RansEvmEpsilonAdjointWallCondition<2> mRansEvmEpsilonAdjointWallCondition2D2N;
    const RansEvmEpsilonAdjointWallCondition<3> mRansEvmEpsilonAdjointWallCondition3D3N;

    const RansEvmVmsMonolithicAdjointWallCondition<2> mRansEvmVmsMonolithicAdjointWallCondition2D2N;
    const RansEvmVmsMonolithicAdjointWallCondition<3> mRansEvmVmsMonolithicAdjointWallCondition3D3N;

    const RansEvmMonolithicKEpsilonVMSAdjointWallCondition<2> mRansEvmMonolithicKEpsilonVMSAdjointWallCondition2D2N;
    const RansEvmMonolithicKEpsilonVMSAdjointWallCondition<3> mRansEvmMonolithicKEpsilonVMSAdjointWallCondition3D3N;
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
    KratosRANSApplication& operator=(KratosRANSApplication const& rOther);

    /// Copy constructor.
    KratosRANSApplication(KratosRANSApplication const& rOther);

    ///@}

}; // Class KratosRANSApplication

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.

#endif // KRATOS_RANS_APPLICATION_H_INCLUDED  defined
