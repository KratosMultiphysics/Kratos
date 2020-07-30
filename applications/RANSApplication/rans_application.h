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
#include "custom_elements/evm_k_epsilon/rans_evm_k_epsilon_epsilon.h"
#include "custom_elements/evm_k_epsilon/rans_evm_k_epsilon_k.h"
#include "custom_elements/evm_k_epsilon/rans_evm_k_epsilon_low_re_epsilon.h"
#include "custom_elements/evm_k_epsilon/rans_evm_k_epsilon_low_re_k.h"

// Condition includes
#include "custom_conditions/evm_k_epsilon/rans_evm_k_epsilon_epsilon_wall.h"
#include "custom_conditions/evm_k_epsilon/rans_evm_k_epsilon_vms_monolithic_wall.h"

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

    /// k-epsilon turbulence model elements
    const RansEvmKEpsilonLowReKElement<2, 3> mRansEvmKEpsilonLowReK2D;
    const RansEvmKEpsilonLowReKElement<3, 4> mRansEvmKEpsilonLowReK3D;

    const RansEvmKEpsilonLowReEpsilonElement<2, 3> mRansEvmKEpsilonLowReEpsilon2D;
    const RansEvmKEpsilonLowReEpsilonElement<3, 4> mRansEvmKEpsilonLowReEpsilon3D;

    const RansEvmKEpsilonKElement<2, 3> mRansEvmKEpsilonK2D;
    const RansEvmKEpsilonKElement<3, 4> mRansEvmKEpsilonK3D;

    const RansEvmKEpsilonEpsilon<2, 3> mRansEvmKEpsilonEpsilon2D;
    const RansEvmKEpsilonEpsilon<3, 4> mRansEvmKEpsilonEpsilon3D;

    /// k-epsilon turbulence model conditions
    const RansEvmKEpsilonEpsilonWall<2> mRansEvmKEpsilonEpsilonWall2D2N;
    const RansEvmKEpsilonEpsilonWall<3> mRansEvmKEpsilonEpsilonWall3D3N;

    const RansEvmKEpsilonVmsMonolithicWall<2> mRansEvmKEpsilonVmsMonolithicWall2D2N;
    const RansEvmKEpsilonVmsMonolithicWall<3> mRansEvmKEpsilonVmsMonolithicWall3D3N;

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
