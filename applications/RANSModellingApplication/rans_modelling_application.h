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

#if !defined(KRATOS_RANS_MODELLING_APPLICATION_H_INCLUDED)
#define KRATOS_RANS_MODELLING_APPLICATION_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/kratos_application.h"

// Primal element includes
#include "custom_elements/evm_k_epsilon/evm_epsilon_element.h"
#include "custom_elements/evm_k_epsilon/evm_k_element.h"
#include "custom_elements/evm_k_epsilon/evm_low_re_epsilon_element.h"
#include "custom_elements/evm_k_epsilon/evm_low_re_k_element.h"

// Primal condition includes
#include "custom_conditions/evm_k_epsilon/evm_epsilon_wall_condition.h"
#include "custom_conditions/evm_k_epsilon/evm_vms_monolithic_wall_condition.h"

// Adjoint element includes
#include "custom_elements/evm_k_epsilon/evm_epsilon_adjoint_element.h"

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
class KratosRANSModellingApplication : public KratosApplication
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of KratosRANSModellingApplication
    KRATOS_CLASS_POINTER_DEFINITION(KratosRANSModellingApplication);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    KratosRANSModellingApplication();

    /// Destructor.
    ~KratosRANSModellingApplication() override
    {
    }

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
        return "KratosRANSModellingApplication";
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
    const EvmLowReKElement<2, 3> mRANSEVMLowReK2D;
    const EvmLowReKElement<3, 4> mRANSEVMLowReK3D;

    const EvmLowReEpsilonElement<2, 3> mRANSEVMLowReEpsilon2D;
    const EvmLowReEpsilonElement<3, 4> mRANSEVMLowReEpsilon3D;

    const EvmKElement<2, 3> mRANSEVMK2D;
    const EvmKElement<3, 4> mRANSEVMK3D;

    const EvmEpsilonElement<2, 3> mRANSEVMEpsilon2D;
    const EvmEpsilonElement<3, 4> mRANSEVMEpsilon3D;

    /// k-epsilon turbulence model conditions
    const EvmEpsilonWallCondition<2> mRansEVMEpsilonWallCondition2D2N;
    const EvmEpsilonWallCondition<3> mRansEVMEpsilonWallCondition3D3N;

    const EVMVMSMonolithicWallCondition<2> mRansEVMVMSMonolithicWallCondition2D2N;
    const EVMVMSMonolithicWallCondition<3> mRansEVMVMSMonolithicWallCondition3D3N;

    // k-epsilon adjoint elements
    const EvmEpsilonAdjointElement<2, 3> mRANSEVMEpsilonAdjoint2D;
    const EvmEpsilonAdjointElement<3, 4> mRANSEVMEpsilonAdjoint3D;


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
    KratosRANSModellingApplication& operator=(KratosRANSModellingApplication const& rOther);

    /// Copy constructor.
    KratosRANSModellingApplication(KratosRANSModellingApplication const& rOther);

    ///@}

}; // Class KratosRANSModellingApplication

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.

#endif // KRATOS_RANS_MODELLING_APPLICATION_H_INCLUDED  defined
