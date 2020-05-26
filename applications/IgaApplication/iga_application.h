//  KRATOS  _____________
//         /  _/ ____/   |
//         / // / __/ /| |
//       _/ // /_/ / ___ |
//      /___/\____/_/  |_| Application

#if !defined(KRATOS_IGA_APPLICATION_H_INCLUDED)
#define  KRATOS_IGA_APPLICATION_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"

#include "custom_elements/shell_3p_element.h"
#include "custom_elements/iga_truss_element.h"
#include "custom_elements/shell_kl_discrete_element.h"

//conditions
#include "custom_conditions/load_condition.h"
#include "custom_conditions/penalty_coupling_condition.h"


namespace Kratos {

///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
*/
class KratosIgaApplication
    : public KratosApplication {
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of KratosIgaApplication
    KRATOS_CLASS_POINTER_DEFINITION(KratosIgaApplication);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    KratosIgaApplication();

    /// Destructor.
    ~KratosIgaApplication() override {}

    ///@}
    ///@name Operations
    ///@{

    void Register() override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information
    std::string Info() const override
    {
        return "KratosIgaApplication";
    }

    /// Print Information
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
        PrintData(rOStream);
    }

    /// Print object's data.
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

    const Shell3pElement mShell3pElement;
    const IgaTrussElement mIgaTrussElement;
    const ShellKLDiscreteElement mShellKLDiscreteElement;

    //Conditions
    const LoadCondition mLoadCondition;
    const PenaltyCouplingCondition mPenaltyCouplingCondition;

    ///@}
    ///@name Private methods
    ///@{

    /// Assignment operator.
    KratosIgaApplication& operator=(KratosIgaApplication const& rOther);

    /// Copy constructor.
    KratosIgaApplication(KratosIgaApplication const& rOther);

    ///@}

}; // class KratosIgaApplication

///@}

} // namespace Kratos

#endif // !defined(KRATOS_IGA_APPLICATION_H_INCLUDED)
