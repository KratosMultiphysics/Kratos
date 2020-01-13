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

//#include "custom_elements/iga_truss_element.h"
#include "custom_elements/shell_kl_discrete_element.h"
#include "custom_elements/iga_shell_3p_element.h"

#include "custom_conditions/iga_check_condition.h"

#include "custom_conditions/coupling_penalty_discrete_condition.h"
#include "custom_conditions/support_penalty_curve_discrete_condition.h"
#include "custom_conditions/support_penalty_point_discrete_condition.h"

#include "custom_conditions/load_point_discrete_condition.h"
#include "custom_conditions/load_surface_discrete_condition.h"

#include "custom_conditions/load_curve_discrete_condition.h"

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

    //const IgaTrussElement mIgaTrussElement;
    const ShellKLDiscreteElement mShellKLDiscreteElement;
    const IgaShell3pElement mIgaShell3pElement;

    const IgaCheckCondition mIgaCheckCondition;

    const CouplingPenaltyDiscreteCondition mCouplingPenaltyDiscreteCondition;

    const SupportPenaltyCurveDiscreteCondition mSupportPenaltyCurveDiscreteCondition;
    const SupportPenaltyPointDiscreteCondition mSupportPenaltyPointDiscreteCondition;

    const LoadSurfaceDiscreteCondition mLoadSurfaceDiscreteCondition;
    const LoadCurveDiscreteCondition mLoadCurveDiscreteCondition;


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
