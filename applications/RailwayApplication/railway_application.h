
//

#pragma once

// System includes
#include <iostream>
#include <string>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"
#include "includes/smart_pointers.h"


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
class KRATOS_API(RAILWAY_APPLICATION) KratosRailwayApplication : public KratosApplication
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of KratosGeoMechanicsApplication
    KRATOS_CLASS_POINTER_DEFINITION(KratosRailwayApplication);

    ///@}
    ///@name Life Cycle
    ///@{

    KratosRailwayApplication();
    ~KratosRailwayApplication() override                                 = default;
    KratosRailwayApplication(const KratosRailwayApplication&)            = delete;
    KratosRailwayApplication& operator=(const KratosRailwayApplication&) = delete;
    KratosRailwayApplication(KratosRailwayApplication&&)                 = delete;
    KratosRailwayApplication& operator=(KratosRailwayApplication&&)      = delete;

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
    std::string Info() const override { return "KratosRailwayApplication"; }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
        PrintData(rOStream);
    }

    ///// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        KRATOS_WATCH("in KratosRailwayApplication")
        KRATOS_WATCH(KratosComponents<VariableData>::GetComponents().size())


    }

    ///@}
    ///@name Friends
    ///@{

    ///@}

private:
    ///@name Static Member Variables
    ///@{

    // static const ApplicationCondition  msApplicationCondition;

    ///@}
    ///@name Member Variables
    ///@{

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


    ///@}

}; // Class KratosRailwayApplication

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos
