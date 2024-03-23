//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         DigitalTwinApplication/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//                   Ihar Antonau
//                   Fabian Meister
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/kratos_application.h"

namespace Kratos {

///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
 */
class KRATOS_API(DIGITAL_TWIN_APPLICATION) KratosDigitalTwinApplication : public KratosApplication {
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of KratosDigitalTwinApplication
    KRATOS_CLASS_POINTER_DEFINITION(KratosDigitalTwinApplication);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    KratosDigitalTwinApplication();

    /// Copy constructor.
    KratosDigitalTwinApplication(KratosDigitalTwinApplication const& rOther) = delete;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    KratosDigitalTwinApplication& operator=(KratosDigitalTwinApplication const& rOther) = delete;

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
        return "KratosDigitalTwinApplication";
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

}; // Class KratosDigitalTwinApplication

///@}

} // namespace Kratos.
