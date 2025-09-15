//    _____  _____ _                  __          __                                                 _ _           _   _
//   / ____|/ ____| |                 \ \        / /                               /\               | (_)         | | (_)
//  | |    | (___ | |__   __ _ _ __ _ _\ \  /\  / / __ __ _ _ __  _ __   ___ _ __ /  \   _ __  _ __ | |_  ___ __ _| |_ _  ___  _ __
//  | |     \___ \| '_ \ / _` | '__| '_ \ \/  \/ / '__/ _` | '_ \| '_ \ / _ \ '__/ /\ \ | '_ \| '_ \| | |/ __/ _` | __| |/ _ \| '_  |
//  | |____ ____) | | | | (_| | |  | |_) \  /\  /| | | (_| | |_) | |_) |  __/ | / ____ \| |_) | |_) | | | (_| (_| | |_| | (_) | | | |
//   \_____|_____/|_| |_|\__,_|_|  | .__/ \/  \/ |_|  \__,_| .__/| .__/ \___|_|/_/    \_\ .__/| .__/|_|_|\___\__,_|\__|_|\___/|_| |_|
//                                 | |                     | |   | |                    | |   | |
//                                 |_|                     |_|   |_|                    |_|   |_|
//
//
//  License: BSD License
//   license: CSharpWrapperApplication/license.txt
//
//  Main authors:    Hubert Balcerzak
//                   Riccardo Rossi
//                   Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_CSHARP_WRAPPER_APPLICATION_H_INCLUDED )
#define  KRATOS_CSHARP_WRAPPER_APPLICATION_H_INCLUDED


// System includes
#include <string>
#include <iostream>

// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"

// Application includes
#include "csharp_wrapper_application_variables.h"

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
 * @class KratosCSharpWrapperApplication
 * @ingroup KratosCSharpWrapperApplication
 * @brief This is a custom application for being able to couple C++ code of Kratos with CSharp
 * @details This applications adds a CSharp wrapper interface
 * @author Hubert Balcerzak
 * @author Riccardo Rossi
 * @author Vicente Mataix Ferrandiz
 */
class KRATOS_API(CSHARP_WRAPPER_APPLICATION) KratosCSharpWrapperApplication
    : public KratosApplication
{

public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of KratosCSharpWrapperApplication
    KRATOS_CLASS_POINTER_DEFINITION(KratosCSharpWrapperApplication);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    KratosCSharpWrapperApplication();

    /// Destructor.
    ~KratosCSharpWrapperApplication() override {}


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
        return "KratosCSharpWrapperApplication";
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
        rOStream << "KratosCSharpWrapperApplication\tNUMBER OF COMPONENTS: " << KratosComponents<VariableData>::GetComponents().size() << std::endl;
        rOStream << "Variables:" << std::endl;
        KratosComponents<VariableData>().PrintData(rOStream);
        rOStream << std::endl;
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

    /// Assignment operator.
    KratosCSharpWrapperApplication& operator=(KratosCSharpWrapperApplication const& rOther);

    /// Copy constructor.
    KratosCSharpWrapperApplication(KratosCSharpWrapperApplication const& rOther);


    ///@}

}; // Class KratosCSharpWrapperApplication

///@}


///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // KRATOS_CSHARP_WRAPPER_APPLICATION_H_INCLUDED  defined
