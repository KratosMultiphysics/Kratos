//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    janosch
//

#if !defined(KRATOS_EXTERNAL_SOLVERS_APPLICATION_H_INCLUDED)
#define KRATOS_EXTERNAL_SOLVERS_APPLICATION_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"
#include "includes/variables.h"

namespace Kratos {

///@name Kratos Globals
///@{

// Variables definition

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

/// registers the linear solvers to kratos
/** registers the linear solvers to kratos
*/
class ExternalSolversApplicationRegisterLinearSolvers
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ExternalSolversApplicationRegisterLinearSolvers
    KRATOS_CLASS_POINTER_DEFINITION(ExternalSolversApplicationRegisterLinearSolvers);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ExternalSolversApplicationRegisterLinearSolvers();

    /// Destructor.
    virtual ~ExternalSolversApplicationRegisterLinearSolvers(){};


    ///@}

private:
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    ExternalSolversApplicationRegisterLinearSolvers& operator=(ExternalSolversApplicationRegisterLinearSolvers const& rOther) = delete;

    /// Copy constructor.
    ExternalSolversApplicationRegisterLinearSolvers(ExternalSolversApplicationRegisterLinearSolvers const& rOther) = delete;

    ///@}

}; // Class ExternalSolversApplicationRegisterLinearSolvers

/// Short class definition.
/** Detail class definition.
*/
class KratosExternalSolversApplication : public KratosApplication {
   public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of KratosExternalSolversApplication
    KRATOS_CLASS_POINTER_DEFINITION(KratosExternalSolversApplication);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    KratosExternalSolversApplication()
        : KratosApplication("ExternalSolversApplication"){};

    /// Destructor.
    ~KratosExternalSolversApplication() override {}

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
    std::string Info() const override {
        return "KratosExternalSolversApplication";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override {
        rOStream << Info();
        PrintData(rOStream);
    }

    ///// Print object's data.
    void PrintData(std::ostream& rOStream) const override {
        KRATOS_WATCH("in KratosExternalSolvers application");
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

    //       static const ApplicationCondition  msApplicationCondition;

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
    KratosExternalSolversApplication& operator=(
        KratosExternalSolversApplication const& rOther);

    /// Copy constructor.
    KratosExternalSolversApplication(
        KratosExternalSolversApplication const& rOther);

    ///@}

};  // Class KratosExternalSolversApplication

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

}  // namespace Kratos.

#endif  // KRATOS_EXTERNAL_SOLVERS_APPLICATION_H_INCLUDED  defined
