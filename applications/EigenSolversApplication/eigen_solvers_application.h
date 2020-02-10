/*
//  KRATOS _______
//        / ____(_)___ ____  ____
//       / __/ / / __ `/ _ \/ __ \
//      / /___/ / /_/ /  __/ / / /
//     /_____/_/\__, /\___/_/ /_/ SolversApplication
//             /____/
//
//  Author: Thomas Oberbichler
*/

#if !defined(KRATOS_EIGENSOLVERS_APPLICATION_H_INCLUDED)
#define KRATOS_EIGENSOLVERS_APPLICATION_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define_python.h"
#include "includes/kratos_application.h"
#include "includes/variables.h"

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

/// registers the linear solvers to kratos
/** registers the linear solvers to kratos
*/
class EigenSolversApplicationRegisterLinearSolvers
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of EigenSolversApplicationRegisterLinearSolvers
    KRATOS_CLASS_POINTER_DEFINITION(EigenSolversApplicationRegisterLinearSolvers);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    EigenSolversApplicationRegisterLinearSolvers();

    /// Destructor.
    virtual ~EigenSolversApplicationRegisterLinearSolvers(){};


    ///@}

private:
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    EigenSolversApplicationRegisterLinearSolvers& operator=(EigenSolversApplicationRegisterLinearSolvers const& rOther) = delete;

    /// Copy constructor.
    EigenSolversApplicationRegisterLinearSolvers(EigenSolversApplicationRegisterLinearSolvers const& rOther) = delete;

    ///@}

}; // Class EigenSolversApplicationRegisterLinearSolvers

class KratosEigenSolversApplication : public KratosApplication
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of KratosEigenSolversApplication
    KRATOS_CLASS_POINTER_DEFINITION(KratosEigenSolversApplication);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    KratosEigenSolversApplication();

    /// Destructor.
    ~KratosEigenSolversApplication() override {}

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
        return "KratosEigenSolversApplication";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream &rOStream) const override
    {
        rOStream << Info();
        PrintData(rOStream);
    }

    /// Print object's data.
    void PrintData(std::ostream &rOStream) const override
    {
        KRATOS_WATCH("in KratosEigenSolversApplication application");
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
    KratosEigenSolversApplication &operator=(KratosEigenSolversApplication const &rOther);

    /// Copy constructor.
    KratosEigenSolversApplication(KratosEigenSolversApplication const &rOther);

    ///@}

}; // class KratosEigenSolversApplication

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos

#endif // defined(KRATOS_EIGENSOLVERS_APPLICATION_H_INCLUDED)
