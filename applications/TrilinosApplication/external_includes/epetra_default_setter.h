//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//  Collaboratos:    Philipp Bucher
//

#if !defined (KRATOS_EPETRA_DEFAULT_SETTER_H_INCLUDED)
#define KRATOS_EPETRA_DEFAULT_SETTER_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "ml_MultiLevelPreconditioner.h"
#include "Teuchos_ParameterList.hpp"

namespace Kratos
{
///@name Kratos Classes
///@{

/// Helper-class to set defaults in the Teuchos::ParameterList
/// using ML_Epetra::SetDefaults
class EpetraDefaultSetter
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of EpetraDefaultSetter
    KRATOS_CLASS_POINTER_DEFINITION(EpetraDefaultSetter);

    ///@}
    ///@name Life Cycle
    ///@{

    EpetraDefaultSetter()
    {
        KRATOS_WARNING("DEPRECATION") << "\"EpetraDefaultSetter\" is deprecated, please use \"TrilinosSolverUtilities::SetEpetraDefaults\"" << std::endl;
    }

    /// Copy constructor.
    EpetraDefaultSetter(EpetraDefaultSetter const& rOther) = delete;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    EpetraDefaultSetter& operator=(EpetraDefaultSetter const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    void SetDefaults(Teuchos::ParameterList& rParameterlist, const std::string& rSettingsName)
    {
        ML_Epetra::SetDefaults(rSettingsName.c_str(), rParameterlist);
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "EpetraDefaultSetter" ;
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "EpetraDefaultSetter";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const {}

    ///@}

}; // Class EpetraDefaultSetter

///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const EpetraDefaultSetter& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}  // namespace Kratos.

#endif // KRATOS_EPETRA_DEFAULT_SETTER_H_INCLUDED defined
