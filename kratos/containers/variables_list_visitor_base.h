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
//


#if !defined(KRATOS_VARIABLES_LIST_VISITOR_BASE_H_INCLUDED )
#define  KRATOS_VARIABLES_LIST_VISITOR_BASE_H_INCLUDED


// System includes
#include <string>
#include <iostream>


// Project includes
#include "includes/define.h"
#include "includes/ublas_interface.h"


namespace Kratos
{

template<class TDataType> class Variable; 

///@addtogroup ApplicationNameApplication
///@{

///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
*/
class VariablesListVisitorBase
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of VariablesListVisitorBase
    KRATOS_CLASS_POINTER_DEFINITION(VariablesListVisitorBase);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    VariablesListVisitorBase(){}

    /// Destructor.
    virtual ~VariablesListVisitorBase(){}

    ///@}
    ///@name Operations
    ///@{

    /// The visit class is called for each variable in the given variables list.
    virtual void Visit(Variable<bool> const& TheVariable){
        ThrowError("bool");
   }

    /// The visit class is called for each variable in the given variables list.
    virtual void Visit(Variable<int> const& TheVariable){
        ThrowError("integer");
    }

    /// The visit class is called for each variable in the given variables list.
    virtual void Visit(Variable<unsigned int> const& TheVariable){
        ThrowError("unsigned integer");
    }

    /// The visit class is called for each variable in the given variables list.
    virtual void Visit(Variable<double> const& TheVariable){
        ThrowError("double");
    }

    /// The visit class is called for each variable in the given variables list.
    virtual void Visit(Variable<array_1d<double,3>> const& TheVariable){
        ThrowError("array_1d3");
    }

    /// The visit class is called for each variable in the given variables list.
    virtual void Visit(Variable<Vector> const& TheVariable){
        ThrowError("Vector");
    }

    /// The visit class is called for each variable in the given variables list.
    virtual void Visit(Variable<Matrix> const& TheVariable){
        ThrowError("Matrix");
    }

    virtual void Visit(VariableData const& TheVariable){ // Fallback method
        KRATOS_ERROR << "Unsupported variable type. The supported varibles types to visits are" 
        << "boo, int, unsinged int, array_1d<double,3>, Vector, and Matrix." << std::endl;
    }

    template<typename TVariablesContainerType> 
    void VisitVariables(TVariablesContainerType const& Variables){
        for (auto p_variable : Variables) {
            p_variable->AcceptVisitor(*this);
        }
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "VariablesListVisitorBase" ;
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "VariablesListVisitorBase";}

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {}

    ///@}
    ///@name Friends
    ///@{


    ///@}

private:

    void ThrowError(std::string VariableType){
         KRATOS_ERROR << "Calling the base class visit method for the " << VariableType << 
         " variable. The main cause of this error is that the visitor class does not support " << VariableType  << std::endl;
    }

     ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    VariablesListVisitorBase& operator=(VariablesListVisitorBase const& rOther) = delete;

    /// Copy constructor.
    VariablesListVisitorBase(VariablesListVisitorBase const& rOther) = delete;

    ///@}

}; // Class VariablesListVisitorBase

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                const VariablesListVisitorBase& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_VARIABLES_LIST_VISITOR_BASE_H_INCLUDED  defined


