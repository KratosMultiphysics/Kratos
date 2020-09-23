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

    virtual void Visit(Variable<bool> const& TheVariable) = 0;
    virtual void Visit(Variable<int> const& TheVariable) = 0;
    virtual void Visit(Variable<unsigned int> const& TheVariable) = 0;
    virtual void Visit(Variable<double> const& TheVariable) = 0;
    virtual void Visit(Variable<array_1d<double,3>> const& TheVariable) = 0;
    virtual void Visit(Variable<Vector> const& TheVariable) = 0;
    virtual void Visit(Variable<Matrix> const& TheVariable) = 0;
    
    virtual void Visit(VariableData const& TheVariable) = 0; // Fallback method

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


