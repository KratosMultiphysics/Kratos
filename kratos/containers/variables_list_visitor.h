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


#if !defined(KRATOS_VARIABLES_LIST_VISITOR_H_INCLUDED )
#define  KRATOS_VARIABLES_LIST_VISITOR_H_INCLUDED


// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/ublas_interface.h"
#include "containers/variable.h"


namespace Kratos
{
///@addtogroup ApplicationNameApplication
///@{

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
template<typename TFunctionType>
class VariablesListVisitor : public VariablesListVisitorBase
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of VariablesListVisitor
    KRATOS_CLASS_POINTER_DEFINITION(VariablesListVisitor);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    VariablesListVisitor(TFunctionType&& TheFunction)
        : mFunction(std::move(TheFunction)) {}

    /// Destructor.
    virtual ~VariablesListVisitor(){}

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    virtual void Visit(Variable<bool> const& TheVariable) override {mFunction(TheVariable);}
    virtual void Visit(Variable<int> const& TheVariable) override {mFunction(TheVariable);}
    virtual void Visit(Variable<unsigned int> const& TheVariable) override {mFunction(TheVariable);}
    virtual void Visit(Variable<double> const& TheVariable) override {mFunction(TheVariable);}
    virtual void Visit(Variable<array_1d<double,3>> const& TheVariable) override {mFunction(TheVariable);}
    virtual void Visit(Variable<Vector> const& TheVariable) override {mFunction(TheVariable);}
    virtual void Visit(Variable<Matrix> const& TheVariable) override {mFunction(TheVariable);}
    
    virtual void Visit(VariableData const& TheVariable) override { // Fallback method
        KRATOS_ERROR << "The variable " << TheVariable << " type is not supported by the variables list visitor." << std::endl;
    } 

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
    virtual std::string Info() const override
    {
std::stringstream buffer;
    buffer << "VariablesListVisitor" ;
    return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override {rOStream << "VariablesListVisitor";}

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override {}

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

        TFunctionType mFunction;


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
    VariablesListVisitor& operator=(VariablesListVisitor const& rOther) = delete;

    /// Copy constructor.
    VariablesListVisitor(VariablesListVisitor const& rOther) = delete;

    ///@}

}; // Class VariablesListVisitor

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// output stream function
template<typename TFunctionType>
inline std::ostream& operator << (std::ostream& rOStream,
                const VariablesListVisitor<TFunctionType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_VARIABLES_LIST_VISITOR_H_INCLUDED  defined


