//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

#if !defined(KRATOS_GET_VALUE_FUNCTOR_H_INCLUDED )
#define  KRATOS_GET_VALUE_FUNCTOR_H_INCLUDED


// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"

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
template< class TVariableType >
class GetValueFunctor
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of GetValueFunctor
    KRATOS_CLASS_POINTER_DEFINITION(GetValueFunctor);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    GetValueFunctor(
            const TVariableType& rVariable 
           )
            : mrVariable(rVariable)
    {}


    /// Destructor.
    virtual ~GetValueFunctor(){}

    ///@}
    ///@name Operators
    ///@{
    typename TVariableType::Type operator()(GlobalPointer< Node >& gp) const
    {
        return gp->GetValue(mrVariable);
    }


    ///@}
    ///@name Operations
    ///@{


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
    virtual std::string Info() const
    {
    std::stringstream buffer;
    buffer << "GetValueFunctor" ;
    return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "GetValueFunctor";}

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {}

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
    const TVariableType& mrVariable;


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

}; // Class GetValueFunctor

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template<class TContainerType>
inline std::istream& operator >> (std::istream& rIStream,
                GetValueFunctor<TContainerType>& rThis)
{return rIStream;}

/// output stream function
template<class TContainerType>
inline std::ostream& operator << (std::ostream& rOStream,
                const GetValueFunctor<TContainerType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_GET_VALUE_FUNCTOR_H_INCLUDED  defined


