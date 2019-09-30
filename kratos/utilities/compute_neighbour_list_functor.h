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

#if !defined(KRATOS_NEIGHBOUR_LIST_FUNCTOR_H_INCLUDED )
#define  KRATOS_NEIGHBOUR_LIST_FUNCTOR_H_INCLUDED


// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "utilities/global_pointer_utilities.h"

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
template< class TContainerType, class TVariableType >
class ComputeNeighbourListFunctor
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ComputeNeighbourListFunctor
    KRATOS_CLASS_POINTER_DEFINITION(ComputeNeighbourListFunctor);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ComputeNeighbourListFunctor(
            const TContainerType& rContainer,
            const TVariableType& rVar
            
            )
            : mrContainer(rContainer), mrVariable(rVar)
    {}

    /// Destructor.
    virtual ~ComputeNeighbourListFunctor(){}

    ///@}
    ///@name Operators
    ///@{
    GlobalPointersVector< typename TVariableType::Type::data_type > operator()(const DataCommunicator& rComm) const
    {
        GlobalPointersVector<typename TVariableType::Type::data_type> gp_list;
        for(auto& item : mrContainer)
        {
            auto& neighbours = item.GetValue(mrVariable);
            for(auto& gp : neighbours.GetContainer())
                gp_list.push_back(gp);
        }
        gp_list.Unique();
        return gp_list;
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
    buffer << "ComputeNeighbourListFunctor" ;
    return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "ComputeNeighbourListFunctor";}

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
    const TContainerType& mrContainer;
    const TVariableType& mrVariable; //note that A REFERENCE is stored


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

}; // Class ComputeNeighbourListFunctor

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template<class TContainerType, class TVariableType>
inline std::istream& operator >> (std::istream& rIStream,
                ComputeNeighbourListFunctor<TContainerType, TVariableType>& rThis)
{return rIStream;}

/// output stream function
template<class TContainerType, class TVariableType>
inline std::ostream& operator << (std::ostream& rOStream,
                const ComputeNeighbourListFunctor<TContainerType, TVariableType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_NEIGHBOUR_LIST_FUNCTOR_H_INCLUDED  defined


