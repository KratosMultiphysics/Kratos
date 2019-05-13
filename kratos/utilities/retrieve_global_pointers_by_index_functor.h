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

#if !defined(KRATOS_RETRIEVE_GLOBAL_POINTERS_BY_INDEX_FUNCTOR_H_INCLUDED )
#define  KRATOS_RETRIEVE_GLOBAL_POINTERS_BY_INDEX_FUNCTOR_H_INCLUDED


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
template< class TContainerType >
class RetrieveGlobalPointersByIndex
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of RetrieveGlobalPointersByIndex
    KRATOS_CLASS_POINTER_DEFINITION(RetrieveGlobalPointersByIndex);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    RetrieveGlobalPointersByIndex(
            const TContainerType& rContainer,
            const std::vector<int>& rIndices
            )
            : mrContainer(rContainer), mrIndices(rIndices)
    {}

    /// Destructor.
    virtual ~RetrieveGlobalPointersByIndex(){}

    ///@}
    ///@name Operators
    ///@{
    GlobalPointersVector<typename TContainerType::value_type> operator()(const DataCommunicator& rComm) const
    {
        return GlobalPointerUtilities::RetrieveGlobalIndexedPointers(mrContainer, mrIndices, rComm );
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
    buffer << "RetrieveGlobalPointersByIndex" ;
    return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "RetrieveGlobalPointersByIndex";}

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
    const std::vector<int>& mrIndices; //note that A REFERENCE is stored


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

}; // Class RetrieveGlobalPointersByIndex

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template<class TContainerType>
inline std::istream& operator >> (std::istream& rIStream,
                RetrieveGlobalPointersByIndex<TContainerType>& rThis)
{return rIStream;}

/// output stream function
template<class TContainerType>
inline std::ostream& operator << (std::ostream& rOStream,
                const RetrieveGlobalPointersByIndex<TContainerType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_RETRIEVE_GLOBAL_POINTERS_BY_INDEX_FUNCTOR_H_INCLUDED  defined


