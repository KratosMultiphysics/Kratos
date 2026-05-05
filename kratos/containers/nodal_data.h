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

#if !defined(KRATOS_NODAL_DATA_H_INCLUDED )
#define  KRATOS_NODAL_DATA_H_INCLUDED

// System includes
#include <iostream>

// External includes

// Project includes
#include "containers/variables_list_data_value_container.h"

namespace Kratos
{
///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

/// Stores all data and dofs which are stored in each elements
/** This class is the container for nodal data storing:
 *  Id : The Id of the node
*/
class KRATOS_API(KRATOS_CORE) NodalData final
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of NodalData
    KRATOS_CLASS_POINTER_DEFINITION(NodalData);

    using IndexType = std::size_t;

    using SizeType = std::size_t;

    using SolutionStepsNodalDataContainerType = VariablesListDataValueContainer;

    using BlockType=VariablesListDataValueContainer::BlockType;

    ///@}
    ///@name Life Cycle
    ///@{

    NodalData(IndexType TheId);

    NodalData(IndexType TheId, VariablesList::Pointer pVariablesList, SizeType NewQueueSize = 1);

    NodalData(IndexType TheId, VariablesList::Pointer pVariablesList, BlockType const * ThisData, SizeType NewQueueSize = 1);

    NodalData(NodalData&&) noexcept = default;

    ///@}
    ///@name Operators
    ///@{

    NodalData& operator=(NodalData const& rOther) = default;

    NodalData& operator=(NodalData&&) noexcept = default;

    ///@}

    ///@}
    ///@name Access
    ///@{

    /// Returns the Id of the Node. Same as GetId to ensure backward compatibility
    IndexType Id() const noexcept
    {
        return mId;
    }

    /// Returns the Id of the Node.
    IndexType GetId() const noexcept
    {
        return mId;
    }

    /// Sets the Id of the Node.
    void SetId(IndexType NewId) noexcept
    {
        mId = NewId;
    }


    void SetSolutionStepData(VariablesListDataValueContainer const& TheData)
    {
        mSolutionStepsNodalData = TheData;
    }

    VariablesListDataValueContainer& GetSolutionStepData() noexcept
    {
        return mSolutionStepsNodalData;
    }

    const VariablesListDataValueContainer& GetSolutionStepData() const noexcept
    {
        return mSolutionStepsNodalData;
    }

    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const;

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const;


    ///@}
    ///@name Friends
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

        IndexType mId;

        SolutionStepsNodalDataContainerType mSolutionStepsNodalData;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    // Only for serializer
    NodalData() : mId(0), mSolutionStepsNodalData(){}

    void save(Serializer& rSerializer) const;

    void load(Serializer& rSerializer);

    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Copy constructor.
    NodalData(NodalData const& rOther) = default;

    ///@}

}; // Class NodalData

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                NodalData& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                const NodalData& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_NODAL_DATA_H_INCLUDED  defined
