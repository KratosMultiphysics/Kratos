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
//

#if !defined(KRATOS_VARIABLES_LIST_DATA_VALUE_CONTAINER_H_INCLUDED )
#define KRATOS_VARIABLES_LIST_DATA_VALUE_CONTAINER_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "containers/variable.h"
#include "containers/variables_list.h"

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

/**
* @class VariablesListDataValueContainer
* @ingroup KratosCore
* @brief A  shared  variable  list gives the position of each variable in the containers sharing it.
* @details The mechanism is very simple. There is an array which stores the local offset for each variable in the container and assigns  the  value−1  for  the  rest  of  the variables
* For more details see P. Dadvand, R. Rossi, E. Oñate: An Object-oriented Environment for Developing Finite Element Codes for Multi-disciplinary Applications. Computational Methods in Engineering. 2010
* @author Pooyan Dadvand
* @author Riccardo Rossi
*/
class KRATOS_API(KRATOS_CORE) VariablesListDataValueContainer final
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of VariablesListDataValueContainer
    KRATOS_CLASS_POINTER_DEFINITION(VariablesListDataValueContainer);

    typedef VariablesList::BlockType BlockType;

    /// Type of the container used for variables
    typedef BlockType* ContainerType;

    typedef std::size_t IndexType;

    typedef std::size_t SizeType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    explicit VariablesListDataValueContainer(const SizeType NewQueueSize = 1);

    /// Copy constructor.
    VariablesListDataValueContainer(VariablesListDataValueContainer const& rOther);

    /// Variables list constructor.
    VariablesListDataValueContainer(VariablesList::Pointer pVariablesList, const SizeType NewQueueSize = 1);

    /// Variables list and data constructor
    VariablesListDataValueContainer(VariablesList::Pointer pVariablesList, BlockType const * ThisData, const SizeType NewQueueSize = 1);

    /// Destructor.
    ~VariablesListDataValueContainer();

    ///@}
    ///@name Operators
    ///@{

    template<class TDataType>
    TDataType& operator()(const Variable<TDataType>& rThisVariable)
    {
        return GetValue(rThisVariable);
    }

    template<class TDataType>
    TDataType& operator()(const Variable<TDataType>& rThisVariable, SizeType QueueIndex)
    {
        return GetValue(rThisVariable, QueueIndex);
    }

    template<class TDataType>
    const TDataType& operator()(const Variable<TDataType>& rThisVariable) const
    {
        return GetValue(rThisVariable);
    }

    template<class TDataType>
    const TDataType& operator()(const Variable<TDataType>& rThisVariable, SizeType QueueIndex) const
    {
        return GetValue(rThisVariable, QueueIndex);
    }

    template<class TDataType> 
    TDataType& operator[](const Variable<TDataType>& rThisVariable)
    {
        return GetValue(rThisVariable);
    }

    template<class TDataType> 
    const TDataType& operator[](const Variable<TDataType>& rThisVariable) const
    {
        return GetValue(rThisVariable);
    }

    /// Assignment operator.
    VariablesListDataValueContainer& operator=(const VariablesListDataValueContainer& rOther);

    ///@}
    ///@name Operations
    ///@{

    template<class TDataType>
    TDataType& GetValue(const Variable<TDataType>& rThisVariable)
    {
        KRATOS_DEBUG_ERROR_IF(!mpVariablesList) << "This container don't have a variables list assigned. A possible reason is creating a node without a model part." << std::endl;
        KRATOS_ERROR_IF_NOT(mpVariablesList->Has(rThisVariable)) << "This container only can store the variables specified in its variables list. The variables list doesn't have this variable:" << rThisVariable << std::endl;
        return *(reinterpret_cast<TDataType*>(Position(rThisVariable)) + rThisVariable.GetComponentIndex());
    }

    template<class TDataType>
    TDataType& GetValue(const Variable<TDataType>& rThisVariable, SizeType QueueIndex)
    {
        KRATOS_DEBUG_ERROR_IF(!mpVariablesList) << "This container don't have a variables list assigned. A possible reason is creating a node without a model part." << std::endl;
        KRATOS_ERROR_IF_NOT(mpVariablesList->Has(rThisVariable)) << "This container only can store the variables specified in its variables list. The variables list doesn't have this variable:" << rThisVariable << std::endl;
        return *(reinterpret_cast<TDataType*>(Position(rThisVariable, QueueIndex)) + rThisVariable.GetComponentIndex());
    }

    template<class TDataType>
    const TDataType& GetValue(const Variable<TDataType>& rThisVariable) const
    {
        KRATOS_DEBUG_ERROR_IF(!mpVariablesList) << "This container don't have a variables list assigned. A possible reason is creating a node without a model part." << std::endl;
        KRATOS_ERROR_IF_NOT(mpVariablesList->Has(rThisVariable)) << "This container only can store the variables specified in its variables list. The variables list doesn't have this variable:" << rThisVariable << std::endl;
        return *(reinterpret_cast<const TDataType*>(Position(rThisVariable)) + rThisVariable.GetComponentIndex());
    }

    template<class TDataType>
    const TDataType& GetValue(const Variable<TDataType>& rThisVariable, SizeType QueueIndex) const
    {
        KRATOS_DEBUG_ERROR_IF(!mpVariablesList) << "This container don't have a variables list assigned. A possible reason is creating a node without a model part." << std::endl;
        KRATOS_ERROR_IF_NOT(mpVariablesList->Has(rThisVariable)) << "This container only can store the variables specified in its variables list. The variables list doesn't have this variable:" << rThisVariable << std::endl;
        return *(reinterpret_cast<const TDataType*>(Position(rThisVariable, QueueIndex)) + rThisVariable.GetComponentIndex());
    }

    template<class TDataType>
    TDataType& FastGetValue(const Variable<TDataType>& rThisVariable)
    {
        return *(reinterpret_cast<TDataType*>(Position(rThisVariable)) + rThisVariable.GetComponentIndex());
    }

    template<class TDataType>
    TDataType* pFastGetValue(const Variable<TDataType>& rThisVariable)
    {
        return (reinterpret_cast<TDataType*>(Position(rThisVariable)) + rThisVariable.GetComponentIndex());
    }

    template<class TDataType>
    TDataType& FastGetValue(const Variable<TDataType>& rThisVariable, SizeType QueueIndex)
    {
        return *(reinterpret_cast<TDataType*>(Position(rThisVariable, QueueIndex)) + rThisVariable.GetComponentIndex());
    }

    template<class TDataType>
    TDataType& FastGetValue(const Variable<TDataType>& rThisVariable, SizeType QueueIndex, SizeType ThisPosition)
    {
        KRATOS_DEBUG_ERROR_IF(!mpVariablesList) << "This container don't have a variables list assigned. A possible reason is creating a node without a model part." << std::endl;
        KRATOS_DEBUG_ERROR_IF_NOT(mpVariablesList->Has(rThisVariable)) << "This container only can store the variables specified in its variables list. The variables list doesn't have this variable:" << rThisVariable << std::endl;
        KRATOS_DEBUG_ERROR_IF((QueueIndex + 1) > mQueueSize) << "Trying to access data from step " << QueueIndex << " but only " << mQueueSize << " steps are stored." << std::endl;
        return *(reinterpret_cast<TDataType*>(Position(QueueIndex) + ThisPosition) + rThisVariable.GetComponentIndex());
    }

    template<class TDataType>
    TDataType& FastGetCurrentValue(const Variable<TDataType>& rThisVariable, SizeType ThisPosition)
    {
        KRATOS_DEBUG_ERROR_IF(!mpVariablesList) << "This container don't have a variables list assigned. A possible reason is creating a node without a model part." << std::endl;
        KRATOS_DEBUG_ERROR_IF_NOT(mpVariablesList->Has(rThisVariable)) << "This container only can store the variables specified in its variables list. The variables list doesn't have this variable:" << rThisVariable << std::endl;
        return *(reinterpret_cast<TDataType*>(mpCurrentPosition + ThisPosition) + rThisVariable.GetComponentIndex());
    }

    template<class TDataType>
    const TDataType& FastGetValue(const Variable<TDataType>& rThisVariable) const
    {
        return *(reinterpret_cast<const TDataType*>(Position(rThisVariable)) + rThisVariable.GetComponentIndex());
    }

    template<class TDataType>
    const TDataType* pFastGetValue(const Variable<TDataType>& rThisVariable) const
    {
        return (reinterpret_cast<const TDataType*>(Position(rThisVariable)) + rThisVariable.GetComponentIndex());
    }

    template<class TDataType>
    const TDataType& FastGetValue(const Variable<TDataType>& rThisVariable, SizeType QueueIndex) const
    {
        return *(reinterpret_cast<const TDataType*>(Position(rThisVariable, QueueIndex)) + rThisVariable.GetComponentIndex());
    }

    template<class TDataType>
    const TDataType& FastGetValue(const Variable<TDataType>& rThisVariable, SizeType QueueIndex, SizeType ThisPosition) const
    {
        return *(reinterpret_cast<const TDataType*>(Position(QueueIndex) + ThisPosition) + rThisVariable.GetComponentIndex());
    }

    template<class TDataType>
    const TDataType& FastGetCurrentValue(const Variable<TDataType>& rThisVariable, SizeType ThisPosition) const
    {
        return *(reinterpret_cast<TDataType*>(mpCurrentPosition + ThisPosition) + rThisVariable.GetComponentIndex());
    }

    SizeType Size() const;

    SizeType QueueSize() const;

    SizeType TotalSize() const;

    template<class TDataType> 
    void SetValue(const Variable<TDataType>& rThisVariable, TDataType const& rValue)
    {
        GetValue(rThisVariable) = rValue;
    }

    template<class TDataType> 
    void SetValue(const Variable<TDataType>& rThisVariable, TDataType const& rValue, SizeType QueueIndex)
    {
        GetValue(rThisVariable, QueueIndex) = rValue;
    }

    void Clear();

    ///@}
    ///@name Access
    ///@{

    VariablesList::Pointer pGetVariablesList();

    const VariablesList::Pointer pGetVariablesList() const;

    VariablesList& GetVariablesList();

    const VariablesList& GetVariablesList() const;

    void SetVariablesList(VariablesList::Pointer pVariablesList);

    void SetVariablesList(VariablesList::Pointer pVariablesList, const SizeType ThisQueueSize);

    void Resize(const SizeType NewSize);

    BlockType* Data();

    const BlockType* Data() const;

    BlockType* Data(const SizeType QueueIndex);

    BlockType* Data(VariableData const & rThisVariable);

    SizeType DataSize();

    SizeType TotalDataSize();

    void AssignData(BlockType* Source, SizeType QueueIndex);

    void CloneFront();

    void PushFront();

    void AssignZero();

    void AssignZero(const SizeType QueueIndex);

    ///@}
    ///@name Inquiry
    ///@{

    template<class TDataType> 
    bool Has(const Variable<TDataType>& rThisVariable) const
    {
        if(!mpVariablesList)
            return false;

        return mpVariablesList->Has(rThisVariable);
    }

    bool IsEmpty();

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const
    {
        return std::string("variables list data value container");
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "variables list data value container";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const;

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

    SizeType mQueueSize;

    BlockType* mpCurrentPosition;

    ContainerType mpData;

    VariablesList::Pointer mpVariablesList;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    inline void Allocate();

    inline void Reallocate();

    void DestructElements(const SizeType ThisIndex);

    void DestructAllElements();

    void AssignData(BlockType* Source, BlockType* Destination);

    inline SizeType LocalOffset(VariableData const & rThisVariable) const;

    inline BlockType* Position(VariableData const & rThisVariable) const
    {
        KRATOS_DEBUG_ERROR_IF(!mpVariablesList) << "This container don't have a variables list assigned. A possible reason is creating a node without a model part." << std::endl;
        KRATOS_DEBUG_ERROR_IF_NOT(mpVariablesList->Has(rThisVariable)) << "This container only can store the variables specified in its variables list. The variables list doesn't have this variable:" << rThisVariable << std::endl;
        return mpCurrentPosition + mpVariablesList->Index(rThisVariable.SourceKey());
    }

    inline BlockType* Position(VariableData const & rThisVariable, const SizeType ThisIndex) const
    {
        KRATOS_DEBUG_ERROR_IF(!mpVariablesList) << "This container don't have a variables list assigned. A possible reason is creating a node without a model part." << std::endl;
        KRATOS_DEBUG_ERROR_IF_NOT(mpVariablesList->Has(rThisVariable)) << "This container only can store the variables specified in its variables list. The variables list doesn't have this variable:" << rThisVariable << std::endl;
        KRATOS_DEBUG_ERROR_IF((ThisIndex + 1) > mQueueSize) << "Trying to access data from step " << ThisIndex << " but only " << mQueueSize << " steps are stored." << std::endl;
        return Position(ThisIndex) + mpVariablesList->Index(rThisVariable.SourceKey());
    }

    inline BlockType* Position() const;

    inline BlockType* Position(const SizeType ThisIndex) const
    {
        KRATOS_DEBUG_ERROR_IF(!mpVariablesList) << "This container don't have a variables list assigned. A possible reason is creating a node without a model part." << std::endl;
        const SizeType total_size = TotalSize();
        BlockType* position = mpCurrentPosition + ThisIndex * mpVariablesList->DataSize();
        return (position < mpData + total_size) ? position : position - total_size;
    }
    
    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

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

    ///@}

}; // Class VariablesListDataValueContainer

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  VariablesListDataValueContainer& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const VariablesListDataValueContainer& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_VARIABLES_LIST_DATA_VALUE_CONTAINER_H_INCLUDED  defined
