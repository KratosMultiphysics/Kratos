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
#include <string>
#include <iostream>
#include <cstddef>
#include <cstring>

// External includes

// Project includes
#include "includes/define.h"
#include "containers/variable.h"
#include "containers/variables_list.h"
#include "includes/global_variables.h"

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
    explicit VariablesListDataValueContainer(SizeType NewQueueSize = 1)
        : mQueueSize(NewQueueSize), mpCurrentPosition(0),
          mpData(0), mpVariablesList(nullptr)
    {
        if(!mpVariablesList)
            return;

        // Allcating memory
        Allocate();

        // Setting the current position at the begining of data
        mpCurrentPosition = mpData;

        const SizeType size = mpVariablesList->DataSize();
        for(VariablesList::const_iterator it_variable = mpVariablesList->begin(); it_variable != mpVariablesList->end() ; it_variable++) {
            BlockType* position = Position(*it_variable);
            for(SizeType i = 0 ; i < mQueueSize ; i++) {
                it_variable->AssignZero(position + i * size);
            }
        }
    }

    /// Copy constructor.
    VariablesListDataValueContainer(VariablesListDataValueContainer const& rOther)
        : mQueueSize(rOther.mQueueSize), mpCurrentPosition(0),
          mpData(0), mpVariablesList(rOther.mpVariablesList)
    {
        if(!mpVariablesList)
            return;

        // Allcating memory
        Allocate();

        // Setting the current position with relative source container offset
        mpCurrentPosition = mpData + (rOther.mpCurrentPosition - rOther.mpData);

        const SizeType size = mpVariablesList->DataSize();
        for(VariablesList::const_iterator it_variable = mpVariablesList->begin(); it_variable != mpVariablesList->end() ; it_variable++) {
            const SizeType offset = LocalOffset(*it_variable);
            for(SizeType i = 0 ; i < mQueueSize ; i++) {
                const SizeType total_offset =  offset + i * size;
                it_variable->Copy(rOther.mpData + total_offset, mpData + total_offset);
            }
        }
    }

    /// Variables list constructor.
    VariablesListDataValueContainer(VariablesList::Pointer  pVariablesList, SizeType NewQueueSize = 1)
        : mQueueSize(NewQueueSize), mpCurrentPosition(0),
          mpData(0), mpVariablesList(pVariablesList)
    {
        if(!mpVariablesList)
            return;

        // Allcating memory
        Allocate();

        // Setting the current position at the begining of data
        mpCurrentPosition = mpData;

        const SizeType size = mpVariablesList->DataSize();
        for(VariablesList::const_iterator it_variable = mpVariablesList->begin(); it_variable != mpVariablesList->end() ; it_variable++) {
            BlockType*  position = Position(*it_variable);
            for(SizeType i = 0 ; i < mQueueSize ; i++) {
                it_variable->AssignZero(position + i * size);
            }
        }
    }

    /// Variables list and data constructor
    VariablesListDataValueContainer(VariablesList::Pointer  pVariablesList, BlockType const * ThisData, SizeType NewQueueSize = 1)
        : mQueueSize(NewQueueSize), mpCurrentPosition(0),
          mpData(0), mpVariablesList(pVariablesList)
    {
        if(!mpVariablesList)
            return;

        // Allcating memory
        Allocate();

        // Setting the current position at the begining of data
        mpCurrentPosition = mpData;

        const SizeType size = mpVariablesList->DataSize();
        for(VariablesList::const_iterator it_variable = mpVariablesList->begin(); it_variable != mpVariablesList->end() ; it_variable++) {
            const SizeType offset = LocalOffset(*it_variable);
            for(SizeType i = 0 ; i < mQueueSize ; i++) {
                const SizeType total_offset =  offset + i * size;
                it_variable->Copy(ThisData + total_offset, mpData + total_offset);
            }
        }
    }

    /// Destructor.
    ~VariablesListDataValueContainer()
    {
        Clear();
    }

    ///@}
    ///@name Operators
    ///@{

//       template<class TDataType> const TDataType& operator()(const VariableData& rThisVariable, SizeType QueueIndex = 0) const
// 	{
// 	  return GetValue(rThisVariable, QueueIndex);
// 	}

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

//       template<class TDataType> TDataType& operator[](const VariableData& rThisVariable)
// 	{
// 	  return GetValue(rThisVariable, 0);
// 	}

//       template<class TDataType> const TDataType& operator[](const VariableData& rThisVariable) const
// 	{
// 	  return GetValue(rThisVariable, 0);
// 	}

    template<class TDataType> TDataType& operator[](const Variable<TDataType>& rThisVariable)
    {
        return GetValue(rThisVariable);
    }

    template<class TDataType> const TDataType& operator[](const Variable<TDataType>& rThisVariable) const
    {
        return GetValue(rThisVariable);
    }

    /// Assignment operator.
    VariablesListDataValueContainer& operator=(const VariablesListDataValueContainer& rOther)
    {
        if(rOther.mpVariablesList == 0) {
            Clear();
        } else if((mpVariablesList == rOther.mpVariablesList) && (mQueueSize == rOther.mQueueSize)) {
            mpCurrentPosition = mpData + (rOther.mpCurrentPosition - rOther.mpData);

            const SizeType size = mpVariablesList->DataSize();
            for(VariablesList::const_iterator it_variable = mpVariablesList->begin(); it_variable != mpVariablesList->end() ; it_variable++) {
                const SizeType offset = LocalOffset(*it_variable);
                for(SizeType i = 0 ; i < mQueueSize ; i++) {
                    const SizeType total_offset =  offset + i * size;
                    it_variable->Assign(rOther.mpData + total_offset, mpData + total_offset);
                }
            }
        } else {
            DestructAllElements();

            mQueueSize = rOther.mQueueSize;
            mpVariablesList = rOther.mpVariablesList;

            Reallocate();

            mpCurrentPosition = mpData + (rOther.mpCurrentPosition - rOther.mpData);

            const SizeType size = mpVariablesList->DataSize();
            for(VariablesList::const_iterator it_variable = mpVariablesList->begin(); it_variable != mpVariablesList->end() ; it_variable++) {
                const SizeType offset = LocalOffset(*it_variable);
                for(SizeType i = 0 ; i < mQueueSize ; i++) {
                    const SizeType total_offset =  offset + i * size;
                    it_variable->Copy(rOther.mpData + total_offset, mpData + total_offset);
                }
            }
        }

        return *this;
    }

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

    SizeType Size() const
    {
        if(!mpVariablesList)
            return 0;

        return mpVariablesList->DataSize();
    }

    SizeType QueueSize() const
    {
        return mQueueSize;
    }

    SizeType TotalSize() const
    {
        if(!mpVariablesList)
            return 0;

        return mQueueSize * mpVariablesList->DataSize();
    }

    template<class TDataType> void SetValue(const Variable<TDataType>& rThisVariable, TDataType const& rValue)
    {
        GetValue(rThisVariable) = rValue;
    }

    template<class TDataType> void SetValue(const Variable<TDataType>& rThisVariable, TDataType const& rValue, SizeType QueueIndex)
    {
        GetValue(rThisVariable, QueueIndex) = rValue;
    }

    void Clear()
    {
        DestructAllElements();
        if(mpData)
            free(mpData);

        mpData = 0;
    }


    ///@}
    ///@name Access
    ///@{

    VariablesList::Pointer pGetVariablesList()
    {
        return mpVariablesList;
    }

    const VariablesList::Pointer pGetVariablesList() const
    {
        return mpVariablesList;
    }

    VariablesList& GetVariablesList()
    {
        return *mpVariablesList;
    }

    const VariablesList& GetVariablesList() const
    {
        return *mpVariablesList;
    }

    void SetVariablesList(VariablesList::Pointer pVariablesList)
    {
        DestructAllElements();

        mpVariablesList = pVariablesList;

        if(!mpVariablesList)
            return;

        Reallocate();

        mpCurrentPosition = mpData;

        const SizeType size = mpVariablesList->DataSize();
        for(VariablesList::const_iterator it_variable = mpVariablesList->begin(); it_variable != mpVariablesList->end() ; it_variable++) {
            BlockType*  position = Position(*it_variable);
            for(SizeType i = 0 ; i < mQueueSize ; i++) {
                it_variable->AssignZero(position + i * size);
            }
        }
    }

    void SetVariablesList(VariablesList::Pointer pVariablesList, SizeType ThisQueueSize)
    {
        DestructAllElements();

        mpVariablesList = pVariablesList;

        mQueueSize = ThisQueueSize;

        if(!mpVariablesList)
            return;

        Reallocate();

        mpCurrentPosition = mpData;

        const SizeType size = mpVariablesList->DataSize();
        for(VariablesList::const_iterator it_variable = mpVariablesList->begin(); it_variable != mpVariablesList->end() ; it_variable++) {
            BlockType*  position = Position(*it_variable);
            for(SizeType i = 0 ; i < mQueueSize ; i++) {
                it_variable->AssignZero(position + i * size);
            }
        }
    }

    void Resize(SizeType NewSize)
    {
        if(mQueueSize == NewSize) // if new size is equal to the provious one
            return; // Do nothing!!

        if(!mpVariablesList)
            return;

         if(mQueueSize > NewSize) { // if new size is smaller
            // Destructing elements out of the new size
            for(SizeType i = NewSize ; i < mQueueSize ; i++)
                DestructElements(i);

            const SizeType size = mpVariablesList->DataSize();

            //allocating memory
            BlockType* temp = (BlockType*)malloc(size * sizeof(BlockType) * NewSize);

            //Copying data to allocated memory
            for(SizeType i = 0 ; i < NewSize ; i++)
                memcpy(temp + i * size, Position(i), size * sizeof(BlockType));

            // Updating the queue size
            mQueueSize = NewSize;

            // freeing the old memory
            free(mpData);

            // Setting data pointer to the allocated memory
            mpData = temp;

            // updating the current position
            mpCurrentPosition = mpData;

        } else {
            // Calculating the difference of new and old sizes
            const SizeType difference = NewSize - mQueueSize;

            //keeping the old queue size
            const SizeType old_size = mQueueSize;

            // Getting the relative offset of current position
            const SizeType current_offset = mpCurrentPosition - mpData;

            //Updating the queue size
            mQueueSize = NewSize;

            // Reallocating for new size
            Reallocate();

            SizeType size = mpVariablesList->DataSize();

            // Updating the mpCurrentPosition
            mpCurrentPosition = mpData + current_offset;

            // moving the region after current position to the end
            const SizeType region_size = old_size * size - current_offset;
            memmove(mpCurrentPosition + difference * size, mpCurrentPosition, region_size * sizeof(BlockType));

            // Assigning zero to the added elements
            for(VariablesList::const_iterator it_variable = mpVariablesList->begin(); it_variable != mpVariablesList->end() ; it_variable++) {
                BlockType*  position = mpCurrentPosition + LocalOffset(*it_variable);
                for(SizeType i = 0 ; i < difference ; i++) {
                    it_variable->AssignZero(position + i * size);
                }
            }

            //moving the current position to the moved place
            mpCurrentPosition +=  difference * size;
        }
    }


// 	void Resize(SizeType NewSize, BlockType* SourceData)
// 	{
// 	    if(mQueueSize == NewSize)
// 		return;

// 	    if(mQueueSize > NewSize)
// 	    {
// 		for(SizeType i = NewSize ; i < mQueueSize ; i++)
// 		    DestructElements(i);

// 		mQueueSize = NewSize;

// 		Reallocate();
// 	    }
// 	    else
// 	    {
// 		SizeType old_size = mQueueSize;
// 		mQueueSize = NewSize;
// 		Reallocate();
// 		SizeType size = mpVariablesList->DataSize();
// 		for(VariablesList::const_iterator it_variable = mpVariablesList->begin(); it_variable != mpVariablesList->end() ; it_variable++) {
// 		    SizeType  offset = LocalOffset(*it_variable);
// 		    for(SizeType i = old_size ; i < mQueueSize ; i++)
// 			it_variable->Copy(SourceData + offset, mpData + offset + i * size);
// 		}
// 	    }
// 	}


    BlockType* Data()
    {
        return mpData;
    }

    const BlockType* Data() const
    {
        return mpData;
    }

    BlockType* Data(SizeType QueueIndex)
    {
        return Position(QueueIndex);
    }

    BlockType* Data(VariableData const & rThisVariable)
    {
        KRATOS_DEBUG_ERROR_IF(!mpVariablesList->Has(rThisVariable)) << "Variable " << rThisVariable.Name() << " is not added to this variables list. Stopping" << std::endl;
        return Position(rThisVariable);
    }

    SizeType DataSize()
    {
        if(!mpVariablesList)
            return 0;

        return mpVariablesList->DataSize() * sizeof(BlockType);
    }

    SizeType TotalDataSize()
    {
        if(!mpVariablesList)
            return 0;

        return mpVariablesList->DataSize() * sizeof(BlockType) * mQueueSize;
    }

    void AssignData(BlockType* Source, SizeType QueueIndex)
    {
        AssignData(Source, Position(QueueIndex));
    }

//	void SetData(BlockType* pThisData)
//	    {
//		if(mpData)
//			free(mpData);
//		mpData = pThisData;
//	    }

    void CloneFront()
    {
        if(mQueueSize == 0) {
            Resize(1);
            return;
        }

        if(mQueueSize == 1)
            return;

        KRATOS_DEBUG_ERROR_IF(!mpVariablesList) << "This container don't have a variables list assigned. A possible reason is creating a node without a model part." << std::endl;

        const SizeType size = mpVariablesList->DataSize();
        BlockType* position = (mpCurrentPosition == mpData) ? mpData + TotalSize() - size :  mpCurrentPosition - size;
        AssignData(mpCurrentPosition, position);
        mpCurrentPosition = position;

    }

    void PushFront()
    {
        if(mQueueSize == 0) {
            Resize(1);
            return;
        }

        if(mQueueSize == 1)
            return;

        KRATOS_DEBUG_ERROR_IF(!mpVariablesList) << "This container don't have a variables list assigned. A possible reason is creating a node without a model part." << std::endl;

        const SizeType size = mpVariablesList->DataSize();
        mpCurrentPosition = (mpCurrentPosition == mpData) ? mpData + TotalSize() - size :  mpCurrentPosition - size;
        AssignZero();

    }

    void AssignZero()
    {
        KRATOS_DEBUG_ERROR_IF(!mpVariablesList) << "This container don't have a variables list assigned. A possible reason is creating a node without a model part." << std::endl;
        for(VariablesList::const_iterator it_variable = mpVariablesList->begin(); it_variable != mpVariablesList->end() ; it_variable++) {
            it_variable->AssignZero(mpCurrentPosition + LocalOffset(*it_variable)); 
        }
    }

    void AssignZero(const SizeType QueueIndex)
    {
        KRATOS_DEBUG_ERROR_IF(!mpVariablesList) << "This container don't have a variables list assigned. A possible reason is creating a node without a model part." << std::endl;
        BlockType* position = Position(QueueIndex);
        for(VariablesList::const_iterator it_variable = mpVariablesList->begin(); it_variable != mpVariablesList->end() ; it_variable++) {
            it_variable->AssignZero(position + LocalOffset(*it_variable));
        }
    }

    ///@}
    ///@name Inquiry
    ///@{

    template<class TDataType> bool Has(const Variable<TDataType>& rThisVariable) const
    {
        if(!mpVariablesList)
            return false;

        return mpVariablesList->Has(rThisVariable);
    }

    bool IsEmpty()
    {
        if(!mpVariablesList)
            return true;

        return mpVariablesList->IsEmpty();
    }

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
    void PrintData(std::ostream& rOStream) const
    {
        if(!mpVariablesList)
            rOStream << "No varaibles list is assigned yet." << std::endl;

        for(VariablesList::const_iterator it_variable = mpVariablesList->begin(); it_variable != mpVariablesList->end() ; it_variable++) {
            rOStream <<"    ";
            for(SizeType i = 0 ; i < mQueueSize ; i++) {
                rOStream << i << ": ";
                it_variable->Print(Position(*it_variable, i), rOStream);
                rOStream << "  ";
            }
            rOStream << std::endl;
        }

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


    inline void Allocate()
    {
        KRATOS_DEBUG_ERROR_IF(!mpVariablesList) << "This container don't have a variables list assigned. A possible reason is creating a node without a model part." << std::endl;
        mpData = (BlockType*)malloc(mpVariablesList->DataSize() * sizeof(BlockType) * mQueueSize);
    }

    inline void Reallocate()
    {
        KRATOS_DEBUG_ERROR_IF(!mpVariablesList) << "This container don't have a variables list assigned. A possible reason is creating a node without a model part." << std::endl;
        mpData = (BlockType*)realloc(mpData, mpVariablesList->DataSize() * sizeof(BlockType) * mQueueSize);
    }

    void DestructElements(SizeType ThisIndex)
    {
        if(!mpVariablesList)
            return;

        if(mpData == 0)
            return;
        BlockType* position = Position(ThisIndex);
        for(VariablesList::const_iterator it_variable = mpVariablesList->begin(); it_variable != mpVariablesList->end() ; it_variable++) {
            it_variable->Destruct(position + LocalOffset(*it_variable));
        }
    }

    void DestructAllElements()
    {
        if(!mpVariablesList)
            return;

        if(mpData == 0)
            return;

        const SizeType size = mpVariablesList->DataSize();
        for(VariablesList::const_iterator it_variable = mpVariablesList->begin(); it_variable != mpVariablesList->end() ; it_variable++) {
            BlockType*  position = mpData + LocalOffset(*it_variable);
            for(SizeType i = 0 ; i < mQueueSize ; i++) {
                it_variable->Destruct(position + i * size);
            }
        }
    }


    void AssignData(BlockType* Source, BlockType* Destination)
    {
        KRATOS_DEBUG_ERROR_IF(!mpVariablesList) << "This container don't have a variables list assigned. A possible reason is creating a node without a model part." << std::endl;
        for(VariablesList::const_iterator it_variable = mpVariablesList->begin(); it_variable != mpVariablesList->end() ; it_variable++) {
            const SizeType offset = LocalOffset(*it_variable);
            it_variable->Assign(Source + offset, Destination + offset);
        }
    }

    inline SizeType LocalOffset(VariableData const & rThisVariable) const
    {
        KRATOS_DEBUG_ERROR_IF(!mpVariablesList) << "This container don't have a variables list assigned. A possible reason is creating a node without a model part." << std::endl;
        return mpVariablesList->Index(rThisVariable.SourceKey());
    }

    inline BlockType* Position(VariableData const & rThisVariable) const
    {
        KRATOS_DEBUG_ERROR_IF(!mpVariablesList) << "This container don't have a variables list assigned. A possible reason is creating a node without a model part." << std::endl;
        KRATOS_DEBUG_ERROR_IF_NOT(mpVariablesList->Has(rThisVariable)) << "This container only can store the variables specified in its variables list. The variables list doesn't have this variable:" << rThisVariable << std::endl;
        return mpCurrentPosition + mpVariablesList->Index(rThisVariable.SourceKey());
    }

    inline BlockType* Position(VariableData const & rThisVariable, SizeType ThisIndex) const
    {
        KRATOS_DEBUG_ERROR_IF(!mpVariablesList) << "This container don't have a variables list assigned. A possible reason is creating a node without a model part." << std::endl;
        KRATOS_DEBUG_ERROR_IF_NOT(mpVariablesList->Has(rThisVariable)) << "This container only can store the variables specified in its variables list. The variables list doesn't have this variable:" << rThisVariable << std::endl;
        KRATOS_DEBUG_ERROR_IF((ThisIndex + 1) > mQueueSize) << "Trying to access data from step " << ThisIndex << " but only " << mQueueSize << " steps are stored." << std::endl;
        return Position(ThisIndex) + mpVariablesList->Index(rThisVariable.SourceKey());
    }

    inline BlockType* Position() const
    {
        return mpCurrentPosition;
    }

    inline BlockType* Position(SizeType ThisIndex) const
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

    void save(Serializer& rSerializer) const
    {
        KRATOS_ERROR_IF(!mpVariablesList) << "Cannot save a container with no variables list assigned" << std::endl;
        KRATOS_ERROR_IF(mpData == 0) << "Cannot save an empty variables list container" << std::endl;

        rSerializer.save("Variables List", mpVariablesList);
        rSerializer.save("QueueSize", mQueueSize);
        if(mpVariablesList->DataSize() != 0 )
            rSerializer.save("QueueIndex", SizeType(mpCurrentPosition-mpData)/mpVariablesList->DataSize());
        else
            rSerializer.save("QueueIndex", SizeType(0));

        const SizeType size = mpVariablesList->DataSize();
        for(VariablesList::const_iterator it_variable = mpVariablesList->begin(); it_variable != mpVariablesList->end() ; it_variable++) {
            BlockType*  position = mpData + LocalOffset(*it_variable);
            for(SizeType i = 0 ; i < mQueueSize ; i++) {
                //rSerializer.save("VariableName", it_variable->Name());
                it_variable->Save(rSerializer, position + i * size);
            }
        }
    }

    void load(Serializer& rSerializer)
    {
        rSerializer.load("Variables List", mpVariablesList);
        rSerializer.load("QueueSize", mQueueSize);
        SizeType queue_index;
        rSerializer.load("QueueIndex", queue_index);
        Allocate();

        // Setting the current position at the begining of data
        if(queue_index > mQueueSize)
            KRATOS_THROW_ERROR(std::invalid_argument, "Invalid Queue index loaded : ", queue_index)
            mpCurrentPosition = mpData + queue_index * mpVariablesList->DataSize();

        std::string name;
        for(SizeType i = 0 ; i < mQueueSize ; i++)
            AssignZero(i);

        const SizeType size = mpVariablesList->DataSize();
        for(VariablesList::const_iterator it_variable = mpVariablesList->begin() ;it_variable != mpVariablesList->end() ; it_variable++) {
            BlockType*  position = mpData + LocalOffset(*it_variable);
            for(SizeType i = 0 ; i < mQueueSize ; i++) {
                //rSerializer.load("VariableName", name);
                it_variable->Load(rSerializer, position + i * size);
            }
        }
    }

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
