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

// System includes
#include <string>
#include <iostream>
#include <cstddef>
#include <cstring>

// External includes

// Project includes
#include "containers/variables_list_data_value_container.h"
#include "includes/global_variables.h"

namespace Kratos
{

VariablesListDataValueContainer::VariablesListDataValueContainer(const SizeType NewQueueSize)
    : mQueueSize(NewQueueSize), mpCurrentPosition(0),
      mpData(0), mpVariablesList(nullptr)
{
    if(!mpVariablesList)
        return;

    // Allocating memory
    Allocate();

    // Setting the current position at the beginning of data
    mpCurrentPosition = mpData;

    const SizeType size = mpVariablesList->DataSize();
    for(VariablesList::const_iterator it_variable = mpVariablesList->begin(); it_variable != mpVariablesList->end() ; it_variable++) {
        BlockType* position = Position(*it_variable);
        for(SizeType i = 0 ; i < mQueueSize ; i++) {
            it_variable->AssignZero(position + i * size);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

VariablesListDataValueContainer::VariablesListDataValueContainer(VariablesListDataValueContainer const& rOther)
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

/***********************************************************************************/
/***********************************************************************************/

VariablesListDataValueContainer::VariablesListDataValueContainer(
    VariablesList::Pointer pVariablesList, 
    const SizeType NewQueueSize
    ) : mQueueSize(NewQueueSize), mpCurrentPosition(0),
        mpData(0), mpVariablesList(pVariablesList)
{
    if(!mpVariablesList)
        return;

    // Allcating memory
    Allocate();

    // Setting the current position at the beginning of data
    mpCurrentPosition = mpData;

    const SizeType size = mpVariablesList->DataSize();
    for(VariablesList::const_iterator it_variable = mpVariablesList->begin(); it_variable != mpVariablesList->end() ; it_variable++) {
        BlockType*  position = Position(*it_variable);
        for(SizeType i = 0 ; i < mQueueSize ; i++) {
            it_variable->AssignZero(position + i * size);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

VariablesListDataValueContainer::VariablesListDataValueContainer(
    VariablesList::Pointer pVariablesList, 
    BlockType const * ThisData, 
    const SizeType NewQueueSize
    ) : mQueueSize(NewQueueSize), mpCurrentPosition(0),
        mpData(0), mpVariablesList(pVariablesList)
{
    if(!mpVariablesList)
        return;

    // Allcating memory
    Allocate();

    // Setting the current position at the beginning of data
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

/***********************************************************************************/
/***********************************************************************************/

VariablesListDataValueContainer::~VariablesListDataValueContainer()
{
    Clear();
}

/***********************************************************************************/
/***********************************************************************************/

VariablesListDataValueContainer& VariablesListDataValueContainer::operator=(const VariablesListDataValueContainer& rOther)
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

/***********************************************************************************/
/***********************************************************************************/

VariablesListDataValueContainer::SizeType VariablesListDataValueContainer::Size() const
{
    if(!mpVariablesList)
        return 0;

    return mpVariablesList->DataSize();
}

/***********************************************************************************/
/***********************************************************************************/

VariablesListDataValueContainer::SizeType VariablesListDataValueContainer::QueueSize() const
{
    return mQueueSize;
}

/***********************************************************************************/
/***********************************************************************************/

VariablesListDataValueContainer::SizeType VariablesListDataValueContainer::TotalSize() const
{
    if(!mpVariablesList)
        return 0;

    return mQueueSize * mpVariablesList->DataSize();
}

/***********************************************************************************/
/***********************************************************************************/

void VariablesListDataValueContainer::Clear()
{
    DestructAllElements();
    if(mpData)
        free(mpData);

    mpData = 0;
}

/***********************************************************************************/
/***********************************************************************************/

VariablesList::Pointer VariablesListDataValueContainer::pGetVariablesList()
{
    return mpVariablesList;
}

/***********************************************************************************/
/***********************************************************************************/

const VariablesList::Pointer VariablesListDataValueContainer::pGetVariablesList() const
{
    return mpVariablesList;
}

/***********************************************************************************/
/***********************************************************************************/

VariablesList& VariablesListDataValueContainer::GetVariablesList()
{
    return *mpVariablesList;
}

/***********************************************************************************/
/***********************************************************************************/

const VariablesList& VariablesListDataValueContainer::GetVariablesList() const
{
    return *mpVariablesList;
}

/***********************************************************************************/
/***********************************************************************************/

void VariablesListDataValueContainer::SetVariablesList(VariablesList::Pointer pVariablesList)
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

/***********************************************************************************/
/***********************************************************************************/

void VariablesListDataValueContainer::SetVariablesList(VariablesList::Pointer pVariablesList, const SizeType ThisQueueSize)
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

/***********************************************************************************/
/***********************************************************************************/

void VariablesListDataValueContainer::Resize(const SizeType NewSize)
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

/***********************************************************************************/
/***********************************************************************************/

VariablesListDataValueContainer::BlockType* VariablesListDataValueContainer::Data()
{
    return mpData;
}

/***********************************************************************************/
/***********************************************************************************/

const VariablesListDataValueContainer::BlockType* VariablesListDataValueContainer::Data() const
{
    return mpData;
}

/***********************************************************************************/
/***********************************************************************************/

VariablesListDataValueContainer::BlockType* VariablesListDataValueContainer::Data(const SizeType QueueIndex)
{
    return Position(QueueIndex);
}

/***********************************************************************************/
/***********************************************************************************/

VariablesListDataValueContainer::BlockType* VariablesListDataValueContainer::Data(VariableData const & rThisVariable)
{
    KRATOS_DEBUG_ERROR_IF(!mpVariablesList->Has(rThisVariable)) << "Variable " << rThisVariable.Name() << " is not added to this variables list. Stopping" << std::endl;
    return Position(rThisVariable);
}

/***********************************************************************************/
/***********************************************************************************/

VariablesListDataValueContainer::SizeType VariablesListDataValueContainer::DataSize()
{
    if(!mpVariablesList)
        return 0;

    return mpVariablesList->DataSize() * sizeof(BlockType);
}

/***********************************************************************************/
/***********************************************************************************/

VariablesListDataValueContainer::SizeType VariablesListDataValueContainer::TotalDataSize()
{
    if(!mpVariablesList)
        return 0;

    return mpVariablesList->DataSize() * sizeof(BlockType) * mQueueSize;
}

/***********************************************************************************/
/***********************************************************************************/

void VariablesListDataValueContainer::AssignData(BlockType* Source, SizeType QueueIndex)
{
    AssignData(Source, Position(QueueIndex));
}

/***********************************************************************************/
/***********************************************************************************/

void VariablesListDataValueContainer::CloneFront()
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

/***********************************************************************************/
/***********************************************************************************/

void VariablesListDataValueContainer::PushFront()
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

/***********************************************************************************/
/***********************************************************************************/

void VariablesListDataValueContainer::AssignZero()
{
    KRATOS_DEBUG_ERROR_IF(!mpVariablesList) << "This container don't have a variables list assigned. A possible reason is creating a node without a model part." << std::endl;
    for(VariablesList::const_iterator it_variable = mpVariablesList->begin(); it_variable != mpVariablesList->end() ; it_variable++) {
        it_variable->AssignZero(mpCurrentPosition + LocalOffset(*it_variable)); 
    }
}

/***********************************************************************************/
/***********************************************************************************/

void VariablesListDataValueContainer::AssignZero(const SizeType QueueIndex)
{
    KRATOS_DEBUG_ERROR_IF(!mpVariablesList) << "This container don't have a variables list assigned. A possible reason is creating a node without a model part." << std::endl;
    BlockType* position = Position(QueueIndex);
    for(VariablesList::const_iterator it_variable = mpVariablesList->begin(); it_variable != mpVariablesList->end() ; it_variable++) {
        it_variable->AssignZero(position + LocalOffset(*it_variable));
    }
}

/***********************************************************************************/
/***********************************************************************************/

bool VariablesListDataValueContainer::IsEmpty()
{
    if(!mpVariablesList)
        return true;

    return mpVariablesList->IsEmpty();
}

/***********************************************************************************/
/***********************************************************************************/

void VariablesListDataValueContainer::PrintData(std::ostream& rOStream) const
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

/***********************************************************************************/
/***********************************************************************************/

inline void VariablesListDataValueContainer::Allocate()
{
    KRATOS_DEBUG_ERROR_IF(!mpVariablesList) << "This container don't have a variables list assigned. A possible reason is creating a node without a model part." << std::endl;
    mpData = (BlockType*)malloc(mpVariablesList->DataSize() * sizeof(BlockType) * mQueueSize);
}

/***********************************************************************************/
/***********************************************************************************/

inline void VariablesListDataValueContainer::Reallocate()
{
    KRATOS_DEBUG_ERROR_IF(!mpVariablesList) << "This container don't have a variables list assigned. A possible reason is creating a node without a model part." << std::endl;
    mpData = (BlockType*)realloc(mpData, mpVariablesList->DataSize() * sizeof(BlockType) * mQueueSize);
}

/***********************************************************************************/
/***********************************************************************************/

void VariablesListDataValueContainer::DestructElements(const SizeType ThisIndex)
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

/***********************************************************************************/
/***********************************************************************************/

void VariablesListDataValueContainer::DestructAllElements()
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

/***********************************************************************************/
/***********************************************************************************/

void VariablesListDataValueContainer::AssignData(BlockType* Source, BlockType* Destination)
{
    KRATOS_DEBUG_ERROR_IF(!mpVariablesList) << "This container don't have a variables list assigned. A possible reason is creating a node without a model part." << std::endl;
    for(VariablesList::const_iterator it_variable = mpVariablesList->begin(); it_variable != mpVariablesList->end() ; it_variable++) {
        const SizeType offset = LocalOffset(*it_variable);
        it_variable->Assign(Source + offset, Destination + offset);
    }
}

/***********************************************************************************/
/***********************************************************************************/

inline VariablesListDataValueContainer::SizeType VariablesListDataValueContainer::LocalOffset(VariableData const & rThisVariable) const
{
    KRATOS_DEBUG_ERROR_IF(!mpVariablesList) << "This container don't have a variables list assigned. A possible reason is creating a node without a model part." << std::endl;
    return mpVariablesList->Index(rThisVariable.SourceKey());
}

/***********************************************************************************/
/***********************************************************************************/

inline VariablesListDataValueContainer::BlockType* VariablesListDataValueContainer::Position() const
{
    return mpCurrentPosition;
}

/***********************************************************************************/
/***********************************************************************************/

void VariablesListDataValueContainer::save(Serializer& rSerializer) const
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

/***********************************************************************************/
/***********************************************************************************/

void VariablesListDataValueContainer::load(Serializer& rSerializer)
{
    rSerializer.load("Variables List", mpVariablesList);
    rSerializer.load("QueueSize", mQueueSize);
    SizeType queue_index;
    rSerializer.load("QueueIndex", queue_index);
    Allocate();

    // Setting the current position at the beginning of data
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

}  // namespace Kratos.
