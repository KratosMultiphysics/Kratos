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


#include "includes/serializer.h"
#include "containers/model.h"
#include "containers/variable.h"
#include "includes/kratos_components.h"


namespace Kratos
{

KRATOS_CREATE_LOCAL_FLAG( Serializer, MPI,  0 );
KRATOS_CREATE_LOCAL_FLAG( Serializer, SHALLOW_GLOBAL_POINTERS_SERIALIZATION, 1 );


Serializer::RegisteredObjectsContainerType Serializer::msRegisteredObjects;

Serializer::RegisteredObjectsNameContainerType Serializer::msRegisteredObjectsName;

VariableData* Serializer::GetVariableData(std::string const & VariableName)
{
    return KratosComponents<VariableData>::pGet(VariableName);
}

/// Sets the Serializer in a state ready to be loaded
/// Note: If the same object is loaded twice before deleting it from memory all its pointers will be duplicated.
void Serializer::SetLoadState() {
    mLoadedPointers.clear();
    SeekBegin();
}

/// Sets the pointer of the stream buffer at the begnining
void Serializer::SeekBegin() {
    mpBuffer->seekg(0, mpBuffer->beg);
}

/// Sets the pointer of the stream buffer at tht end 
void Serializer::SeekEnd() {
    mpBuffer->seekg(0, mpBuffer->end);
}

void Serializer::load(std::string const & rTag, ModelPart*& pValue)
{
    PointerType pointer_type = SP_INVALID_POINTER;
    void* p_pointer;
    read(pointer_type);

    if(pointer_type != SP_INVALID_POINTER)
    {
        read(p_pointer);
        LoadedPointersContainerType::iterator i_pointer = mLoadedPointers.find(p_pointer);
        if(i_pointer == mLoadedPointers.end())
        {
            if(pointer_type == SP_BASE_CLASS_POINTER)
            {
                KRATOS_ERROR_IF(!pValue) << "an already constructed modelpart must be passed to load a ModelPart" <<std::endl;
            }
            else if(pointer_type == SP_DERIVED_CLASS_POINTER)
            {
                KRATOS_ERROR << "should not find SP_DERIVED_CLASS_POINTER for ModelPart load" << std::endl;
            }

            // Load the pointer address before loading the content
            mLoadedPointers[p_pointer]=&pValue;
            load(rTag, *pValue);
        }
        else
        {
            KRATOS_ERROR << "modelpart has already been serialized - should not be done twice!" << std::endl;
        }
    }
}

void Serializer::load(std::string const & rTag, Kratos::unique_ptr<ModelPart>& pValue)
{
    PointerType pointer_type = SP_INVALID_POINTER;
    void* p_pointer;
    read(pointer_type);

    if(pointer_type != SP_INVALID_POINTER)
    {
        read(p_pointer);
        LoadedPointersContainerType::iterator i_pointer = mLoadedPointers.find(p_pointer);
        if(i_pointer == mLoadedPointers.end())
        {
            if(pointer_type == SP_BASE_CLASS_POINTER)
            {
                KRATOS_ERROR_IF(!pValue) << "an already constructed modelpart must be passed to load a ModelPart" <<std::endl;
            }
            else if(pointer_type == SP_DERIVED_CLASS_POINTER)
            {
                KRATOS_ERROR << "should not find SP_DERIVED_CLASS_POINTER for ModelPart load" << std::endl;
            }

            // Load the pointer address before loading the content
            mLoadedPointers[p_pointer]=&pValue;
            load(rTag, *pValue);
        }
        else
        {
            KRATOS_ERROR << "modelpart has already been serialized - should not be done twice!" << std::endl;
        }
    }
}

void Serializer::load(std::string const & rTag, Kratos::shared_ptr<ModelPart>& pValue)
{
    PointerType pointer_type = SP_INVALID_POINTER;
    void* p_pointer;
    read(pointer_type);

    if(pointer_type != SP_INVALID_POINTER)
    {
        read(p_pointer);
        LoadedPointersContainerType::iterator i_pointer = mLoadedPointers.find(p_pointer);
        if(i_pointer == mLoadedPointers.end())
        {
            if(pointer_type == SP_BASE_CLASS_POINTER)
            {
                KRATOS_ERROR_IF(!pValue) << "an already constructed modelpart must be passed to load a ModelPart" <<std::endl;
            }
            else if(pointer_type == SP_DERIVED_CLASS_POINTER)
            {
                KRATOS_ERROR << "should not find SP_DERIVED_CLASS_POINTER for ModelPart load" << std::endl;
            }

            // Load the pointer address before loading the content
            mLoadedPointers[p_pointer]=&pValue;
            load(rTag, *pValue);
        }
        else
        {
            KRATOS_ERROR << "modelpart has already been serialized - should not be done twice!" << std::endl;
        }
    }
}


}

