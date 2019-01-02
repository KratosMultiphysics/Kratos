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

Serializer::RegisteredObjectsContainerType Serializer::msRegisteredObjects;

Serializer::RegisteredObjectsNameContainerType Serializer::msRegisteredObjectsName;

VariableData* Serializer::GetVariableData(std::string const & VariableName)
{
    return KratosComponents<VariableData>::pGet(VariableName);
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

                load(rTag, *pValue);
            }
            else if(pointer_type == SP_DERIVED_CLASS_POINTER)
            {
                KRATOS_ERROR << "should not find SP_DERIVED_CLASS_POINTER for ModelPart load" << std::endl;
            }
            mLoadedPointers[p_pointer]=&pValue;
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

                load(rTag, *pValue);
            }
            else if(pointer_type == SP_DERIVED_CLASS_POINTER)
            {
                KRATOS_ERROR << "should not find SP_DERIVED_CLASS_POINTER for ModelPart load" << std::endl;
            }
            mLoadedPointers[p_pointer]=&pValue;
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

                load(rTag, *pValue);
            }
            else if(pointer_type == SP_DERIVED_CLASS_POINTER)
            {
                KRATOS_ERROR << "should not find SP_DERIVED_CLASS_POINTER for ModelPart load" << std::endl;
            }
            mLoadedPointers[p_pointer]=&pValue;
        }
        else
        {
            KRATOS_ERROR << "modelpart has already been serialized - should not be done twice!" << std::endl;
        }
    }
}


}

