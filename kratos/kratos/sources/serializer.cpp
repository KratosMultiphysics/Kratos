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



}

