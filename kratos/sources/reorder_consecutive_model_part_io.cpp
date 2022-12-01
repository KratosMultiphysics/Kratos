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

// Project includes
#include "includes/reorder_consecutive_model_part_io.h"


namespace Kratos
{

ModelPartIO::SizeType ReorderConsecutiveModelPartIO::ReorderedNodeId(ModelPartIO::SizeType NodeId)
{
    IdMapType::iterator i = mNodeIdMap.find(NodeId);
    if(i != mNodeIdMap.end())
        return i->second;

    mNodeIdMap.insert(IdMapType::value_type(NodeId, ++mNumberOfNodes));
    return mNumberOfNodes;
}

ModelPartIO::SizeType ReorderConsecutiveModelPartIO::ReorderedElementId(ModelPartIO::SizeType ElementId)
{
    IdMapType::iterator i = mElementIdMap.find(ElementId);
    if(i != mElementIdMap.end())
        return i->second;

    mElementIdMap.insert(IdMapType::value_type(ElementId, ++mNumberOfElements));
    return mNumberOfElements;
}

ModelPartIO::SizeType ReorderConsecutiveModelPartIO::ReorderedConditionId(ModelPartIO::SizeType ConditionId)
{
    IdMapType::iterator i = mConditionIdMap.find(ConditionId);
    if(i != mConditionIdMap.end())
        return i->second;

    mConditionIdMap.insert(IdMapType::value_type(ConditionId, ++mNumberOfConditions));
    return mNumberOfConditions;
}

} // namespace Kratos.
