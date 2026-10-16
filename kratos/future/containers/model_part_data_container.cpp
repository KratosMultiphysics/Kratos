//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//  Collaborator:    Vicente Mataix Ferrandiz
//
//

// System includes

// External includes

// Project includes
#include "future/containers/model_part_data_container.h"

namespace Kratos::Future
{

namespace
{

// Chunk size heuristic: the requested size, or the current entity count (at least one slot)
std::size_t InitialChunkSize(std::size_t RequestedChunkSize, std::size_t NumberOfEntities)
{
    return (RequestedChunkSize > 0) ? RequestedChunkSize : std::max<std::size_t>(NumberOfEntities, 1);
}

} // anonymous namespace

ModelPartDataContainer::ModelPartDataContainer(ModelPart& rModelPart, std::size_t ChunkSize)
    : mrModelPart(rModelPart),
      mNodes(rModelPart.GetBufferSize(), InitialChunkSize(ChunkSize, rModelPart.NumberOfNodes())),
      mElements(rModelPart.GetBufferSize(), InitialChunkSize(ChunkSize, rModelPart.NumberOfElements())),
      mConditions(rModelPart.GetBufferSize(), InitialChunkSize(ChunkSize, rModelPart.NumberOfConditions())),
      mGeometries(rModelPart.GetBufferSize(), InitialChunkSize(ChunkSize, rModelPart.NumberOfGeometries())),
      mMasterSlaveConstraints(rModelPart.GetBufferSize(), InitialChunkSize(ChunkSize, rModelPart.NumberOfMasterSlaveConstraints()))
{
    KRATOS_TRY

    mNodes.RegisterEntities(mrModelPart.NodesBegin(), mrModelPart.NodesEnd());
    mElements.RegisterEntities(mrModelPart.ElementsBegin(), mrModelPart.ElementsEnd());
    mConditions.RegisterEntities(mrModelPart.ConditionsBegin(), mrModelPart.ConditionsEnd());
    mGeometries.RegisterEntities(mrModelPart.GeometriesBegin(), mrModelPart.GeometriesEnd());
    mMasterSlaveConstraints.RegisterEntities(mrModelPart.MasterSlaveConstraintsBegin(), mrModelPart.MasterSlaveConstraintsEnd());

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void ModelPartDataContainer::Update()
{
    KRATOS_TRY

    mNodes.RegisterEntities(mrModelPart.NodesBegin(), mrModelPart.NodesEnd());
    mElements.RegisterEntities(mrModelPart.ElementsBegin(), mrModelPart.ElementsEnd());
    mConditions.RegisterEntities(mrModelPart.ConditionsBegin(), mrModelPart.ConditionsEnd());
    mGeometries.RegisterEntities(mrModelPart.GeometriesBegin(), mrModelPart.GeometriesEnd());
    mMasterSlaveConstraints.RegisterEntities(mrModelPart.MasterSlaveConstraintsBegin(), mrModelPart.MasterSlaveConstraintsEnd());

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void ModelPartDataContainer::CloneStepData(StepCategory Category)
{
    KRATOS_TRY

    mNodes.CloneStepData(Category);
    mElements.CloneStepData(Category);
    mConditions.CloneStepData(Category);
    mGeometries.CloneStepData(Category);
    mMasterSlaveConstraints.CloneStepData(Category);

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

EntityDataContainer& ModelPartDataContainer::Nodes()
{
    return mNodes;
}

/***********************************************************************************/
/***********************************************************************************/

EntityDataContainer& ModelPartDataContainer::Elements()
{
    return mElements;
}

/***********************************************************************************/
/***********************************************************************************/

EntityDataContainer& ModelPartDataContainer::Conditions()
{
    return mConditions;
}

/***********************************************************************************/
/***********************************************************************************/

EntityDataContainer& ModelPartDataContainer::Geometries()
{
    return mGeometries;
}

/***********************************************************************************/
/***********************************************************************************/

EntityDataContainer& ModelPartDataContainer::MasterSlaveConstraints()
{
    return mMasterSlaveConstraints;
}

/***********************************************************************************/
/***********************************************************************************/

ModelPart& ModelPartDataContainer::GetModelPart()
{
    return mrModelPart;
}

/***********************************************************************************/
/***********************************************************************************/

std::string ModelPartDataContainer::Info() const
{
    return "ModelPartDataContainer of " + mrModelPart.Name();
}

/***********************************************************************************/
/***********************************************************************************/

void ModelPartDataContainer::PrintInfo(std::ostream& rOStream) const
{
    rOStream << Info();
}

/***********************************************************************************/
/***********************************************************************************/

void ModelPartDataContainer::PrintData(std::ostream& rOStream) const
{
    rOStream << "Nodes:" << std::endl;
    mNodes.PrintData(rOStream);
    rOStream << "Elements:" << std::endl;
    mElements.PrintData(rOStream);
    rOStream << "Conditions:" << std::endl;
    mConditions.PrintData(rOStream);
    rOStream << "Geometries:" << std::endl;
    mGeometries.PrintData(rOStream);
    rOStream << "MasterSlaveConstraints:" << std::endl;
    mMasterSlaveConstraints.PrintData(rOStream);
}

} // namespace Kratos::Future
