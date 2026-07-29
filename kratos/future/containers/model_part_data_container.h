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

#pragma once

// System includes

// External includes

// Project includes
#include "includes/model_part.h"
#include "future/containers/entity_data_container.h"

namespace Kratos::Future
{
///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

/**
 * @class ModelPartDataContainer
 * @ingroup KratosCore
 * @brief Experimental DataContainer-backed storage for the entities of a ModelPart.
 * @details Aggregates one EntityDataContainer per entity kind of a ModelPart (nodes; elements, conditions and geometries; master-slave constraints), snapshotting the entities present at construction: each kind's container is created with the ModelPart's buffer size and a chunk size matching the entity count, and every entity is registered in the container's sorted-by-Id iteration order.
 *
 * This is a PARALLEL storage path (Phase II): the ModelPart and its entities are never modified; their existing storage (solution-step data, DataValueContainer) is completely independent from the data stored here.
 *
 * Semantics and limitations:
 * - Holds a reference to the ModelPart: the manager must not outlive the owning Model.
 * - Entities added to the ModelPart after construction are picked up by Update(); removals are NOT tracked (unregister manually, leaving a slot hole — see EntityDataContainer).
 * - The buffer size is snapshotted at construction; a later ModelPart::SetBufferSize is not reflected.
 * - Element/Condition storage is keyed by the element/condition Id: unlike the legacy path (where GetValue reaches the shared Geometry data), two elements sharing a geometry get INDEPENDENT values here. This is intentional.
 * - Advance historical data by calling CloneStepData(StepCategory::TimeStep) right after ModelPart::CloneTimeStep.
 * @see EntityDataContainer, DataContainer
 * @author Pooyan Dadvand
 */
class KRATOS_API(KRATOS_CORE) ModelPartDataContainer
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ModelPartDataContainer
    KRATOS_CLASS_POINTER_DEFINITION(ModelPartDataContainer);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Constructor.
     * @details Snapshots the entities of the given ModelPart: registers every current entity of each supported kind, using the ModelPart's buffer size for historical variables.
     * @param rModelPart The model part whose entities are managed. Must outlive this object.
     * @param ChunkSize Initial slot capacity per kind; 0 (default) sizes each kind to max(current entity count, 1).
     */
    explicit ModelPartDataContainer(ModelPart& rModelPart, std::size_t ChunkSize = 0);

    /// Destructor.
    ~ModelPartDataContainer() = default;

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Register the entities added to the ModelPart since construction (or the last Update).
     * @details Registration is idempotent, so already known entities are untouched. Removals from the ModelPart are NOT tracked.
     */
    void Update();

    /**
     * @brief Clone the current step onto the next step slot for all entity kinds.
     * @details Call right after ModelPart::CloneTimeStep to advance the DataContainer-backed historical variables in lockstep with the legacy solution-step data.
     * @param Category The step category being cloned.
     */
    void CloneStepData(StepCategory Category);

    ///@}
    ///@name Access
    ///@{

    /**
     * @brief Get the storage of the nodes.
     * @return Reference to the nodes' EntityDataContainer.
     */
    EntityDataContainer& Nodes();

    /**
     * @brief Get the storage of the elements.
     * @details Keyed by element Id — independent per element, unlike the legacy GetValue path reaching the shared Geometry data.
     * @return Reference to the elements' EntityDataContainer.
     */
    EntityDataContainer& Elements();

    /**
     * @brief Get the storage of the conditions.
     * @details Keyed by condition Id — independent per condition, unlike the legacy GetValue path reaching the shared Geometry data.
     * @return Reference to the conditions' EntityDataContainer.
     */
    EntityDataContainer& Conditions();

    /**
     * @brief Get the storage of the geometries.
     * @details Covers the geometries stored in the ModelPart geometry container. Standalone geometries can use a raw EntityDataContainer, provided they carry a meaningful unique Id.
     * @return Reference to the geometries' EntityDataContainer.
     */
    EntityDataContainer& Geometries();

    /**
     * @brief Get the storage of the master-slave constraints.
     * @return Reference to the constraints' EntityDataContainer.
     */
    EntityDataContainer& MasterSlaveConstraints();

    /**
     * @brief Get the managed model part.
     * @return Reference to the model part.
     */
    ModelPart& GetModelPart();

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

private:
    ///@name Member Variables
    ///@{

    ModelPart& mrModelPart;             /// The managed model part (non-owning; must outlive this object)

    EntityDataContainer mNodes;         /// DataContainer-backed storage of the nodes

    EntityDataContainer mElements;      /// DataContainer-backed storage of the elements (keyed by element Id)

    EntityDataContainer mConditions;    /// DataContainer-backed storage of the conditions (keyed by condition Id)

    EntityDataContainer mGeometries;    /// DataContainer-backed storage of the geometries

    EntityDataContainer mMasterSlaveConstraints; /// DataContainer-backed storage of the master-slave constraints

    ///@}

}; // Class ModelPartDataContainer

///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator << (
    std::ostream& rOStream,
    const ModelPartDataContainer& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

///@} addtogroup block

} // namespace Kratos::Future
