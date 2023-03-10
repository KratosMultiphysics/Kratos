//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//                   Philipp Bucher (https://github.com/philbucher)
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/model_part.h"
#include "utilities/parallel_utilities.h"

namespace Kratos
{
/**
 * @namespace EntitiesUtilities
 * @ingroup KratosCore
 * @brief This namespace includes several utilities necessaries for the computation of entities functions in a efficient way
 * @author Vicente Mataix Ferrandiz
 * @author Philipp Bucher
 */
namespace EntitiesUtilities
{
    /**
     * @brief This method initializes all the active entities (conditions, elements, constraints)
     * @param rModelPart The model part of the problem to solve
     */
    void KRATOS_API(KRATOS_CORE) InitializeAllEntities(ModelPart& rModelPart);

    /**
     * @brief This method calls InitializeSolution for all the entities (conditions, elements, constraints)
     * @param rModelPart The model part of the problem to solve
     */
    void KRATOS_API(KRATOS_CORE) InitializeSolutionStepAllEntities(ModelPart& rModelPart);

    /**
     * @brief This method calls FinalizeSolutionStep for all the entities (conditions, elements, constraints)
     * @param rModelPart The model part of the problem to solve
     */
    void KRATOS_API(KRATOS_CORE) FinalizeSolutionStepAllEntities(ModelPart& rModelPart);

    /**
     * @brief This method calls InitializeNonLinearIteration for all the entities (conditions, elements, constraints)
     * @param rModelPart The model part of the problem to solve
     */
    void KRATOS_API(KRATOS_CORE) InitializeNonLinearIterationAllEntities(ModelPart& rModelPart);

    /**
     * @brief This method calls FinalizeNonLinearIteration for all the entities (conditions, elements, constraints)
     * @param rModelPart The model part of the problem to solve
     */
    void KRATOS_API(KRATOS_CORE) FinalizeNonLinearIterationAllEntities(ModelPart& rModelPart);

    /**
     * @brief This method returns the appropriate TEntitytype container (elements, conditions, and nodes) from model part
     * @param rModelPart The model part of the problem to solve
     */
    template<class TEntityType>
    KRATOS_API(KRATOS_CORE) PointerVectorSet<TEntityType, IndexedObject>& GetEntities(ModelPart& rModelPart);

    /**
     * @brief This method initializes all the active entities
     * @param rModelPart The model part of the problem to solve
     */
    template<class TEntityType>
    void InitializeEntities(ModelPart& rModelPart)
    {
        KRATOS_TRY

        // Array of entities
        auto& r_entities_array = GetEntities<TEntityType>(rModelPart);

        // The current process info
        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

        // Initialize
        block_for_each(
            r_entities_array,
            [&r_current_process_info](TEntityType& rEntity) {
                // If the entity is active
                if (rEntity.IsActive()) {
                    rEntity.Initialize(r_current_process_info);
                }
            }
        );

        KRATOS_CATCH("")
    }

    /**
     * @brief This method calls InitializeSolutionStep for all the entities
     * @param rModelPart The model part of the problem to solve
     */
    template<class TEntityType>
    void InitializeSolutionStepEntities(ModelPart& rModelPart)
    {
        KRATOS_TRY

        // The current process info
        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

        // Initialize
        block_for_each(
            GetEntities<TEntityType>(rModelPart),
            [&r_current_process_info](TEntityType& rEntity){
                rEntity.InitializeSolutionStep(r_current_process_info);
            }
        );

        KRATOS_CATCH("")
    }

    /**
     * @brief This method calls FinalizeSolutionStep for all the entities
     * @param rModelPart The model part of the problem to solve
     */
    template<class TEntityType>
    void FinalizeSolutionStepEntities(ModelPart& rModelPart)
    {
        KRATOS_TRY

        // The current process info
        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

        // Initialize
        block_for_each(
            GetEntities<TEntityType>(rModelPart),
            [&r_current_process_info](TEntityType& rEntity){
                rEntity.FinalizeSolutionStep(r_current_process_info);
            }
        );

        KRATOS_CATCH("")
    }

    /**
     * @brief This method calls InitializeNonLinearIteration for all the entities
     * @param rModelPart The model part of the problem to solve
     */
    template<class TEntityType>
    void InitializeNonLinearIterationEntities(ModelPart& rModelPart)
    {
        KRATOS_TRY

        // The current process info
        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

        // Initialize
        block_for_each(
            GetEntities<TEntityType>(rModelPart),
            [&r_current_process_info](TEntityType& rEntity){
                rEntity.InitializeNonLinearIteration(r_current_process_info);
            }
        );

        KRATOS_CATCH("")
    }

    /**
     * @brief This method calls FinalizeNonLinearIteration for all the entities
     * @param rModelPart The model part of the problem to solve
     */
    template<class TEntityType>
    void FinalizeNonLinearIterationEntities(ModelPart& rModelPart)
    {
        KRATOS_TRY

        // The current process info
        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

        // Initialize
        block_for_each(
            GetEntities<TEntityType>(rModelPart),
            [&r_current_process_info](TEntityType& rEntity){
                rEntity.FinalizeNonLinearIteration(r_current_process_info);
            }
        );

        KRATOS_CATCH("")
    }

    ///@}

}; // namespace EntitiesUtilities
}  // namespace Kratos
