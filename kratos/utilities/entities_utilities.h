//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_ENTITIES_UTILITIES)
#define KRATOS_ENTITIES_UTILITIES

// System includes

// External includes

// Project includes
#include "includes/model_part.h"

namespace Kratos
{
/**
 * @namespace EntitiesUtilities
 * @ingroup KratosCore
 * @brief This namespace includes several utilities necessaries for the computation of entities functions in a efficient way
 * @author Vicente Mataix Ferrandiz
 */
namespace EntitiesUtilities
{
    /**
     * @brief This method initializes all the conditions
     * @param rModelPart The model of the problem to solve
     */
    template<class TEntityType>
    PointerVectorSet<TEntityType, IndexedObject>& KRATOS_API(KRATOS_CORE) GetEntities(ModelPart& rModelPart);

    /**
     * @brief This method initializes all the conditions
     * @param rModelPart The model of the problem to solve
     */
    template<class TEntityType>
    void InitializeEntities(ModelPart& rModelPart)
    {
        KRATOS_TRY

        // The number of entities
        auto& r_entities_array = GetEntities<TEntityType>(rModelPart);
        const int number_of_entities = static_cast<int>(r_entities_array.size());
        const auto it_ent_begin = r_entities_array.begin();

        // The current process info
        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

        // Initialize
        #pragma omp parallel
        {
            #pragma omp for schedule(guided, 512)
            for (int i_ent = 0; i_ent < number_of_entities; ++i_ent) {
                auto it_ent = it_ent_begin + i_ent;

                // Detect if the entition is active or not. If the user did not make any choice the entition
                // It is active by default
                bool entition_is_active = true;
                if (it_ent->IsDefined(ACTIVE))
                    entition_is_active = it_ent->Is(ACTIVE);

                if (entition_is_active) {
                    it_ent->Initialize(r_current_process_info);
                }
            }
        }

        KRATOS_CATCH("")
    }

    /**
     * @brief This method initializes all the entities (conditions, elements, constraints)
     * @param rModelPart The model of the problem to solve
     */
    void InitializeAllEntities(ModelPart& rModelPart)
    {
        KRATOS_TRY

        InitializeEntities<Condition>(rModelPart);
        InitializeEntities<Element>(rModelPart);
        InitializeEntities<MasterSlaveConstraint>(rModelPart);

        KRATOS_CATCH("")
    }

    ///@}
    ///@name Private  Acces
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

}; // namespace EntitiesUtilities
}  // namespace Kratos
#endif /* KRATOS_ENTITIES_UTILITIES defined */
