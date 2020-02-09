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
 * @namespace EntitiesUtilities
 * @ingroup KratosCore
 * @brief This namespace includes several utilities necessaries for the computation of entities functions in a efficient way
 * @author Vicente Mataix Ferrandiz
 */
namespace EntitiesUtilities
{
    /**
     * @brief This method initializes all the entities (conditions, elements, constraints)
     * @param rModelPart The model of the problem to solve
     */
    void KRATOS_API(KRATOS_CORE) InitializeEntities(ModelPart& rModelPart);

    /**
     * @brief This method initializes all the conditions
     * @param rModelPart The model of the problem to solve
     */
    void KRATOS_API(KRATOS_CORE) InitializeConditions(ModelPart& rModelPart);

    /**
     * @brief This method initializes all the elements
     * @param rModelPart The model of the problem to solve
     */
    void KRATOS_API(KRATOS_CORE) InitializeElements(ModelPart& rModelPart);

    /**
     * @brief This method initializes all the constraints
     * @param rModelPart The model of the problem to solve
     */
    void KRATOS_API(KRATOS_CORE) InitializeMasterSlaveConstraints(ModelPart& rModelPart);

}; // namespace EntitiesUtilities
}  // namespace Kratos
#endif /* KRATOS_ENTITIES_UTILITIES defined */
