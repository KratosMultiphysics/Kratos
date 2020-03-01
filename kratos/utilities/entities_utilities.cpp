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

// System includes

// External includes

// Project includes
#include "utilities/entities_utilities.h"

namespace Kratos
{
void EntitiesUtilities::InitializeAllEntities(ModelPart& rModelPart)
{
    KRATOS_TRY

    EntitiesUtilities::InitializeEntities<Condition>(rModelPart);
    EntitiesUtilities::InitializeEntities<Element>(rModelPart);
    EntitiesUtilities::InitializeEntities<MasterSlaveConstraint>(rModelPart);

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template<>
ModelPart::ElementsContainerType& EntitiesUtilities::GetEntities<Element>(ModelPart& rModelPart)
{
    return rModelPart.Elements();
}

/***********************************************************************************/
/***********************************************************************************/

template<>
ModelPart::ConditionsContainerType& EntitiesUtilities::GetEntities<Condition>(ModelPart& rModelPart)
{
    return rModelPart.Conditions();
}

/***********************************************************************************/
/***********************************************************************************/

template<>
ModelPart::MasterSlaveConstraintContainerType& EntitiesUtilities::GetEntities<MasterSlaveConstraint>(ModelPart& rModelPart)
{
    return rModelPart.MasterSlaveConstraints();
}

} // namespace Kratos
