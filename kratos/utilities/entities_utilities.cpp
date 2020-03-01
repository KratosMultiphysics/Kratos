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
namespace EntitiesUtilities
{
template<>
PointerVectorSet<Element, IndexedObject>& GetEntities<Element>(ModelPart& rModelPart)
{
    return rModelPart.Elements();
}

/***********************************************************************************/
/***********************************************************************************/

template<>
PointerVectorSet<Condition, IndexedObject>& GetEntities<Condition>(ModelPart& rModelPart)
{
    return rModelPart.Conditions();
}

/***********************************************************************************/
/***********************************************************************************/

template<>
PointerVectorSet<MasterSlaveConstraint, IndexedObject>& GetEntities<MasterSlaveConstraint>(ModelPart& rModelPart)
{
    return rModelPart.MasterSlaveConstraints();
}

} // namespace EntitiesUtilities
} // namespace Kratos
