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
//                   Philipp Bucher (https://github.com/philbucher)
//

// System includes

// External includes

// Project includes
#include "utilities/entities_utilities.h"

namespace Kratos
{
namespace EntitiesUtilities
{

template<> KRATOS_API(KRATOS_CORE) PointerVectorSet<Element, IndexedObject>& GetEntities<Element>(ModelPart& rModelPart);
template<> KRATOS_API(KRATOS_CORE) PointerVectorSet<Condition, IndexedObject>& GetEntities<Condition>(ModelPart& rModelPart);
template<> KRATOS_API(KRATOS_CORE) PointerVectorSet<MasterSlaveConstraint, IndexedObject>& GetEntities<MasterSlaveConstraint>(ModelPart& rModelPart);

void InitializeAllEntities(ModelPart& rModelPart)
{
    KRATOS_TRY

    InitializeEntities<Element>(rModelPart);
    InitializeEntities<Condition>(rModelPart);
    InitializeEntities<MasterSlaveConstraint>(rModelPart);

    KRATOS_CATCH("")
}

void InitializeSolutionStepAllEntities(ModelPart& rModelPart)
{
    KRATOS_TRY

    InitializeSolutionStepEntities<Element>(rModelPart);
    InitializeSolutionStepEntities<Condition>(rModelPart);
    InitializeSolutionStepEntities<MasterSlaveConstraint>(rModelPart);

    KRATOS_CATCH("")
}

void FinalizeSolutionStepAllEntities(ModelPart& rModelPart)
{
    KRATOS_TRY

    FinalizeSolutionStepEntities<Element>(rModelPart);
    FinalizeSolutionStepEntities<Condition>(rModelPart);
    FinalizeSolutionStepEntities<MasterSlaveConstraint>(rModelPart);

    KRATOS_CATCH("")
}

void InitializeNonLinearIterationAllEntities(ModelPart& rModelPart)
{
    KRATOS_TRY

    InitializeNonLinearIterationEntities<Element>(rModelPart);
    InitializeNonLinearIterationEntities<Condition>(rModelPart);
    InitializeNonLinearIterationEntities<MasterSlaveConstraint>(rModelPart);

    KRATOS_CATCH("")
}

void FinalizeNonLinearIterationAllEntities(ModelPart& rModelPart)
{
    KRATOS_TRY

    FinalizeNonLinearIterationEntities<Element>(rModelPart);
    FinalizeNonLinearIterationEntities<Condition>(rModelPart);
    FinalizeNonLinearIterationEntities<MasterSlaveConstraint>(rModelPart);

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

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

/***********************************************************************************/
/***********************************************************************************/

} // namespace EntitiesUtilities
} // namespace Kratos
