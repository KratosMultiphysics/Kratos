//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//                   Máté Kelemen
//

// Project includes
#include "containers/container_expression/expression_utilities.h"


namespace Kratos {


template <class TContainer>
void ExpressionUtilities::Clone(const ContainerExpression<TContainer>& rSource,
                                ContainerExpression<TContainer>& rTarget)
{
    rTarget.SetExpression(rSource.pGetExpression());
}


#define KRATOS_INSTANTIATE_EXPRESSION_UTILITY(CONTAINER_TYPE)                                               \
    template void ExpressionUtilities::Clone<CONTAINER_TYPE>(const ContainerExpression<CONTAINER_TYPE>&,    \
                                                             ContainerExpression<CONTAINER_TYPE>&)


KRATOS_INSTANTIATE_EXPRESSION_UTILITY(ModelPart::NodesContainerType);


KRATOS_INSTANTIATE_EXPRESSION_UTILITY(ModelPart::ElementsContainerType);


KRATOS_INSTANTIATE_EXPRESSION_UTILITY(ModelPart::ConditionsContainerType);


#undef KRATOS_INSTANTIATE_EXPRESSION_UTILITY


} // namespace Kratos
