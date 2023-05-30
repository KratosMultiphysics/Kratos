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
#include "containers/container_expression/expressions/literal/literal_expression.h"
#include "containers/container_expression/expressions/literal/literal_flat_expression.h"
#include "containers/container_expression/expressions/binary/binary_expression.h"


namespace Kratos {


template <class TContainer>
void ExpressionUtilities::Clone(const ContainerExpression<TContainer>& rSource,
                                ContainerExpression<TContainer>& rTarget)
{
    rTarget.SetExpression(rSource.pGetExpression());
}


template <class TContainer>
void ExpressionUtilities::Pow(ContainerExpression<TContainer>& rBase, double Exponent)
{
    rBase.SetExpression(BinaryExpression<BinaryOperations::Power>::Create(
        rBase.pGetExpression(),
        LiteralExpression<double>::Create(Exponent, rBase.GetContainer().size())
    ));
}


template <class TContainer>
void ExpressionUtilities::Pow(ContainerExpression<TContainer>& rBase,
                              const ContainerExpression<TContainer>& rExponent)
{
    KRATOS_ERROR_IF(rBase.GetContainer().size() != rExponent.GetContainer().size())
        << "Size mismatch in power operation.\n"
        << "\tLeft operand : " << rBase << '\n'
        << "\tRight operand: " << rExponent;
    rBase.SetExpression(BinaryExpression<BinaryOperations::Power>::Create(
        rBase.pGetExpression(),
        rExponent.pGetExpression()
    ));
}


#define KRATOS_INSTANTIATE_EXPRESSION_UTILITY(CONTAINER_TYPE)                                               \
    template void ExpressionUtilities::Clone<CONTAINER_TYPE>(const ContainerExpression<CONTAINER_TYPE>&,    \
                                                             ContainerExpression<CONTAINER_TYPE>&);         \
    template void ExpressionUtilities::Pow<CONTAINER_TYPE>(ContainerExpression<CONTAINER_TYPE>&,            \
                                                           double);                                         \
    template void ExpressionUtilities::Pow<CONTAINER_TYPE>(ContainerExpression<CONTAINER_TYPE>&,            \
                                                           const ContainerExpression<CONTAINER_TYPE>&)


KRATOS_INSTANTIATE_EXPRESSION_UTILITY(ModelPart::NodesContainerType);


KRATOS_INSTANTIATE_EXPRESSION_UTILITY(ModelPart::ElementsContainerType);


KRATOS_INSTANTIATE_EXPRESSION_UTILITY(ModelPart::ConditionsContainerType);


#undef KRATOS_INSTANTIATE_EXPRESSION_UTILITY


} // namespace Kratos
