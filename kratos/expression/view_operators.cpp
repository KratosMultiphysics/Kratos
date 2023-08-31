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
#include "expression/view_operators.h"
#include "expression/unary_slice_expression.h"
#include "expression/unary_reshape_expression.h"
#include "expression/unary_combine_expression.h"


namespace Kratos {


Expression::Pointer Slice(const Expression::ConstPointer& rpExpression,
                          std::size_t Offset,
                          std::size_t Stride)
{
    return UnarySliceExpression::Create(
        rpExpression,
        Offset,
        Stride
    );
}


Expression::Pointer Reshape(const Expression::ConstPointer& rpExpression,
                            const std::vector<std::size_t>& rNewShape)
{
    return Reshape(rpExpression, rNewShape.begin(), rNewShape.end());
}


Expression::Pointer Comb(const std::vector<Expression::ConstPointer>& rExpressions)
{
    return Comb(rExpressions.begin(), rExpressions.end());
}


} // namespace Kratos
