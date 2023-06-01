//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Máté Kelemen
//

#pragma once

// Project includes
#include "containers/container_expression/expressions/expression.h"


namespace Kratos {


class ExpressionIO
{
public:
    virtual ~ExpressionIO() = default;

    virtual Expression::Pointer Read() = 0;

    virtual void Write(const Expression& rExpression) = 0;

protected:
    double EvaluateExpression(const Expression& rExpression,
                              Expression::IndexType EntityIndex,
                              Expression::IndexType EntityDataBeginIndex,
                              Expression::IndexType ComponentIndex) const
    {
        return rExpression.Evaluate(EntityIndex, EntityDataBeginIndex, ComponentIndex);
    }
}; // class ExpressionIO


} // namespace Kratos
