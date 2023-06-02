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
#include "includes/define.h"


namespace Kratos {


class KRATOS_API(KRATOS_CORE) ExpressionInput
{
public:
    virtual ~ExpressionInput() = default;

    virtual Expression::Pointer Execute() const = 0;

    Expression::Pointer operator()() const
    {
        return this->Execute();
    }

protected:
    double EvaluateExpression(const Expression& rExpression,
                              Expression::IndexType EntityIndex,
                              Expression::IndexType EntityDataBeginIndex,
                              Expression::IndexType ComponentIndex) const
    {
        return rExpression.Evaluate(EntityIndex, EntityDataBeginIndex, ComponentIndex);
    }
}; // class ExpressionInput


class KRATOS_API(KRATOS_CORE) ExpressionOutput
{
public:
    virtual ~ExpressionOutput() = default;

    virtual void Execute(const Expression& rExpression) = 0;

    Expression::Pointer operator()(const Expression& rExpression)
    {
        this->Execute(rExpression);
    }

protected:
    double EvaluateExpression(const Expression& rExpression,
                              Expression::IndexType EntityIndex,
                              Expression::IndexType EntityDataBeginIndex,
                              Expression::IndexType ComponentIndex) const
    {
        return rExpression.Evaluate(EntityIndex, EntityDataBeginIndex, ComponentIndex);
    }
}; // class ExpressionOutput


} // namespace Kratos
