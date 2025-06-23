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
#include "expression/expression.h"
#include "includes/define.h"
#include "includes/smart_pointers.h"


namespace Kratos {


class KRATOS_API(KRATOS_CORE) ExpressionInput
{
public:
    /// @name Member Aliases
    /// @{

    KRATOS_CLASS_POINTER_DEFINITION(ExpressionInput);

    /// @}
    /// @name  Life Cycle
    /// @{

    virtual ~ExpressionInput() = default;

    /// @}
    /// @name Operations
    /// @{

    virtual Expression::Pointer Execute() const = 0;

    Expression::Pointer operator()() const
    {
        return this->Execute();
    }

    /// @}


protected:
    /// @name Protected Operations
    /// @{

    double EvaluateExpression(const Expression& rExpression,
                              Expression::IndexType EntityIndex,
                              Expression::IndexType EntityDataBeginIndex,
                              Expression::IndexType ComponentIndex) const
    {
        return rExpression.Evaluate(EntityIndex, EntityDataBeginIndex, ComponentIndex);
    }

    /// @}
}; // class ExpressionInput


class KRATOS_API(KRATOS_CORE) ExpressionOutput
{
public:
    /// @name Member Aliases
    /// @{

    KRATOS_CLASS_POINTER_DEFINITION(ExpressionOutput);

    /// @}
    /// @name  Life Cycle
    /// @{

    virtual ~ExpressionOutput() = default;

    /// @}
    /// @name Operations
    /// @{

    virtual void Execute(const Expression& rExpression) = 0;

    void operator()(const Expression& rExpression)
    {
        this->Execute(rExpression);
    }

    /// @}

protected:
    /// @name Protected Operations
    /// @{

    double EvaluateExpression(const Expression& rExpression,
                              Expression::IndexType EntityIndex,
                              Expression::IndexType EntityDataBeginIndex,
                              Expression::IndexType ComponentIndex) const
    {
        return rExpression.Evaluate(EntityIndex, EntityDataBeginIndex, ComponentIndex);
    }

    /// @}
}; // class ExpressionOutput


} // namespace Kratos
