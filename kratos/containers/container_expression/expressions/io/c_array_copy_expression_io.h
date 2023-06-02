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
//

#pragma once

// System includes
#include <variant>

// Project includes
#include "containers/container_expression/expressions/expression.h"
#include "includes/define.h"
#include "expression_io.h"


namespace Kratos {


class KRATOS_API(KRATOS_CORE) CArrayCopyExpressionInput: public ExpressionInput
{
public:
    ///@name Type definitions
    ///@{

    using IndexType = std::size_t;

    using RawArrayType = std::variant<int const*, double const*>;

    KRATOS_CLASS_POINTER_DEFINITION(CArrayCopyExpressionInput);

    ///@}
    ///@name  Life Cycle
    ///@{

    template<class TRawDataType>
    CArrayCopyExpressionInput(
        TRawDataType const* pBegin,
        const int NumberOfEntities,
        int const* pShapeBegin,
        const int ShapeSize);

    ~CArrayCopyExpressionInput() override = default;

    ///@}
    ///@name Operations
    ///@{

    Expression::Pointer Execute() const override;

    ///@}

private:
    ///@name Private member variables
    ///@{

    RawArrayType mpCArray;

    const int mNumberOfEntities;

    int const * const mpShapeBegin;

    const int mShapeSize;

    ///@}
}; // class CArrayCopyExpressionInput


class KRATOS_API(KRATOS_CORE) CArrayCopyExpressionOutput: public ExpressionOutput
{
public:
    ///@name Type definitions
    ///@{

    using IndexType = std::size_t;

    using RawArrayType = std::variant<int*, double*>;

    KRATOS_CLASS_POINTER_DEFINITION(CArrayCopyExpressionOutput);

    ///@}
    ///@name  Life Cycle
    ///@{

    template<class TRawDataType>
    CArrayCopyExpressionOutput(
        TRawDataType* pBegin,
        const int NumberOfEntities,
        int const* pShapeBegin,
        const int ShapeSize);

    ~CArrayCopyExpressionOutput() override = default;

    ///@}
    ///@name Operations
    ///@{

    void Execute(const Expression& rExpression) override;

    ///@}

private:
    ///@name Private member variables
    ///@{

    RawArrayType mpCArray;

    const int mNumberOfEntities;

    int const * const mpShapeBegin;

    const int mShapeSize;

    ///@}
}; // class CArrayCopyExpressionOutput


} // namespace Kratos
