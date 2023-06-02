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
#include "expression_io.h"
#include "includes/define.h"

namespace Kratos {

class KRATOS_API(KRATOS_CORE) CArrayMoveExpressionInput : public ExpressionInput
{
public:
    ///@name Type definitions
    ///@{

    using IndexType = std::size_t;

    using RawArrayType = std::variant<int*, double*>;

    KRATOS_CLASS_POINTER_DEFINITION(CArrayMoveExpressionInput);

    ///@}
    ///@name  Life Cycle
    ///@{

    template <class TRawDataType>
    CArrayMoveExpressionInput(
        TRawDataType* pBegin,
        const int NumberOfEntities,
        int const* pShapeBegin,
        const int ShapeSize);

    ~CArrayMoveExpressionInput() override = default;

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

    int const* const mpShapeBegin;

    const int mShapeSize;

    ///@}
}; // class CArrayMoveExpressionInput

} // namespace Kratos
