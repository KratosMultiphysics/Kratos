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
#include "containers/container_expression/container_expression.h"
#include "containers/container_expression/expressions/expression.h"
#include "containers/container_expression/traits.h"
#include "includes/define.h"
#include "expression_io.h"

namespace Kratos {

class KRATOS_API(KRATOS_CORE) CArrayExpressionIO
{
public:
    ///@name Public classes
    ///@{

    class KRATOS_API(KRATOS_CORE) CArrayExpressionInput: public ExpressionInput
    {
    public:
        ///@name Type definitions
        ///@{

        using IndexType = std::size_t;

        using RawArrayType = std::variant<int const*, double const*>;

        KRATOS_CLASS_POINTER_DEFINITION(CArrayExpressionInput);

        ///@}
        ///@name  Life Cycle
        ///@{

        template<class TRawDataType>
        KRATOS_API(KRATOS_CORE) CArrayExpressionInput(
            TRawDataType const* pBegin,
            const int NumberOfEntities,
            int const* pShapeBegin,
            const int ShapeSize);

        ~CArrayExpressionInput() override = default;

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

        std::vector<IndexType> mShape;

        ///@}
    }; // class CArrayExpressionInput


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

    class KRATOS_API(KRATOS_CORE) CArrayExpressionOutput: public ExpressionOutput
    {
    public:
        ///@name Type definitions
        ///@{

        using IndexType = std::size_t;

        using RawArrayType = std::variant<int*, double*>;

        KRATOS_CLASS_POINTER_DEFINITION(CArrayExpressionOutput);

        ///@}
        ///@name  Life Cycle
        ///@{

        template<class TRawDataType>
        KRATOS_API(KRATOS_CORE) CArrayExpressionOutput(
            TRawDataType* pBegin,
            const int mSize);

        ~CArrayExpressionOutput() override = default;

        ///@}
        ///@name Operations
        ///@{

        void Execute(const Expression& rExpression) override;

        ///@}

    private:
        ///@name Private member variables
        ///@{

        RawArrayType mpCArray;

        const int mSize;

        ///@}
    }; // class CArrayExpressionOutput

    ///@}
    ///@name Public static operations
    ///@{

    template<class TRawDataType, class TContainerType, MeshType TMeshType>
    static void Read(
        ContainerExpression<TContainerType, TMeshType>& rContainerExpression,
        TRawDataType const* pBegin,
        const int NumberOfEntities,
        int const* pShapeBegin,
        const int ShapeSize);

    template<class TRawDataType, class TContainerType, MeshType TMeshType>
    static void Move(
        ContainerExpression<TContainerType, TMeshType>& rContainerExpression,
        TRawDataType* pBegin,
        const int NumberOfEntities,
        int const* pShapeBegin,
        const int ShapeSize);

    template<class TRawDataType, class TContainerType, MeshType TMeshType>
    static void Write(
        const ContainerExpression<TContainerType, TMeshType>& rContainerExpression,
        TRawDataType* pBegin,
        const int mSize);

    ///@}
};

} // namespace Kratos
