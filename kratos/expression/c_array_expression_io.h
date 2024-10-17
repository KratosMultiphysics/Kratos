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
#include <vector>

// Project includes
#include "expression/container_expression.h"
#include "expression/expression.h"
#include "expression/traits.h"
#include "includes/define.h"
#include "expression_io.h"

namespace Kratos {

class KRATOS_API(KRATOS_CORE) CArrayExpressionIO
{
public:
    ///@name Public classes
    ///@{

    class KRATOS_API(KRATOS_CORE) Input: public ExpressionInput
    {
    public:
        ///@name Type definitions
        ///@{

        using IndexType = std::size_t;

        using RawArrayType = std::variant<int const*, double const*>;

        KRATOS_CLASS_POINTER_DEFINITION(Input);

        ///@}
        ///@name  Life Cycle
        ///@{

        template<class TRawDataType>
        KRATOS_API(KRATOS_CORE) Input(
            TRawDataType const* pBegin,
            const int NumberOfEntities,
            int const* pShapeBegin,
            const int ShapeSize);

        ~Input() override = default;

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
    }; // class Input


    class KRATOS_API(KRATOS_CORE) MoveInput : public ExpressionInput
    {
    public:
        ///@name Type definitions
        ///@{

        using IndexType = std::size_t;

        using RawArrayType = std::variant<int*, double*>;

        KRATOS_CLASS_POINTER_DEFINITION(MoveInput);

        ///@}
        ///@name  Life Cycle
        ///@{

        template <class TRawDataType>
        MoveInput(
            TRawDataType* pBegin,
            const int NumberOfEntities,
            int const* pShapeBegin,
            const int ShapeSize);

        ~MoveInput() override = default;

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
    }; // class MoveInput

    class KRATOS_API(KRATOS_CORE) Output: public ExpressionOutput
    {
    public:
        ///@name Type definitions
        ///@{

        using IndexType = std::size_t;

        using RawArrayType = std::variant<int*, double*>;

        KRATOS_CLASS_POINTER_DEFINITION(Output);

        ///@}
        ///@name  Life Cycle
        ///@{

        template<class TRawDataType>
        KRATOS_API(KRATOS_CORE) Output(
            TRawDataType* pBegin,
            const int mSize);

        ~Output() override = default;

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
    }; // class Output

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

    template<class TContainerType, MeshType TMeshType>
    static void Read(
        ContainerExpression<TContainerType, TMeshType>& rContainerExpression,
        const Vector& rValues,
        const std::vector<IndexType>& rShape);

    template<class TRawDataType, class TContainerType, MeshType TMeshType>
    static void Move(
        ContainerExpression<TContainerType, TMeshType>& rContainerExpression,
        TRawDataType* pBegin,
        const int NumberOfEntities,
        int const* pShapeBegin,
        const int ShapeSize);

    template<class TContainerType, MeshType TMeshType>
    static void Move(
        ContainerExpression<TContainerType, TMeshType>& rContainerExpression,
        Vector& rValues,
        const std::vector<IndexType>& rShape);

    template<class TRawDataType, class TContainerType, MeshType TMeshType>
    static void Write(
        const ContainerExpression<TContainerType, TMeshType>& rContainerExpression,
        TRawDataType* pBegin,
        const int mSize);

    template<class TContainerType, MeshType TMeshType>
    static void Write(
        const ContainerExpression<TContainerType, TMeshType>& rContainerExpression,
        Vector& rValues);

    ///@}
};

} // namespace Kratos
