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

// Project includes
#include "includes/define.h"
#include "containers/container_expression/expressions/expression.h"
#include "containers/container_expression/expressions/literal/literal_flat_expression.h"

namespace Kratos {

///@name Kratos Classes
///@{

/**
 * @brief Iterator class for expressions generic lazy expressions
 *
 */
template<class TExpressionType>
class ExpressionIterator {
public:
    ///@name Type definitions
    ///@{

    using IndexType = std::size_t;

    using DataType = double;

    using ConstIteratorType = ExpressionIterator<Expression>;

    KRATOS_CLASS_POINTER_DEFINITION(ExpressionIterator);

    ///@}
    ///@name Life cycle
    ///@{

    ExpressionIterator(Expression::Pointer pExpression)
        : mpExpression(pExpression),
          mEntityIndex(0),
          mEntityDataBeginIndex(0),
          mComponentIndex(0),
          mFlattenedShapeSize(pExpression->GetFlattenedShapeSize())
    {}

    ///@}
    ///@name Public operations
    ///@{

    Expression::Pointer GetExpression() { return mpExpression; }

    ConstIteratorType cbegin() { return ConstIteratorType(mpExpression); }

    ConstIteratorType cend()
    {
        ConstIteratorType result(mpExpression);
        result.mEntityIndex = mpExpression->NumberOfEntities();
        result.mEntityDataBeginIndex = result.mEntityIndex * result.mFlattenedShapeSize;
        return result;
    }

    ///@}
    ///@name Public operators
    ///@{

    double operator*() const
    {
        return mpExpression->Evaluate(mEntityIndex, mEntityDataBeginIndex, mComponentIndex);
    }

    bool operator==(const ConstIteratorType& rOther) const
    {
        return (&*mpExpression == &*rOther.mpExpression) && (mEntityIndex == rOther.mEntityIndex) && (mComponentIndex == rOther.mComponentIndex);
    }

    bool operator!=(const ConstIteratorType& rOther) const
    {
        return !this->operator==(rOther);
    }

    ConstIteratorType operator+(const IndexType n) const
    {
        ConstIteratorType result(mpExpression);
        result.mComponentIndex = (mComponentIndex + n) % mFlattenedShapeSize;
        result.mEntityIndex += n / mFlattenedShapeSize;
        result.mEntityDataBeginIndex = result.mEntityIndex * mFlattenedShapeSize;
        return result;
    }

    ConstIteratorType& operator++()
    {
        ++mComponentIndex;
        if (mComponentIndex == mFlattenedShapeSize) {
            mComponentIndex = 0;
            ++mEntityIndex;
            mEntityDataBeginIndex = mEntityIndex * mFlattenedShapeSize;
        }
        return *this;
    }

    ConstIteratorType operator++(int)
    {
        ConstIteratorType temp = *this;
        ++*this;
        return temp;
    }

    ///@}

private:
    ///@name Private member variables
    ///@{

    Expression::Pointer mpExpression;

    IndexType mEntityIndex;

    IndexType mEntityDataBeginIndex;

    IndexType mComponentIndex;

    IndexType mFlattenedShapeSize;

    ///@}
};

/**
 * @brief Iterator class for expressions with concrete data arrays
 *
 */
template<>
class ExpressionIterator<LiteralFlatExpression<char>> {
public:
    ///@name Type definitions
    ///@{

    using ExpressionType = LiteralFlatExpression<char>;

    using DataType = typename ExpressionType::DataType;

    using IteratorType = DataType*;

    using ConstIteratorType = DataType const*;

    KRATOS_CLASS_POINTER_DEFINITION(ExpressionIterator);

    ///@}
    ///@name Life cycle
    ///@{

    ExpressionIterator(typename ExpressionType::Pointer pExpression): mpExpression(pExpression) {}

    ///@}
    ///@name Public operations
    ///@{

    typename ExpressionType::Pointer GetExpression() { return mpExpression; }

    IteratorType begin() { return mpExpression->begin(); }

    IteratorType end() { return mpExpression->end(); }

    ConstIteratorType cbegin() { return mpExpression->cbegin(); }

    ConstIteratorType cend() { return mpExpression->cend(); }

    ///@}

private:
    ///@name Private member variables
    ///@{

    typename ExpressionType::Pointer mpExpression;

    ///@}
};

template<>
class ExpressionIterator<LiteralFlatExpression<int>> {
public:
    ///@name Type definitions
    ///@{

    using ExpressionType = LiteralFlatExpression<int>;

    using DataType = typename ExpressionType::DataType;

    using IteratorType = DataType*;

    using ConstIteratorType = DataType const*;

    KRATOS_CLASS_POINTER_DEFINITION(ExpressionIterator);

    ///@}
    ///@name Life cycle
    ///@{

    ExpressionIterator(typename ExpressionType::Pointer pExpression): mpExpression(pExpression) {}

    ///@}
    ///@name Public operations
    ///@{

    typename ExpressionType::Pointer GetExpression() { return mpExpression; }

    IteratorType begin() { return mpExpression->begin(); }

    IteratorType end() { return mpExpression->end(); }

    ConstIteratorType cbegin() { return mpExpression->cbegin(); }

    ConstIteratorType cend() { return mpExpression->cend(); }

    ///@}

private:
    ///@name Private member variables
    ///@{

    typename ExpressionType::Pointer mpExpression;

    ///@}
};

template<>
class ExpressionIterator<LiteralFlatExpression<double>> {
public:
    ///@name Type definitions
    ///@{

    using ExpressionType = LiteralFlatExpression<double>;

    using DataType = typename ExpressionType::DataType;

    using IteratorType = DataType*;

    using ConstIteratorType = DataType const*;

    KRATOS_CLASS_POINTER_DEFINITION(ExpressionIterator);

    ///@}
    ///@name Life cycle
    ///@{

    ExpressionIterator(typename ExpressionType::Pointer pExpression): mpExpression(pExpression) {}

    ///@}
    ///@name Public operations
    ///@{

    typename ExpressionType::Pointer GetExpression() { return mpExpression; }

    IteratorType begin() { return mpExpression->begin(); }

    IteratorType end() { return mpExpression->end(); }

    ConstIteratorType cbegin() { return mpExpression->cbegin(); }

    ConstIteratorType cend() { return mpExpression->cend(); }

    ///@}

private:
    ///@name Private member variables
    ///@{

    typename ExpressionType::Pointer mpExpression;

    ///@}
};

} // namespace Kratos