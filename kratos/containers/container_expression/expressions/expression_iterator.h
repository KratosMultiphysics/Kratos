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
#include <iterator>

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

    using value_type = double;

    using IndexType = std::size_t;

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

    /**
     * @brief Copy constructor
     *
     * @param rOther
     */
    ExpressionIterator(const ExpressionIterator& rOther)
        : mpExpression(rOther.mpExpression),
          mEntityIndex(rOther.mEntityIndex),
          mEntityDataBeginIndex(rOther.mEntityDataBeginIndex),
          mComponentIndex(rOther.mComponentIndex),
          mFlattenedShapeSize(rOther.mFlattenedShapeSize)
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

    ConstIteratorType& operator=(const ExpressionIterator& rOther)
    {
        mpExpression = rOther.mpExpression;
        mEntityIndex = rOther.mEntityIndex;
        mEntityDataBeginIndex = rOther.mEntityDataBeginIndex;
        mComponentIndex = rOther.mComponentIndex;
        mFlattenedShapeSize = rOther.mFlattenedShapeSize;
        return *this;
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

#define KRATOS_EXPRESSION_ITERATOR_SPECIALIZED_DEFINITION(DATA_TYPE)\
template<>                                                          \
class ExpressionIterator<LiteralFlatExpression<DATA_TYPE>> {        \
public:                                                             \
                                                                    \
    using ExpressionType = LiteralFlatExpression<DATA_TYPE>;        \
                                                                    \
    using value_type = typename ExpressionType::DataType;           \
                                                                    \
    using IteratorType = value_type*;                               \
                                                                    \
    using ConstIteratorType = value_type const*;                    \
                                                                    \
    KRATOS_CLASS_POINTER_DEFINITION(ExpressionIterator);            \
                                                                    \
    ExpressionIterator(typename ExpressionType::Pointer pExpression)\
        : mpExpression(pExpression) {}                              \
                                                                    \
    typename ExpressionType::Pointer GetExpression() {              \
        return mpExpression;                                        \
    }                                                               \
                                                                    \
    IteratorType begin() { return mpExpression->begin(); }          \
                                                                    \
    IteratorType end() { return mpExpression->end(); }              \
                                                                    \
    ConstIteratorType cbegin() { return mpExpression->cbegin(); }   \
                                                                    \
    ConstIteratorType cend() { return mpExpression->cend(); }       \
                                                                    \
private:                                                            \
                                                                    \
    typename ExpressionType::Pointer mpExpression;                  \
                                                                    \
};

KRATOS_EXPRESSION_ITERATOR_SPECIALIZED_DEFINITION(char)
KRATOS_EXPRESSION_ITERATOR_SPECIALIZED_DEFINITION(int)
KRATOS_EXPRESSION_ITERATOR_SPECIALIZED_DEFINITION(double)

#undef KRATOS_EXPRESSION_ITERATOR_SPECIALIZED_DEFINITION

} // namespace Kratos