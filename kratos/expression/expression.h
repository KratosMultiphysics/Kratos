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

#pragma once

// System includes
#include <atomic>
#include <ostream>
#include <string>
#include <vector>
#include <iterator>

// Project includes
#include "includes/define.h"
#include "intrusive_ptr/intrusive_ptr.hpp"

namespace Kratos {

///@name Kratos Classes
///@{

/**
 * @brief Base class or all the expression types.
 *
 * This is used to represent an expression for arithmetic evaluation
 *
 */
class KRATOS_API(KRATOS_CORE) Expression {
public:

    ///@name Type definitions
    ///@{

    using Pointer = Kratos::intrusive_ptr<Expression>;

    using ConstPointer = Kratos::intrusive_ptr<const Expression>;

    using IndexType = std::size_t;

    ///@}
    ///@name Public classes
    ///@{

    class ExpressionIterator {
    public:
        ///@name Type definitions
        ///@{

        using difference_type = IndexType;

        using value_type = double;

        using pointer = double*;

        using reference = double&;

        using iterator_category = std::input_iterator_tag;

        KRATOS_CLASS_POINTER_DEFINITION(ExpressionIterator);

        ///@}
        ///@name Life cycle
        ///@{

        ExpressionIterator();

        ExpressionIterator(Expression::ConstPointer pExpression);

        /**
         * @brief Copy constructor
         *
         * @param rOther
         */
        ExpressionIterator(const ExpressionIterator& rOther);

        ///@}
        ///@name Public operations
        ///@{

        Expression::ConstPointer GetExpression() const;

        ///@}
        ///@name Public operators
        ///@{

        double operator*() const;

        bool operator==(const ExpressionIterator& rOther) const;

        bool operator!=(const ExpressionIterator& rOther) const;

        ExpressionIterator& operator=(const ExpressionIterator& rOther);

        ExpressionIterator& operator++();

        ExpressionIterator operator++(int);

        ///@}

    private:
        ///@name Private member variables
        ///@{

        Expression::ConstPointer mpExpression;

        IndexType mEntityIndex;

        IndexType mEntityDataBeginIndex;

        IndexType mItemComponentIndex;

        IndexType mItemComponentCount;

        ///@}
        ///@name Friend classes
        ///@{

        friend class Expression;

        ///@}
    };

    ///@}
    ///@name Iterator type definitions
    ///@{

    using value_type = double;

    using size_type = IndexType;

    using const_iterator = ExpressionIterator;

    ///@}
    ///@name Life cycle
    ///@{

    Expression(const IndexType NumberOfEntities) : mNumberOfEntities(NumberOfEntities) {}

    virtual ~Expression() = default;

    ///@}
    ///@name Public operations
    ///@{

    /**
     * @brief Evalute the expression for the given entity data start index and component index and returns the value
     *
     * @param EntityIndex           Index of the entity.
     * @param EntityDataBeginIndex  Index at which entity data starts.
     * @param ComponentIndex        Component index.
     * @return double               Evaluated expression.
     * @todo Move to private.
     */
    virtual double Evaluate(
        const IndexType EntityIndex,
        const IndexType EntityDataBeginIndex,
        const IndexType ComponentIndex) const = 0;

    /**
     * @brief Get the Shape of the expression
     *
     * @return const std::vector<IndexType>     Size of each dimension is in the vector elements.
     */
    virtual const std::vector<IndexType> GetItemShape() const = 0;

    /**
     * @brief Get the maximum number of entities allowed for this expression.
     *
     * @return IndexType
     */
    inline IndexType NumberOfEntities() const { return mNumberOfEntities; };

    /**
     * @brief Get the Local Size of the expression
     *
     * @return IndexType
     */
    IndexType GetItemComponentCount() const;

    /**
     * @brief Get the Max Depth of the lazy expression tree.
     *
     * Returns the maximum depth of the lazy expression tree.
     *
     * @warning This is a recursive computation, hence this should not
     *          be done repeatedly unless necessary.
     *
     * @return IndexType Max depth of the lazy expression tree.
     */
    virtual IndexType GetMaxDepth() const = 0;

    ///@}
    ///@name Input and output
    ///@{

    IndexType size() const;

    const_iterator begin() const;

    const_iterator end() const;

    const_iterator cbegin() const;

    const_iterator cend() const;

    virtual std::string Info() const = 0;

    ///@}

private:
    friend class ExpressionInput;

    friend class ExpressionOutput;

    ///@name Private member variables
    ///@{

    const IndexType mNumberOfEntities;

    ///@}
    ///@name Private operations
    ///@{

    //*********************************************
    // this block is needed for refcounting in the @ref intrusive ptr
    mutable std::atomic<int> mReferenceCounter{0};

    friend void intrusive_ptr_add_ref(const Expression* x)
    {
        x->mReferenceCounter.fetch_add(1, std::memory_order_relaxed);
    }

    friend void intrusive_ptr_release(const Expression* x)
    {
        if (x->mReferenceCounter.fetch_sub(1, std::memory_order_release) == 1) {
            std::atomic_thread_fence(std::memory_order_acquire);
            delete x;
        }
    }

    //*********************************************

    ///@}
};

/// @}
/// output stream functions
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const Expression& rThis)
{
    return rOStream << rThis.Info();
}

} // namespace Kratos