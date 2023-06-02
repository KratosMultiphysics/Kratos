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

    using Pointer = Kratos::intrusive_ptr<const Expression>;

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

        ExpressionIterator(Expression::Pointer pExpression);

        /**
         * @brief Copy constructor
         *
         * @param rOther
         */
        ExpressionIterator(const ExpressionIterator& rOther);

        ///@}
        ///@name Public operations
        ///@{

        Expression::Pointer GetExpression() const;

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

        Expression::Pointer mpExpression;

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


/// @name View Operators
/// @{

/** @brief Construct an expression containing a subset of the components of all items.
 *
 *  @details Slicing is done on each entitiy's data array, and not on the flattened
 *           expression. For example:
 *
 *           Assume an @ref Expression of shape [5] and 2 entities with
 *           the following data in the flattened representation:
 *           @code
 *           data = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
 *                   <---- 1 ----> <----- 2 ----->
 *           @endcode
 *           Data for entity 1 is represented with <--1-->.
 *
 *           Let @code Offset = 1 @code and @code Stride = 3 @endcode. The resulting sliced expression
 *           then represents the following data:
 *
 *           output_data = [2, 3, 4, 7, 8, 9]
 *           output container shape = [3] = equal to Stride.
 *
 *           Slicing will always create a one dimensional array even if the input expression is multidimensional.
 *           @see Reshape to reshape the one dimensional array to the desired shape if required.
 *
 *           This creates a lazy expression, hence it has a constant cost complexity irrespective of the data size.
 *
 *  @param Offset Index of the first component to begin slicing at.
 *  @param Stride Number of components from the offset in the sliced item.
 */
Expression::Pointer Slice(std::size_t Offset, std::size_t Stride);

/** @brief Construct an expression with identical data but interpreted with a new item shape.
 *
 *  @details Reshaping is done on each entitiy's data array, and not on the flattened
 *           expression. For example:
 *
 *           Assume an @ref Expression of shape [2, 3] and 2 entities with
 *           following data in the flattened representation:
 *           @code
 *           data = [[[1, 2, 3], [4, 5, 6]], [[7, 8, 9], [10, 11, 12]]]
 *                   <-------- 1 --------->  <----------- 2 ----------->
 *           @endcode
 *           The underlying data of the reshaped expression is interpreted as follows:
 *           @code
 *           output_data = [[[1, 2], [3, 4], [5, 6]], [[7, 8], [9, 10], [11, 12]]]
 *           output container shape = [3, 2]
 *           @endcode
 *           This creates a lazy expression, hence it has a constant cost complexity irrespective of the data size.
 *
 *  @param rExpression Expression to reshape.
 *  @param NewShapeBegin Iterator pointing to the first component of the new shape.
 *  @param NewShapeEnd Iterator past the last component of the new shape.
 */
template <class TIterator, std::enable_if_t<std::is_same_v<typename std::iterator_traits<TIterator>::value_type,Expression::Pointer>,bool> = true>
Expression::Pointer Reshape(TIterator NewShapeBegin,
                            TIterator NewShapeEnd);

/** @brief Construct an expression with identical data but interpreted with a new item shape.
 *
 *  @details Reshaping is done on each entitiy's data array, and not on the flattened
 *           expression. For example:
 *
 *           Assume an @ref Expression of shape [2, 3] and 2 entities with
 *           following data in the flattened representation:
 *           @code
 *           data = [[[1, 2, 3], [4, 5, 6]], [[7, 8, 9], [10, 11, 12]]]
 *                   <-------- 1 --------->  <----------- 2 ----------->
 *           @endcode
 *           The underlying data of the reshaped expression is interpreted as follows:
 *           @code
 *           output_data = [[[1, 2], [3, 4], [5, 6]], [[7, 8], [9, 10], [11, 12]]]
 *           output container shape = [3, 2]
 *           @endcode
 *           This creates a lazy expression, hence it has a constant cost complexity irrespective of the data size.
 *
 *  @param rExpression Expression to reshape.
 *  @param rNewShape New shape to used to reshape the existing expression.
 */
Expression::Pointer Reshape(const std::vector<std::size_t>& rNewShape);

/** @brief Append the components of an expression to the current expression's components.
 *
 *  @details Each item in the combed expression contains all components of both the current
 *           and the provided expression in order.
 *
 *           For example, assume the current expression has the following data with item shape [2]
 *           and 3 items in total:
 *           @code
 *           data = [1, 2, 3, 4, 5, 6]
 *                   ----  ----  ----
 *           @endcode
 *           Let @a rpOther contain the following data:
 *           @code
 *           rpOther = data{7, 8, 9} with 3 items, and item shape = []
 *                          -  -  -
 *           @endcode
 *
 *           The resulting expression has item shape [5] with 3 items:
 *           @code
 *           output_data = [1, 2, 7, 3, 4, 8, 5, 6, 9]
 *                          -------  -------  -------
 *           @endcode
 *
 *           This creates a lazy expression, hence it has a constant cost complexity irrespective of the
 *           data size (The expression won't be evaluated unless @ref Expression::Evaluate is called).
 *
 *  @param rTarget Expression to comb components of other expressions into.
 *  @param rpOther Expressions to comb components from.
 */
Expression::Pointer Comb(const Expression::Pointer& rpOther);

/** @brief Append the components of a set of expressions to the current expression's components.
 *
 *  @details This method combines a set of expressions into the current one as explained in the following example:
 *           All provided expressions in @a rOthers must have the same number of
 *           items as the current expression. Combing is done in the following order:
 *               1. Components of the current expression
 *               2. Components of @a rOthers in the order of the input array
 *
 *           For example, assume the current expression has the following data with item shape [2]
 *           and 3 items in total:
 *           @code
 *           data = [1, 2, 3, 4, 5, 6]
 *                   ----  ----  ----
 *           @endcode
 *           Let @a rOthers contain the following expressions:
 *           @code
 *           rOthers[0] = data{7, 8, 9} with 3 items, and item shape = []
 *                             -  -  -
 *           @endcode
 *           @code
 *           rOthers[1] = data{10, 11, 12, 13, 14, 15} with 3 items, and item shape = [2]
 *                             ------  ------  ------
 *           @endcode
 *
 *           The resulting expression has item shape [5] with 3 items:
 *           @code
 *           output_data = [1, 2, 7, 10, 11, 3, 4, 8, 12, 13, 5, 6, 9, 14, 15]
 *                          ---------------  ---------------  ---------------
 *           @endcode
 *
 *           This creates a lazy expression, hence it has a constant cost complexity irrespective of the
 *           data size (The expression won't be evaluated unless @ref Expression::Evaluate is called).
 *
 *  @param rTarget Expression to comb components of other expressions into.
 *  @param rOthers Expressions to comb components from.
 */
Expression::Pointer Comb(const std::vector<Expression::Pointer>& rOthers);

/** @brief Append the components of a set of expressions to the current expression's components.
 *
 *  @details This method combines a set of expressions into the current one as explained in the following example:
 *           All provided expressions in @a rOthers must have the same number of
 *           items as the current expression. Combing is done in the following order:
 *               1. Components of the current expression
 *               2. Components of @a rOthers in the order of the input array
 *
 *           For example, assume the current expression has the following data with item shape [2]
 *           and 3 items in total:
 *           @code
 *           data = [1, 2, 3, 4, 5, 6]
 *                   ----  ----  ----
 *           @endcode
 *           Let @a rOthers contain the following expressions:
 *           @code
 *           rOthers[0] = data{7, 8, 9} with 3 items, and item shape = []
 *                             -  -  -
 *           @endcode
 *           @code
 *           rOthers[1] = data{10, 11, 12, 13, 14, 15} with 3 items, and item shape = [2]
 *                             ------  ------  ------
 *           @endcode
 *
 *           The resulting expression has item shape [5] with 3 items:
 *           @code
 *           output_data = [1, 2, 7, 10, 11, 3, 4, 8, 12, 13, 5, 6, 9, 14, 15]
 *                          ---------------  ---------------  ---------------
 *           @endcode
 *
 *           This creates a lazy expression, hence it has a constant cost complexity irrespective of the
 *           data size (The expression won't be evaluated unless @ref Expression::Evaluate is called).
 *
 *  @param rTarget Expression to comb components of other expressions into.
 *  @param OthersBegin Iterator pointing to the first expression to comb components from.
 *  @param OthersEnd Iterator past the last expression to comb components from.
 */
template <class TIterator>
Expression::Pointer Comb(TIterator OthersBegin, TIterator OthersEnd);

/// @}
/// output stream functions
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const Expression& rThis)
{
    return rOStream << rThis.Info();
}

} // namespace Kratos