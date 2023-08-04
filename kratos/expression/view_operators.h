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

// Project includes
#include "expression/expression.h"
#include "expression/unary_reshape_expression.h"
#include "expression/unary_combine_expression.h"
#include "includes/define.h"


namespace Kratos {


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
KRATOS_API(KRATOS_CORE) Expression::Pointer Slice(const Expression::ConstPointer& rpExpression,
                                                  std::size_t Offset,
                                                  std::size_t Stride);


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
template <class TIterator>
Expression::Pointer Reshape(const Expression::ConstPointer& rpExpression,
                            TIterator NewShapeBegin,
                            TIterator NewShapeEnd)
{
    static_assert(std::is_same_v<typename std::iterator_traits<TIterator>::value_type,std::size_t>);
    return UnaryReshapeExpression::Create(rpExpression, NewShapeBegin, NewShapeEnd);
}


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
Expression::Pointer Reshape(const Expression::ConstPointer& rpExpression,
                            const std::vector<std::size_t>& rNewShape);


/** @brief Append the components of a set of expressions to the current expression's components.
 *
 *  @details This method combines a set of expressions into the current one as explained in the following example:
 *           All provided expressions in @a rExpressions must have the same number of items.
 *
 *           For example, let @a rExpressions contain the following expressions:
 *           @code
 *           rExpressions[0] = data{1, 2, 3} with 3 items, and item shape = []
 *                             -  -  -
 *           @endcode
 *           @code
 *           rExpressions[1] = data{4, 5, 6, 7, 8, 9} with 3 items, and item shape = [2]
 *                             ----  ----  ----
 *           @endcode
 *
 *           The resulting expression has item shape [3] with 3 items:
 *           @code
 *           output_data = [1, 4, 5, 2, 6, 7, 3, 8, 9]
 *                          -------  -------  -------
 *           @endcode
 *
 *           This creates a lazy expression, hence it has a constant cost complexity irrespective of the
 *           data size (The expression won't be evaluated unless @ref Expression::Evaluate is called).
 *
 *  @param rExpressions Expressions to comb components from.
 */
Expression::Pointer Comb(const std::vector<Expression::ConstPointer>& rExpressions);


/** @brief Append the components of a set of expressions to the current expression's components.
 *
 *  @details This method combines a set of expressions into the current one as explained in the following example:
 *           All provided expressions in the { @a ExpressionBegin, @a ExpressionEnd } range must have the same number of items.
 *
 *           For example, let the { @a ExpressionBegin, @a ExpressionEnd } range contain the following expressions:
 *           @code
 *           ExpressionBegin[0] = data{1, 2, 3} with 3 items, and item shape = []
 *                             -  -  -
 *           @endcode
 *           @code
 *           ExpressionBegin[1] = data{4, 5, 6, 7, 8, 9} with 3 items, and item shape = [2]
 *                             ----  ----  ----
 *           @endcode
 *
 *           The resulting expression has item shape [3] with 3 items:
 *           @code
 *           output_data = [1, 4, 5, 2, 6, 7, 3, 8, 9]
 *                          -------  -------  -------
 *           @endcode
 *
 *           This creates a lazy expression, hence it has a constant cost complexity irrespective of the
 *           data size (The expression won't be evaluated unless @ref Expression::Evaluate is called).
 *
 *  @param ExpressionBegin Iterator pointing to the first expression to comb components from.
 *  @param ExpressionEnd Iterator past the last expression to comb components from.
 */
template <class TIterator>
Expression::Pointer Comb(TIterator ExpressionBegin, TIterator ExpressionEnd)
{
    static_assert(std::is_same_v<typename std::iterator_traits<TIterator>::value_type,Expression::ConstPointer>
                  || std::is_same_v<typename std::iterator_traits<TIterator>::value_type,Expression::Pointer>);
    return UnaryCombineExpression::Create(ExpressionBegin, ExpressionEnd);
}


/// @}


} // namespace Kratos
