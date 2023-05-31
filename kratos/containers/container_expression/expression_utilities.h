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
#include "containers/container_expression/expressions/unary/unary_combine_expression.h"
#include "containers/container_expression/expressions/unary/unary_reshape_expression.h"
#include "containers/container_expression/container_expression.h"
#include "includes/define.h"
#include <iterator>


namespace Kratos {


struct KRATOS_API(KRATOS_CORE) ExpressionUtilities
{
    /// @name Copy operations
    /// @{

    /** @brief Clones the existing data container.
     *  @param rSource Source expression to copy from.
     *  @param rTarget Target expression to copy to.
     *  @details Clone an existing expression. This is light weight operation
     *           since this just clones the expression pointer. No underlying
     *           data is copied.
     */
    template <class TContainer>
    static void Clone(const ContainerExpression<TContainer>& rSource,
                      ContainerExpression<TContainer>& rTarget);

    /// @}
    /// @name Views
    /// @{

    /** @brief Returns a slice of the provided expression. Slicing is based on the entity values.
     *
     *  @details Slicing is done on each entitiy's data array, and not on the flattened
     *           expression. For example:
     *
     *           Assume a @ref SpecializedContainerExpression with an expression of shape [5] and 2 entities with
     *           the following data values in the flattened representation.
     *
     *           data = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
     *                   <---- 1 ----> <----- 2 ----->
     *
     *           Data for entity 1 is represented with <--1-->.
     *
     *           If the Offset is 1, and Stride is 3 then, the @ref SpecializedContainerExpression output
     *           of this method will produce the lazy expression which will give following data when
     *           SpecializedContainerExpression::Evaluate is called.
     *
     *           output_data = [2, 3, 4, 7, 8, 9]
     *           output container shape = [3] = equal to Stride.
     *
     *           Slicing will always create a one dimensional vector even if the input has more than one dimensions.
     *           @see Reshape to reshape the one dimensional vector to the desired shape if required.
     *
     *           This creates a lazy expression, hence it has a constant cost complexity irrespective of the data size.
     *
     *  @param Offset Offset of the component to start slicing at.
     *  @param Stride Number of components from the offset in the sliced entity.
     */
    template <class TContainer>
    static void Slice(ContainerExpression<TContainer>& rExpression,
                      std::size_t Offset,
                      std::size_t Stride);

    /** @brief Returns current expression reshaped to specified shape.
     *
     *  @details Reshaping is done on each entitiy's data array, and not on the flattened
     *           expression. For example:
     *
     *           Assume a SpecializedContainerExpression with an expression of shape [2, 3] and 2 entities with
     *           following data values in the flattened representation.
     *
     *           data = [[[1, 2, 3], [4, 5, 6]], [[7, 8, 9], [10, 11, 12]]]
     *                   <-------- 1 --------->  <----------- 2 ----------->
     *
     *           If the rShape = [3, 2] then the returned @ref SpecializedContainerExpression will have a lazy
     *           expression which will return the following output data when Evaluate is called.
     *
     *           output_data = [[[1, 2], [3, 4], [5, 6]], [[7, 8], [9, 10], [11, 12]]]
     *           output container shape = [3, 2]
     *
     *           This creates a lazy expression, hence it has a constant cost complexity irrespective of the data size.
     *
     *  @param rExpression Expression to reshape.
     *  @param rNewShape New shape to used to reshape the existing expression.
     */
    template <class TContainer>
    static void Reshape(ContainerExpression<TContainer>& rExpression,
                        const std::vector<std::size_t>& rNewShape);

    /** @brief Returns current expression reshaped to specified shape.
     *
     *  @details Reshaping is done on each entitiy's data array, and not on the flattened
     *           expression. For example:
     *
     *           Assume a SpecializedContainerExpression with an expression of shape [2, 3] and 2 entities with
     *           following data values in the flattened representation.
     *
     *           data = [[[1, 2, 3], [4, 5, 6]], [[7, 8, 9], [10, 11, 12]]]
     *                   <-------- 1 --------->  <----------- 2 ----------->
     *
     *           If the rShape = [3, 2] then the returned @ref SpecializedContainerExpression will have a lazy
     *           expression which will return the following output data when Evaluate is called.
     *
     *           output_data = [[[1, 2], [3, 4], [5, 6]], [[7, 8], [9, 10], [11, 12]]]
     *           output container shape = [3, 2]
     *
     *           This creates a lazy expression, hence it has a constant cost complexity irrespective of the data size.
     *
     *  @param rExpression Expression to reshape.
     *  @param NewShapeBegin Iterator to the first component of the array defining the new shape.
     *  @param NewShapeEnd Iterator past the last component of the new shape array.
     */
    template <class TContainer, class TIterator>
    static void Reshape(ContainerExpression<TContainer>& rExpression,
                        TIterator NewShapeBegin,
                        TIterator NewShapeEnd)
    {
        rExpression.SetExpression(UnaryReshapeExpression::Create(
            rExpression.pGetExpression(),
            NewShapeBegin,
            NewShapeEnd
        ));
    }

    /** @brief Append the components of a set of expressions to a target expression's components.
     *
     *  @details This method provides a combined container expression as explained in the following example.
     *           All provided container expressions in @a rOthers must have the same number of
     *           items as the @a rTarget container expression. Combing is done in the following order:
     *               1. Components of @a rTarget
     *               2. Components of @a rOthers in the order of the input array
     *
     *           For example, assume the current expression has the following data with item shape [2]
     *           and 3 items in total:
     *           @code
     *           data = [1, 2, 3, 4, 5, 6]
     *                   ----  ----  ----
     *           @endcode
     *           Let @a rOther contain the following data:
     *           rOther = data{7, 8, 9} with 3 items, and item shape = []
     *                         -  -  -
     *
     *           The resulting expression has item shape [3] with 3 items:
     *           @code
     *           output_data = [1, 2, 7, 3, 4, 8, 5, 6, 9]
     *                          -------  -------  -------
     *           @endcode
     *
     *           This creates a lazy expression, hence it has a constant cost complexity irrespective of the
     *           data size (The expression won't be evaluated unless Evaluate is called).
     *
     *  @param rTarget Expression to comb components of other expressions into.
     *  @param rOther Expression to comb components from.
     */
    template <class TContainer>
    static void Comb(ContainerExpression<TContainer>& rTarget,
                     const ContainerExpression<TContainer>& rOther);

    /** @brief Append the components of a set of expressions to a target expression's components.
     *
     *  @details This method provides a combined container expression as explained in the following example.
     *           All provided container expressions in @a rOthers must have the same number of
     *           items as the @a rTarget container expression. Combing is done in the following order:
     *               1. Components of @a rTarget
     *               2. Components of @a rOthers in the order of the input array
     *
     *           For example, assume the current expression has the following data with item shape [2]
     *           and 3 items in total:
     *           @code
     *           data = [1, 2, 3, 4, 5, 6]
     *                   ----  ----  ----
     *           @endcode
     *           Let @a rOthers contain the following expressions:
     *           rOthers[0] = data{7, 8, 9} with 3 items, and item shape = []
     *                             -  -  -
     *           rOthers[1] = data{10, 11, 12, 13, 14, 15} with 3 items, and item shape = [2]
     *                             ------  ------  ------
     *
     *           The resulting expression has item shape [5] with 3 items:
     *           @code
     *           output_data = [1, 2, 7, 10, 11, 3, 4, 8, 12, 13, 5, 6, 9, 14, 15]
     *                          ---------------  ---------------  ---------------
     *           @endcode
     *
     *           This creates a lazy expression, hence it has a constant cost complexity irrespective of the
     *           data size (The expression won't be evaluated unless Evaluate is called).
     *
     *  @param rTarget Expression to comb components of other expressions into.
     *  @param rListOfOthers Expressions to comb components from.
     */
    template <class TContainer>
    static void Comb(ContainerExpression<TContainer>& rTarget,
                     const std::vector<typename ContainerExpression<TContainer>::Pointer>& rOthers);

    /** @brief Append the components of a set of expressions to a target expression's components.
     *
     *  @details This method provides a combined container expression as explained in the following example.
     *           All provided container expressions in @a OthersBegin to @a OthersEnd must have the same number of
     *           items as the @a rTarget container expression. Combing is done in the following order:
     *               1. Components of @a rTarget
     *               2. Components of the rest of the expressions defined in the @a OthersBegin : @a OthersEnd range
     *
     *           For example, assume the current expression has the following data with item shape [2]
     *           and 3 items in total:
     *           @code
     *           data = [1, 2, 3, 4, 5, 6]
     *                   ----  ----  ----
     *           @endcode
     *           Let @a OthersBegin : @a OthersEnd contain the following expressions:
     *           OthersBegin[0] = data{7, 8, 9} with 3 items, and item shape = []
     *                                 -  -  -
     *           OthersBegin[1] = data{10, 11, 12, 13, 14, 15} with 3 items, and item shape = [2]
     *                                 ------  ------  ------
     *
     *           The resulting expression has item shape [5] with 3 items:
     *           @code
     *           output_data = [1, 2, 7, 10, 11, 3, 4, 8, 12, 13, 5, 6, 9, 14, 15]
     *                          ---------------  ---------------  ---------------
     *           @endcode
     *
     *           This creates a lazy expression, hence it has a constant cost complexity irrespective of the
     *           data size (The expression won't be evaluated unless Evaluate is called).
     *
     *  @param rTarget Expression to comb components of other expressions into.
     *  @param OthersBegin Iterator pointing to the first expression to comb components from.
     *  @param OthersEnd Iterator past the last expression to comb components from.
     */
    template <class TContainer, class TIterator>
    static void Comb(ContainerExpression<TContainer>& rTarget,
                     TIterator OthersBegin,
                     TIterator OthersEnd)
    {
        static_assert(std::is_same_v<
            typename std::iterator_traits<TIterator>::value_type,
            typename ContainerExpression<TContainer>::Pointer
        >);
        std::vector<Expression::Pointer> expressions {rTarget.pGetExpression()};
        std::transform(OthersBegin,
                       OthersEnd,
                       std::back_inserter(expressions),
                       [](const auto& rpExpression){return rpExpression->pGetExpression();});
        rTarget.SetExpression(UnaryCombineExpression::Create(
            expressions.begin(),
            expressions.end()
        ));
    }

    /// @}
    /// @name Arithmetic Operations
    /// @{

    /// @brief Raise each component of an expression to the power of an exponent.
    template <class TContainer>
    static void Pow(ContainerExpression<TContainer>& rBase, double Exponent);

    /// @brief Raise each component of a base expression to the power of an exponent expression component-wise.
    template <class TContainer>
    static void Pow(ContainerExpression<TContainer>& rBase,
                    const ContainerExpression<TContainer>& rExponent);

    /// @}
}; // struct ExpressionUtilities


} // namespace Kratos
