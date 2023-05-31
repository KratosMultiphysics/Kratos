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
#include "containers/container_expression/expressions/unary/unary_reshape_expression.h"
#include "containers/container_expression/container_expression.h"
#include "includes/define.h"


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

    /// @}
    /// @name Arithmetic Operations
    /// @{

    template <class TContainer>
    static void Pow(ContainerExpression<TContainer>& rBase, double Exponent);

    template <class TContainer>
    static void Pow(ContainerExpression<TContainer>& rBase,
                    const ContainerExpression<TContainer>& rExponent);

    /// @}
}; // struct ExpressionUtilities


} // namespace Kratos
