//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main author:     Suneth Warnakulasuriya
//                   Máté Kelemen
//

#pragma once

// System includes
#include <type_traits>
#include <iterator>

// Project includes
#include "includes/define.h"
#include "includes/data_communicator.h"
#include "expression/container_expression.h"
#include "expression/expression.h"
#include "expression/unary_reshape_expression.h"
#include "expression/unary_combine_expression.h"

// Application includes

namespace Kratos
{

///@name Kratos Classes
///@{

class KRATOS_API(KRATOS_CORE) ExpressionUtils
{
public:
    ///@name Type definitions
    ///@{

    using IndexType = std::size_t;

    ///@}
    ///@name Static entity value operations
    ///@{

    /**
     * @brief Collapse the lazy expression tree structure in the expression.
     *
     * @details This method returns an expression which is created by collapsing the
     *          tree structure in the @ref rpExpression. This is useful in cases, when the
     *          expression tree becomes large and memory intensive, so the tree can be
     *          collapsed to a one leaf expression releasing memory. This has to evaluate
     *          the expression for each entity, hence this is a computationly expensive task.
     *
     * This method is optimized and compatible with OpenMP and MPI.
     *
     * @param rpExpression Expression to collapse the lazy expression tree.
     * @return Expression::ConstPointer Collapsed expression.
     */
    static Expression::ConstPointer Collapse(const Expression::ConstPointer& rpExpression);

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
    static Expression::ConstPointer Slice(
        const Expression::ConstPointer& rpExpression,
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
     *  @param rpExpression Expression to reshape.
     *  @param NewShapeBegin Iterator pointing to the first component of the new shape.
     *  @param NewShapeEnd Iterator past the last component of the new shape.
     */
    template <class TIterator>
    static Expression::ConstPointer Reshape(
        const Expression::ConstPointer& rpExpression,
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
     *  @param rpExpression Expression to reshape.
     *  @param rNewShape New shape to used to reshape the existing expression.
     */
    static Expression::ConstPointer Reshape(
        const Expression::ConstPointer& rpExpression,
        const std::vector<IndexType>& rNewShape);

    /** @brief Append the components of a set of expressions to the current expression's components.
     *
     *  @details This method combines a set of expressions into the current one as explained in the following example:
     *           All provided expressions in @a rpExpressions must have the same number of items.
     *
     *           For example, let @a rpExpressions contain the following expressions:
     *           @code
     *           rpExpressions[0] = data{1, 2, 3} with 3 items, and item shape = []
     *                             -  -  -
     *           @endcode
     *           @code
     *           rpExpressions[1] = data{4, 5, 6, 7, 8, 9} with 3 items, and item shape = [2]
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
     *  @param rpExpressions Expressions to comb components from.
     */
    static Expression::ConstPointer Comb(const std::vector<Expression::ConstPointer>& rpExpressions);


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
    static Expression::ConstPointer Comb(TIterator ExpressionBegin, TIterator ExpressionEnd)
    {
        static_assert(std::is_same_v<typename std::iterator_traits<TIterator>::value_type,Expression::ConstPointer>
                    || std::is_same_v<typename std::iterator_traits<TIterator>::value_type,Expression::Pointer>);
        return UnaryCombineExpression::Create(ExpressionBegin, ExpressionEnd);
    }

    /**
     * @brief Returns an expression which represents the absolute value of the given expression.
     *
     */
    static Expression::ConstPointer Abs(const Expression::ConstPointer& rpExpression);

    /**
     * @brief Returns an expression which raises the input expression to geven power.
     *
     * @return Expression::ConstPointer
     */
    static Expression::ConstPointer Pow(
        const Expression::ConstPointer& rpExpression,
        const double Power);

    /**
     * @brief Returns an expression which raises the input expression's each component to the specific powers in rPowerpExpression.
     *
     */
    static Expression::ConstPointer Pow(
        const Expression::ConstPointer& rpExpression,
        const Expression::ConstPointer& rpPowerpExpression);

    /**
     * @brief Scales the given expression by a constant value.
     *
     */
    static Expression::ConstPointer Scale(
        const Expression::ConstPointer& rpExpression,
        const double Scale);

    /**
     * @brief Scales the given expression by specific weights for each component given in @ref rScaleExpression.
     *
     */
    static Expression::ConstPointer Scale(
        const Expression::ConstPointer& rpExpression,
        const Expression::ConstPointer& rpScaleExpression);

    /**
     * @brief Returns an expression having min value from all the components for each entity.
     *
     */
    static Expression::ConstPointer EntityMin(const Expression::ConstPointer& rpExpression);

    /**
     * @brief Returns an expression having max value from all the components for each entity.
     *
     */
    static Expression::ConstPointer EntityMax(const Expression::ConstPointer& rpExpression);

    /**
     * @brief Returns an expression having sum of component values for each entity.
     *
     */
    static Expression::ConstPointer EntitySum(const Expression::ConstPointer& rpExpression);

    ///@}
    ///@name Static scalar reduction operations
    ///@{

    /**
     * @brief Returns the sum of the expression assuming it is a flat vector [Shape is not considered].
     *
     */
    static double Sum(
        const Expression::ConstPointer& rpExpression,
        const DataCommunicator& rDataCommunicator);

    /**
     * @brief Returns the inf norm of the expression assuming it is a flat vector [Shape is not considered].
     *
     */
    static double NormInf(
        const Expression::ConstPointer& rpExpression,
        const DataCommunicator& rDataCommunicator);

    /**
     * @brief Returns the L2 norm of the expression assuming it is a flat vector [Shape is not considered].
     *
     */
    static double NormL2(
        const Expression::ConstPointer& rpExpression,
        const DataCommunicator& rDataCommunicator);

    /**
     * @brief Returns the P norm of the expression assuming it is a flat vector [Shape is not considered].
     *
     */
    static double NormP(
        const Expression::ConstPointer& rpExpression,
        const double P,
        const DataCommunicator& rDataCommunicator);

    /**
     * @brief Returns the inner product between two expressions.
     *
     * @details This does not consider shapes of the expressions. They should have same size of flattened vectors.
     *
     * @throws If the flattend size mismatch.
     * @throws If the number of entities mismatch.
     *
     */
    static double InnerProduct(
        const Expression::ConstPointer& rpExpression1,
        const Expression::ConstPointer& rpExpression2,
        const DataCommunicator& rDataCommunicator);

    ///@}
    ///@name Static Container expression operations
    ///@{

    #ifndef KRATOS_EXPRESSION_UTILS_CEXP_METHOD_1
    #define KRATOS_EXPRESSION_UTILS_CEXP_METHOD_1(METHOD_NAME)                      \
        template <class TContainerType>                                             \
        static ContainerExpression<TContainerType> METHOD_NAME(                     \
            const ContainerExpression<TContainerType>& rContainerEXpression)        \
        {                                                                           \
            auto copy = rContainerEXpression;                                       \
            copy.SetExpression(METHOD_NAME(rContainerEXpression.pGetExpression())); \
            return copy;                                                            \
        }
    #endif

    #ifndef KRATOS_EXPRESSION_UTILS_CEXP_METHOD_2
    #define KRATOS_EXPRESSION_UTILS_CEXP_METHOD_2(METHOD_NAME)                                     \
        template <class TContainerType>                                                            \
        static double METHOD_NAME(const ContainerExpression<TContainerType>& rContainerExpression) \
        {                                                                                          \
            return METHOD_NAME(                                                                    \
                rContainerExpression.pGetExpression(),                                             \
                rContainerExpression.GetModelPart().GetCommunicator().GetDataCommunicator());      \
        }
    #endif

    #ifndef KRATOS_EXPRESSION_UTILS_CEXP_METHOD_3
    #define KRATOS_EXPRESSION_UTILS_CEXP_METHOD_3(METHOD_NAME)                                   \
        template <class TContainerType>                                                          \
        static ContainerExpression<TContainerType> METHOD_NAME(                                  \
            const ContainerExpression<TContainerType>& rContainerExpression, const double Value) \
        {                                                                                        \
            auto copy = rContainerExpression;                                                    \
            copy.SetExpression(METHOD_NAME(rContainerExpression.pGetExpression(), Value));       \
            return copy;                                                                         \
        }                                                                                        \
        template <class TContainerType>                                                          \
        static ContainerExpression<TContainerType> METHOD_NAME(                                  \
            const ContainerExpression<TContainerType>& rContainerExpression1,                    \
            const ContainerExpression<TContainerType>& rContainerExpression2)                    \
        {                                                                                        \
            auto copy = rContainerExpression1;                                                   \
            copy.SetExpression(METHOD_NAME(rContainerExpression1.pGetExpression(),               \
                                        rContainerExpression2.pGetExpression()));             \
            return copy;                                                                         \
        }
    #endif

    KRATOS_EXPRESSION_UTILS_CEXP_METHOD_1(Collapse)
    KRATOS_EXPRESSION_UTILS_CEXP_METHOD_1(Abs)
    KRATOS_EXPRESSION_UTILS_CEXP_METHOD_1(EntityMin)
    KRATOS_EXPRESSION_UTILS_CEXP_METHOD_1(EntityMax)
    KRATOS_EXPRESSION_UTILS_CEXP_METHOD_1(EntitySum)
    KRATOS_EXPRESSION_UTILS_CEXP_METHOD_2(Sum)
    KRATOS_EXPRESSION_UTILS_CEXP_METHOD_2(NormInf)
    KRATOS_EXPRESSION_UTILS_CEXP_METHOD_2(NormL2)
    KRATOS_EXPRESSION_UTILS_CEXP_METHOD_3(Pow)
    KRATOS_EXPRESSION_UTILS_CEXP_METHOD_3(Scale)

    #undef KRATOS_EXPRESSION_UTILS_CEXP_METHOD_1
    #undef KRATOS_EXPRESSION_UTILS_CEXP_METHOD_2
    #undef KRATOS_EXPRESSION_UTILS_CEXP_METHOD_3

    template<class TContainerType>
    static ContainerExpression<TContainerType> Slice(
        const ContainerExpression<TContainerType>& rContainerExpression,
        std::size_t Offset,
        std::size_t Stride)
    {
        auto copy = rContainerExpression;
        copy.SetExpression(Slice(rContainerExpression.pGetExpression(), Offset, Stride));
        return copy;
    }

    template <class TContainerType, class TIterator>
    static ContainerExpression<TContainerType> Reshape(
        const ContainerExpression<TContainerType>& rContainerExpression,
        TIterator NewShapeBegin,
        TIterator NewShapeEnd)
    {
        auto copy = rContainerExpression;
        copy.SetExpression(Reshape(rContainerExpression.pGetExpression(), NewShapeBegin, NewShapeEnd));
        return copy;
    }

    template<class TContainerType>
    static ContainerExpression<TContainerType> Reshape(
        const ContainerExpression<TContainerType>& rContainerExpression,
        const std::vector<IndexType>& rNewShape)
    {
        return Reshape(rContainerExpression, rNewShape.begin(), rNewShape.end());
    }

    template <class TContainerType>
    static ContainerExpression<TContainerType> Comb(const std::vector<typename ContainerExpression<TContainerType>::Pointer>& rpContainerExpressions)
    {
        std::vector<Expression::ConstPointer> exps;
        std::transform(rpContainerExpressions.begin(), rpContainerExpressions.end(), std::back_inserter(exps), [](const auto& V) { return V->pGetExpression(); });
        KRATOS_ERROR_IF(rpContainerExpressions.empty())
            << "Empty container list provided for combing.";
        auto copy = *rpContainerExpressions.front();
        copy.SetExpression(UnaryCombineExpression::Create(exps.begin(), exps.end()));
        return copy;
    }

    template <class TContainerType>
    static double NormP(
        const ContainerExpression<TContainerType>& rContainerExpression,
        const double P)
    {
        return NormP(rContainerExpression.pGetExpression(), P, rContainerExpression.GetModelPart().GetCommunicator().GetDataCommunicator());
    }

    template<class TContainerType>
    static double InnerProduct(
        const ContainerExpression<TContainerType>& rContainerExpression1,
        const ContainerExpression<TContainerType>& rContainerExpression2)
    {
        return InnerProduct(rContainerExpression1.pGetExpression(), rContainerExpression2.pGetExpression(), rContainerExpression1.GetModelPart().GetCommunicator().GetDataCommunicator());
    }

    ///@}
};

///@}
}