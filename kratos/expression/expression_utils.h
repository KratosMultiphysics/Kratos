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
     * @brief Returns an expression which represents the component wise absolute value of the given expression.
     *
     * @details If the input @a rpExpression is \f$\underline{\mathbb{u}}\f$, where \f$u_{ij}\f$ representss \f$j^{th}\f$ component of the flattened entity data for \f$i^{th}\f$ entity
     *          , then the returned expression \f$\left|\underline{\mathbb{u}}\right|\f$
     *
     *          \f[
     *              \left|\underline{\mathbb{u}}\right| = \left|u_{ij}\right|
     *          \f]
     *
     * @return Expression::ConstPointer Expression which computes component wise absolute value.
     */
    static Expression::ConstPointer Abs(const Expression::ConstPointer& rpExpression);

    /**
     * @brief Returns an expression which raises each component to the given power.
     *
     * @details If the input @a rpExpression is \f$\underline{\mathbb{u}}\f$, where \f$u_{ij}\f$ representss \f$j^{th}\f$ component of the flattened entity data for \f$i^{th}\f$ entity
     *          , where P is specified by @a Power.
     *
     *          \f[
     *              Pow(\underline{\mathbb{u}}, P) = u_{ij}^P
     *          \f]
     *
     * @return Expression::ConstPointer Expression which raises each component to specified power.
     */
    static Expression::ConstPointer Pow(
        const Expression::ConstPointer& rpExpression,
        const double Power);

    /**
     * @brief Returns an expression which raises each component to the given power from another expression.
     *
     * @details If the input @a rpExpression is \f$\underline{\mathbb{u}}\f$, where \f$u_{ij}\f$ representss \f$j^{th}\f$ component of the flattened entity data for \f$i^{th}\f$ entity
     *          and the @a rpPowerpExpression is \f$\underline{\mathbb{P}}\f$, where \f$p_{ij}\f$ representss \f$j^{th}\f$ component of the flattened entity data for \f$i^{th}\f$ entity
     *          , then the returned expression can be illustrated as below.
     *
     *          \f[
     *              Pow(\underline{\mathbb{u}}, \underline{\mathbb{p}}) = u_{ij}^{p_{ij}}
     *          \f]
     *
     *          If the given @a rpPowerpExpression is a scalar expression then the following will be returned
     *
     *          \f[
     *              Pow(\underline{\mathbb{u}}, \underline{\mathbb{p}}) = u_{ij}^{p_{i}}
     *          \f]
     * @throws If the number of entities mismatch between @a rpExpression and @a rpPowerpExpression
     * @throws If the shape of the @a rpExpression and @a rpPowerpExpression does not match or @a rpPowerpExpression is not representing a scalar expression.
     *
     * @return Expression::ConstPointer Expression which raises each component to specified power given by the @a rpPowerpExpression.
     */
    static Expression::ConstPointer Pow(
        const Expression::ConstPointer& rpExpression,
        const Expression::ConstPointer& rpPowerpExpression);

    /**
     * @brief Returns an expression which scales each component to the specified value
     *
     * @details If the input @a rpExpression is \f$\underline{\mathbb{u}}\f$, where \f$u_{ij}\f$ representss \f$j^{th}\f$ component of the flattened entity data for \f$i^{th}\f$ entity
     *          , where Scale is specified by @a Scale.
     *
     *          \f[
     *              Scale(\underline{\mathbb{u}}, S) = u_{ij}^S
     *          \f]
     *
     * @return Expression::ConstPointer Expression which scales the given expression by specified scale @a Scale
     */
    static Expression::ConstPointer Scale(
        const Expression::ConstPointer& rpExpression,
        const double Scale);

    /**
     * @brief Returns an expression which scales each component by a value from another expression.
     *
     * @details If the input @a rpExpression is \f$\underline{\mathbb{u}}\f$, where \f$u_{ij}\f$ representss \f$j^{th}\f$ component of the flattened entity data for \f$i^{th}\f$ entity
     *          and the @a rpScaleExpression is \f$\underline{\mathbb{s}}\f$, where \f$s_{ij}\f$ representss \f$j^{th}\f$ component of the flattened entity data for \f$i^{th}\f$ entity
     *          , then the returned expression can be illustrated as below.
     *
     *          \f[
     *              Scale(\underline{\mathbb{u}}, \underline{\mathbb{s}}) = u_{ij}^{s_{ij}}
     *          \f]
     *
     *          If the given @a rpScaleExpression is a scalar expression then the following will be returned
     *
     *          \f[
     *              Scale(\underline{\mathbb{u}}, \underline{\mathbb{s}}) = u_{ij}^{s_{i}}
     *          \f]
     * @throws If the number of entities mismatch between @a rpExpression and @a rpScaleExpression
     * @throws If the shape of the @a rpExpression and @a rpScaleExpression does not match or @a rpScaleExpression is not representing a scalar expression.
     *
     * @return Expression::ConstPointer Expression which scales each component to specified value given by the @a rpScaleExpression.
     */
    static Expression::ConstPointer Scale(
        const Expression::ConstPointer& rpExpression,
        const Expression::ConstPointer& rpScaleExpression);

    /**
     * @brief Returns an expression having min value from all the components for each entity.
     *
     * @details If the input @a rpExpression is \f$\underline{\mathbb{u}}\f$, where \f$u_{ij}\f$ representss \f$j^{th}\f$ component of the flattened (having \f$N\f$ total components) entity data for \f$i^{th}\f$ entity
     *          , Following illustrates the returned expression which is always a scalar expression having \f$m_i\f$ representing the \f$i^{th}\f$ entity data.
     *
     *          \f[
     *              EntityMin(\underline{\mathbb{u}}) = m_{i}
     *          \f]
     *
     *          Where,
     *          \f[
     *              m_{i} = \min_{j\in \left[0, N\right)} u_{ij}
     *          \f]
     *
     * @return Expression::ConstPointer Scalar expression having the minimum of all components for each entity in the input expression.
     */
    static Expression::ConstPointer EntityMin(const Expression::ConstPointer& rpExpression);

    /**
     * @brief Returns an expression having max value from all the components for each entity.
     *
     * @details If the input @a rpExpression is \f$\underline{\mathbb{u}}\f$, where \f$u_{ij}\f$ representss \f$j^{th}\f$ component of the flattened (having \f$N\f$ total components) entity data for \f$i^{th}\f$ entity
     *          , Following illustrates the returned expression which is always a scalar expression having \f$m_i\f$ representing the \f$i^{th}\f$ entity data.
     *
     *          \f[
     *              EntityMax(\underline{\mathbb{u}}) = m_{i}
     *          \f]
     *
     *          Where,
     *          \f[
     *              m_{i} = \max_{j\in \left[0, N\right)} u_{ij}
     *          \f]
     *
     * @return Expression::ConstPointer Scalar expression having the maximum of all components for each entity in the input expression.
     */
    static Expression::ConstPointer EntityMax(const Expression::ConstPointer& rpExpression);

    /**
     * @brief Returns an expression having sum of component values for each entity.
     *
     * @details If the input @a rpExpression is \f$\underline{\mathbb{u}}\f$, where \f$u_{ij}\f$ representss \f$j^{th}\f$ component of the flattened (having \f$N\f$ total components) entity data for \f$i^{th}\f$ entity
     *          , Following illustrates the returned expression which is always a scalar expression having \f$m_i\f$ representing the \f$i^{th}\f$ entity data.
     *
     *          \f[
     *              EntitySum(\underline{\mathbb{u}}) = m_{i}
     *          \f]
     *
     *          Where,
     *          \f[
     *              m_{i} = \sum_{j\in \left[0, N\right)} u_{ij}
     *          \f]
     *
     * @return Expression::ConstPointer Scalar expression having the sum of all components for each entity in the input expression.
     */
    static Expression::ConstPointer EntitySum(const Expression::ConstPointer& rpExpression);

    ///@}
    ///@name Static scalar reduction operations
    ///@{

    /**
     * @brief Returns the sum of the expression assuming it is a flat vector [Shape is not considered].
     *
     * @details If the input @a rpExpression is \f$\underline{\mathbb{u}}\f$, where \f$u_{ij}\f$ representss \f$j^{th}\f$ component of the flattened (having \f$N\f$ total components) entity data for \f$i^{th}\f$ entity (having \f$M\f$ entities)
     *          , Following illustrates the returned value where the entity data is flattened.
     *
     *          \f[
     *              Sum(\underline{\mathbb{u}}) = \sum_{i\in \left[0, M\right)} \sum_{j\in \left[0, N\right)} u_{ij}
     *          \f]
     *
     *          This method is compatible with shared and distributed memory parallelized runs.
     *
     * @param rpExpression          Expressions to be summed.
     * @param rDataCommunicator     Data communicator for MPI communication.
     * @return double               The sum of all the components in all the entities.
     */
    static double Sum(
        const Expression::ConstPointer& rpExpression,
        const DataCommunicator& rDataCommunicator);

    /**
     * @brief Returns the infinity norm of the expression assuming it is a flat vector [Shape is not considered].
     *
     * @details If the input @a rpExpression is \f$\underline{\mathbb{u}}\f$, where \f$u_{ij}\f$ representss \f$j^{th}\f$ component of the flattened (having \f$N\f$ total components) entity data for \f$i^{th}\f$ entity (having \f$M\f$ entities)
     *          , Following illustrates the returned value where the entity data is flattened.
     *
     *          \f[
     *              NormInf(\underline{\mathbb{u}}) = \max_{(i,j)\in \left[0, M\right)\times \left[0, N\right)} \left|u_{ij}\right|
     *          \f]
     *
     *          This method is compatible with shared and distributed memory parallelized runs.
     *
     * @param rpExpression          Expressions to be summed.
     * @param rDataCommunicator     Data communicator for MPI communication.
     * @return double               The infinity norm of all the components in all the entities.
     */
    static double NormInf(
        const Expression::ConstPointer& rpExpression,
        const DataCommunicator& rDataCommunicator);

    /**
     * @brief Returns the L2 norm of the expression assuming it is a flat vector [Shape is not considered].
     *
     * @details If the input @a rpExpression is \f$\underline{\mathbb{u}}\f$, where \f$u_{ij}\f$ representss \f$j^{th}\f$ component of the flattened (having \f$N\f$ total components) entity data for \f$i^{th}\f$ entity (having \f$M\f$ entities)
     *          , Following illustrates the returned value where the entity data is flattened.
     *
     *          \f[
     *              NormL2(\underline{\mathbb{u}}) = \sqrt {\sum_{i\in \left[0, M\right)} \sum_{j\in \left[0, N\right)} u_{ij}^2}
     *          \f]
     *
     *          This method is compatible with shared and distributed memory parallelized runs.
     *
     * @param rpExpression          Expressions to be summed.
     * @param rDataCommunicator     Data communicator for MPI communication.
     * @return double               The L2 norm of all the components in all the entities.
     */
    static double NormL2(
        const Expression::ConstPointer& rpExpression,
        const DataCommunicator& rDataCommunicator);

    /**
     * @brief Returns the P norm of the expression assuming it is a flat vector [Shape is not considered].
     *
     * @details If the input @a rpExpression is \f$\underline{\mathbb{u}}\f$, where \f$u_{ij}\f$ representss \f$j^{th}\f$ component of the flattened (having \f$N\f$ total components) entity data for \f$i^{th}\f$ entity (having \f$M\f$ entities)
     *          , Following illustrates the returned value where the entity data is flattened.
     *
     *          \f[
     *              NormP(\underline{\mathbb{u}}, P) = \left(\sum_{i\in \left[0, M\right)} \sum_{j\in \left[0, N\right)} \left|u_{ij}\right|^P\right)^{1/P}
     *          \f]
     *
     *          This method is compatible with shared and distributed memory parallelized runs.
     *
     * @param rpExpression          Expressions to be summed.
     * @param rDataCommunicator     Data communicator for MPI communication.
     * @param P                     P norm coefficent.
     * @return double               The P norm of all the components in all the entities.
     */
    static double NormP(
        const Expression::ConstPointer& rpExpression,
        const double P,
        const DataCommunicator& rDataCommunicator);

    /**
     * @brief Returns the inner product between two expressions.
     *
     * @details If the input @a rpExpression1 is \f$\underline{\mathbb{u}}\f$, where \f$u_{ij}\f$ representss \f$j^{th}\f$ component of the flattened (having \f$N\f$ total components) entity data for \f$i^{th}\f$ entity (having \f$M\f$ entities)
     *          and the input @a rpExpression2 is \f$\underline{\mathbb{v}}\f$, where \f$v_{ij}\f$ representss \f$j^{th}\f$ component of the flattened (having \f$N\f$ total components) entity data for \f$i^{th}\f$ entity (having \f$M\f$ entities),
     *          Following illustrates the returned value where the entity data is flattened. This does not consider shapes of the expressions. They should have same size of flattened vectors.
     *
     *          \f[
     *              InnerProduct(\underline{\mathbb{u}}, \underline{\mathbb{v}}) = \sum_{i\in \left[0, M\right)} \sum_{j\in \left[0, N\right)} u_{ij} \times v_{ij}
     *          \f]
     *
     *          This method is compatible with shared and distributed memory parallelized runs.
     *
     * @throws If the flattend size mismatch.
     * @throws If the number of entities mismatch.
     *
     * @param rpExpression1         Expression1.
     * @param rpExpression1         Expression2.
     * @param rDataCommunicator     Data communicator for MPI communication.
     * @return double               The inner product of @a rpExpression1 and @a rpExpression2
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
        /**                                                                         \
        @brief Returns a container expression having METHOD_NAME evaluated on the given @a rContainerExpression's expression. \
        @see METHOD_NAME.                                           \
        @tparam TContainerType Container type.                                     \
        @param ContainerExpression<TContainerType> Container expression to apply METHOD_NAME.                    \
        @return ContainerExpression<TContainerType> Resulting container expression. **/                          \
        template <class TContainerType>                                             \
        static ContainerExpression<TContainerType> METHOD_NAME(                     \
            const ContainerExpression<TContainerType>& rContainerExpression)        \
        {                                                                           \
            auto copy = rContainerExpression;                                       \
            copy.SetExpression(METHOD_NAME(rContainerExpression.pGetExpression())); \
            return copy;                                                            \
        }
    #endif

    #ifndef KRATOS_EXPRESSION_UTILS_CEXP_METHOD_2
    #define KRATOS_EXPRESSION_UTILS_CEXP_METHOD_2(METHOD_NAME)                                     \
        /**                                                                                        \
        @brief Returns a reduced value by evaluating METHOD_NAME on the given @a rContainerExpression's expression.                                                                         \
        @see METHOD_NAME.                                                          \
        @tparam TContainerType Container type.                                                    \
        @param ContainerExpression<TContainerType> Container expression to apply METHOD_NAME.                                                                                  \
        @return double Resulting scalar value after evaluating METHOD_NAME on @a rContainerExpression. **/                                                                     \
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
        /**                                                                                      \
        @brief Returns a container expression by evaluating METHOD_NAME on the given @a rContainerExpression's expression using @a Value                                                      \
        @see METHOD_NAME.                                                        \
        @tparam TContainerType Container type.                                                  \
        @param ContainerExpression<TContainerType> Container expression to apply METHOD_NAME.                                                                                \
        @return ContainerExpression<TContainerType> Resulting container expression. **/                                                                                         \
        template <class TContainerType>                                                          \
        static ContainerExpression<TContainerType> METHOD_NAME(                                  \
            const ContainerExpression<TContainerType>& rContainerExpression, const double Value) \
        {                                                                                        \
            auto copy = rContainerExpression;                                                    \
            copy.SetExpression(METHOD_NAME(rContainerExpression.pGetExpression(), Value));       \
            return copy;                                                                         \
        }                                                                                        \
        /**                                                                                      \
         @brief Returns a container expression by evaluating METHOD_NAME on the given @a rContainerExpression's expression using @a rContainerExpression2's expression                                      \
        @see METHOD_NAME.                                                        \
        @tparam TContainerType Container type.                                                  \
        @param ContainerExpression<TContainerType> Container expression1 to apply METHOD_NAME.                                                                                \
        @param ContainerExpression<TContainerType> Container expression1 to apply METHOD_NAME.                                                                                \
        @return ContainerExpression<TContainerType> Resulting container expression. **/                                                                                         \
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

    /**
     * @brief Slicing given container expression's expression.
     * @see Slice.
     * @tparam TContainerType Container type
     * @param rContainerExpression  Container expression to be sliced.
     * @param Offset Index of the first component to begin slicing at.
     * @param Stride Number of components from the offset in the sliced item.
     * @return ContainerExpression<TContainerType>      Resulting sliced container expression.
     */
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

    /**
     * @brief Reshape the data in the given container expression's expression.
     * @see Reshape.
     * @tparam TContainerType Container type
     * @tparam TIterator Iterator type for the shape.
     * @param rContainerExpression  Container expression to be sliced.
     * @param NewShapeBegin Starting iterator for the new shape.
     * @param NewShapeEnd Ending iterator for the new shape
     * @return ContainerExpression<TContainerType> Resulting reshaped container expression.
     */
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

    /**
     * @brief Reshape the data in the given container expression's expression.
     * @see Reshape.
     * @tparam TContainerType Container type
     * @param rContainerExpression  Container expression to be sliced.
     * @param rNewShape New shape.
     * @return ContainerExpression<TContainerType> Resulting reshaped container expression.
     */
    template<class TContainerType>
    static ContainerExpression<TContainerType> Reshape(
        const ContainerExpression<TContainerType>& rContainerExpression,
        const std::vector<IndexType>& rNewShape)
    {
        return Reshape(rContainerExpression, rNewShape.begin(), rNewShape.end());
    }

    /** @brief Append the components of a set of container expressions' expression to the current container expression's expression components.
     * @see Comb
     * @throws If the @a rpContainerExpressions is empty.
     * @tparam TContainerType Type of the data container
     * @param rpContainerExpressions  List of container expressions to comb through.
     * @return ContainerExpression<TContainerType>  Resulting container expression.
     */
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

    /**
     * @brief Computes the P norm of the given container expression's expression.
     * @see NormP
     * @tparam TContainerType   Container type.
     * @param rContainerExpression  Input container expression.
     * @param P Norm coefficient.
     * @return double   Resulting p norm.
     */
    template <class TContainerType>
    static double NormP(
        const ContainerExpression<TContainerType>& rContainerExpression,
        const double P)
    {
        return NormP(rContainerExpression.pGetExpression(), P, rContainerExpression.GetModelPart().GetCommunicator().GetDataCommunicator());
    }

    /**
     * @brief Computes inner product between two container expressions's expressions.
     * @see InnerProduct
     * @tparam TContainerType Container type.
     * @param rContainerExpression1    Container expression 1.
     * @param rContainerExpression2    Container expression 2.
     * @return double                  Inner product.
     */
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
