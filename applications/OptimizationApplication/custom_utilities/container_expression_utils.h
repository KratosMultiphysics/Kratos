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
//

#pragma once

// System includes

// Project includes
#include "includes/define.h"
#include "spaces/ublas_space.h"
#include "spatial_containers/spatial_containers.h"
#include "expression/container_expression.h"
#include "expression/expression_utils.h"
#include "collective_expression.h"

// Application includes

namespace Kratos
{

///@name Kratos Classes
///@{

class KRATOS_API(OPTIMIZATION_APPLICATION) ContainerExpressionUtils
{
public:
    ///@name Type definitions
    ///@{

    using IndexType = std::size_t;

    using SparseSpaceType = UblasSpace<double, CompressedMatrix, Vector>;

    using SparseMatrixType = SparseSpaceType::MatrixType;

    ///@}
    ///@name Static operations
    ///@{

    /**
     * @brief Calculate max L2 norm of the evaluated expression for each entitiy.
     *
     * This calculates L2 norm of the each entity expression by evaluating, and then returns
     * the maximum of those. This is also an expensive operation.
     *
     * This method is optimized and compatible with OpenMP and MPI.
     *
     * @tparam TContainerType
     * @param rContainer                Container data
     * @return double                   Max L2 norm
     */
    template<class TContainerType>
    static double EntityMaxNormL2(const ContainerExpression<TContainerType>& rContainer);

    /**
     * @brief Computes inner product between two container expressions by evaluating
     * both expressions for each entity in their containers, hence this is an
     * expensive operation. Both container expressions should have the same types of containe expressions,
     * therefore they may have the same container expressions lists.
     *
     * This method is optimized and compatible with OpenMP and MPI.
     *
     * @param rContainer1               Container expressions 1
     * @param rContainer2               Container expressions 2
     * @return double                   Output of the inner product.
     */
    static double InnerProduct(
        const CollectiveExpression& rContainer1,
        const CollectiveExpression& rContainer2);

    /**
     * @brief Calculates matrix vector product for container variable data
     *
     * This computes matrix and vector product between rMatrix and the contaienr variable data rInput.
     * This is an expensive operation since this involve evaluating the expression for each
     * entity and doing matrix vector multiplication.
     *
     * This is only compatible with OpenMP.
     *
     * @tparam TContainerType
     * @param rOutput                       Output container containing matrix vector multiplication.
     * @param rMatrix                       Sparse matrix
     * @param rInput                        Input container values to be used with matrix product.
     */
    template<class TContainerType>
    static void ProductWithEntityMatrix(
        ContainerExpression<TContainerType>& rOutput,
        const SparseMatrixType& rMatrix,
        const ContainerExpression<TContainerType>& rInput);

    /**
     * @brief Calculates matrix vector product for container variable data
     *
     * This computes matrix and vector product between rMatrix and the contaienr variable data rInput.
     * This is an expensive operation since this involve evaluating the expression for each
     * entity and doing matrix vector multiplication.
     *
     * This is only compatible with OpenMP.
     *
     * @tparam TContainerType
     * @param rOutput                       Output container containing matrix vector multiplication.
     * @param rMatrix                       dense matrix
     * @param rInput                        Input container values to be used with matrix product.
     */
    template<class TContainerType>
    static void ProductWithEntityMatrix(
        ContainerExpression<TContainerType>& rOutput,
        const Matrix& rMatrix,
        const ContainerExpression<TContainerType>& rInput);

    /**
     * @brief Transposes a dense matrix.
     *
     * This method is optimized for OpenMP.
     *
     * @param rOutput
     * @param rInput
     */
    static void Transpose(
        Matrix& rOutput,
        const Matrix& rInput);

    /**
     * @brief Transposes a sparse matrix.
     *
     * This method is optimized for OpenMP.
     *
     * @param rOutput
     * @param rInput
     */
    static void Transpose(
        SparseMatrixType& rOutput,
        const SparseMatrixType& rInput);

    /**
     * @brief Computes number ofneighbour entities for the container
     *
     * This method computes number of entities present around nodes in the
     * given TContainerType in the model parts provided in the rOutput.
     *
     * This method is optimized and compatible with OpenMP and MPI.
     *
     * @tparam TContainerType
     * @param rOutput                       Output containing number of entities around each node.
     */
    template<class TContainerType>
    static void ComputeNumberOfNeighbourEntities(
        ContainerExpression<ModelPart::NodesContainerType>& rOutput);

    /**
     * @brief Maps container data to nodal data
     *
     * This method mapps containr data given in rInput to nodal data (rOutput). Mapping is done using
     * rNeighbourEntities (which should consist of number of neighbour entities surrounding each node.).
     *
     * All the model parts in rOutput, rInput and rNeighbourEntities should be the same.
     *
     * This method is optimized and compatible with OpenMP and MPI.
     *
     * @tparam TContainerType
     * @param rOutput               Output containing mapped nodal values.
     * @param rInput                Input containing values in the TContainerType which needs to be mapped to nodes.
     * @param rNeighbourEntities    Number of neighbour entities present around each node.
     */
    template<class TContainerType>
    static void MapContainerVariableToNodalVariable(
        ContainerExpression<ModelPart::NodesContainerType>& rOutput,
        const ContainerExpression<TContainerType>& rInput,
        const ContainerExpression<ModelPart::NodesContainerType>& rNeighbourEntities);

    /**
     * @brief Maps nodal values to container variable data.
     *
     * This method maps nodal data to given container type.
     *
     * rOutput and rInput should have the same model part.
     *
     * This method is optimized and compatible with OpenMP and MPI.
     *
     * @tparam TContainerType
     * @param rOutput               Output containing mapped data.
     * @param rInput                Nodal data to be mapped to the other type of container.
     */
    template<class TContainerType>
    static void MapNodalVariableToContainerVariable(
        ContainerExpression<TContainerType>& rOutput,
        const ContainerExpression<ModelPart::NodesContainerType>& rInput);

    /**
     * @brief Computes nodal and entity wise matrix multiplication.
     *
     * This method first distributes rNodalValues to nodes. Then it calculates
     * the matrix from the rEntities using rMatrixVariable. Then nodal values of each
     * entitiy in rEntities is computed. Finally the matrix (from rMatrixVariable in each entity)
     * and vector (from nodal values for each entity taken from rNodalValues) multiplication is carried out.
     * The solution of this multiplication is stored in rOutput.
     *
     * This method is optimized and compatible with OpenMP and MPI.
     *
     * @tparam TContainerType
     * @param rOutput                   Output containing the matrix vectormultiplication
     * @param rNodalValues              Nodal values used as the vector.
     * @param rMatrixVariable           The variable used to obtain matrix from each variable.
     * @param rEntities                 Entities to compute the matrix.
     */
    template<class TContainerType>
    static void ComputeNodalVariableProductWithEntityMatrix(
        ContainerExpression<ModelPart::NodesContainerType>& rOutput,
        const ContainerExpression<ModelPart::NodesContainerType>& rNodalValues,
        const Variable<Matrix>& rMatrixVariable,
        TContainerType& rEntities);

    ///@}
    ///@name Static operations derrived from Kratos::ExpressionUtils
    ///@{

    #ifndef KRATOS_OPTAPP_EXPRESSION_UTILS_CEXP_METHOD_1
    #define KRATOS_OPTAPP_EXPRESSION_UTILS_CEXP_METHOD_1(METHOD_NAME)                              \
        static CollectiveExpression METHOD_NAME(const CollectiveExpression& rCollectiveExpression) \
        {                                                                                          \
            KRATOS_TRY                                                                             \
            auto result = rCollectiveExpression;                                                   \
            auto r_list_of_container_expressions = result.GetContainerExpressions();               \
            for (IndexType i = 0; i < r_list_of_container_expressions.size(); ++i) {               \
                std::visit(                                                                        \
                    [](auto& pResult) {                                                            \
                        *pResult = ExpressionUtils::METHOD_NAME(*pResult);                         \
                    },                                                                             \
                    r_list_of_container_expressions[i]);                                           \
            }                                                                                      \
            return result;                                                                         \
            KRATOS_CATCH("");                                                                      \
        }
    #endif

    #ifndef KRATOS_OPTAPP_EXPRESSION_UTILS_CEXP_METHOD_2
    #define KRATOS_OPTAPP_EXPRESSION_UTILS_CEXP_METHOD_2(METHOD_NAME)                \
        static double METHOD_NAME(const CollectiveExpression& rCollectiveExpression) \
        {                                                                            \
            KRATOS_TRY                                                               \
            double value = 0.0;                                                      \
            auto r_list_of_container_expressions =                                   \
                rCollectiveExpression.GetContainerExpressions();                     \
            for (IndexType i = 0; i < r_list_of_container_expressions.size(); ++i) { \
                value += std::visit(                                                 \
                    [](const auto& pResult) {                                        \
                        return ExpressionUtils::METHOD_NAME(*pResult);               \
                    },                                                               \
                    r_list_of_container_expressions[i]);                             \
            }                                                                        \
            return value;                                                            \
            KRATOS_CATCH("");                                                        \
        }
    #endif

    #ifndef KRATOS_OPTAPP_EXPRESSION_UTILS_CEXP_METHOD_3
    #define KRATOS_OPTAPP_EXPRESSION_UTILS_CEXP_METHOD_3(METHOD_NAME)                            \
        static CollectiveExpression METHOD_NAME(                                                 \
            const CollectiveExpression& rCollectiveExpression, const double V)                   \
        {                                                                                        \
            KRATOS_TRY                                                                           \
            auto result = rCollectiveExpression;                                                 \
            auto r_list_of_container_expressions = result.GetContainerExpressions();             \
            for (IndexType i = 0; i < r_list_of_container_expressions.size(); ++i) {             \
                std::visit(                                                                      \
                    [V](auto& pResult) {                                                         \
                        *pResult = ExpressionUtils::METHOD_NAME(*pResult, V);                    \
                    },                                                                           \
                    r_list_of_container_expressions[i]);                                         \
            }                                                                                    \
            return result;                                                                       \
            KRATOS_CATCH("");                                                                    \
        }                                                                                        \
        static CollectiveExpression METHOD_NAME(                                                 \
            const CollectiveExpression& rCollectiveExpression1,                                  \
            const CollectiveExpression& rCollectiveExpression2)                                  \
        {                                                                                        \
            KRATOS_TRY                                                                           \
                                                                                                \
            KRATOS_ERROR_IF_NOT(rCollectiveExpression1.IsCompatibleWith(rCollectiveExpression2)) \
                << "Unsupported collective variable data holders provided for "                  \
                "\""                                                                          \
                << #METHOD_NAME << "\"."                                                         \
                << "\nLeft operand : " << rCollectiveExpression1                                 \
                << "\nRight operand: " << rCollectiveExpression2 << std::endl;                   \
                                                                                                \
            auto result = rCollectiveExpression1;                                                \
            auto r_list_of_container_expressions = result.GetContainerExpressions();             \
            const auto& r_right_container_expressions =                                          \
                rCollectiveExpression2.GetContainerExpressions();                                \
            for (IndexType i = 0; i < r_list_of_container_expressions.size(); ++i) {             \
                std::visit(                                                                      \
                    [&r_right_container_expressions, i](auto& pResult) {                         \
                        auto p_right = std::get<std::decay_t<decltype(pResult)>>(                \
                            r_right_container_expressions[i]);                                   \
                        *pResult = ExpressionUtils::METHOD_NAME(*pResult, *p_right);             \
                    },                                                                           \
                    r_list_of_container_expressions[i]);                                         \
            }                                                                                    \
            return result;                                                                       \
                                                                                                \
            KRATOS_CATCH("");                                                                    \
        }
    #endif

    KRATOS_OPTAPP_EXPRESSION_UTILS_CEXP_METHOD_1(Collapse)
    KRATOS_OPTAPP_EXPRESSION_UTILS_CEXP_METHOD_1(Abs)
    KRATOS_OPTAPP_EXPRESSION_UTILS_CEXP_METHOD_1(EntityMin)
    KRATOS_OPTAPP_EXPRESSION_UTILS_CEXP_METHOD_1(EntityMax)
    KRATOS_OPTAPP_EXPRESSION_UTILS_CEXP_METHOD_1(EntitySum)
    KRATOS_OPTAPP_EXPRESSION_UTILS_CEXP_METHOD_2(Sum)
    KRATOS_OPTAPP_EXPRESSION_UTILS_CEXP_METHOD_2(NormInf)
    KRATOS_OPTAPP_EXPRESSION_UTILS_CEXP_METHOD_2(NormL2)
    KRATOS_OPTAPP_EXPRESSION_UTILS_CEXP_METHOD_3(Pow)
    KRATOS_OPTAPP_EXPRESSION_UTILS_CEXP_METHOD_3(Scale)

    #undef KRATOS_OPTAPP_EXPRESSION_UTILS_CEXP_METHOD_1
    #undef KRATOS_OPTAPP_EXPRESSION_UTILS_CEXP_METHOD_2
    #undef KRATOS_OPTAPP_EXPRESSION_UTILS_CEXP_METHOD_3

    ///@}
};

///@}
}