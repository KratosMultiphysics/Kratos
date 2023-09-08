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
     * @brief Calculate infinity norm of the evaluated expression for each entitiy.
     *
     * This method calculates the infinity norm of the expression by evaluating
     * the expression for each entity in the container, hence this is an expensive
     * operation.
     *
     * This method is optimized and compatible with OpenMP and MPI.
     *
     * @tparam TContainerType
     * @param rContainer                Container data
     * @return double                   Infinity norm
     */
    template<class TContainerType>
    static double NormInf(const ContainerExpression<TContainerType>& rContainer);

    /**
     * @brief Calculate infinity norm of the evaluated expressions for each entitiy.
     *
     * This method calculates the infinity norm of the expressions by evaluating
     * the expression for each entity in the container. Then the infinity norm of all the
     * collective expressions are computed, hence this is an expensive
     * operation.
     *
     * This method is optimized and compatible with OpenMP and MPI.
     *
     * @param rContainer                Collective expressions.
     * @return double                   Infinity norm.
     */
    static double NormInf(const CollectiveExpression& rContainer);

    /**
     * @brief Calculate L2 norm of the evaluated expression for each entitiy.
     *
     * This method calculates the L2 norm of the expression by evaluating
     * the expression for each entity in the container, hence this is an expensive
     * operation.
     *
     * This method is optimized and compatible with OpenMP and MPI.
     *
     * @tparam TContainerType
     * @param rContainer                Container data
     * @return double                   L2 norm
     */
    template<class TContainerType>
    static double NormL2(const ContainerExpression<TContainerType>& rContainer);

    /**
     * @brief Calculate L2 norm of the evaluated expressions for each entitiy.
     *
     * This method calculates the L2 norm of the expression by evaluating
     * the expression for each entity in the container. Then the L2 norm of all the
     * collective expressions are computed, hence this is an expensive
     * operation.
     *
     * This method is optimized and compatible with OpenMP and MPI.
     *
     * @param rContainer                Collective expressions.
     * @return double                   L2 norm
     */
    static double NormL2(const CollectiveExpression& rContainer);

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
     * expensive operation. Both containers should have the same model part,
     * therefore they will have the same containers.
     *
     * This method is optimized and compatible with OpenMP and MPI.
     *
     * @tparam TContainerType
     * @param rContainer1               Container variable data 1
     * @param rContainer2               Container variable data 2
     * @return double                   Output of the inner product.
     */
    template<class TContainerType>
    static double InnerProduct(
        const ContainerExpression<TContainerType>& rContainer1,
        const ContainerExpression<TContainerType>& rContainer2);

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
};

///@}
}