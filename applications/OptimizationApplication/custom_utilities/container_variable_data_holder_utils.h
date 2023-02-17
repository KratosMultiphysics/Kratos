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

// Application includes
#include "custom_utilities/container_variable_data_holder/container_variable_data_holder_base.h"

namespace Kratos
{

///@name Kratos Classes
///@{

class KRATOS_API(OPTIMIZATION_APPLICATION) ContainerVariableDataHolderUtils
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

    template<class TContainerType>
    static double NormInf(const ContainerVariableDataHolderBase<TContainerType>& rContainer);

    template<class TContainerType>
    static double NormL2(const ContainerVariableDataHolderBase<TContainerType>& rContainer);

    template<class TContainerType>
    static double EntityMaxNormL2(const ContainerVariableDataHolderBase<TContainerType>& rContainer);

    template<class TContainerType>
    static double InnerProduct(
        const ContainerVariableDataHolderBase<TContainerType>& rContainer1,
        const ContainerVariableDataHolderBase<TContainerType>& rContainer2);

    template<class TContainerType>
    static void ProductWithEntityMatrix(
        ContainerVariableDataHolderBase<TContainerType>& rOutput,
        const SparseMatrixType& rMatrix,
        const ContainerVariableDataHolderBase<TContainerType>& rInput);

    template<class TContainerType>
    static void ProductWithEntityMatrix(
        ContainerVariableDataHolderBase<TContainerType>& rOutput,
        const Matrix& rMatrix,
        const ContainerVariableDataHolderBase<TContainerType>& rInput);

    static void Transpose(
        Matrix& rOutput,
        const Matrix& rInput);

    static void Transpose(
        SparseMatrixType& rOutput,
        const SparseMatrixType& rInput);

    template<class TContainerType>
    static void ComputeNumberOfNeighbourEntities(
        ContainerVariableDataHolderBase<ModelPart::NodesContainerType>& rOutput,
        const ContainerVariableDataHolderBase<TContainerType>& rInput);

    template<class TContainerType>
    static void MapContainerVariableDataHolderToNodalVariableDataHolder(
        ContainerVariableDataHolderBase<ModelPart::NodesContainerType>& rOutput,
        const ContainerVariableDataHolderBase<TContainerType>& rInput,
        const ContainerVariableDataHolderBase<ModelPart::NodesContainerType>& rNeighbourEntities);

    template<class TContainerType>
    static void MapNodalVariableDataHolderToContainerVariableDataHolder(
        ContainerVariableDataHolderBase<TContainerType>& rOutput,
        const ContainerVariableDataHolderBase<ModelPart::NodesContainerType>& rInput);

    template<class TContainerType>
    static void ComputeVariableDataHolderProductWithEntityMatrix(
        ContainerVariableDataHolderBase<ModelPart::NodesContainerType>& rOutput,
        const ContainerVariableDataHolderBase<ModelPart::NodesContainerType>& rNodalValues,
        const Variable<Matrix>& rMatrixVariable,
        TContainerType& rEntities);

    ///@}
};

///@}
}