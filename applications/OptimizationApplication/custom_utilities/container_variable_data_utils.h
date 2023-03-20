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
#include "custom_utilities/container_variable_data/container_variable_data.h"

namespace Kratos
{

///@name Kratos Classes
///@{

class KRATOS_API(OPTIMIZATION_APPLICATION) ContainerVariableDataUtils
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
    static void Pow(
        ContainerVariableData<TContainerType>& rOutputContainer,
        const ContainerVariableData<TContainerType>& rInputContainer1,
        const ContainerVariableData<TContainerType>& rInputContainer2);

    template<class TContainerType>
    static void Pow(
        ContainerVariableData<TContainerType>& rOutputContainer,
        const ContainerVariableData<TContainerType>& rInputContainer1,
        const double Value);

    template<class TContainerType>
    static double NormInf(const ContainerVariableData<TContainerType>& rContainer);

    template<class TContainerType>
    static double NormL2(const ContainerVariableData<TContainerType>& rContainer);

    template<class TContainerType>
    static double EntityMaxNormL2(const ContainerVariableData<TContainerType>& rContainer);

    template<class TContainerType>
    static double InnerProduct(
        const ContainerVariableData<TContainerType>& rContainer1,
        const ContainerVariableData<TContainerType>& rContainer2);

    template<class TContainerType>
    static void ProductWithEntityMatrix(
        ContainerVariableData<TContainerType>& rOutput,
        const SparseMatrixType& rMatrix,
        const ContainerVariableData<TContainerType>& rInput);

    template<class TContainerType>
    static void ProductWithEntityMatrix(
        ContainerVariableData<TContainerType>& rOutput,
        const Matrix& rMatrix,
        const ContainerVariableData<TContainerType>& rInput);

    static void Transpose(
        Matrix& rOutput,
        const Matrix& rInput);

    static void Transpose(
        SparseMatrixType& rOutput,
        const SparseMatrixType& rInput);

    template<class TContainerType>
    static void ComputeNumberOfNeighbourEntities(
        ContainerVariableData<ModelPart::NodesContainerType>& rOutput);

    template<class TContainerType>
    static void MapContainerVariableDataToNodalVariableData(
        ContainerVariableData<ModelPart::NodesContainerType>& rOutput,
        const ContainerVariableData<TContainerType>& rInput,
        const ContainerVariableData<ModelPart::NodesContainerType>& rNeighbourEntities);

    template<class TContainerType>
    static void MapNodalVariableDataToContainerVariableData(
        ContainerVariableData<TContainerType>& rOutput,
        const ContainerVariableData<ModelPart::NodesContainerType>& rInput);

    template<class TContainerType>
    static void ComputeVariableDataProductWithEntityMatrix(
        ContainerVariableData<ModelPart::NodesContainerType>& rOutput,
        const ContainerVariableData<ModelPart::NodesContainerType>& rNodalValues,
        const Variable<Matrix>& rMatrixVariable,
        TContainerType& rEntities);

    ///@}
};

///@}
}