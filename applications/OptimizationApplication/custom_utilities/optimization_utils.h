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
#include <string>
#include <vector>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/data_communicator.h"
#include "includes/kratos_parameters.h"
#include "containers/model.h"
#include "tensor_adaptors/tensor_adaptor.h"

// Application includes

namespace Kratos
{

///@name Kratos Classes
///@{

class KRATOS_API(OPTIMIZATION_APPLICATION) OptimizationUtils
{
public:
    ///@name Type definitions
    ///@{

    using IndexType = std::size_t;

    ///@}
    ///@name Public classes
    ///@{

    template<class TEntityPointType>
    struct KDTreeThreadLocalStorage
    {
        explicit KDTreeThreadLocalStorage(const IndexType MaxNumberOfNeighbors, const IndexType Stride)
        {
            mNeighbourEntityPoints.resize(MaxNumberOfNeighbors);
            mResultingSquaredDistances.resize(MaxNumberOfNeighbors);
            mListOfWeights.resize(MaxNumberOfNeighbors);
            mListOfDampedWeights.resize(Stride, std::vector<double>(MaxNumberOfNeighbors));
        }

        std::vector<TEntityPointType> mNeighbourEntityPoints;
        std::vector<double> mResultingSquaredDistances;
        std::vector<double> mListOfWeights;
        std::vector<std::vector<double>> mListOfDampedWeights;
    };

    ///@}
    ///@name Static operations
    ///@{

    template<class TContainerType>
    static GeometryData::KratosGeometryType GetContainerEntityGeometryType(
        const TContainerType& rContainer,
        const DataCommunicator& rDataCommunicator);

    template<class TContainerType, class TDataType>
    static bool IsVariableExistsInAllContainerProperties(
        const TContainerType& rContainer,
        const Variable<TDataType>& rVariable,
        const DataCommunicator& rDataCommunicator);

    template<class TContainerType, class TDataType>
    static bool IsVariableExistsInAtLeastOneContainerProperties(
        const TContainerType& rContainer,
        const Variable<TDataType>& rVariable,
        const DataCommunicator& rDataCommunicator);

    template<class TContainerType>
    static void CreateEntitySpecificPropertiesForContainer(
        ModelPart& rModelPart,
        TContainerType& rContainer,
        const bool IsRecursive);

    template<class TContainerType, class TDataType>
    static void UpdatePropertiesVariableWithRootValueRecursively(
        TContainerType& rContainer,
        const Variable<TDataType>& rVariable);

    template<class TDataType>
    static IndexType GetVariableDimension(
        const Variable<TDataType>& rVariable,
        const IndexType DomainSize);

    static void SetSolutionStepVariablesList(
        ModelPart& rDestinationModelPart,
        const ModelPart& rOriginModelPart);

    static bool IsSolutionStepVariablesListASubSet(
        const ModelPart& rMainSetModelPart,
        const ModelPart& rSubSetModelPart);

    static std::vector<std::string> GetSolutionStepVariableNamesList(const ModelPart& rModelPart);

    static std::vector<std::vector<ModelPart*>> GetComponentWiseModelParts(
        Model& rModel,
        Parameters Settings);

    /**
     * @brief Maps container data to nodal data
     *
     * This method returns a mapped container data given in @p rInput to nodal given by @p pNodes.
     *
     * @throws std::runtime_error if the container in @p rInput is not of type conditions or elements.
     * @throws std::runtime_error if the entities in @p rInput 's container contain nodes which are not found in the @p pNodes.
     *
     * @param rInput    Input containing values in which needs to be mapped to nodes (container should be conditions or elements).
     * @param pNodes    List of nodes to which the values will be mapped out.
     * @returns         Returns a tensor adaptor having the container @p pNodes with the mapped values.
     */
    static TensorAdaptor<double>::Pointer MapContainerDataToNodalData(
        const TensorAdaptor<double>& rInput,
        ModelPart::NodesContainerType::Pointer pNodes);

    /**
     * @brief Maps nodal data to a specific container data
     *
     * This method returns a mapped nodal data given in @p rInput to container given by @p pEntities.
     *
     * @throws std::runtime_error if the container in @p rInput is not a nodal container.
     * @throws std::runtime_error if the nodal container in @p rInput is not same as nodal container in @p rNeighbourCount .
     * @throws std::runtime_error if the entities in @p rInput 's container contain nodes which are not found in the nodal container of @p rInput .
     *
     * @param rInput            Input containing values in which needs to be mapped to conditions / elements (container should be nodal).
     * @param pEntities         List of entities to which the values will be mapped out.
     * @param rNeighbourCount   Number of neighbour entities around each node.
     * @returns                 Returns a tensor adaptor having the container @p pEntities with the mapped values.
     */
    template<class TContainerPointerType>
    static TensorAdaptor<double>::Pointer MapNodalDataToContainerData(
        const TensorAdaptor<double>& rInput,
        TContainerPointerType pEntities,
        const TensorAdaptor<int>& rNeighbourCount);

    ///@}
};

///@}
}