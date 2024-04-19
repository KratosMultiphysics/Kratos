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
#include "includes/model_part.h"
#include "includes/data_communicator.h"
#include "includes/kratos_parameters.h"
#include "containers/model.h"

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
        TContainerType& rContainer);

    template<class TDataType>
    static IndexType GetVariableDimension(
        const Variable<TDataType>& rVariable,
        const IndexType DomainSize);

    static void CopySolutionStepVariablesList(
        ModelPart& rDestinationModelPart,
        const ModelPart& rOriginModelPart);

    static bool IsSolutionStepVariablesListASubSet(
        const ModelPart& rMainSetModelPart,
        const ModelPart& rSubSetModelPart);

    static std::vector<std::vector<ModelPart*>> GetComponentWiseModelParts(
        Model& rModel,
        Parameters Settings);

    ///@}
};

///@}
}