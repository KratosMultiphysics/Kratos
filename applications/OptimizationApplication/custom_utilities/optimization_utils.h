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
#include <algorithm>

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
    ///@name Static operations
    ///@{

    template<class TEntity>
    inline static const array_1d<double, 3> GetEntityPosition(const TEntity& rEntity)
    {
        if constexpr(std::is_same_v<TEntity, ModelPart::NodeType>) {
            return rEntity.Coordinates();
        } else if constexpr(std::is_same_v<TEntity, ModelPart::ConditionType>) {
            return rEntity.GetGeometry().Center();
        } else if constexpr(std::is_same_v<TEntity, ModelPart::ElementType>) {
            return rEntity.GetGeometry().Center();
        } else {
            static_assert(!std::is_same_v<TEntity, TEntity>, "Unsupported entity type.");
            return array_1d<double, 3>(0.0);
        }
    }

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

    ///@}
};

///@}
}