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
#include <vector>
#include <variant>
#include <unordered_map>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/data_communicator.h"

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

    using SensitivityFieldVariableTypes = std::variant<const Variable<double>*, const Variable<array_1d<double, 3>>*>;

    using SensitivityModelPartVariablesListMap = std::unordered_map<ModelPart*, std::vector<SensitivityFieldVariableTypes>>;

    using SensitivityVariableModelPartsListMap = std::unordered_map<SensitivityFieldVariableTypes, std::vector<ModelPart*>>;

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

    template<class TContainerType>
    static IndexType GetNumberOfContainerItemsWithFlag(
        const TContainerType& rContainer,
        const DataCommunicator& rDataCommunicator,
        const Flags& rFlag,
        const bool FlagValue = true);

    static void ReverseSensitivityModelPartVariablesListMap(
        SensitivityVariableModelPartsListMap& rOutput,
        const SensitivityModelPartVariablesListMap& rSensitivityModelPartVariableInfo);

    ///@}
};

///@}
}