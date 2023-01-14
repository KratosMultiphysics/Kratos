//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 license: HDF5Application/license.txt
//
//  Main author:     Suneth Warnakulasuriya
//

#pragma once

// System includes
#include <vector>

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

    ///@}
    ///@name Static operations
    ///@{

    template<class TContainerType>
    static void GetContainerIds(
        const TContainerType& rContainer,
        std::vector<IndexType>& rOutput);

    template<class TContainerType, class TDataType>
    static void GetContainerVariableToVector(
        const TContainerType& rContainer,
        const Variable<TDataType>& rVariable,
        const IndexType DomainSize,
        Vector& rOutput);

    template<class TContainerType>
    static GeometryData::KratosGeometryType GetContainerEntityGeometryType(
        const TContainerType& rContainer,
        const DataCommunicator& rDataCommunicator);

    template<class TContainerType>
    static void GetContainerPropertiesVariableToVector(
        const TContainerType& rContainer,
        const Variable<double>& rVariable,
        Vector& rOutput);

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

    template<class TDataType>
    static IndexType GetLocalSize(
        const IndexType DomainSize);

    static double CalculateVectorL2Norm(const Vector& rInput);

    template<class TContainerType>
    static void AssignVectorToContainerProperties(
        TContainerType& rContainer,
        const Variable<double>& rPropertiesVariable,
        const Vector& rValues);

    template<class TContainerType, class TDataType>
    static void AssignVectorToContainer(
        TContainerType& rContainer,
        const IndexType DomainSize,
        const Variable<TDataType>& rVariable,
        const Vector& rValues);

    template<class TContainerType>
    static void CreateEntitySpecificPropertiesForContainer(
        ModelPart& rModelPart,
        TContainerType& rContainer);

    ///@}
private:
    ///@name Private operations
    ///@{

    template<class TDataType>
    static inline void AssignValueToVector(
        const TDataType& rValue,
        const IndexType ValueComponentIndex,
        const IndexType VectoStartingIndex,
        Vector& rOutput);

    template<class TDataType>
    static inline void AssignValueFromVector(
        TDataType& rOutput,
        const IndexType ValueComponentIndex,
        const IndexType VectoStartingIndex,
        const Vector& rInput);

    ///@}
};

///@}
}