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
        std::vector<IndexType>& rOutput,
        const TContainerType& rContainer);

    template<class TContainerType, class TDataType>
    static void GetContainerVariableToVector(
        Vector& rOutput,
        const TContainerType& rContainer,
        const Variable<TDataType>& rVariable,
        const IndexType DomainSize);

    template<class TDataType>
    static void GetHistoricalContainerVariableToVector(
        Vector& rOutput,
        const ModelPart& rModelPart,
        const Variable<TDataType>& rVariable,
        const IndexType DomainSize);

    template<class TContainerType>
    static GeometryData::KratosGeometryType GetContainerEntityGeometryType(
        const TContainerType& rContainer,
        const DataCommunicator& rDataCommunicator);

    template<class TContainerType, class TDataType>
    static void GetContainerPropertiesVariableToVector(
        Vector& rOutput,
        const TContainerType& rContainer,
        const Variable<TDataType>& rVariable,
        const IndexType DomainSize);

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

    template<class TContainerType, class TDataType>
    static void AssignVectorToContainerPropertiesVariable(
        TContainerType& rContainer,
        const Variable<TDataType>& rPropertiesVariable,
        const IndexType DomainSize,
        const Vector& rValues);

    template<class TContainerType, class TDataType>
    static void AssignVectorToContainerVariable(
        TContainerType& rContainer,
        const Variable<TDataType>& rVariable,
        const IndexType DomainSize,
        const Vector& rValues);

    template<class TDataType>
    static void AssignVectorToContainerHistoricalVariable(
        ModelPart& rModelPart,
        const Variable<TDataType>& rVariable,
        const IndexType DomainSize,
        const Vector& rValues);

    template<class TContainerType>
    static void CreateEntitySpecificPropertiesForContainer(
        ModelPart& rModelPart,
        TContainerType& rContainer);

    static void AddVectors(
        Vector& rOutput,
        const Vector& rA,
        const Vector& rB);

    static void SubstractVectors(
        Vector& rOutput,
        const Vector& rA,
        const Vector& rB);

    static void MultiplyVector(
        Vector& rOutput,
        const Vector& rA,
        const double Multiplier);

    static void DivideVector(
        Vector& rOutput,
        const Vector& rA,
        const double Divisor);

    static double NormInf(
        const Vector& rOutput);

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