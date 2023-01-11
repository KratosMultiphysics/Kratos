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

class KRATOS_API(OPTIMIZATION_APPLICATION) OptimizationVariableUtils
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

    template<class TContainerType>
    static void GetContainerVariableToMatrix(
        const TContainerType& rContainer,
        const Variable<array_1d<double, 3>>& rVariable,
        Matrix& rOutput);

    template<class TContainerType>
    static void GetContainerVariableToVector(
        const TContainerType& rContainer,
        const Variable<double>& rVariable,
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

    ///@}
private:
    ///@name Private operations
    ///@{


    ///@}
};

///@}
}