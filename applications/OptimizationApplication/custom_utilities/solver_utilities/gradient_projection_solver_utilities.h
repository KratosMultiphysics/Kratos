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

// Application includes

namespace Kratos
{

///@name Kratos Classes
///@{

class KRATOS_API(OPTIMIZATION_APPLICATION) GradientProjectionSolverUtilities
{
public:
    ///@name Type definitions
    ///@{

    using IndexType = std::size_t;

    ///@}
    ///@name Static operations
    ///@{

    template<class TContainerType, class TDataType>
    static void CalculateProjectedSearchDirectionAndCorrection(
        TContainerType& rContainer,
        const IndexType DomainSize,
        const Variable<TDataType>& rSearchDirectionVariable,
        const Variable<TDataType>& rSearchDirectionCorrectionVariable,
        const Vector& rConstraintValues,
        const Vector& rObjectiveGradient,
        const std::vector<Vector>& rConstraintsGradients);

    template<class TContainerType, class TDataType>
    static void CalculateControlChange(
        TContainerType& rContainer,
        const DataCommunicator& rDataCommunicator,
        const Variable<TDataType>& rSearchDirectionVariable,
        const Variable<TDataType>& rSearchDirectionCorrectionVariable,
        const Variable<TDataType>& rControlChangeVariable,
        const double StepSize,
        const double MaxCorrectionShare);

    ///@}
private:
    ///@name Private operations
    ///@{

    template<class DataType>
    static IndexType inline GetLocalSize(const IndexType DomainSize);

    template<class TDataType>
    static void inline AddValue(
        TDataType& rVariable,
        const IndexType VariableComponentIndex,
        const IndexType DofStartingIndex,
        const Vector& rValues);

    template<class TDataType>
    static double inline CalculateValueNormSquare(
        const TDataType& rValue);

    ///@}
};

///@}
}