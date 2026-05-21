//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

#pragma once

// System includes
#include <vector>

// Project includes
#include "includes/define.h"

// Application includes
#include "collective_expression.h"

namespace Kratos {

class KRATOS_API(OPTIMIZATION_APPLICATION) CollectiveExpressionIO
{
public:
    ///@name Type definitions
    ///@{

    using VariableType = std::variant<
                                const Variable<int>*,
                                const Variable<double>*,
                                const Variable<array_1d<double, 3>>*,
                                const Variable<array_1d<double, 4>>*,
                                const Variable<array_1d<double, 6>>*,
                                const Variable<array_1d<double, 9>>*,
                                const Variable<Vector>*,
                                const Variable<Matrix>*>;

    ///@}
    ///@name Public classes
    ///@{

    class HistoricalVariable
    {
    public:
        ///@name Life cycle
        ///@{

        KRATOS_CLASS_POINTER_DEFINITION(HistoricalVariable);

        HistoricalVariable(const VariableType& rVariable): mVariable(rVariable) {}

        ///@}
        ///@name Public member variables
        ///@{

        const VariableType mVariable;

        ///@}
    };

    class NonHistoricalVariable
    {
    public:
        ///@name Life cycle
        ///@{

        KRATOS_CLASS_POINTER_DEFINITION(NonHistoricalVariable);

        NonHistoricalVariable(const VariableType& rVariable): mVariable(rVariable) {}

        ///@}
        ///@name Public member variables
        ///@{

        const VariableType mVariable;

        ///@}
    };

    class PropertiesVariable
    {
    public:
        ///@name Life cycle
        ///@{

        KRATOS_CLASS_POINTER_DEFINITION(PropertiesVariable);

        PropertiesVariable(const VariableType& rVariable): mVariable(rVariable) {}

        ///@}
        ///@name Public member variables
        ///@{

        const VariableType mVariable;

        ///@}
    };

    using ContainerVariableType = std::variant<
                                        HistoricalVariable::Pointer,
                                        NonHistoricalVariable::Pointer,
                                        PropertiesVariable::Pointer>;

    ///@}
    ///@name Public static operations
    ///@{

    template<class TRawDataType>
    static void Read(
        CollectiveExpression& rCollectiveExpression,
        TRawDataType const* pBegin,
        int const* NumberOfEntities,
        int const** pListShapeBegin,
        int const* ShapeSizes,
        const int NumberOfContainers);

    template<class TRawDataType>
    static void Move(
        CollectiveExpression& rCollectiveExpression,
        TRawDataType* pBegin,
        int const* NumberOfEntities,
        int const** pListShapeBegin,
        int const* ShapeSizes,
        const int NumberOfContainers);

    template<class TRawDataType>
    static void Write(
        const CollectiveExpression& rCollectiveExpression,
        TRawDataType* pBegin,
        const int Size);

    static void Read(
        CollectiveExpression& rCollectiveExpression,
        const ContainerVariableType& rContainerVariable);

    static void Read(
        CollectiveExpression& rCollectiveExpression,
        const std::vector<ContainerVariableType>& rContainerVariables);

    static void Write(
        const CollectiveExpression& rCollectiveExpression,
        const ContainerVariableType& rContainerVariable);

    static void Write(
        const CollectiveExpression& rCollectiveExpression,
        const std::vector<ContainerVariableType>& rContainerVariables);

    ///@}
};

} // namespace Kratos
