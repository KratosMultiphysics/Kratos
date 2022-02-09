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

#if !defined(KRATOS_RANS_TEST_UTILITIES_H_INCLUDED)
#define KRATOS_RANS_TEST_UTILITIES_H_INCLUDED

// System includes
#include <functional>

// External includes

// Project includes
#include "includes/model_part.h"

namespace Kratos
{
namespace RansApplicationTestUtilities
{
using IndexType = std::size_t;

template <class TDataType>
KRATOS_API(RANS_APPLICATION) void AssignRandomValues(
    TDataType& rValue,
    const std::string& rSeed,
    const double MinValue = 0.0,
    const double MaxValue = 1.0);

KRATOS_API(RANS_APPLICATION) ModelPart& CreateTestModelPart(
    Model& rModel,
    const std::string& rElementName,
    const std::string& rConditionName,
    const std::function<void(ModelPart& rModelPart)>& rAddNodalSolutionStepVariablesFuncion,
    const std::function<void(ModelPart::NodeType&)>& rAddDofsFunction,
    const std::function<void(Properties&)>& rSetProperties,
    const int BufferSize = 2);

KRATOS_API(RANS_APPLICATION) ModelPart& CreateScalarVariableTestModelPart(
    Model& rModel,
    const std::string& rElementName,
    const std::string& rConditionName,
    const std::function<void(ModelPart& rModelPart)>& rAddNodalSolutionStepVariablesFuncion,
    const std::function<void(Properties&)>& rSetProperties,
    const Variable<double>& rDofVariable,
    const int BufferSize = 2,
    const bool DoInitializeElements = true,
    const bool DoInitializeConditions = true);

template <class TDataType>
KRATOS_API(RANS_APPLICATION) void RandomFillNodalHistoricalVariable(
    ModelPart& rModelPart,
    const Variable<TDataType>& rVariable,
    const double MinValue = 0.0,
    const double MaxValue = 1.0,
    const int Step = 0);

template <class TContainerType, class TDataType>
KRATOS_API(RANS_APPLICATION) void RandomFillContainerVariable(
    ModelPart& rModelPart,
    const Variable<TDataType>& rVariable,
    const double MinValue = 0.0,
    const double MaxValue = 1.0);

template <class TContainerType>
KRATOS_API(RANS_APPLICATION) void TestEquationIdVector(
    ModelPart& rModelPart);

template <class TContainerType>
KRATOS_API(RANS_APPLICATION) void TestGetDofList(
    ModelPart& rModelPart,
    const Variable<double>& rVariable);

} // namespace RansApplicationTestUtilities
} // namespace Kratos

#endif