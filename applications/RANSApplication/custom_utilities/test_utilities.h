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
#include <functional>

// External includes

// Project includes
#include "includes/model_part.h"

namespace Kratos
{
namespace RansApplicationTestUtilities
{
using IndexType = std::size_t;

KRATOS_API(RANS_APPLICATION) ModelPart& CreateScalarVariableTestModelPart(
    Model& rModel,
    const std::string& rElementName,
    const std::string& rConditionName,
    const std::function<void(Properties&)>& rSetElementProperties,
    const std::function<void(Properties&)>& rSetConditionProperties,
    const std::function<void(ModelPart& rModelPart)>& rAddNodalSolutionStepVariablesFuncion,
    const Variable<double>& rDofVariable,
    const int BufferSize = 2,
    const bool DoInitializeElements = true,
    const bool DoInitializeConditions = true);

KRATOS_API(RANS_APPLICATION) void CheckElementsAndConditions(
    const ModelPart& rModelPart);

} // namespace RansApplicationTestUtilities
} // namespace Kratos