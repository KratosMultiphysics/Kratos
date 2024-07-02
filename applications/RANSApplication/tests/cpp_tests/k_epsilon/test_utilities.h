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

#if !defined(KRATOS_K_EPSILON_HIGH_RE_TEST_UTILITIES_H_INCLUDED)
#define KRATOS_K_EPSILON_HIGH_RE_TEST_UTILITIES_H_INCLUDED

// System includes
#include <string>

// External includes

// Project includes
#include "containers/model.h"
#include "includes/model_part.h"

// Application includes
#include "custom_utilities/test_utilities.h"

namespace Kratos
{
namespace KEpsilonTestUtilities
{
ModelPart& RansKEpsilonK2D3NSetUp(
    Model& rModel,
    const std::string& rElementName);

ModelPart& RansKEpsilonEpsilon2D3NSetUp(
    Model& rModel,
    const std::string& rElementName);

ModelPart& RansKEpsilonEpsilon2D2NSetUp(
    Model& rModel,
    const std::string& rConditionName);

} // namespace KEpsilonTestUtilities
} // namespace Kratos

#endif // KRATOS_K_EPSILON_HIGH_RE_TEST_UTILITIES_H_INCLUDED