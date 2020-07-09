//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//
//

#if !defined(KRATOS_EVM_K_EPSILON_HIGH_RE_TEST_UTILITIES_H_INCLUDED)
#define KRATOS_EVM_K_EPSILON_HIGH_RE_TEST_UTILITIES_H_INCLUDED

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
namespace EvmKEpsilonTestUtilities
{
ModelPart& RansEvmKEpsilonK2D3N_SetUp(Model& rModel, const std::string& rElementName);

ModelPart& RansEvmKEpsilonEpsilon2D3N_SetUp(Model& rModel, const std::string& rElementName);

ModelPart& RansEvmKEpsilonEpsilon2D2N_SetUp(Model& rModel, const std::string& rConditionName);
} // namespace EvmKEpsilonTestUtilities
} // namespace Kratos

#endif // KRATOS_EVM_K_EPSILON_HIGH_RE_TEST_UTILITIES_H_INCLUDED