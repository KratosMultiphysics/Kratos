//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Eduard Gómez
//

// System includes

/* External includes */

/* Project includes */
#include "tests/cpp_tests/strategies/strategies/explicit_strategy_testing_utils.h"
#include "solving_strategies/strategies/explicit_solving_strategy_bfecc.h"

namespace Kratos
{
namespace Testing
{
namespace BFECCTestInternals
{

using namespace ExplicitStrategyTestingUtils;
using ExplicitSolvingStrategyBFECCType = ExplicitSolvingStrategyBFECC<SparseSpaceType, LocalSpaceType>;

} // namespace BFECCTestInternals

KRATOS_TEST_CASE_IN_SUITE(ExplicitSolvingStrategyBFECC, KratosCoreFastSuite)
{
    using namespace BFECCTestInternals;
    RunTest<ExplicitSolvingStrategyBFECCType>(1e-2);
}

KRATOS_TEST_CASE_IN_SUITE(ExplicitSolvingStrategyBFECCConvergence, KratosCoreFastSuite)
{
    using namespace BFECCTestInternals;
    ConvergenceTest<ExplicitSolvingStrategyBFECCType>(2);
}

} // namespace Testing
} // namespace Kratos
