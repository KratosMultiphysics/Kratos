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
#include "testing/testing.h"
#include "spaces/ublas_space.h"
#include "solving_strategies/strategies/implicit_solving_strategy.h"
#include "linear_solvers/linear_solver.h"
#include "custom_strategies/rom_builder_and_solver.h"

namespace Kratos {
namespace Testing {
namespace ROMBuilderAndSolverTests {

using SparseSpaceType = UblasSpace<double, CompressedMatrix, boost::numeric::ublas::vector<double>>;
using LocalSpaceType = UblasSpace<double, Matrix, Vector>;
using LinearSolverType = LinearSolver<SparseSpaceType, LocalSpaceType >;
using BuilderAndSolverType = BuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType >;
using ROMBuilderAndSolverType = ROMBuilderAndSolver<SparseSpaceType, LocalSpaceType, LinearSolverType>;

KRATOS_TEST_CASE_IN_SUITE(ROMBuilderAndSolver, RomApplicationFastSuite)
{
    KRATOS_CHECK_NEAR(1.0, 1.0, 1e-8);
}



}
}
}