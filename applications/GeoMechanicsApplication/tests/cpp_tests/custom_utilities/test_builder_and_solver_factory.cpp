// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Richard Faasse
//

#include "custom_utilities/builder_and_solver_factory.hpp"
#include "includes/expect.h"
#include "linear_solvers/linear_solver.h"
#include "spaces/ublas_space.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite_without_kernel.h"

using namespace Kratos;
using SparseSpaceType  = UblasSpace<double, CompressedMatrix, Vector>;
using LocalSpaceType   = UblasSpace<double, Matrix, Vector>;
using LinearSolverType = LinearSolver<SparseSpaceType, LocalSpaceType>;
using BuilderAndSolverFactoryType = BuilderAndSolverFactory<SparseSpaceType, LocalSpaceType, LinearSolverType>;

namespace Kratos::Testing
{

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, CreateBuilderAndSolver_ReturnsCorrectType_WhenBlockBuilderIsTrue)
{
    const std::string validBlockParameters = R"(
    {
        "block_builder": true
    }
    )";

    const auto solver = std::make_shared<LinearSolverType>();

    const auto builder_and_solver =
        BuilderAndSolverFactoryType::Create(Parameters{validBlockParameters}, solver);
    const auto block_builder =
        dynamic_cast<const ResidualBasedBlockBuilderAndSolver<SparseSpaceType, LocalSpaceType, LinearSolverType>*>(
            builder_and_solver.get());
    EXPECT_NE(block_builder, nullptr);
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, CreateBuilderAndSolver_ReturnsCorrectType_WhenBlockBuilderIsFalse)
{
    const std::string validBlockParameters = R"(
    {
        "block_builder": false
    }
    )";

    const auto solver = std::make_shared<LinearSolverType>();
    const auto builder_and_solver =
        BuilderAndSolverFactoryType::Create(Parameters{validBlockParameters}, solver);
    const auto elimination_builder =
        dynamic_cast<const ResidualBasedEliminationBuilderAndSolver<SparseSpaceType, LocalSpaceType, LinearSolverType>*>(
            builder_and_solver.get());
    EXPECT_NE(elimination_builder, nullptr);
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, CreateBuilderAndSolver_Throws_WhenBlockBuilderIsUndefined)
{
    const auto solver = std::make_shared<LinearSolverType>();
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        BuilderAndSolverFactoryType::Create(Parameters{}, solver),
        "the block_builder parameter is not defined, aborting BuilderAndSolverCreation")
}

} // namespace Kratos::Testing
