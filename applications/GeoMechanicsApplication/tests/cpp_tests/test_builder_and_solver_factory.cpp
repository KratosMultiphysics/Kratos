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

#include "testing/testing.h"
#include "custom_utilities/builder_and_solver_factory.hpp"
#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"

using namespace Kratos;
using SparseSpaceType = UblasSpace<double, CompressedMatrix, Vector>;
using LocalSpaceType = UblasSpace<double, Matrix, Vector>;
using LinearSolverType = LinearSolver<SparseSpaceType, LocalSpaceType>;
using BuilderAndSolverFactoryType = BuilderAndSolverFactory<SparseSpaceType, LocalSpaceType, LinearSolverType>;
namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(CreateBuilderAndSolver_ReturnsCorrectType_WhenBlockBuilderIsTrue, KratosGeoMechanicsFastSuite)
{
    const std::string validBlockParameters = R"(
    {
        "block_builder": true
    }
    )";

    const auto solver = std::make_shared<LinearSolverType>();

    const auto builder_and_solver = BuilderAndSolverFactoryType::Create(Parameters{validBlockParameters}, solver);
    const auto block_builder = dynamic_cast<const ResidualBasedBlockBuilderAndSolver<SparseSpaceType, LocalSpaceType, LinearSolverType>*>(builder_and_solver.get());
    KRATOS_EXPECT_NE(block_builder, nullptr);
}

KRATOS_TEST_CASE_IN_SUITE(CreateBuilderAndSolver_ReturnsCorrectType_WhenBlockBuilderIsFalse, KratosGeoMechanicsFastSuite)
{
    const std::string validBlockParameters = R"(
    {
        "block_builder": false
    }
    )";

    const auto solver = std::make_shared<LinearSolverType>();
    const auto builder_and_solver = BuilderAndSolverFactoryType::Create(Parameters{validBlockParameters}, solver);
    const auto elimination_builder = dynamic_cast<const ResidualBasedEliminationBuilderAndSolver<SparseSpaceType, LocalSpaceType, LinearSolverType>*>(builder_and_solver.get());
    KRATOS_EXPECT_NE(elimination_builder, nullptr);
}

KRATOS_TEST_CASE_IN_SUITE(CreateBuilderAndSolver_Throws_WhenBlockBuilderIsUndefined, KratosGeoMechanicsFastSuite)
{
    const auto solver = std::make_shared<LinearSolverType>();
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
            BuilderAndSolverFactoryType::Create(Parameters{}, solver),
            "the block_builder parameter is not defined, aborting BuilderAndSolverCreation")
}

}
