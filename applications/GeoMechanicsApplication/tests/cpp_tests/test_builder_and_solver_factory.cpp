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

namespace Kratos::Testing
{

template<class TSparseSpaceType, class TDenseSpaceType>
class DummyLinearSolver : public LinearSolver<TSparseSpaceType, TDenseSpaceType>
{};

KRATOS_TEST_CASE_IN_SUITE(CreateBuildAndSolver_ReturnsCorrectType_WhenBlockBuilderIsTrue, KratosGeoMechanicsFastSuite)
{
    const std::string validBlockParameters = R"(
    {
        "block_builder": true
    }
    )";


    auto solver = std::make_shared<DummyLinearSolver<SparseSpaceType, LocalSpaceType>>();

    auto builder_and_solver = BuilderAndSolverFactory<SparseSpaceType, LocalSpaceType, LinearSolver<SparseSpaceType, LocalSpaceType>>::Create(Parameters{validBlockParameters}, solver);
    auto block_builder = dynamic_cast<ResidualBasedBlockBuilderAndSolver<SparseSpaceType, LocalSpaceType, LinearSolver<SparseSpaceType, LocalSpaceType>>*>(builder_and_solver.get());
    KRATOS_EXPECT_NE(block_builder, nullptr);
}





}
