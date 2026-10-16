//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//

// System includes
#include <array>
#include <vector>

// External includes
#include <amgcl/relaxation/ilu0_chow_patel.hpp>
#include <amgcl/value_type/static_matrix.hpp>

// Project includes
#include "testing/testing.h"

namespace Kratos::Testing {
namespace {

using BlockType = amgcl::static_matrix<double, 2, 2>;
using RhsType = amgcl::static_matrix<double, 2, 1>;
using BackendType = amgcl::backend::builtin<BlockType>;
using MatrixType = typename BackendType::matrix;
using RelaxationType = amgcl::relaxation::ilu0_chow_patel<BackendType>;

BlockType MakeBlock(
    const double A00,
    const double A01,
    const double A10,
    const double A11)
{
    BlockType result;
    result(0, 0) = A00;
    result(0, 1) = A01;
    result(1, 0) = A10;
    result(1, 1) = A11;
    return result;
}

RhsType MakeRhs(const double Value0, const double Value1)
{
    RhsType result;
    result(0) = Value0;
    result(1) = Value1;
    return result;
}

void TestBlockFactorization(const bool SymmetricScaling)
{
    const std::array<ptrdiff_t, 3> row_ptr{{0, 2, 4}};
    const std::array<ptrdiff_t, 4> columns{{0, 1, 0, 1}};
    const std::array<BlockType, 4> values{{
        MakeBlock(4.0, 1.0, 2.0, 3.0),
        MakeBlock(1.0, 2.0, 0.0, 1.0),
        MakeBlock(2.0, 0.0, 1.0, 1.0),
        MakeBlock(5.0, 1.0, 1.0, 4.0)}};
    const MatrixType matrix(2, 2, row_ptr, columns, values);

    const std::vector<RhsType> expected_solution{
        MakeRhs(1.0, -2.0),
        MakeRhs(0.5, 3.0)};
    const std::vector<RhsType> rhs{
        values[0] * expected_solution[0] + values[1] * expected_solution[1],
        values[2] * expected_solution[0] + values[3] * expected_solution[1]};

    RelaxationType::params parameters;
    parameters.sweeps = 1;
    parameters.omega = 1.0;
    parameters.symmetric_scaling = SymmetricScaling;
    parameters.solve.serial = true;

    for (int repetition = 0; repetition < 3; ++repetition) {
        RelaxationType relaxation(matrix, parameters, BackendType::params());
        std::vector<RhsType> solution(2);
        relaxation.apply(matrix, rhs, solution);

        for (std::size_t i = 0; i < solution.size(); ++i) {
            KRATOS_EXPECT_NEAR(solution[i](0), expected_solution[i](0), 1.0e-12);
            KRATOS_EXPECT_NEAR(solution[i](1), expected_solution[i](1), 1.0e-12);
        }
    }
}

} // namespace

KRATOS_TEST_CASE_IN_SUITE(AmgclIlu0ChowPatelBlockRowScaling, KratosCoreFastSuite)
{
    TestBlockFactorization(false);
}

KRATOS_TEST_CASE_IN_SUITE(AmgclIlu0ChowPatelBlockSymmetricScaling, KratosCoreFastSuite)
{
    TestBlockFactorization(true);
}

KRATOS_TEST_CASE_IN_SUITE(AmgclIlu0ChowPatelPreservesSmallMatrixScale, KratosCoreFastSuite)
{
    using ScalarBackendType = amgcl::backend::builtin<double>;
    using ScalarRelaxationType = amgcl::relaxation::ilu0_chow_patel<ScalarBackendType>;

    const std::array<ptrdiff_t, 3> row_ptr{{0, 2, 4}};
    const std::array<ptrdiff_t, 4> columns{{0, 1, 0, 1}};
    const std::array<double, 4> values{{4.0e-12, 1.0e-12, 2.0e-12, 3.0e-12}};
    const typename ScalarBackendType::matrix matrix(2, 2, row_ptr, columns, values);
    const std::vector<double> expected_solution{2.0, -1.0};
    const std::vector<double> rhs{
        values[0] * expected_solution[0] + values[1] * expected_solution[1],
        values[2] * expected_solution[0] + values[3] * expected_solution[1]};

    ScalarRelaxationType::params parameters;
    parameters.sweeps = 1;
    parameters.omega = 1.0;
    parameters.solve.serial = true;

    ScalarRelaxationType relaxation(matrix, parameters, ScalarBackendType::params());
    std::vector<double> solution(2);
    relaxation.apply(matrix, rhs, solution);

    KRATOS_EXPECT_NEAR(solution[0], expected_solution[0], 1.0e-12);
    KRATOS_EXPECT_NEAR(solution[1], expected_solution[1], 1.0e-12);
}

} // namespace Kratos::Testing
