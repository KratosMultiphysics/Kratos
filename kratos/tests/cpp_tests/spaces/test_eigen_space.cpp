//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// The Eigen-backed types only exist under the Eigen backend
#ifdef KRATOS_USE_EIGEN_BACKEND

// System includes
#include <cmath>
#include <vector>

// External includes

// Project includes
#include "includes/define.h"
#include "testing/testing.h"
#include "spaces/eigen_space.h"
#include "containers/model.h"

namespace Kratos::Testing {

namespace {

/// Fills an Eigen CSR matrix with the same tri-diagonal pattern used by the
/// UblasSpace tests: diagonal 4.5, sub-diagonal -1.123, super-diagonal 2.336.
void FillTriDiagonalMatrix(EigenCompressedMatrix<double>& rMatrix, const std::size_t Size)
{
    const std::size_t nnz = 3 * Size - 2;
    rMatrix = EigenCompressedMatrix<double>(Size, Size, nnz);

    auto row_indices = rMatrix.index1_data();
    auto col_indices = rMatrix.index2_data();
    auto values = rMatrix.value_data();

    std::size_t counter = 0;
    row_indices[0] = 0;
    for (std::size_t i = 0; i < Size; ++i) {
        if (i >= 1) {
            col_indices[counter] = i - 1;
            values[counter] = -1.123;
            ++counter;
        }
        col_indices[counter] = i;
        values[counter] = 4.5;
        ++counter;
        if (i + 1 < Size) {
            col_indices[counter] = i + 1;
            values[counter] = 2.336;
            ++counter;
        }
        row_indices[i + 1] = counter;
    }
    rMatrix.set_filled(Size + 1, nnz);
}

} // namespace

KRATOS_TEST_CASE_IN_SUITE(EigenSpaceNormSparseMatrix, KratosCoreFastSuite)
{
    using SparseSpaceType = TEigenSparseSpace<double>;

    const std::size_t size = 10;
    SparseSpaceType::MatrixType mat;
    FillTriDiagonalMatrix(mat, size);

    // Same reference values as the UblasSpace test
    KRATOS_EXPECT_NEAR(16.216110045260546, SparseSpaceType::TwoNorm(mat), 1e-12);
    KRATOS_EXPECT_NEAR(31.131, SparseSpaceType::JacobiNorm(mat), 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(EigenSpaceNormDenseMatrix, KratosCoreFastSuite)
{
    using LocalSpaceType = TEigenDenseSpace<double>;

    const std::size_t size = 10;
    LocalSpaceType::MatrixType mat(size, size);
    mat.setZero();

    for (std::size_t i = 0; i < mat.size1(); ++i) {
        mat(i, i) = 4.5;
        if (i >= 1) {mat(i, i - 1) = -1.123;}
        if (i + 1 < mat.size2()) {mat(i, i + 1) = 2.336;}
    }

    // Same reference values as the UblasSpace test
    KRATOS_EXPECT_NEAR(16.216110045260546, LocalSpaceType::TwoNorm(mat), 1e-12);
    KRATOS_EXPECT_NEAR(31.131, LocalSpaceType::JacobiNorm(mat), 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(EigenSpaceCheckAndCorrectZeroDiagonalValues, KratosCoreFastSuite)
{
    using SparseSpaceType = TEigenSparseSpace<double>;
    using SparseMatrixType = typename SparseSpaceType::MatrixType;

    Model current_model;
    ModelPart& r_model_part = current_model.CreateModelPart("Main");
    auto& r_process_info = r_model_part.GetProcessInfo();
    r_process_info.SetValue(BUILD_SCALE_FACTOR, 1.0);

    // Diagonal matrix with entries 0, 1, 2, ..., 11 (first diagonal is zero)
    SparseMatrixType matrix12x12(12, 12, 12);
    auto row_indices = matrix12x12.index1_data();
    auto col_indices = matrix12x12.index2_data();
    auto values = matrix12x12.value_data();
    for (std::size_t i = 0; i < 12; ++i) {
        row_indices[i] = i;
        col_indices[i] = i;
        values[i] = static_cast<double>(i);
    }
    row_indices[12] = 12;
    matrix12x12.set_filled(13, 12);

    SparseSpaceType::VectorType vector12(12);
    vector12.setZero();

    const double norm = SparseSpaceType::CheckAndCorrectZeroDiagonalValues(r_process_info, matrix12x12, vector12, SCALING_DIAGONAL::NO_SCALING);
    KRATOS_EXPECT_DOUBLE_EQ(norm, 1.0);
    KRATOS_EXPECT_DOUBLE_EQ(matrix12x12(0, 0), 1.0);
}

KRATOS_TEST_CASE_IN_SUITE(EigenSpaceGetScaleNorm, KratosCoreFastSuite)
{
    using SparseSpaceType = TEigenSparseSpace<double>;
    using SparseMatrixType = typename SparseSpaceType::MatrixType;

    Model current_model;
    ModelPart& r_model_part = current_model.CreateModelPart("Main");
    auto& r_process_info = r_model_part.GetProcessInfo();
    r_process_info.SetValue(BUILD_SCALE_FACTOR, 3.0);

    // Diagonal matrix with entries 1, 2, ..., 12
    SparseMatrixType matrix12x12(12, 12, 12);
    auto row_indices = matrix12x12.index1_data();
    auto col_indices = matrix12x12.index2_data();
    auto values = matrix12x12.value_data();
    for (std::size_t i = 0; i < 12; ++i) {
        row_indices[i] = i;
        col_indices[i] = i;
        values[i] = static_cast<double>(i + 1);
    }
    row_indices[12] = 12;
    matrix12x12.set_filled(13, 12);

    // Same reference values as the UblasSpace test
    double norm = SparseSpaceType::GetScaleNorm(r_process_info, matrix12x12, SCALING_DIAGONAL::NO_SCALING);
    KRATOS_EXPECT_DOUBLE_EQ(norm, 1.0);
    norm = SparseSpaceType::GetScaleNorm(r_process_info, matrix12x12, SCALING_DIAGONAL::CONSIDER_PRESCRIBED_DIAGONAL);
    KRATOS_EXPECT_DOUBLE_EQ(norm, 3.0);
    norm = SparseSpaceType::GetScaleNorm(r_process_info, matrix12x12, SCALING_DIAGONAL::CONSIDER_NORM_DIAGONAL);
    KRATOS_EXPECT_NEAR(norm, 2.124591464, 1.0e-6);
    norm = SparseSpaceType::GetScaleNorm(r_process_info, matrix12x12, SCALING_DIAGONAL::CONSIDER_MAX_DIAGONAL);
    KRATOS_EXPECT_DOUBLE_EQ(norm, 12.0);
    norm = SparseSpaceType::GetAveragevalueDiagonal(matrix12x12);
    KRATOS_EXPECT_DOUBLE_EQ(norm, 6.5);
    norm = SparseSpaceType::GetMinDiagonal(matrix12x12);
    KRATOS_EXPECT_DOUBLE_EQ(norm, 1.0);
}

// -------------------------------------------------------------------------
// Reference checks: run the space operations on known data and compare the
// results against straightforward hand-written computations.
// -------------------------------------------------------------------------

KRATOS_TEST_CASE_IN_SUITE(EigenSpaceVectorOpsReference, KratosCoreFastSuite)
{
    using SpaceType = TEigenSparseSpace<double>;

    const std::size_t size = 7;
    SpaceType::VectorType x(size), y(size), z(size);
    std::vector<double> rx(size), ry(size), rz(size);
    for (std::size_t i = 0; i < size; ++i) {
        x[i] = rx[i] = 0.5 * static_cast<double>(i) - 1.0;
        y[i] = ry[i] = 2.0 - static_cast<double>(i % 3);
    }

    // Dot / TwoNorm
    double reference_dot = 0.0, reference_norm = 0.0;
    for (std::size_t i = 0; i < size; ++i) {
        reference_dot += rx[i] * ry[i];
        reference_norm += rx[i] * rx[i];
    }
    KRATOS_EXPECT_NEAR(SpaceType::Dot(x, y), reference_dot, 1e-12);
    KRATOS_EXPECT_NEAR(SpaceType::TwoNorm(x), std::sqrt(reference_norm), 1e-12);

    // ScaleAndAdd (both overloads)
    SpaceType::ScaleAndAdd(1.5, x, -0.5, y, z);
    for (std::size_t i = 0; i < size; ++i) {
        rz[i] = 1.5 * rx[i] - 0.5 * ry[i];
        KRATOS_EXPECT_NEAR(z[i], rz[i], 1e-12);
    }

    SpaceType::ScaleAndAdd(0.25, x, 3.0, z);
    for (std::size_t i = 0; i < size; ++i) {
        rz[i] = 0.25 * rx[i] + 3.0 * rz[i];
        KRATOS_EXPECT_NEAR(z[i], rz[i], 1e-12);
    }

    // InplaceMult / UnaliasedAdd / Assign / Set / SetToZero
    SpaceType::InplaceMult(z, -2.0);
    SpaceType::UnaliasedAdd(z, 0.75, x);
    SpaceType::Assign(y, -1.0, z);
    for (std::size_t i = 0; i < size; ++i) {
        rz[i] = -2.0 * rz[i] + 0.75 * rx[i];
        ry[i] = -rz[i];
        KRATOS_EXPECT_NEAR(y[i], ry[i], 1e-12);
    }

    SpaceType::Set(y, 3.14);
    KRATOS_EXPECT_NEAR(y[size-1], 3.14, 1e-12);

    SpaceType::SetToZero(y);
    KRATOS_EXPECT_NEAR(y[0], 0.0, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(EigenSpaceSpMVReference, KratosCoreFastSuite)
{
    using SpaceType = TEigenSparseSpace<double>;

    const std::size_t size = 10;
    SpaceType::MatrixType mat;
    FillTriDiagonalMatrix(mat, size);
    const double diagonal = 4.5, sub_diagonal = -1.123, super_diagonal = 2.336;
    auto entry = [&](const std::size_t i, const std::size_t j) {
        if (i == j) return diagonal;
        if (j + 1 == i) return sub_diagonal;
        if (i + 1 == j) return super_diagonal;
        return 0.0;
    };

    SpaceType::VectorType x(size), y(size);
    for (std::size_t i = 0; i < size; ++i) {
        x[i] = 1.0 + 0.1 * static_cast<double>(i);
    }

    // Mult
    SpaceType::Mult(mat, x, y);
    for (std::size_t i = 0; i < size; ++i) {
        double reference = 0.0;
        for (std::size_t j = 0; j < size; ++j) reference += entry(i, j) * x[j];
        KRATOS_EXPECT_NEAR(y[i], reference, 1e-12);
    }

    // TransposeMult
    SpaceType::TransposeMult(mat, x, y);
    for (std::size_t i = 0; i < size; ++i) {
        double reference = 0.0;
        for (std::size_t j = 0; j < size; ++j) reference += entry(j, i) * x[j];
        KRATOS_EXPECT_NEAR(y[i], reference, 1e-12);
    }

    // Norms and diagonal queries
    double frobenius = 0.0, jacobi = 0.0;
    for (std::size_t i = 0; i < size; ++i) {
        for (std::size_t j = 0; j < size; ++j) {
            frobenius += entry(i, j) * entry(i, j);
            if (i != j) jacobi += std::abs(entry(i, j));
        }
    }
    KRATOS_EXPECT_NEAR(SpaceType::TwoNorm(mat), std::sqrt(frobenius), 1e-12);
    KRATOS_EXPECT_NEAR(SpaceType::JacobiNorm(mat), jacobi, 1e-12);
    KRATOS_EXPECT_NEAR(SpaceType::GetDiagonalNorm(mat), std::sqrt(static_cast<double>(size)) * diagonal, 1e-12);
    KRATOS_EXPECT_NEAR(SpaceType::GetMaxDiagonal(mat), diagonal, 1e-12);
    KRATOS_EXPECT_NEAR(SpaceType::GetMinDiagonal(mat), diagonal, 1e-12);

    // SetToZero keeps the graph but zeroes the values
    const std::size_t nnz_before = mat.nnz();
    SpaceType::SetToZero(mat);
    KRATOS_EXPECT_EQ(mat.nnz(), nnz_before);
    SpaceType::Mult(mat, x, y);
    for (std::size_t i = 0; i < size; ++i) {
        KRATOS_EXPECT_NEAR(y[i], 0.0, 1e-12);
    }
}

KRATOS_TEST_CASE_IN_SUITE(EigenSpaceElementInsertion, KratosCoreFastSuite)
{
    // Element insertion through operator() leaves the Eigen matrix in
    // uncompressed mode; the space operations must normalize it transparently
    using SparseSpaceType = TEigenSparseSpace<double>;

    SparseSpaceType::MatrixType mat(3, 3);
    mat(0, 0) = 2.0;
    mat(1, 1) = 3.0;
    mat(2, 0) = -1.0;
    mat(2, 2) = 4.0;

    SparseSpaceType::VectorType x(3), y(3);
    for (std::size_t i = 0; i < 3; ++i) {
        x[i] = 1.0;
    }

    SparseSpaceType::Mult(mat, x, y);
    KRATOS_EXPECT_NEAR(y[0], 2.0, 1e-12);
    KRATOS_EXPECT_NEAR(y[1], 3.0, 1e-12);
    KRATOS_EXPECT_NEAR(y[2], 3.0, 1e-12);

    KRATOS_EXPECT_NEAR(SparseSpaceType::TwoNorm(mat), std::sqrt(4.0 + 9.0 + 1.0 + 16.0), 1e-12);
    KRATOS_EXPECT_NEAR(SparseSpaceType::GetMaxDiagonal(mat), 4.0, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(EigenSpaceUblasParityGraph, KratosCoreFastSuite)
{
    using EigenSpaceType = TEigenSparseSpace<double>;

    const std::size_t size = 6;
    EigenSpaceType::MatrixType e_mat;
    FillTriDiagonalMatrix(e_mat, size);

    // Graph queries against the known tri-diagonal pattern
    KRATOS_EXPECT_EQ(EigenSpaceType::GraphDegree(0, e_mat), 2);
    KRATOS_EXPECT_EQ(EigenSpaceType::GraphDegree(1, e_mat), 3);
    KRATOS_EXPECT_EQ(EigenSpaceType::GraphDegree(size - 1, e_mat), 2);

    std::vector<std::size_t> neighbors;
    EigenSpaceType::GraphNeighbors(1, e_mat, neighbors);
    KRATOS_EXPECT_EQ(neighbors.size(), 3);
    KRATOS_EXPECT_EQ(neighbors[0], 0);
    KRATOS_EXPECT_EQ(neighbors[1], 1);
    KRATOS_EXPECT_EQ(neighbors[2], 2);
}

} // namespace Kratos::Testing

#endif // KRATOS_USE_EIGEN_BACKEND
