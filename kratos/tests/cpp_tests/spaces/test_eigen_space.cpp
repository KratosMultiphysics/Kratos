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

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "testing/testing.h"
#include "spaces/eigen_space.h"
#include "spaces/ublas_space.h"
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

/// The same tri-diagonal matrix as a ublas CompressedMatrix, for parity checks.
void FillTriDiagonalMatrix(CompressedMatrix& rMatrix, const std::size_t Size)
{
    rMatrix = CompressedMatrix(Size, Size);
    for (std::size_t i = 0; i < Size; ++i) {
        if (i >= 1) {rMatrix.push_back(i, i - 1, -1.123);}
        rMatrix.push_back(i, i, 4.5);
        if (i + 1 < Size) {rMatrix.push_back(i, i + 1, 2.336);}
    }
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
// Cross-backend parity: run the same operation through TUblasSparseSpace and
// TEigenSparseSpace on identical data and compare the results.
// -------------------------------------------------------------------------

KRATOS_TEST_CASE_IN_SUITE(EigenSpaceUblasParityVectorOps, KratosCoreFastSuite)
{
    using EigenSpaceType = TEigenSparseSpace<double>;
    using UblasSpaceType = TUblasSparseSpace<double>;

    const std::size_t size = 7;
    EigenSpaceType::VectorType ex(size), ey(size), ez(size);
    UblasSpaceType::VectorType ux(size), uy(size), uz(size);
    for (std::size_t i = 0; i < size; ++i) {
        ex[i] = ux[i] = 0.5 * static_cast<double>(i) - 1.0;
        ey[i] = uy[i] = 2.0 - static_cast<double>(i % 3);
    }

    // Dot / TwoNorm
    KRATOS_EXPECT_NEAR(EigenSpaceType::Dot(ex, ey), UblasSpaceType::Dot(ux, uy), 1e-12);
    KRATOS_EXPECT_NEAR(EigenSpaceType::TwoNorm(ex), UblasSpaceType::TwoNorm(ux), 1e-12);

    // ScaleAndAdd (both overloads)
    EigenSpaceType::ScaleAndAdd(1.5, ex, -0.5, ey, ez);
    UblasSpaceType::ScaleAndAdd(1.5, ux, -0.5, uy, uz);
    for (std::size_t i = 0; i < size; ++i) {
        KRATOS_EXPECT_NEAR(ez[i], uz[i], 1e-12);
    }

    EigenSpaceType::ScaleAndAdd(0.25, ex, 3.0, ez);
    UblasSpaceType::ScaleAndAdd(0.25, ux, 3.0, uz);
    for (std::size_t i = 0; i < size; ++i) {
        KRATOS_EXPECT_NEAR(ez[i], uz[i], 1e-12);
    }

    // InplaceMult / UnaliasedAdd / Assign / Set / SetToZero
    EigenSpaceType::InplaceMult(ez, -2.0);
    UblasSpaceType::InplaceMult(uz, -2.0);
    EigenSpaceType::UnaliasedAdd(ez, 0.75, ex);
    UblasSpaceType::UnaliasedAdd(uz, 0.75, ux);
    EigenSpaceType::Assign(ey, -1.0, ez);
    UblasSpaceType::Assign(uy, -1.0, uz);
    for (std::size_t i = 0; i < size; ++i) {
        KRATOS_EXPECT_NEAR(ey[i], uy[i], 1e-12);
    }

    EigenSpaceType::Set(ey, 3.14);
    UblasSpaceType::Set(uy, 3.14);
    KRATOS_EXPECT_NEAR(ey[size-1], uy[size-1], 1e-12);

    EigenSpaceType::SetToZero(ey);
    UblasSpaceType::SetToZero(uy);
    KRATOS_EXPECT_NEAR(ey[0], uy[0], 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(EigenSpaceUblasParitySpMV, KratosCoreFastSuite)
{
    using EigenSpaceType = TEigenSparseSpace<double>;
    using UblasSpaceType = TUblasSparseSpace<double>;

    const std::size_t size = 10;
    EigenSpaceType::MatrixType e_mat;
    UblasSpaceType::MatrixType u_mat;
    FillTriDiagonalMatrix(e_mat, size);
    FillTriDiagonalMatrix(u_mat, size);

    EigenSpaceType::VectorType ex(size), ey(size);
    UblasSpaceType::VectorType ux(size), uy(size);
    for (std::size_t i = 0; i < size; ++i) {
        ex[i] = ux[i] = 1.0 + 0.1 * static_cast<double>(i);
    }

    // Mult
    EigenSpaceType::Mult(e_mat, ex, ey);
    UblasSpaceType::Mult(u_mat, ux, uy);
    for (std::size_t i = 0; i < size; ++i) {
        KRATOS_EXPECT_NEAR(ey[i], uy[i], 1e-12);
    }

    // TransposeMult
    EigenSpaceType::TransposeMult(e_mat, ex, ey);
    UblasSpaceType::TransposeMult(u_mat, ux, uy);
    for (std::size_t i = 0; i < size; ++i) {
        KRATOS_EXPECT_NEAR(ey[i], uy[i], 1e-12);
    }

    // Norms and diagonal queries
    KRATOS_EXPECT_NEAR(EigenSpaceType::TwoNorm(e_mat), UblasSpaceType::TwoNorm(u_mat), 1e-12);
    KRATOS_EXPECT_NEAR(EigenSpaceType::JacobiNorm(e_mat), UblasSpaceType::JacobiNorm(u_mat), 1e-12);
    KRATOS_EXPECT_NEAR(EigenSpaceType::GetDiagonalNorm(e_mat), UblasSpaceType::GetDiagonalNorm(u_mat), 1e-12);
    KRATOS_EXPECT_NEAR(EigenSpaceType::GetMaxDiagonal(e_mat), UblasSpaceType::GetMaxDiagonal(u_mat), 1e-12);
    KRATOS_EXPECT_NEAR(EigenSpaceType::GetMinDiagonal(e_mat), UblasSpaceType::GetMinDiagonal(u_mat), 1e-12);

    // SetToZero keeps the graph but zeroes the values
    EigenSpaceType::SetToZero(e_mat);
    UblasSpaceType::SetToZero(u_mat);
    KRATOS_EXPECT_EQ(e_mat.nnz(), static_cast<std::size_t>(u_mat.nnz()));
    EigenSpaceType::Mult(e_mat, ex, ey);
    UblasSpaceType::Mult(u_mat, ux, uy);
    for (std::size_t i = 0; i < size; ++i) {
        KRATOS_EXPECT_NEAR(ey[i], uy[i], 1e-12);
        KRATOS_EXPECT_NEAR(ey[i], 0.0, 1e-12);
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
