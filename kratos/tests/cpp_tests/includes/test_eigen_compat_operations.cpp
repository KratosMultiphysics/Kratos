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
#include "testing/testing.h"
// NOTE: this translation unit deliberately includes BOTH the uBLAS and the
// Eigen interfaces and calls the same free functions (prod, inner_prod,
// noalias, ...) on both families of types. It is the permanent regression
// test that the Eigen compatibility overloads never collide with the ublas
// overloads injected by the using-directive in ublas_interface.h.
#include "includes/ublas_interface.h"
#include "includes/kratos_eigen_interface.h"

namespace Kratos::Testing {

namespace {

// 2x2 test data used by all cases, filled identically in both backends.
template<class TMatrixType>
void FillMatrix(TMatrixType& rM)
{
    rM(0, 0) = 1.0; rM(0, 1) = 2.0;
    rM(1, 0) = 3.0; rM(1, 1) = 4.0;
}

template<class TVectorType>
void FillVector(TVectorType& rV)
{
    rV[0] = 5.0;
    rV[1] = -2.0;
}

} // namespace

KRATOS_TEST_CASE_IN_SUITE(EigenCompatProd, KratosCoreFastSuite)
{
    Matrix ublas_A(2, 2), ublas_B(2, 2);
    EigenMatrix<double> eigen_A(2, 2), eigen_B(2, 2);
    FillMatrix(ublas_A); FillMatrix(ublas_B);
    FillMatrix(eigen_A); FillMatrix(eigen_B);

    // dense x dense
    const Matrix ublas_C = prod(ublas_A, ublas_B);
    const EigenMatrix<double> eigen_C = prod(eigen_A, eigen_B);
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t j = 0; j < 2; ++j)
            KRATOS_EXPECT_NEAR(eigen_C(i, j), ublas_C(i, j), 1e-12);

    // dense x vector
    Vector ublas_x(2);
    EigenVector<double> eigen_x(2);
    FillVector(ublas_x); FillVector(eigen_x);
    const Vector ublas_y = prod(ublas_A, ublas_x);
    const EigenVector<double> eigen_y = prod(eigen_A, eigen_x);
    KRATOS_EXPECT_NEAR(eigen_y[0], ublas_y[0], 1e-12);
    KRATOS_EXPECT_NEAR(eigen_y[1], ublas_y[1], 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(EigenCompatSparseProd, KratosCoreFastSuite)
{
    // ublas sparse
    CompressedMatrix ublas_A(2, 2);
    ublas_A(0, 0) = 2.0; ublas_A(0, 1) = -1.0; ublas_A(1, 1) = 3.0;
    Vector ublas_x(2);
    FillVector(ublas_x);
    const Vector ublas_y = prod(ublas_A, ublas_x);

    // eigen sparse (same pattern and values)
    EigenCompressedMatrix<double> eigen_A(2, 2, 3);
    auto* row_indices = eigen_A.index1_data().begin();
    auto* col_indices = eigen_A.index2_data().begin();
    auto* values = eigen_A.value_data().begin();
    row_indices[0] = 0; row_indices[1] = 2; row_indices[2] = 3;
    col_indices[0] = 0; col_indices[1] = 1; col_indices[2] = 1;
    values[0] = 2.0; values[1] = -1.0; values[2] = 3.0;
    eigen_A.set_filled(3, 3);

    EigenVector<double> eigen_x(2);
    FillVector(eigen_x);
    const EigenVector<double> eigen_y = prod(eigen_A, eigen_x);

    KRATOS_EXPECT_NEAR(eigen_y[0], ublas_y[0], 1e-12);
    KRATOS_EXPECT_NEAR(eigen_y[1], ublas_y[1], 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(EigenCompatInnerOuterProd, KratosCoreFastSuite)
{
    Vector ublas_x(2), ublas_y(2);
    EigenVector<double> eigen_x(2), eigen_y(2);
    FillVector(ublas_x); FillVector(eigen_x);
    ublas_y[0] = 0.5; ublas_y[1] = 4.0;
    eigen_y[0] = 0.5; eigen_y[1] = 4.0;

    KRATOS_EXPECT_NEAR(inner_prod(eigen_x, eigen_y), inner_prod(ublas_x, ublas_y), 1e-12);

    const Matrix ublas_op = outer_prod(ublas_x, ublas_y);
    const EigenMatrix<double> eigen_op = outer_prod(eigen_x, eigen_y);
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t j = 0; j < 2; ++j)
            KRATOS_EXPECT_NEAR(eigen_op(i, j), ublas_op(i, j), 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(EigenCompatTransAndNorms, KratosCoreFastSuite)
{
    Matrix ublas_A(2, 2);
    EigenMatrix<double> eigen_A(2, 2);
    FillMatrix(ublas_A); FillMatrix(eigen_A);

    const Matrix ublas_At = trans(ublas_A);
    const EigenMatrix<double> eigen_At = trans(eigen_A);
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t j = 0; j < 2; ++j)
            KRATOS_EXPECT_NEAR(eigen_At(i, j), ublas_At(i, j), 1e-12);

    Vector ublas_x(2);
    EigenVector<double> eigen_x(2);
    FillVector(ublas_x); FillVector(eigen_x);
    KRATOS_EXPECT_NEAR(norm_1(eigen_x), norm_1(ublas_x), 1e-12);
    KRATOS_EXPECT_NEAR(norm_2(eigen_x), norm_2(ublas_x), 1e-12);
    KRATOS_EXPECT_NEAR(norm_inf(eigen_x), norm_inf(ublas_x), 1e-12);
    KRATOS_EXPECT_NEAR(sum(eigen_x), sum(ublas_x), 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(EigenCompatNoalias, KratosCoreFastSuite)
{
    Matrix ublas_A(2, 2);
    EigenMatrix<double> eigen_A(2, 2);
    FillMatrix(ublas_A); FillMatrix(eigen_A);
    Vector ublas_x(2), ublas_y(2);
    EigenVector<double> eigen_x(2), eigen_y(2);
    FillVector(ublas_x); FillVector(eigen_x);

    noalias(ublas_y) = prod(ublas_A, ublas_x);
    noalias(eigen_y) = prod(eigen_A, eigen_x);
    KRATOS_EXPECT_NEAR(eigen_y[0], ublas_y[0], 1e-12);
    KRATOS_EXPECT_NEAR(eigen_y[1], ublas_y[1], 1e-12);

    noalias(ublas_y) += prod(ublas_A, ublas_x);
    noalias(eigen_y) += prod(eigen_A, eigen_x);
    KRATOS_EXPECT_NEAR(eigen_y[0], ublas_y[0], 1e-12);
    KRATOS_EXPECT_NEAR(eigen_y[1], ublas_y[1], 1e-12);

    noalias(ublas_y) -= prod(ublas_A, ublas_x);
    noalias(eigen_y) -= prod(eigen_A, eigen_x);
    KRATOS_EXPECT_NEAR(eigen_y[0], ublas_y[0], 1e-12);
    KRATOS_EXPECT_NEAR(eigen_y[1], ublas_y[1], 1e-12);

    // matrix target as well
    Matrix ublas_C(2, 2);
    EigenMatrix<double> eigen_C(2, 2);
    noalias(ublas_C) = prod(ublas_A, ublas_A);
    noalias(eigen_C) = prod(eigen_A, eigen_A);
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t j = 0; j < 2; ++j)
            KRATOS_EXPECT_NEAR(eigen_C(i, j), ublas_C(i, j), 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(EigenCompatRowColumnSubrange, KratosCoreFastSuite)
{
    Matrix ublas_A(2, 2);
    EigenMatrix<double> eigen_A(2, 2);
    FillMatrix(ublas_A); FillMatrix(eigen_A);

    // reads
    KRATOS_EXPECT_NEAR(row(eigen_A, 1)[0], row(ublas_A, 1)[0], 1e-12);
    KRATOS_EXPECT_NEAR(column(eigen_A, 1)[0], column(ublas_A, 1)(0), 1e-12);

    // write-through
    row(ublas_A, 0) *= 2.0;
    row(eigen_A, 0) *= 2.0;
    KRATOS_EXPECT_NEAR(eigen_A(0, 1), ublas_A(0, 1), 1e-12);

    Vector ublas_x(4);
    EigenVector<double> eigen_x(4);
    for (std::size_t i = 0; i < 4; ++i) {
        ublas_x[i] = static_cast<double>(i) + 1.0;
        eigen_x[i] = static_cast<double>(i) + 1.0;
    }
    KRATOS_EXPECT_NEAR(subrange(eigen_x, 1, 3)[1], subrange(ublas_x, 1, 3)[1], 1e-12);

    subrange(ublas_x, 1, 3) *= 3.0;
    subrange(eigen_x, 1, 3) *= 3.0;
    KRATOS_EXPECT_NEAR(eigen_x[2], ublas_x[2], 1e-12);
    KRATOS_EXPECT_NEAR(eigen_x[3], ublas_x[3], 1e-12);
}

} // namespace Kratos::Testing
