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
#include <boost/numeric/ublas/matrix.hpp>

// Project includes
#include "testing/testing.h"
#include "includes/kratos_eigen_interface.h"

namespace Kratos::Testing {

KRATOS_TEST_CASE_IN_SUITE(EigenCompressedMatrixCSRConstruction, KratosCoreFastSuite)
{
    // ublas-style (rows, cols, nnz) constructor must yield a compressed
    // matrix with pre-allocated, writable CSR arrays.
    EigenCompressedMatrix<double> matrix(3, 3, 5);

    KRATOS_EXPECT_TRUE(matrix.isCompressed());
    KRATOS_EXPECT_EQ(matrix.size1(), 3);
    KRATOS_EXPECT_EQ(matrix.size2(), 3);
    // As in ublas: nnz() is the filled count (0 before the row pointers are
    // written), while the storage proxies span the allocated capacity.
    KRATOS_EXPECT_EQ(matrix.nnz(), 0);
    KRATOS_EXPECT_EQ(matrix.index1_data().size(), 4);
    KRATOS_EXPECT_EQ(matrix.index2_data().size(), 5);
    KRATOS_EXPECT_EQ(matrix.value_data().size(), 5);

    // Fill the CSR arrays by hand, exactly as the block builder-and-solver
    // does: row pointers, sorted column indices, values, then set_filled.
    // Pattern: [ (0,0) (0,2) ; (1,1) ; (2,0) (2,2) ]
    auto* row_indices = matrix.index1_data().begin();
    auto* col_indices = matrix.index2_data().begin();
    auto* values = matrix.value_data().begin();

    row_indices[0] = 0;
    row_indices[1] = 2;
    row_indices[2] = 3;
    row_indices[3] = 5;

    col_indices[0] = 0;
    col_indices[1] = 2;
    col_indices[2] = 1;
    col_indices[3] = 0;
    col_indices[4] = 2;

    values[0] = 1.0;
    values[1] = 2.0;
    values[2] = 3.0;
    values[3] = 4.0;
    values[4] = 5.0;

    matrix.set_filled(4, 5);
    KRATOS_EXPECT_EQ(matrix.nnz(), 5);

    // Read back through ublas-style element access
    KRATOS_EXPECT_NEAR(matrix(0, 0), 1.0, 1e-12);
    KRATOS_EXPECT_NEAR(matrix(0, 2), 2.0, 1e-12);
    KRATOS_EXPECT_NEAR(matrix(1, 1), 3.0, 1e-12);
    KRATOS_EXPECT_NEAR(matrix(2, 0), 4.0, 1e-12);
    KRATOS_EXPECT_NEAR(matrix(2, 2), 5.0, 1e-12);

    // Entries outside the pattern read as zero through const access
    const auto& r_const_matrix = matrix;
    KRATOS_EXPECT_NEAR(r_const_matrix(0, 1), 0.0, 1e-12);
    KRATOS_EXPECT_NEAR(r_const_matrix(1, 0), 0.0, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(EigenCompressedMatrixResize, KratosCoreFastSuite)
{
    EigenCompressedMatrix<double> matrix(3, 3, 3);
    auto* row_indices = matrix.index1_data().begin();
    row_indices[0] = 0; row_indices[1] = 1; row_indices[2] = 2; row_indices[3] = 3;
    auto* col_indices = matrix.index2_data().begin();
    col_indices[0] = 0; col_indices[1] = 1; col_indices[2] = 2;
    auto* values = matrix.value_data().begin();
    values[0] = 1.0; values[1] = 2.0; values[2] = 3.0;

    // ublas-style non-preserving resize drops values and structure
    matrix.resize(5, 5, false);
    KRATOS_EXPECT_EQ(matrix.size1(), 5);
    KRATOS_EXPECT_EQ(matrix.size2(), 5);
    KRATOS_EXPECT_EQ(matrix.nnz(), 0);
}

KRATOS_TEST_CASE_IN_SUITE(EigenMatrixUblasMemberSurface, KratosCoreFastSuite)
{
    EigenMatrix<double> matrix(2, 3);
    KRATOS_EXPECT_EQ(matrix.size1(), 2);
    KRATOS_EXPECT_EQ(matrix.size2(), 3);

    for (std::size_t i = 0; i < matrix.size1(); ++i)
        for (std::size_t j = 0; j < matrix.size2(); ++j)
            matrix(i, j) = static_cast<double>(i * 3 + j);

    // Preserving resize keeps existing values (ublas resize default semantics)
    matrix.resize(3, 3, true);
    KRATOS_EXPECT_EQ(matrix.size1(), 3);
    KRATOS_EXPECT_NEAR(matrix(1, 2), 5.0, 1e-12);

    // Assignment from an Eigen expression works
    EigenMatrix<double> other(3, 3);
    other = matrix.transpose();
    KRATOS_EXPECT_NEAR(other(2, 1), 5.0, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(EigenVectorUblasMemberSurface, KratosCoreFastSuite)
{
    EigenVector<double> vector(4, 1.5);
    KRATOS_EXPECT_EQ(static_cast<std::size_t>(vector.size()), 4);
    for (std::size_t i = 0; i < 4; ++i) {
        KRATOS_EXPECT_NEAR(vector[i], 1.5, 1e-12);
    }

    vector[2] = 3.0;
    vector.resize(6, true);
    KRATOS_EXPECT_EQ(static_cast<std::size_t>(vector.size()), 6);
    KRATOS_EXPECT_NEAR(vector[2], 3.0, 1e-12);

    vector.resize(2, false);
    KRATOS_EXPECT_EQ(static_cast<std::size_t>(vector.size()), 2);
}

KRATOS_TEST_CASE_IN_SUITE(EigenCompressedMatrixFromUblasDense, KratosCoreFastSuite)
{
    // Construction from a dense uBLAS matrix expression gathers the nonzero
    // entries into compressed storage, mirroring ublas::compressed_matrix's
    // own dense constructor.
    boost::numeric::ublas::matrix<double> dense(3, 3, 0.0);
    dense(0, 0) = 1.0;
    dense(0, 2) = 2.0;
    dense(1, 1) = 3.0;
    dense(2, 0) = 4.0;
    dense(2, 2) = 5.0;

    const EigenCompressedMatrix<double> matrix(dense);

    KRATOS_EXPECT_TRUE(matrix.isCompressed());
    KRATOS_EXPECT_EQ(matrix.size1(), 3);
    KRATOS_EXPECT_EQ(matrix.size2(), 3);
    KRATOS_EXPECT_EQ(matrix.nnz(), 5);

    const auto row_indices = matrix.index1_data();
    const auto col_indices = matrix.index2_data();
    const auto values = matrix.value_data();

    KRATOS_EXPECT_EQ(row_indices[0], 0);
    KRATOS_EXPECT_EQ(row_indices[1], 2);
    KRATOS_EXPECT_EQ(row_indices[2], 3);
    KRATOS_EXPECT_EQ(row_indices[3], 5);

    KRATOS_EXPECT_EQ(col_indices[0], 0);
    KRATOS_EXPECT_EQ(col_indices[1], 2);
    KRATOS_EXPECT_EQ(col_indices[2], 1);
    KRATOS_EXPECT_EQ(col_indices[3], 0);
    KRATOS_EXPECT_EQ(col_indices[4], 2);

    KRATOS_EXPECT_NEAR(values[0], 1.0, 1e-12);
    KRATOS_EXPECT_NEAR(values[1], 2.0, 1e-12);
    KRATOS_EXPECT_NEAR(values[2], 3.0, 1e-12);
    KRATOS_EXPECT_NEAR(values[3], 4.0, 1e-12);
    KRATOS_EXPECT_NEAR(values[4], 5.0, 1e-12);

    // A uBLAS matrix *expression* (not just a container) is accepted too
    const EigenCompressedMatrix<double> scaled(2.0 * dense);
    KRATOS_EXPECT_EQ(scaled.nnz(), 5);
    KRATOS_EXPECT_NEAR(scaled.value_data()[4], 10.0, 1e-12);
}

} // namespace Kratos::Testing
