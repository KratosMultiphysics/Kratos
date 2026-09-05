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

// The mixed uBLAS/Eigen idiom layer only exists under the Eigen backend.
#ifdef KRATOS_USE_EIGEN_BACKEND

// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "containers/array_1d.h"
#include "includes/ublas_interface.h"

namespace Kratos::Testing {

namespace {

BoundedMatrix<double, 3, 3> DiagonalBounded(const double Value)
{
    BoundedMatrix<double, 3, 3> result;
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j)
            result(i, j) = (i == j) ? Value : 0.0;
    return result;
}

} // namespace

KRATOS_TEST_CASE_IN_SUITE(EigenUblasNoAliasTargets, KratosCoreFastSuite)
{
    const auto R = DiagonalBounded(2.0);
    const array_1d<double, 3> x{1.0, 2.0, 3.0};
    Matrix M = ZeroMatrix(3, 3);
    M(0, 0) = 1.0; M(1, 1) = 2.0; M(2, 2) = 3.0;
    Vector v(3);
    v[0] = 1.0; v[1] = 1.0; v[2] = 1.0;

    // Eigen target, Eigen RHS
    array_1d<double, 3> a;
    noalias(a) = prod(R, x);
    KRATOS_EXPECT_EQ(a[2], 6.0);
    noalias(a) += prod(R, x);
    KRATOS_EXPECT_EQ(a[2], 12.0);
    noalias(a) -= prod(R, x);
    KRATOS_EXPECT_EQ(a[2], 6.0);

    // Eigen target, uBLAS RHS
    noalias(a) = ZeroVector(3);
    KRATOS_EXPECT_EQ(a[0], 0.0);
    noalias(a) += v;
    KRATOS_EXPECT_EQ(a[1], 1.0);
    noalias(a) -= 2.0 * v;
    KRATOS_EXPECT_EQ(a[1], -1.0);

    // Eigen bounded-matrix target, uBLAS and Eigen RHS
    BoundedMatrix<double, 3, 3> B;
    noalias(B) = ZeroMatrix(3, 3);
    KRATOS_EXPECT_EQ(B(1, 1), 0.0);
    noalias(B) = prod(R, R);
    KRATOS_EXPECT_EQ(B(1, 1), 4.0);
    noalias(B) += IdentityMatrix(3);
    KRATOS_EXPECT_EQ(B(1, 1), 5.0);
    KRATOS_EXPECT_EQ(B(0, 1), 0.0);

    // uBLAS vector target, Eigen RHS (Map path)
    Vector w(3);
    noalias(w) = prod(R, x);
    KRATOS_EXPECT_EQ(w[2], 6.0);
    noalias(w) += 2.0 * x;
    KRATOS_EXPECT_EQ(w[2], 12.0);
    noalias(w) -= x;
    KRATOS_EXPECT_EQ(w[2], 9.0);

    // uBLAS vector target, uBLAS RHS (must keep the exact ublas semantics)
    noalias(w) = prod(M, v);
    KRATOS_EXPECT_EQ(w[2], 3.0);

    // uBLAS matrix target, Eigen RHS and uBLAS RHS
    Matrix K(3, 3);
    noalias(K) = prod(R, R);
    KRATOS_EXPECT_EQ(K(2, 2), 4.0);
    noalias(K) += IdentityMatrix(3);
    KRATOS_EXPECT_EQ(K(2, 2), 5.0);
    noalias(K) = prod(M, M);
    KRATOS_EXPECT_EQ(K(2, 2), 9.0);
}

KRATOS_TEST_CASE_IN_SUITE(EigenUblasMixedProducts, KratosCoreFastSuite)
{
    const auto R = DiagonalBounded(2.0);
    const array_1d<double, 3> x{1.0, 2.0, 3.0};
    Matrix M = ZeroMatrix(3, 3);
    M(0, 0) = 1.0; M(1, 1) = 2.0; M(2, 2) = 3.0;
    Vector v(3);
    v[0] = 1.0; v[1] = 2.0; v[2] = 3.0;

    // Eigen matrix x uBLAS vector
    const Vector r1 = prod(R, v);
    KRATOS_EXPECT_EQ(r1[2], 6.0);

    // uBLAS matrix x Eigen vector
    const Vector r2 = prod(M, x);
    KRATOS_EXPECT_EQ(r2[2], 9.0);

    // also absorbable into an array_1d
    const array_1d<double, 3> r3 = prod(M, x);
    KRATOS_EXPECT_EQ(r3[1], 4.0);

    // Eigen matrix x uBLAS matrix and reverse
    const Matrix m1 = prod(R, M);
    KRATOS_EXPECT_EQ(m1(2, 2), 6.0);
    const Matrix m2 = prod(M, R);
    KRATOS_EXPECT_EQ(m2(1, 1), 4.0);
    const BoundedMatrix<double, 3, 3> m3 = prod(M, R);
    KRATOS_EXPECT_EQ(m3(0, 0), 2.0);

    // mixed inner products, both orders
    KRATOS_EXPECT_EQ(inner_prod(x, v), 14.0);
    KRATOS_EXPECT_EQ(inner_prod(v, x), 14.0);

    // mixed outer products, both orders
    const Matrix o1 = outer_prod(x, v);
    KRATOS_EXPECT_EQ(o1(2, 2), 9.0);
    const Matrix o2 = outer_prod(v, x);
    KRATOS_EXPECT_EQ(o2(2, 1), 6.0);

    // trans on the Eigen side composes with the mixed product
    const Vector r4 = prod(trans(R), v);
    KRATOS_EXPECT_EQ(r4[0], 2.0);
}

KRATOS_TEST_CASE_IN_SUITE(EigenUblasMixedAdditive, KratosCoreFastSuite)
{
    const array_1d<double, 3> x{1.0, 2.0, 3.0};
    Vector v(3);
    v[0] = 10.0; v[1] = 20.0; v[2] = 30.0;

    const array_1d<double, 3> sum_xv = x + v;
    KRATOS_EXPECT_EQ(sum_xv[2], 33.0);
    const array_1d<double, 3> sum_vx = v + x;
    KRATOS_EXPECT_EQ(sum_vx[2], 33.0);
    const array_1d<double, 3> diff_xv = x - v;
    KRATOS_EXPECT_EQ(diff_xv[2], -27.0);
    const array_1d<double, 3> diff_vx = v - x;
    KRATOS_EXPECT_EQ(diff_vx[2], 27.0);

    // compound assignment of Eigen expressions into the uBLAS containers
    v += 2.0 * x;
    KRATOS_EXPECT_EQ(v[2], 36.0);
    v -= x;
    KRATOS_EXPECT_EQ(v[2], 33.0);

    Matrix K = ZeroMatrix(3, 3);
    const auto R = DiagonalBounded(2.0);
    K += R;
    KRATOS_EXPECT_EQ(K(1, 1), 2.0);
    K -= 0.5 * R;
    KRATOS_EXPECT_EQ(K(1, 1), 1.0);

    // matrix mixed +/-
    const BoundedMatrix<double, 3, 3> sum_matrix = R + IdentityMatrix(3);
    KRATOS_EXPECT_EQ(sum_matrix(1, 1), 3.0);
    KRATOS_EXPECT_EQ(sum_matrix(0, 1), 0.0);
    const BoundedMatrix<double, 3, 3> diff_matrix = IdentityMatrix(3) - R;
    KRATOS_EXPECT_EQ(diff_matrix(1, 1), -1.0);
}

KRATOS_TEST_CASE_IN_SUITE(EigenUblasProxyAccessors, KratosCoreFastSuite)
{
    BoundedMatrix<double, 3, 3> R;
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j)
            R(i, j) = static_cast<double>(3 * i + j);

    // row read into an array_1d
    const array_1d<double, 3> second_row = row(R, 1);
    KRATOS_EXPECT_EQ(second_row[0], 3.0);
    KRATOS_EXPECT_EQ(second_row[2], 5.0);

    // row write-through from a column-vector shaped object
    const array_1d<double, 3> replacement{-1.0, -2.0, -3.0};
    row(R, 0) = replacement;
    KRATOS_EXPECT_EQ(R(0, 1), -2.0);

    // column read
    const array_1d<double, 3> second_column = column(R, 1);
    KRATOS_EXPECT_EQ(second_column[0], -2.0);
    KRATOS_EXPECT_EQ(second_column[2], 7.0);

    // vector subrange read and write
    array_1d<double, 6> long_array{0.0, 1.0, 2.0, 3.0, 4.0, 5.0};
    const array_1d<double, 3> middle = subrange(long_array, 2, 5);
    KRATOS_EXPECT_EQ(middle[0], 2.0);
    KRATOS_EXPECT_EQ(middle[2], 4.0);
    subrange(long_array, 0, 3) = middle;
    KRATOS_EXPECT_EQ(long_array[1], 3.0);

    // matrix subrange
    const BoundedMatrix<double, 2, 2> corner = subrange(R, 1, 3, 1, 3);
    KRATOS_EXPECT_EQ(corner(0, 0), 4.0);
    KRATOS_EXPECT_EQ(corner(1, 1), 8.0);
}

KRATOS_TEST_CASE_IN_SUITE(EigenUblasElementWiseAndNorms, KratosCoreFastSuite)
{
    const array_1d<double, 3> x{1.0, 2.0, 3.0};
    const array_1d<double, 3> y{2.0, 3.0, 4.0};
    Vector v(3);
    v[0] = 2.0; v[1] = 2.0; v[2] = 2.0;

    // pure Eigen element-wise
    const array_1d<double, 3> ep = element_prod(x, y);
    KRATOS_EXPECT_EQ(ep[2], 12.0);
    const array_1d<double, 3> ed = element_div(y, x);
    KRATOS_EXPECT_EQ(ed[2], 4.0 / 3.0);

    // mixed element-wise
    const array_1d<double, 3> epm = element_prod(x, v);
    KRATOS_EXPECT_EQ(epm[2], 6.0);
    const array_1d<double, 3> edm = element_div(v, x);
    KRATOS_EXPECT_EQ(edm[2], 2.0 / 3.0);

    // norms
    const auto R = DiagonalBounded(2.0);
    KRATOS_EXPECT_NEAR(norm_frobenius(R), std::sqrt(12.0), 1e-12);
    KRATOS_EXPECT_NEAR(norm_2(x), std::sqrt(14.0), 1e-12);
    KRATOS_EXPECT_EQ(norm_1(x), 6.0);
    KRATOS_EXPECT_EQ(norm_inf(x), 3.0);
    KRATOS_EXPECT_EQ(sum(x), 6.0);
}

KRATOS_TEST_CASE_IN_SUITE(EigenBoundedTypesSemantics, KratosCoreFastSuite)
{
    // (size, value) partial-fill constructor of the bounded vector
    const BoundedVector<double, 3> filled(3, 1.5);
    KRATOS_EXPECT_EQ(filled[2], 1.5);

    // construction from uBLAS expressions
    const BoundedVector<double, 3> zeroed = ZeroVector(3);
    KRATOS_EXPECT_EQ(zeroed[1], 0.0);
    const BoundedMatrix<double, 2, 2> identity = IdentityMatrix(2);
    KRATOS_EXPECT_EQ(identity(0, 0), 1.0);
    KRATOS_EXPECT_EQ(identity(0, 1), 0.0);

    // resize to the static size is accepted; a different size throws
    BoundedVector<double, 3> v(3, 2.0);
    v.resize(3);
    KRATOS_EXPECT_EQ(v[0], 2.0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(v.resize(2), "cannot be resized");
    BoundedMatrix<double, 2, 3> m;
    m.resize(2, 3);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(m.resize(3, 3), "cannot be resized");

    // swap
    BoundedVector<double, 3> v1(3, 1.0);
    BoundedVector<double, 3> v2(3, 2.0);
    v1.swap(v2);
    KRATOS_EXPECT_EQ(v1[0], 2.0);
    KRATOS_EXPECT_EQ(v2[0], 1.0);

    // comparison
    KRATOS_EXPECT_FALSE(v1 == v2);
    v2 = v1;
    KRATOS_EXPECT_TRUE(v1 == v2);
}

} // namespace Kratos::Testing

#endif // KRATOS_USE_EIGEN_BACKEND
