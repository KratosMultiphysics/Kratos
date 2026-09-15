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
#include <cstring>
#include <numeric>
#include <sstream>
#include <vector>

// External includes

// Project includes
#include "testing/testing.h"
// NOTE: the uBLAS idiom (prod, inner_prod, noalias, trans, row, project, ...)
// exercised below must compile and give the same values under both
// linear-algebra backends; the Eigen-only parts are guarded.
#include "includes/default_interface.h"
#include "containers/array_1d.h"
#include "containers/sparse_graph.h"
#include "geometries/point.h"
#include "utilities/math_utils.h"
#include "spaces/ublas_space.h"

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

// Classes deriving from the dense container types, calling noalias(...)
// unqualified from their member functions. Under the Eigen backend the
// inherited member Eigen::MatrixBase::noalias() must not hide the free form.
class NoaliasTestPoint : public Point
{
public:
    using Point::Point;

    void AssignCoordinates(const array_1d<double, 3>& rOther)
    {
        noalias(this->Coordinates()) = rOther;
    }

    void AddAndSubtract(const array_1d<double, 3>& rAdd, const array_1d<double, 3>& rSub)
    {
        noalias(*this) += rAdd;
        noalias(Coordinates()) -= rSub;
    }

    void WriteInto(array_1d<double, 3>& rTarget) const
    {
        noalias(rTarget) = this->Coordinates();
    }
};

class NoaliasTestArray : public array_1d<double, 3>
{
public:
    void AssignScaled(const array_1d<double, 3>& rOther)
    {
        noalias(*this) = 2.0 * rOther;
    }
};

} // namespace

KRATOS_TEST_CASE_IN_SUITE(EigenCompatProd, KratosCoreFastSuite)
{
    Matrix ublas_A(2, 2), ublas_B(2, 2);
    Matrix eigen_A(2, 2), eigen_B(2, 2);
    FillMatrix(ublas_A); FillMatrix(ublas_B);
    FillMatrix(eigen_A); FillMatrix(eigen_B);

    // dense x dense
    const Matrix ublas_C = prod(ublas_A, ublas_B);
    const Matrix eigen_C = prod(eigen_A, eigen_B);
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t j = 0; j < 2; ++j)
            KRATOS_EXPECT_NEAR(eigen_C(i, j), ublas_C(i, j), 1e-12);

    // dense x vector
    Vector ublas_x(2);
    Vector eigen_x(2);
    FillVector(ublas_x); FillVector(eigen_x);
    const Vector ublas_y = prod(ublas_A, ublas_x);
    const Vector eigen_y = prod(eigen_A, eigen_x);
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
    CompressedMatrix eigen_A(2, 2, 3);
    auto* row_indices = eigen_A.index1_data().begin();
    auto* col_indices = eigen_A.index2_data().begin();
    auto* values = eigen_A.value_data().begin();
    row_indices[0] = 0; row_indices[1] = 2; row_indices[2] = 3;
    col_indices[0] = 0; col_indices[1] = 1; col_indices[2] = 1;
    values[0] = 2.0; values[1] = -1.0; values[2] = 3.0;
    eigen_A.set_filled(3, 3);

    Vector eigen_x(2);
    FillVector(eigen_x);
    const Vector eigen_y = prod(eigen_A, eigen_x);

    KRATOS_EXPECT_NEAR(eigen_y[0], ublas_y[0], 1e-12);
    KRATOS_EXPECT_NEAR(eigen_y[1], ublas_y[1], 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(EigenCompatInnerOuterProd, KratosCoreFastSuite)
{
    Vector ublas_x(2), ublas_y(2);
    Vector eigen_x(2), eigen_y(2);
    FillVector(ublas_x); FillVector(eigen_x);
    ublas_y[0] = 0.5; ublas_y[1] = 4.0;
    eigen_y[0] = 0.5; eigen_y[1] = 4.0;

    KRATOS_EXPECT_NEAR(inner_prod(eigen_x, eigen_y), inner_prod(ublas_x, ublas_y), 1e-12);

    const Matrix ublas_op = outer_prod(ublas_x, ublas_y);
    const Matrix eigen_op = outer_prod(eigen_x, eigen_y);
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t j = 0; j < 2; ++j)
            KRATOS_EXPECT_NEAR(eigen_op(i, j), ublas_op(i, j), 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(EigenCompatTransAndNorms, KratosCoreFastSuite)
{
    Matrix ublas_A(2, 2);
    Matrix eigen_A(2, 2);
    FillMatrix(ublas_A); FillMatrix(eigen_A);

    const Matrix ublas_At = trans(ublas_A);
    const Matrix eigen_At = trans(eigen_A);
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t j = 0; j < 2; ++j)
            KRATOS_EXPECT_NEAR(eigen_At(i, j), ublas_At(i, j), 1e-12);

    // The references are the uBLAS containers explicitly (DenseVector/DenseMatrix
    // stay uBLAS in both backends; Vector/Matrix are Eigen under the Eigen one),
    // so the norms are checked against the uBLAS semantics on either backend.
    DenseVector<double> ublas_x(2);
    Vector eigen_x(2);
    FillVector(ublas_x); FillVector(eigen_x);
    KRATOS_EXPECT_NEAR(norm_1(eigen_x), norm_1(ublas_x), 1e-12);
    KRATOS_EXPECT_NEAR(norm_2(eigen_x), norm_2(ublas_x), 1e-12);
    KRATOS_EXPECT_NEAR(norm_inf(eigen_x), norm_inf(ublas_x), 1e-12);
    KRATOS_EXPECT_NEAR(sum(eigen_x), sum(ublas_x), 1e-12);

    // Matrix norms: ublas norm_1 is the maximum absolute column sum and norm_inf
    // the maximum absolute row sum (not the coefficient-wise lp norms). A
    // non-square matrix with mixed signs so the two, and the Frobenius norm,
    // all differ.
    DenseMatrix<double> ublas_M(2, 3);
    Matrix eigen_M(2, 3);
    for (std::size_t i = 0; i < 2; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            const double v = (i == 0 ? 1.0 : -1.0) * static_cast<double>(i * 3 + j + 1);
            ublas_M(i, j) = v;
            eigen_M(i, j) = v;
        }
    }
    KRATOS_EXPECT_NEAR(norm_1(eigen_M), norm_1(ublas_M), 1e-12);
    KRATOS_EXPECT_NEAR(norm_1(eigen_M), 9.0, 1e-12);   // |3| + |-6|
    KRATOS_EXPECT_NEAR(norm_inf(eigen_M), norm_inf(ublas_M), 1e-12);
    KRATOS_EXPECT_NEAR(norm_inf(eigen_M), 15.0, 1e-12); // |-4| + |-5| + |-6|
    KRATOS_EXPECT_NEAR(norm_frobenius(eigen_M), norm_frobenius(ublas_M), 1e-12);

    // An empty matrix gives 0 in ublas and must not trip Eigen's empty-reduction assertion
    const Matrix eigen_empty(0, 0);
    KRATOS_EXPECT_EQ(norm_1(eigen_empty), 0.0);
    KRATOS_EXPECT_EQ(norm_inf(eigen_empty), 0.0);
}

KRATOS_TEST_CASE_IN_SUITE(EigenCompatNoalias, KratosCoreFastSuite)
{
    Matrix ublas_A(2, 2);
    Matrix eigen_A(2, 2);
    FillMatrix(ublas_A); FillMatrix(eigen_A);
    Vector ublas_x(2), ublas_y(2);
    Vector eigen_x(2), eigen_y(2);
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
    Matrix eigen_C(2, 2);
    noalias(ublas_C) = prod(ublas_A, ublas_A);
    noalias(eigen_C) = prod(eigen_A, eigen_A);
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t j = 0; j < 2; ++j)
            KRATOS_EXPECT_NEAR(eigen_C(i, j), ublas_C(i, j), 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(EigenCompatNoaliasInsideDerivedClass, KratosCoreFastSuite)
{
    array_1d<double, 3> a;
    a[0] = 1.0; a[1] = 2.0; a[2] = 3.0;
    array_1d<double, 3> b;
    b[0] = 0.5; b[1] = -1.0; b[2] = 4.0;

    NoaliasTestPoint point(0.0, 0.0, 0.0);
    point.AssignCoordinates(a);
    KRATOS_EXPECT_VECTOR_NEAR(point.Coordinates(), a, 1e-12);

    point.AddAndSubtract(a, b);
    array_1d<double, 3> expected;
    for (std::size_t i = 0; i < 3; ++i) expected[i] = 2.0 * a[i] - b[i];
    KRATOS_EXPECT_VECTOR_NEAR(point.Coordinates(), expected, 1e-12);

    array_1d<double, 3> copy;
    point.WriteInto(copy);
    KRATOS_EXPECT_VECTOR_NEAR(copy, expected, 1e-12);

    NoaliasTestArray array;
    array.AssignScaled(a);
    for (std::size_t i = 0; i < 3; ++i) KRATOS_EXPECT_NEAR(array[i], 2.0 * a[i], 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(EigenCompatRowColumnSubrange, KratosCoreFastSuite)
{
    Matrix ublas_A(2, 2);
    Matrix eigen_A(2, 2);
    FillMatrix(ublas_A); FillMatrix(eigen_A);

    // reads
    KRATOS_EXPECT_NEAR(row(eigen_A, 1)[0], row(ublas_A, 1)[0], 1e-12);
    KRATOS_EXPECT_NEAR(column(eigen_A, 1)[0], column(ublas_A, 1)(0), 1e-12);

    // write-through
    row(ublas_A, 0) *= 2.0;
    row(eigen_A, 0) *= 2.0;
    KRATOS_EXPECT_NEAR(eigen_A(0, 1), ublas_A(0, 1), 1e-12);

    Vector ublas_x(4);
    Vector eigen_x(4);
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


KRATOS_TEST_CASE_IN_SUITE(EigenCompatLazyFactories, KratosCoreFastSuite)
{
    // The uBLAS factories are lazy in both backends and assign to every dense type
    Matrix A = ZeroMatrix(3, 3);
    KRATOS_EXPECT_EQ(A.size1(), 3);
    KRATOS_EXPECT_EQ(A(1, 1), 0.0);
    noalias(A) = IdentityMatrix(3);
    KRATOS_EXPECT_EQ(A(1, 1), 1.0);
    KRATOS_EXPECT_EQ(A(0, 1), 0.0);
    A += IdentityMatrix(3);
    KRATOS_EXPECT_EQ(A(2, 2), 2.0);
    const Matrix S = ScalarMatrix(2, 2, 4.0);
    KRATOS_EXPECT_EQ(S(1, 0), 4.0);
    const Vector z = ZeroVector(4);
    KRATOS_EXPECT_EQ(z.size(), 4);
    KRATOS_EXPECT_EQ(z[3], 0.0);
    const Vector sv = ScalarVector(3, 2.0);
    KRATOS_EXPECT_EQ(sv[2], 2.0);
    const Vector uv = UnitVector(3, 1);
    KRATOS_EXPECT_EQ(uv[0], 0.0);
    KRATOS_EXPECT_EQ(uv[1], 1.0);
    BoundedMatrix<double, 3, 3> B = ZeroMatrix(3, 3);
    KRATOS_EXPECT_EQ(B(2, 2), 0.0);
    B = IdentityMatrix(3);
    KRATOS_EXPECT_EQ(B(2, 2), 1.0);
    const BoundedVector<double, 3> bv = ZeroVector(3);
    KRATOS_EXPECT_EQ(bv[0], 0.0);
    const array_1d<double, 3> az = ZeroVector(3);
    KRATOS_EXPECT_EQ(az[2], 0.0);
    const Matrix Pz = prod(ZeroMatrix(3, 3), A);
    KRATOS_EXPECT_EQ(Pz(0, 0), 0.0);
    CompressedMatrix Sp = ZeroMatrix(4, 4);
    KRATOS_EXPECT_EQ(Sp.size1(), 4);
    KRATOS_EXPECT_EQ(Sp.nnz(), 0);
    // element insertion keeps the CSR arrays packed (row pointers usable right away)
    Sp(1, 2) = 3.0;
    Sp(0, 0) = 1.0;
    Sp(3, 1) = 2.0;
    KRATOS_EXPECT_EQ(Sp.nnz(), 3);
    KRATOS_EXPECT_EQ(static_cast<std::size_t>(Sp.index1_data()[4]), 3);
    KRATOS_EXPECT_EQ(static_cast<std::size_t>(Sp.index1_data()[2]), 2);
    KRATOS_EXPECT_EQ(Sp(1, 2), 3.0);
    Matrix Inv(3, 3);
    Inv.assign(IdentityMatrix(3));
    KRATOS_EXPECT_EQ(Inv(1, 1), 1.0);
}

KRATOS_TEST_CASE_IN_SUITE(EigenCompatVectorTransposeSemantics, KratosCoreFastSuite)
{
    // uBLAS: trans() of a vector is the identity, prod(v, M) is v^T M
    Vector v(3); v[0] = 1; v[1] = 2; v[2] = 3;
    Vector w(2); w[0] = 10; w[1] = 20;
    const Matrix O = outer_prod(v, trans(w));
    KRATOS_EXPECT_EQ(O.size1(), 3);
    KRATOS_EXPECT_EQ(O.size2(), 2);
    KRATOS_EXPECT_EQ(O(2, 1), 60.0);
    Matrix M(2, 3);
    M(0,0) = 1; M(0,1) = 2; M(0,2) = 3; M(1,0) = 4; M(1,1) = 5; M(1,2) = 6;
    const Vector mv = prod(M, v);
    KRATOS_EXPECT_EQ(mv.size(), 2);
    KRATOS_EXPECT_EQ(mv[0], 14.0);
    KRATOS_EXPECT_EQ(mv[1], 32.0);
    const Vector vm = prod(w, M);
    KRATOS_EXPECT_EQ(vm.size(), 3);
    KRATOS_EXPECT_EQ(vm[0], 90.0);
    KRATOS_EXPECT_EQ(vm[2], 150.0);
    const Vector vm2 = prod(trans(w), M);
    KRATOS_EXPECT_EQ(vm2[1], 120.0);
    const Matrix MtM = prod(trans(M), M);
    KRATOS_EXPECT_EQ(MtM.size1(), 3);
    KRATOS_EXPECT_EQ(MtM(0, 0), 17.0);
    M = trans(M); // alias-safe self transpose
    KRATOS_EXPECT_EQ(M.size1(), 3);
    KRATOS_EXPECT_EQ(M(2, 1), 6.0);
#ifdef KRATOS_USE_EIGEN_BACKEND
    const Matrix Mt = trans(M);
    const Matrix pinned = prod<Matrix>(M, Mt); // result-type pinning (ambiguous in uBLAS itself)
    KRATOS_EXPECT_EQ(pinned.size1(), 3);
#endif
    array_1d<double, 3> a3; a3[0] = 1; a3[1] = 1; a3[2] = 1;
    KRATOS_EXPECT_EQ(inner_prod(v, a3), 6.0);
    // an N x 1 bounded matrix is a matrix by type: its transpose is 1 x N
    BoundedMatrix<double, 3, 1> column; column(0, 0) = 1; column(1, 0) = 2; column(2, 0) = 3;
    BoundedMatrix<double, 1, 3> row_matrix;
    noalias(row_matrix) = trans(column);
    KRATOS_EXPECT_EQ(row_matrix(0, 2), 3.0);
    const BoundedMatrix<double, 3, 3> outer = prod(column, row_matrix);
    KRATOS_EXPECT_EQ(outer(2, 1), 6.0);
#ifdef KRATOS_USE_EIGEN_BACKEND
    KRATOS_EXPECT_EQ(&trans(v), &v);
    KRATOS_EXPECT_EQ(&trans(a3), &a3);
#endif
}

KRATOS_TEST_CASE_IN_SUITE(EigenCompatProxies, KratosCoreFastSuite)
{
    Vector v(3); v[0] = 1; v[1] = 2; v[2] = 3;
    Matrix R(3, 3); R.clear();
    row(R, 1) = v;
    KRATOS_EXPECT_EQ(R(1, 2), 3.0);
    const MatrixRow<Matrix> r1 = row(R, 1);
    KRATOS_EXPECT_EQ(r1.size(), 3);
    for (std::size_t i = 0; i < r1.size(); ++i) {
        KRATOS_EXPECT_EQ(r1[i], v[i]);
    }
    const Vector r1v = row(R, 1);
    KRATOS_EXPECT_EQ(r1v[2], 3.0);
    const Matrix Ro = outer_prod(row(R, 1), trans(row(R, 1)));
    KRATOS_EXPECT_EQ(Ro(2, 2), 9.0);
    KRATOS_EXPECT_EQ(inner_prod(row(R, 1), v), 14.0);
    column(R, 0) = v;
    KRATOS_EXPECT_EQ(R(2, 0), 3.0);
    const MatrixColumn c0 = column(R, 0);
    KRATOS_EXPECT_EQ(c0.size(), 3);
    KRATOS_EXPECT_EQ(c0[2], 3.0);
    noalias(row(R, 2)) = v;
    KRATOS_EXPECT_EQ(R(2, 2), 3.0);
    Vector long_v(6);
    for (std::size_t i = 0; i < 6; ++i) long_v[i] = static_cast<double>(i);
    Vector w(2); w[0] = 10; w[1] = 20;
    subrange(long_v, 2, 4) = w;
    KRATOS_EXPECT_EQ(long_v[2], 10.0);
    KRATOS_EXPECT_EQ(long_v[3], 20.0);
    const Vector sr = subrange(long_v, 1, 3);
    KRATOS_EXPECT_EQ(sr.size(), 2);
    const VectorRange vr = project(long_v, range(4, 6));
    KRATOS_EXPECT_EQ(vr.size(), 2);
    KRATOS_EXPECT_EQ(vr[1], 5.0);
    project(long_v, range(0, 2)) *= 2.0;
    KRATOS_EXPECT_EQ(long_v[1], 2.0);
    Matrix Big(4, 4); Big.clear();
    subrange(Big, 1, 3, 1, 3) = IdentityMatrix(2);
    KRATOS_EXPECT_EQ(Big(2, 2), 1.0);
    const MatrixRange mr = project(Big, range(1, 3), range::all());
    KRATOS_EXPECT_EQ(mr.size1(), 2);
    KRATOS_EXPECT_EQ(mr.size2(), 4);
    noalias(subrange(Big, 0, 2, 0, 2)) += IdentityMatrix(2);
    KRATOS_EXPECT_EQ(Big(0, 0), 1.0);
    const Matrix& r_big = Big;
    const Vector crow = row(r_big, 1);
    KRATOS_EXPECT_EQ(crow[1], 2.0);
    VectorSlice vsl = project(long_v, slice(0, 2, 3));
    KRATOS_EXPECT_EQ(vsl.size(), 3);
    KRATOS_EXPECT_EQ(vsl[1], 10.0);
    vsl[0] = 100.0;
    KRATOS_EXPECT_EQ(long_v[0], 100.0);
    BoundedMatrix<double, 3, 3> B = IdentityMatrix(3);
    const array_1d<double, 3> arow = row(B, 2);
    KRATOS_EXPECT_EQ(arow[2], 1.0);
    // r1 is a live view: column(R, 0) = v above changed its first entry to 2
    const double sip = std::inner_product(r1.begin(), r1.end(), v.begin(), 0.0);
    KRATOS_EXPECT_EQ(sip, 15.0);
}

KRATOS_TEST_CASE_IN_SUITE(EigenCompatStorageAccessAndResize, KratosCoreFastSuite)
{
    Vector v(3); v[0] = 1; v[1] = 2; v[2] = 3;
    KRATOS_EXPECT_EQ(v.data().size(), 3);
    KRATOS_EXPECT_EQ(*v.data().begin(), 1.0);
    KRATOS_EXPECT_EQ(&v.data()[0], &v[0]);
    Vector w(2); w[0] = 10; w[1] = 20;
    std::memcpy(&v.data()[0], &w.data()[0], 2 * sizeof(double));
    KRATOS_EXPECT_EQ(v[0], 10.0);
    Matrix Mm(2, 2, 3.0);
    KRATOS_EXPECT_EQ(Mm.data().size(), 4);
    KRATOS_EXPECT_EQ(Mm.data()[3], 3.0);
    array_1d<double, 3> a3; a3[0] = 7;
    KRATOS_EXPECT_EQ(a3.data()[0], 7.0);

    // resize keeps the entries by default (uBLAS semantics)
    Vector gv(2); gv[0] = 1; gv[1] = 2;
    gv.resize(3);
    KRATOS_EXPECT_EQ(gv.size(), 3);
    KRATOS_EXPECT_EQ(gv[1], 2.0);
    gv.resize(1, false);
    KRATOS_EXPECT_EQ(gv.size(), 1);
    Matrix gm(2, 2, 5.0);
    gm.resize(3, 3);
    KRATOS_EXPECT_EQ(gm(1, 1), 5.0);
    gm.resize(2, 2, false);
    KRATOS_EXPECT_EQ(gm.size1(), 2);

    // non-arithmetic scalars as plain containers
    DenseVector<Matrix> vm(2);
    vm[0] = Matrix(2, 2, 1.0);
    vm.resize(4);
    KRATOS_EXPECT_EQ(vm[0](1, 1), 1.0);
    DenseVector<bool> vb(3, false); vb[1] = true;
    KRATOS_EXPECT_TRUE(vb[1]);
    DenseMatrix<unsigned int> mu(2, 3, 7u);
    KRATOS_EXPECT_EQ(mu(1, 2), 7u);

    // move-only scalars (a container of sparse graphs): the preserving resize
    // must not require copy assignment
    DenseVector<SparseGraph<>> graphs(2);
    graphs[0].AddEntry(1, 2);
    graphs.resize(4);
    KRATOS_EXPECT_EQ(graphs.size(), 4);
    KRATOS_EXPECT_TRUE(graphs[0].Has(1, 2));
    graphs.resize(1, false);
    KRATOS_EXPECT_EQ(graphs.size(), 1);

    // uBLAS text format round trip
    std::stringstream ss;
    ss << w;
    KRATOS_EXPECT_EQ(ss.str(), "[2](10,20)");
    Vector back;
    ss >> back;
    KRATOS_EXPECT_EQ(back.size(), 2);
    KRATOS_EXPECT_EQ(back[1], 20.0);
    std::stringstream ss2;
    ss2 << Mm;
    KRATOS_EXPECT_EQ(ss2.str(), "[2,2]((3,3),(3,3))");
    Matrix mback;
    ss2 >> mback;
    KRATOS_EXPECT_EQ(mback(1, 1), 3.0);

#ifdef KRATOS_USE_EIGEN_BACKEND
    // the plugins: resize(..., preserve) and size1()/size2() on raw Eigen objects
    Eigen::MatrixXd raw(2, 2); raw.setOnes();
    raw.resize(3, 3, true);
    KRATOS_EXPECT_EQ(raw(1, 1), 1.0);
    KRATOS_EXPECT_EQ(raw.size1(), 3);
    KRATOS_EXPECT_EQ(prod(Mm, Mm).size1(), 2);
    double* p_data = v.data();
    KRATOS_EXPECT_EQ(p_data, &v[0]);
#endif
}

KRATOS_TEST_CASE_IN_SUITE(EigenCompatLUFactorization, KratosCoreFastSuite)
{
    Matrix L(3, 3);
    L(0,0) = 0; L(0,1) = 2; L(0,2) = 1; L(1,0) = 1; L(1,1) = 1; L(1,2) = 0; L(2,0) = 3; L(2,1) = 0; L(2,2) = 1; // det = -5
    Matrix factors(L);
    permutation_matrix<std::size_t> pm(3);
    const std::size_t singular = lu_factorize(factors, pm);
    KRATOS_EXPECT_EQ(singular, 0);
    double det = 1.0;
    for (std::size_t i = 0; i < 3; ++i) det *= (pm(i) == i) ? factors(i, i) : -factors(i, i);
    KRATOS_EXPECT_NEAR(det, -5.0, 1e-12);
    Matrix inverse(3, 3);
    inverse.assign(IdentityMatrix(3));
    lu_substitute(factors, pm, inverse);
    const Matrix identity = prod(L, inverse);
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j)
            KRATOS_EXPECT_NEAR(identity(i, j), i == j ? 1.0 : 0.0, 1e-12);
    Vector rhs(3); rhs[0] = 1; rhs[1] = 2; rhs[2] = 3;
    Vector solution(rhs);
    lu_substitute(factors, pm, solution);
    const Vector check = prod(L, solution);
    for (std::size_t i = 0; i < 3; ++i) KRATOS_EXPECT_NEAR(check[i], rhs[i], 1e-12);
    Matrix singular_matrix(2, 2);
    singular_matrix(0,0) = 1; singular_matrix(0,1) = 2; singular_matrix(1,0) = 2; singular_matrix(1,1) = 4;
    permutation_matrix<std::size_t> pm2(2);
    KRATOS_EXPECT_NE(lu_factorize(singular_matrix, pm2), 0);
    // through MathUtils on the same data
    Matrix math_inverse(3, 3);
    double math_det;
    MathUtils<double>::InvertMatrix(L, math_inverse, math_det);
    KRATOS_EXPECT_NEAR(math_det, -5.0, 1e-12);
    KRATOS_EXPECT_NEAR(math_inverse(0, 1), inverse(0, 1), 1e-12);
    Matrix M4 = IdentityMatrix(4);
    M4(0, 3) = 2.0;
    KRATOS_EXPECT_NEAR(MathUtils<double>::Cofactor(M4, 0, 0), 1.0, 1e-12);
    KRATOS_EXPECT_NEAR(MathUtils<double>::Cofactor(M4, 3, 0), -2.0, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(EigenCompatCompressedMatrixIterators, KratosCoreFastSuite)
{
    // The ublas row/entry iterator idiom over a CSR matrix, as the mesh
    // refinement and cutting utilities of MeshingApplication use it.
    compressed_matrix<int> coord(4, 4);
    coord.push_back(0, 1, -2);
    coord.push_back(1, 2, -1);
    coord.push_back(2, 3, -2);

    using i1_t = compressed_matrix<int>::iterator1;
    using i2_t = compressed_matrix<int>::iterator2;

    std::size_t marked_count = 0;
    for (i1_t i1 = coord.begin1(); i1 != coord.end1(); ++i1) {
        for (i2_t i2 = i1.begin(); i2 != i1.end(); ++i2) {
            if (coord(i2.index1(), i2.index2()) == -2) {
                ++marked_count;
            }
        }
    }
    KRATOS_EXPECT_EQ(marked_count, 2);

    // index1()/index2() report the position of the entry, and the entries are
    // writable both through operator() and through the iterator itself
    std::vector<std::pair<std::size_t, std::size_t>> visited;
    for (i1_t i1 = coord.begin1(); i1 != coord.end1(); ++i1) {
        for (i2_t i2 = i1.begin(); i2 != i1.end(); ++i2) {
            visited.emplace_back(i2.index1(), i2.index2());
            if (coord(i2.index1(), i2.index2()) == -2) {
                coord(i2.index1(), i2.index2()) = 7;
            } else {
                *i2 = 5;
            }
        }
    }
    KRATOS_EXPECT_EQ(visited.size(), 3);
    KRATOS_EXPECT_EQ(visited[0].first, 0);
    KRATOS_EXPECT_EQ(visited[0].second, 1);
    KRATOS_EXPECT_EQ(visited[2].first, 2);
    KRATOS_EXPECT_EQ(visited[2].second, 3);
    KRATOS_EXPECT_EQ(coord(0, 1), 7);
    KRATOS_EXPECT_EQ(coord(1, 2), 5);
    KRATOS_EXPECT_EQ(coord(2, 3), 7);

    // const traversal visits the same entries
    const compressed_matrix<int>& r_const_coord = coord;
    std::size_t const_count = 0;
    for (auto i1 = r_const_coord.begin1(); i1 != r_const_coord.end1(); ++i1) {
        KRATOS_EXPECT_TRUE(i1.index1() < 4);
        for (auto i2 = i1.begin(); i2 != i1.end(); ++i2) {
            KRATOS_EXPECT_TRUE(*i2 == 7 || *i2 == 5);
            ++const_count;
        }
    }
    KRATOS_EXPECT_EQ(const_count, 3);
}

KRATOS_TEST_CASE_IN_SUITE(EigenCompatUblasSpaceNames, KratosCoreFastSuite)
{
    // The ublas space spellings must keep naming the backend's spaces, so code
    // written against them (Mapping's MPI definitions, FSI, Dam) compiles in
    // both backends.
    using DenseSpaceType = TUblasDenseSpace<double>;
    using LocalSpaceType = UblasSpace<double, Matrix, Vector>;
    using SparseSpaceType = UblasSpace<double, CompressedMatrix, Vector>;

    Matrix a(2, 2);
    a(0, 0) = 1.0; a(0, 1) = 2.0;
    a(1, 0) = 3.0; a(1, 1) = 4.0;
    Vector x(2);
    x[0] = 1.0; x[1] = 1.0;
    Vector y(2);

    LocalSpaceType::Mult(a, x, y);
    KRATOS_EXPECT_DOUBLE_EQ(y[0], 3.0);
    KRATOS_EXPECT_DOUBLE_EQ(y[1], 7.0);

    DenseSpaceType::Mult(a, x, y);
    KRATOS_EXPECT_DOUBLE_EQ(y[0], 3.0);
    KRATOS_EXPECT_EQ(TUblasSparseSpace<double>::Size(x), 2);

    CompressedMatrix k(2, 2, 2);
    k.push_back(0, 0, 3.0);
    k.push_back(1, 1, 4.0);
    k.complete_index1_data();
    KRATOS_EXPECT_DOUBLE_EQ(SparseSpaceType::TwoNorm(k), 5.0);
}

} // namespace Kratos::Testing
