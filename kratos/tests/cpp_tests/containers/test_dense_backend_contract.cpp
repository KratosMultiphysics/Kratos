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
#include <sstream>
#include <type_traits>

// External includes

// Project includes
#include "testing/testing.h"
#include "containers/array_1d.h"
#include "includes/ublas_interface.h"

namespace Kratos::Testing {

// This suite pins down the contract that BOTH linear-algebra backends must
// honour for the fixed-size dense types (array_1d, BoundedVector,
// BoundedMatrix). It is the executable specification the Eigen-backed
// implementations are written against; every assertion here also holds for
// the historical uBLAS implementations.

KRATOS_TEST_CASE_IN_SUITE(DenseBackendArray1dLayout, KratosCoreFastSuite)
{
    // Variable component access (variable.h GetValueByIndex), DataTypeTraits
    // contiguity and the MPI message buffers all reinterpret an
    // array_1d<double, N> as a plain double[N] starting at the object address.
    static_assert(sizeof(array_1d<double, 3>) == 3 * sizeof(double));
    static_assert(sizeof(array_1d<double, 4>) == 4 * sizeof(double));
    static_assert(sizeof(array_1d<double, 6>) == 6 * sizeof(double));
    static_assert(sizeof(array_1d<double, 9>) == 9 * sizeof(double));
    static_assert(std::is_standard_layout_v<array_1d<double, 3>>);
    static_assert(alignof(array_1d<double, 4>) == alignof(double));

    array_1d<double, 4> a4;
    KRATOS_EXPECT_EQ(reinterpret_cast<double*>(&a4), &a4[0]);
    KRATOS_EXPECT_EQ(&a4[0] + 1, &a4[1]);
    KRATOS_EXPECT_EQ(&a4[0] + 3, &a4[3]);

    // The component-variable access pattern itself
    array_1d<double, 3> a3;
    a3[0] = 1.0; a3[1] = 2.0; a3[2] = 3.0;
    double* p_component = static_cast<double*>(static_cast<void*>(&a3));
    KRATOS_EXPECT_EQ(p_component[0], 1.0);
    KRATOS_EXPECT_EQ(p_component[1], 2.0);
    KRATOS_EXPECT_EQ(p_component[2], 3.0);
}

KRATOS_TEST_CASE_IN_SUITE(DenseBackendArray1dConstructors, KratosCoreFastSuite)
{
    // (size, value) fills the first `size` entries
    const array_1d<double, 3> filled(3, 1.5);
    for (std::size_t i = 0; i < 3; ++i) {
        KRATOS_EXPECT_EQ(filled[i], 1.5);
    }

    // initializer list
    const array_1d<double, 3> listed{1.0, 2.0, 3.0};
    KRATOS_EXPECT_EQ(listed[0], 1.0);
    KRATOS_EXPECT_EQ(listed[1], 2.0);
    KRATOS_EXPECT_EQ(listed[2], 3.0);

    // copy + comparison
    const array_1d<double, 3> copied(listed);
    KRATOS_EXPECT_TRUE(copied == listed);
    array_1d<double, 3> other(listed);
    other[2] = -3.0;
    KRATOS_EXPECT_FALSE(other == listed);

    // both access operators
    KRATOS_EXPECT_EQ(listed(1), listed[1]);

    // size is the static size
    KRATOS_EXPECT_EQ(listed.size(), 3);
}

KRATOS_TEST_CASE_IN_SUITE(DenseBackendArray1dResizeAndClear, KratosCoreFastSuite)
{
    array_1d<double, 3> a{1.0, 2.0, 3.0};

    // resize to the static size preserves the values (generic code shared
    // with dynamic vectors calls resize unconditionally)
    a.resize(3);
    KRATOS_EXPECT_EQ(a[0], 1.0);
    KRATOS_EXPECT_EQ(a[2], 3.0);

    // resize without preservation value-initializes
    a.resize(3, false);
    KRATOS_EXPECT_EQ(a[0], 0.0);
    KRATOS_EXPECT_EQ(a[2], 0.0);

    array_1d<double, 3> b{1.0, 2.0, 3.0};
    b.clear();
    KRATOS_EXPECT_EQ(b[0], 0.0);
    KRATOS_EXPECT_EQ(b[2], 0.0);
}

KRATOS_TEST_CASE_IN_SUITE(DenseBackendArray1dStreamFormat, KratosCoreFastSuite)
{
    // The uBLAS text format is the cross-backend interchange format: python
    // __str__, error messages and regression references all rely on it.
    const array_1d<double, 3> a{1.0, 2.5, -3.0};
    std::stringstream buffer;
    buffer << a;
    KRATOS_EXPECT_EQ(buffer.str(), "[3](1,2.5,-3)");
}

KRATOS_TEST_CASE_IN_SUITE(DenseBackendArray1dHash, KratosCoreFastSuite)
{
    const array_1d<double, 3> a{1.0, 2.0, 3.0};
    const array_1d<double, 3> b{1.0, 2.0, 3.0};
    const array_1d<double, 3> c{3.0, 2.0, 1.0};

    std::hash<array_1d<double, 3>> hasher;
    KRATOS_EXPECT_EQ(hasher(a), hasher(b));
    KRATOS_EXPECT_NE(hasher(a), hasher(c));
}

KRATOS_TEST_CASE_IN_SUITE(DenseBackendArray1dIterators, KratosCoreFastSuite)
{
    array_1d<double, 3> a{1.0, 2.0, 3.0};

    double sum = 0.0;
    for (const auto value : a) {
        sum += value;
    }
    KRATOS_EXPECT_EQ(sum, 6.0);

    for (auto& r_value : a) {
        r_value *= 2.0;
    }
    KRATOS_EXPECT_EQ(a[0], 2.0);
    KRATOS_EXPECT_EQ(a[2], 6.0);

    KRATOS_EXPECT_EQ(static_cast<std::size_t>(std::distance(a.begin(), a.end())), a.size());
}

KRATOS_TEST_CASE_IN_SUITE(DenseBackendBoundedTypesBasics, KratosCoreFastSuite)
{
    BoundedMatrix<double, 2, 3> m;
    KRATOS_EXPECT_EQ(m.size1(), 2);
    KRATOS_EXPECT_EQ(m.size2(), 3);

    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t j = 0; j < 3; ++j)
            m(i, j) = static_cast<double>(3 * i + j);

    // row-major contiguous storage (matches ublas bounded_matrix row_major)
    KRATOS_EXPECT_EQ(&m(0, 0) + 1, &m(0, 1));
    KRATOS_EXPECT_EQ(&m(0, 0) + 3, &m(1, 0));

    BoundedVector<double, 3> v;
    v[0] = 1.0; v(1) = 2.0; v[2] = 3.0;
    KRATOS_EXPECT_EQ(v.size(), 3);
    KRATOS_EXPECT_EQ(&v[0] + 2, &v[2]);

    // full-size resize is accepted and preserves
    m.resize(2, 3);
    KRATOS_EXPECT_EQ(m(1, 2), 5.0);
    v.resize(3);
    KRATOS_EXPECT_EQ(v[2], 3.0);
}

KRATOS_TEST_CASE_IN_SUITE(DenseBackendBoundedStreamFormat, KratosCoreFastSuite)
{
    BoundedVector<double, 3> v;
    v[0] = 1.0; v[1] = 2.5; v[2] = -3.0;
    std::stringstream vector_buffer;
    vector_buffer << v;
    KRATOS_EXPECT_EQ(vector_buffer.str(), "[3](1,2.5,-3)");

    BoundedMatrix<double, 2, 2> m;
    m(0, 0) = 1.0; m(0, 1) = 2.0;
    m(1, 0) = 3.0; m(1, 1) = 4.0;
    std::stringstream matrix_buffer;
    matrix_buffer << m;
    KRATOS_EXPECT_EQ(matrix_buffer.str(), "[2,2]((1,2),(3,4))");
}

KRATOS_TEST_CASE_IN_SUITE(DenseBackendArray1dNonArithmeticScalars, KratosCoreFastSuite)
{
    // array_1d must work as a plain container for non-arithmetic element
    // types under both backends (used in-tree with Vector, nested array_1d,
    // bool, int and size 0).
    array_1d<Vector, 2> vector_pair;
    vector_pair[0] = Vector(3, 1.0);
    vector_pair[1] = Vector(5, 2.0);
    KRATOS_EXPECT_EQ(vector_pair[0].size(), 3);
    KRATOS_EXPECT_EQ(vector_pair[1].size(), 5);
    array_1d<Vector, 2> vector_pair_copy(vector_pair);
    KRATOS_EXPECT_EQ(vector_pair_copy[1][0], 2.0);

    array_1d<array_1d<double, 3>, 4> nested;
    nested[2][1] = 42.0;
    KRATOS_EXPECT_EQ(nested[2][1], 42.0);
    const array_1d<array_1d<double, 3>, 4> nested_copy(nested);
    KRATOS_EXPECT_EQ(nested_copy[2][1], 42.0);

    array_1d<bool, 3> flags;
    flags[0] = true; flags[1] = false; flags[2] = true;
    KRATOS_EXPECT_TRUE(flags[0]);
    KRATOS_EXPECT_FALSE(flags[1]);

    array_1d<int, 4> integers(4, 7);
    KRATOS_EXPECT_EQ(integers[3], 7);

    array_1d<std::size_t, 3> sizes;
    sizes[0] = 1; sizes[1] = 2; sizes[2] = 3;
    KRATOS_EXPECT_EQ(sizes[2], 3);

    // size 0 must at least instantiate and iterate empty
    array_1d<double, 0> empty;
    KRATOS_EXPECT_EQ(empty.size(), 0);
    std::size_t count = 0;
    for (auto it = empty.begin(); it != empty.end(); ++it) ++count;
    KRATOS_EXPECT_EQ(count, 0);
}

KRATOS_TEST_CASE_IN_SUITE(DenseBackendAliasSafeAssignment, KratosCoreFastSuite)
{
    // operator= must be alias-safe (uBLAS evaluates through a protective
    // temporary; the in-tree VMS adjoint elements rely on it with
    // M = trans(M)). The temporary-free fast path is noalias().
    BoundedMatrix<double, 3, 3> M;
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j)
            M(i, j) = static_cast<double>(3 * i + j);

    M = trans(M);
    KRATOS_EXPECT_EQ(M(0, 1), 3.0);
    KRATOS_EXPECT_EQ(M(1, 0), 1.0);
    KRATOS_EXPECT_EQ(M(2, 0), 2.0);

    M = -trans(M);
    KRATOS_EXPECT_EQ(M(1, 0), -3.0);

    // self-referencing vector assignments
    array_1d<double, 3> a{1.0, 2.0, 3.0};
    a = a + a;
    KRATOS_EXPECT_EQ(a[2], 6.0);
    a -= 0.5 * a;
    KRATOS_EXPECT_EQ(a[2], 3.0);
}

KRATOS_TEST_CASE_IN_SUITE(DenseBackendMixedIdioms, KratosCoreFastSuite)
{
    // The ublas idiom surface that element/condition code uses on the
    // fixed-size types, including mixes with the (always-uBLAS) dynamic
    // Matrix/Vector. Every line must compile and give the same values under
    // both backends.
    BoundedMatrix<double, 3, 3> R;
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j)
            R(i, j) = (i == j) ? 2.0 : 0.0;

    array_1d<double, 3> x{1.0, 2.0, 3.0};

    // fixed x fixed
    array_1d<double, 3> y;
    noalias(y) = prod(R, x);
    KRATOS_EXPECT_EQ(y[0], 2.0);
    KRATOS_EXPECT_EQ(y[2], 6.0);

    y = prod(trans(R), x);
    KRATOS_EXPECT_EQ(y[1], 4.0);

    KRATOS_EXPECT_EQ(inner_prod(x, y), 1.0 * 2.0 + 2.0 * 4.0 + 3.0 * 6.0);
    KRATOS_EXPECT_NEAR(norm_2(x), std::sqrt(14.0), 1e-12);

    // dynamic x fixed mixes
    Matrix M = ZeroMatrix(3, 3);
    M(0, 0) = 1.0; M(1, 1) = 2.0; M(2, 2) = 3.0;

    array_1d<double, 3> z;
    noalias(z) = prod(M, x);
    KRATOS_EXPECT_EQ(z[0], 1.0);
    KRATOS_EXPECT_EQ(z[1], 4.0);
    KRATOS_EXPECT_EQ(z[2], 9.0);

    Vector v = prod(M, x);
    KRATOS_EXPECT_EQ(v.size(), 3);
    KRATOS_EXPECT_EQ(v[2], 9.0);

    Vector w(3);
    w[0] = 1.0; w[1] = 1.0; w[2] = 1.0;
    KRATOS_EXPECT_EQ(inner_prod(x, w), 6.0);
    KRATOS_EXPECT_EQ(inner_prod(w, x), 6.0);

    noalias(w) = prod(R, x);
    KRATOS_EXPECT_EQ(w[2], 6.0);

#ifdef KRATOS_USE_EIGEN_BACKEND
    // size1()/size2() are available on lazy Eigen expressions too (MatrixBase
    // plugin), matching the uBLAS expression surface.
    KRATOS_EXPECT_EQ(prod(M, M).size1(), 3);
    KRATOS_EXPECT_EQ(trans(M).size2(), 3);
#endif

    // row extraction from a dynamic matrix into an array_1d
    M(0, 1) = -1.0;
    array_1d<double, 3> first_row = row(M, 0);
    KRATOS_EXPECT_EQ(first_row[0], 1.0);
    KRATOS_EXPECT_EQ(first_row[1], -1.0);
    KRATOS_EXPECT_EQ(first_row[2], 0.0);

    // ZeroVector / ZeroMatrix absorption
    array_1d<double, 3> zeroed = ZeroVector(3);
    KRATOS_EXPECT_EQ(zeroed[0], 0.0);
    BoundedMatrix<double, 3, 3> zero_matrix = ZeroMatrix(3, 3);
    KRATOS_EXPECT_EQ(zero_matrix(1, 1), 0.0);
    BoundedMatrix<double, 3, 3> identity = IdentityMatrix(3);
    KRATOS_EXPECT_EQ(identity(1, 1), 1.0);
    KRATOS_EXPECT_EQ(identity(0, 1), 0.0);

    // arithmetic mixes
    array_1d<double, 3> combined = x + w;
    KRATOS_EXPECT_EQ(combined[2], 9.0);
    combined = x - w;
    KRATOS_EXPECT_EQ(combined[2], -3.0);
    combined = 2.0 * x;
    KRATOS_EXPECT_EQ(combined[2], 6.0);

    // outer product
    Matrix outer = outer_prod(x, w);
    KRATOS_EXPECT_EQ(outer.size1(), 3);
    KRATOS_EXPECT_EQ(outer.size2(), 3);
    KRATOS_EXPECT_EQ(outer(2, 2), 3.0 * 6.0);
}

} // namespace Kratos::Testing
