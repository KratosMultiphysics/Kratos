#ifndef AMGCL_BACKEND_INTERFACE_HPP
#define AMGCL_BACKEND_INTERFACE_HPP

/*
The MIT License

Copyright (c) 2012-2016 Denis Demidov <dennis.demidov@gmail.com>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/

/**
 * \file   amgcl/backend/interface.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Backend interface required for AMG.
 */

#include <cmath>

#include <boost/static_assert.hpp>
#include <boost/type_traits.hpp>

#include <amgcl/value_type/interface.hpp>
#include <amgcl/util.hpp>

namespace amgcl {

/// Provided backends.
namespace backend {

/**
 * \defgroup backends Provided backends
 * \brief Backends implemented in AMGCL.
 *
 * A backend in AMGCL is a class that defines matrix and vector types together
 * with several operations on them, such as creation, matrix-vector products,
 * vector sums, inner products etc.  The AMG hierarchy is moved to the
 * specified backend upon construction. The solution phase then uses types and
 * operations defined in the backend. This enables transparent acceleration of
 * the solution phase with OpenMP, OpenCL, CUDA, or any other technologies.
 */

/**
 * \defgroup backend_interface Backend interface
 * \brief Backend interface specification.
 *
 * One has to specify these templates in order to define a new backend.
 */

/** \addtogroup backend_interface
 * @{
 */

/// Metafunction that returns value type of a matrix or a vector type.
template <class T, class Enable = void>
struct value_type {
    typedef typename T::value_type type;
};

/// Implementation for function returning the number of rows in a matrix.
/** \note Used in rows() */
template <class Matrix, class Enable = void>
struct rows_impl {
    typedef typename Matrix::ROWS_NOT_IMPLEMENTED type;
};

/// Implementation for function returning the number of columns in a matrix.
/** \note Used in cols() */
template <class Matrix, class Enable = void>
struct cols_impl {
    typedef typename Matrix::COLS_NOT_IMPLEMENTED type;
};

template <class Matrix, class Enable = void>
struct ptr_data_impl {
    typedef typename Matrix::PTR_DATA_NOT_IMPLEMENTED type;
};

template <class Matrix, class Enable = void>
struct col_data_impl {
    typedef typename Matrix::COL_DATA_NOT_IMPLEMENTED type;
};

template <class Matrix, class Enable = void>
struct val_data_impl {
    typedef typename Matrix::VAL_DATA_NOT_IMPLEMENTED type;
};

/// Implementation for function returning the number of nonzeros in a matrix.
/** \note Used in nonzeros() */
template <class Matrix, class Enable = void>
struct nonzeros_impl {
    typedef typename Matrix::NONZEROS_NOT_IMPLEMENTED type;
};

/// Implementation for function returning the number of nonzeros in a matrix row.
/** \note Used in row_nonzeros() */
template <class Matrix, class Enable = void>
struct row_nonzeros_impl {
    typedef typename Matrix::ROW_NONZEROS_NOT_IMPLEMENTED type;
};

/// Metafunction returning the row iterator type for a matrix type.
/**
 * \note This only has to be implemented in the backend if support for serial
 * smoothers (Gauss-Seidel or ILU0) is required.
 */
template <class Matrix, class Enable = void>
struct row_iterator {
    typedef typename Matrix::ROW_ITERATOR_NOT_IMPLEMENTED type;
};

/// Implementation for function returning row iterator for a matrix.
/**
 * \note This only has to be implemented in the backend if support for serial
 * smoothers (Gauss-Seidel or ILU0) is required.
 * \note Used in row_begin()
 */
template <class Matrix, class Enable = void>
struct row_begin_impl {
    typedef typename Matrix::ROW_BEGIN_NOT_IMPLEMENTED type;
};

/// Implementation for matrix-vector product.
/** \note Used in spmv() */
template <class Alpha, class Matrix, class Vector1, class Beta, class Vector2, class Enable = void>
struct spmv_impl {
    typedef typename Matrix::SPMV_NOT_IMPLEMENTED type;
};

/// Implementation for residual error compuatation.
/** \note Used in residual() */
template <class Matrix, class Vector1, class Vector2, class Vector3, class Enable = void>
struct residual_impl {
    typedef typename Matrix::RESIDUAL_NOT_IMPLEMENTED type;
};

/// Implementation for zeroing out a vector.
/** \note Used in clear() */
template <class Vector, class Enable = void>
struct clear_impl {
    typedef typename Vector::CLEAR_NOT_IMPLEMENTED type;
};

/// Implementation for vector copy.
/** \note Used in copy() */
template <class Vector1, class Vector2, class Enable = void>
struct copy_impl {
    typedef typename Vector1::COPY_NOT_IMPLEMENTED type;
};

/// Implementation for copying data to backend.
/** \note Used in copy_to_backend() */
template <class Vector, class Enable = void>
struct copy_to_backend_impl {
    typedef typename Vector::COPY_TO_BACKEND_NOT_IMPLEMENTED type;
};

/// Implementation for inner product.
/** \note Used in inner_product() */
template <class Vector1, class Vector2, class Enable = void>
struct inner_product_impl {
    typedef typename Vector1::INNER_PRODUCT_NOT_IMPLEMENTED type;
};

/// Implementation for linear combination of two vectors.
/** \note Used in axpby() */
template <class A, class Vector1, class B, class Vector2, class Enable = void>
struct axpby_impl {
    typedef typename Vector1::AXPBY_NOT_IMPLEMENTED type;
};

/// Implementation for linear combination of three vectors.
/** \note Used in axpbypcz() */
template <class A, class Vector1, class B, class Vector2, class C, class Vector3, class Enable = void>
struct axpbypcz_impl {
    typedef typename Vector1::AXPBYPCZ_NOT_IMPLEMENTED type;
};

/// Implementation for element-wize vector product.
/** \note Used in vmul() */
template <class Alpha, class Vector1, class Vector2, class Beta, class Vector3, class Enable = void>
struct vmul_impl {
    typedef typename Vector1::VMUL_NOT_IMPLEMENTED type;
};

/** @} */

/// Returns the number of rows in a matrix.
template <class Matrix>
size_t rows(const Matrix &matrix) {
    return rows_impl<Matrix>::get(matrix);
}

/// Returns the number of columns in a matrix.
template <class Matrix>
size_t cols(const Matrix &matrix) {
    return cols_impl<Matrix>::get(matrix);
}

template <class Matrix>
typename ptr_data_impl<Matrix>::type
ptr_data(const Matrix &matrix) {
    return ptr_data_impl<Matrix>::get(matrix);
}

template <class Matrix>
typename col_data_impl<Matrix>::type
col_data(const Matrix &matrix) {
    return col_data_impl<Matrix>::get(matrix);
}

template <class Matrix>
typename val_data_impl<Matrix>::type
val_data(const Matrix &matrix) {
    return val_data_impl<Matrix>::get(matrix);
}

/// Returns the number of nonzeros in a matrix.
template <class Matrix>
size_t nonzeros(const Matrix &matrix) {
    return nonzeros_impl<Matrix>::get(matrix);
}

/// Returns row iterator for a matrix.
template <class Matrix>
typename row_iterator<Matrix>::type
row_begin(const Matrix &matrix, size_t row) {
    return row_begin_impl<Matrix>::get(matrix, row);
}

/// Returns number of nonzeros in a matrix row.
template <class Matrix>
size_t row_nonzeros(const Matrix &A, size_t row) {
    return row_nonzeros_impl<Matrix>::get(A, row);
}

/// Performs matrix-vector product.
/**
 * \f[y = \alpha A x + \beta y.\f]
 */
template <class Alpha, class Matrix, class Vector1, class Beta, class Vector2>
void spmv(
        Alpha alpha,
        const Matrix &A,
        const Vector1 &x,
        Beta beta,
        Vector2 &y)
{
    spmv_impl<Alpha, Matrix, Vector1, Beta, Vector2>::apply(alpha, A, x, beta, y);
}

/// Computes residual error.
/**
 * \f[r = rhs - Ax.\f]
 */
template <class Matrix, class Vector1, class Vector2, class Vector3>
void residual(const Vector1 &rhs, const Matrix &A, const Vector2 &x, Vector3 &r)
{
    residual_impl<Matrix, Vector1, Vector2, Vector3>::apply(rhs, A, x, r);
}

/// Zeros out a vector.
template <class Vector>
void clear(Vector &x)
{
    clear_impl<Vector>::apply(x);
}

/// Vector copy.
template <class Vector1, class Vector2>
void copy(const Vector1 &x, Vector2 &y)
{
    copy_impl<Vector1, Vector2>::apply(x, y);
}

/// Copy data to backend.
template <class Vector>
void copy_to_backend(const std::vector<typename value_type<Vector>::type> &data, Vector &x)
{
    copy_to_backend_impl<Vector>::apply(data, x);
}

/// Computes inner product of two vectors.
template <class Vector1, class Vector2>
typename math::inner_product_impl<
    typename value_type<Vector1>::type
    >::return_type
inner_product(const Vector1 &x, const Vector2 &y)
{
    return inner_product_impl<Vector1, Vector2>::get(x, y);
}

/// Computes linear combination of two vectors.
/**
 * \f[y = ax + by.\f]
 */
template <class A, class Vector1, class B, class Vector2>
void axpby(A a, Vector1 const &x, B b, Vector2 &y) {
    axpby_impl<A, Vector1, B, Vector2>::apply(a, x, b, y);
}

/// Computes linear combination of three vectors.
/**
 * \f[z = ax + by + cz.\f]
 */
template <class A, class Vector1, class B, class Vector2, class C, class Vector3>
void axpbypcz(A a, Vector1 const &x, B b, Vector2 const &y, C c, Vector3 &z) {
    axpbypcz_impl<A, Vector1, B, Vector2, C, Vector3>::apply(a, x, b, y, c, z);
}

/// Computes element-wize vector product.
/**
 * \f[z = \alpha xy + \beta z.\f]
 */
template <class Alpha, class Vector1, class Vector2, class Beta, class Vector3>
void vmul(Alpha alpha, const Vector1 &x, const Vector2 &y, Beta beta, Vector3 &z)
{
    vmul_impl<Alpha, Vector1, Vector2, Beta, Vector3>::apply(alpha, x, y, beta, z);
}

/// Is the relaxation supported by the backend?
template <class Backend, template <class> class Relaxation, class Enable = void>
struct relaxation_is_supported : boost::true_type {};

/// Is the coarsening supported by the backend?
template <class Backend, class Coarsening, class Enable = void>
struct coarsening_is_supported : boost::true_type {};

} // namespace backend
} // namespace amgcl


#endif
