#ifndef AMGCL_ADAPTER_EIGEN_HPP
#define AMGCL_ADAPTER_EIGEN_HPP

/*
The MIT License

Copyright (c) 2012-2019 Denis Demidov <dennis.demidov@gmail.com>

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
\file    amgcl/adapter/eigen.hpp
\author  Denis Demidov <dennis.demidov@gmail.com>
\brief   Adapters for Eigen types to be used with builtin backend.
\ingroup adapters
*/

#include <type_traits>
#include <Eigen/SparseCore>
#include <amgcl/util.hpp>
#include <amgcl/backend/builtin.hpp>

namespace amgcl {
namespace backend {

//---------------------------------------------------------------------------
// Backend interface specialization for Eigen types
//---------------------------------------------------------------------------
template <class T, class Enable = void>
struct is_eigen_sparse_matrix : std::false_type {};

template <class T, class Enable = void>
struct is_eigen_type : std::false_type {};

template <typename Scalar, int Flags, typename Storage>
struct is_eigen_sparse_matrix<
    Eigen::MappedSparseMatrix<Scalar, Flags, Storage>
    > : std::true_type
{};

template <typename Scalar, int Flags, typename Storage>
struct is_eigen_sparse_matrix<
    Eigen::SparseMatrix<Scalar, Flags, Storage>
    > : std::true_type
{};

template <class T>
struct is_eigen_type<
    T,
    typename std::enable_if<
        std::is_arithmetic<typename T::Scalar>::value &&
        std::is_base_of<Eigen::EigenBase<T>, T>::value
        >::type
    > : std::true_type
{};

template <class T>
struct value_type<
    T,
    typename std::enable_if<is_eigen_type<T>::value>::type
    >
{
    typedef typename T::Scalar type;
};

template <class T>
struct rows_impl<
    T,
    typename std::enable_if<is_eigen_sparse_matrix<T>::value>::type
    >
{
    static size_t get(const T &matrix) {
        return matrix.rows();
    }
};

template <class T>
struct cols_impl<
    T,
    typename std::enable_if<is_eigen_sparse_matrix<T>::value>::type
    >
{
    static size_t get(const T &matrix) {
        return matrix.cols();
    }
};

template <class T>
struct nonzeros_impl<
    T,
    typename std::enable_if<is_eigen_sparse_matrix<T>::value>::type
    >
{
    static size_t get(const T &matrix) {
        return matrix.nonZeros();
    }
};

template <class T>
struct row_iterator <
    T,
    typename std::enable_if<is_eigen_sparse_matrix<T>::value>::type
    >
{
    typedef typename T::InnerIterator type;
};

template <class T>
struct row_begin_impl <
    T,
    typename std::enable_if<is_eigen_sparse_matrix<T>::value>::type
    >
{
    typedef typename row_iterator<T>::type iterator;
    static iterator get(const T &matrix, size_t row) {
        return iterator(matrix, row);
    }
};

} // namespace backend
} // namespace amgcl

#define AMGCL_USE_EIGEN_VECTORS_WITH_BUILTIN_BACKEND()                         \
    namespace amgcl { namespace backend {                                      \
        template <class T >                                                    \
        struct is_builtin_vector< Eigen::Matrix<T, Eigen::Dynamic, 1> >        \
          : std::true_type {};                                                 \
    } }

#endif
