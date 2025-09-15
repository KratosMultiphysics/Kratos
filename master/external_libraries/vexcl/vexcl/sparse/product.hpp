#ifndef VEXCL_SPARSE_PRODUCT_HPP
#define VEXCL_SPARSE_PRODUCT_HPP

/*
The MIT License

Copyright (c) 2012-2018 Denis Demidov <dennis.demidov@gmail.com>

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
 * \file   vexcl/sparse/product.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Sparse matrix-vector product terminal.
 */

#include <vexcl/operations.hpp>

namespace vex {
namespace sparse {

struct matrix_vector_product_terminal {};

typedef vector_expression<
        typename boost::proto::terminal< matrix_vector_product_terminal >::type
        > matrix_vector_product_expression;

template <class Matrix, class Vector>
struct matrix_vector_product : matrix_vector_product_expression
{
    const Matrix &A;
    typename boost::proto::result_of::as_child<const Vector, vector_domain>::type x;

    matrix_vector_product(const Matrix &A, const Vector &x)
        : A(A), x(boost::proto::as_child(x)) {}
};

} // namespace sparse

namespace traits {

template <> struct is_vector_expr_terminal< sparse::matrix_vector_product_terminal >
    : std::true_type {};

template <class Matrix, class Vector>
struct terminal_is_value< sparse::matrix_vector_product<Matrix, Vector> >
    : std::true_type {};

template <class Matrix, class Vector>
struct terminal_preamble< sparse::matrix_vector_product<Matrix, Vector> > {
    static void get(backend::source_generator &src,
            const sparse::matrix_vector_product<Matrix, Vector> &term,
            const backend::command_queue &q, const std::string &prm_name,
            detail::kernel_generator_state_ptr state)
    {
        term.A.terminal_preamble(term.x, src, q, prm_name, state);
    }
};

template <class Matrix, class Vector>
struct local_terminal_init< sparse::matrix_vector_product<Matrix, Vector> > {
    static void get(backend::source_generator &src,
            const sparse::matrix_vector_product<Matrix, Vector> &term,
            const backend::command_queue &q, const std::string &prm_name,
            detail::kernel_generator_state_ptr state)
    {
        term.A.local_terminal_init(term.x, src, q, prm_name, state);
    }
};

template <class Matrix, class Vector>
struct kernel_param_declaration< sparse::matrix_vector_product<Matrix, Vector> > {
    static void get(backend::source_generator &src,
            const sparse::matrix_vector_product<Matrix, Vector> &term,
            const backend::command_queue &q, const std::string &prm_name,
            detail::kernel_generator_state_ptr state)
    {
        term.A.kernel_param_declaration(term.x, src, q, prm_name, state);
    }
};

template <class Matrix, class Vector>
struct partial_vector_expr< sparse::matrix_vector_product<Matrix, Vector> > {
    static void get(backend::source_generator &src,
            const sparse::matrix_vector_product<Matrix, Vector> &term,
            const backend::command_queue &q, const std::string &prm_name,
            detail::kernel_generator_state_ptr state)
    {
        term.A.partial_vector_expr(term.x, src, q, prm_name, state);
    }
};

template <class Matrix, class Vector>
struct kernel_arg_setter< sparse::matrix_vector_product<Matrix, Vector> > {
    static void set(const sparse::matrix_vector_product<Matrix, Vector> &term,
            backend::kernel &kernel, unsigned part, size_t index_offset,
            detail::kernel_generator_state_ptr state)
    {
        term.A.kernel_arg_setter(term.x, kernel, part, index_offset, state);
    }
};

template <class Matrix, class Vector>
struct expression_properties< sparse::matrix_vector_product<Matrix, Vector> > {
    static void get(const sparse::matrix_vector_product<Matrix, Vector> &term,
            std::vector<backend::command_queue> &queue_list,
            std::vector<size_t> &partition,
            size_t &size
            )
    {
        term.A.expression_properties(term.x, queue_list, partition, size);
    }
};

} // namespace traits

} // namespace vex

#endif
