#ifndef VEXCL_SPMAT_INLINE_SPMV_HPP
#define VEXCL_SPMAT_INLINE_SPMV_HPP

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
 * \file   vexcl/spmat/inline_spmv.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Inline SpMV operation.
 */

namespace vex {

struct inline_spmv_terminal {};

typedef vector_expression<
    typename boost::proto::terminal< inline_spmv_terminal >::type
    > inline_spmv_terminal_expression;

template <class M, class V>
struct inline_spmv : inline_spmv_terminal_expression {
    typedef typename V::value_type value_type;

    const M &A;
    const V &x;

    inline_spmv(const M &A, const V &x) : A(A), x(x) {}
};

/// Inlines a sparse matrix - vector product.
/**
 * When applied to a matrix-vector product, the product becomes inlineable.
 * That is, it may be used in any vector expression (not just additive
 * expressions). The user has to guarantee the function is only used in single-device expressions.
 *
 * Example:
 \code
 // Get maximum residual value:
 eps = sum( fabs(f - vex::make_inline(A * x)) );
 \endcode
 */
#ifdef DOXYGEN
template <class MVProdExpr>
auto make_inline(const MVProdExpr &expr);
#else
template <typename val_t, typename col_t, typename idx_t>
inline_spmv< SpMat<val_t, col_t, idx_t>, vector<val_t> >
make_inline(const additive_operator< SpMat<val_t, col_t, idx_t>, vector<val_t> > &base)
#endif
{
    precondition(base.x.nparts() == 1, "Can not inline multi-device SpMV operation.");

    return inline_spmv< SpMat<val_t, col_t, idx_t>, vector<val_t> >(base.A, base.x);
}

#ifdef VEXCL_MULTIVECTOR_HPP
struct mv_inline_spmv_terminal {};

typedef multivector_expression<
    typename boost::proto::terminal< mv_inline_spmv_terminal >::type
    > mv_inline_spmv_terminal_expression;

template <class M, class V>
struct mv_inline_spmv : mv_inline_spmv_terminal_expression {
    const M &A;
    const V &x;

    mv_inline_spmv(const M &A, const V &x) : A(A), x(x) {}
};

/// Inlines a sparse matrix - multivector product.
/**
 * When applied to a matrix-multivector product, the product becomes
 * inlineable.  That is, it may be used in any multivector expression (not just
 * additive expression). This is only possible in single-device contexts, so
 * user has to guarantee that.
 *
 * Example:
 \code
 // Get maximum residual value:
 eps = sum( fabs(f - vex::make_inline(A * x)) );
 \endcode
 */
template <typename val_t, typename col_t, typename idx_t, class V>
mv_inline_spmv<SpMat<val_t, col_t, idx_t>, V>
make_inline(const multiadditive_operator<SpMat<val_t, col_t, idx_t>, V> &base) {
    precondition(base.x(0).nparts() == 1, "Can not inline multi-device SpMV operation.");

    return mv_inline_spmv<SpMat<val_t, col_t, idx_t>, V>(base.A, base.x);
}
#endif

// Allow inline products to participate in vector expressions:
namespace traits {

template <>
struct is_vector_expr_terminal< inline_spmv_terminal > : std::true_type {};

template <>
struct proto_terminal_is_value< inline_spmv_terminal > : std::true_type {};

#ifdef VEXCL_MULTIVECTOR_HPP
template <>
struct is_multivector_expr_terminal< mv_inline_spmv_terminal >
    : std::true_type
{ };

template <>
struct proto_terminal_is_value< mv_inline_spmv_terminal >
    : std::true_type
{ };

template <size_t I, class M, class V>
struct component< I, mv_inline_spmv<M, V> > {
    typedef inline_spmv< M, vector<typename V::sub_value_type> > type;
};
#endif

template <class M, class V>
struct terminal_preamble< inline_spmv<M, V> > {
    static void get(backend::source_generator &src,
            const inline_spmv<M, V>&,
            const backend::command_queue &queue, const std::string &prm_name,
            detail::kernel_generator_state_ptr state)
    {
        M::inline_preamble(src, queue, prm_name, state);
    }
};

template <class M, class V>
struct kernel_param_declaration< inline_spmv<M, V> > {
    static void get(backend::source_generator &src,
            const inline_spmv<M, V>&,
            const backend::command_queue &queue, const std::string &prm_name,
            detail::kernel_generator_state_ptr state)
    {
        M::inline_parameters(src, queue, prm_name, state);
    }
};

template <class M, class V>
struct partial_vector_expr< inline_spmv<M, V> > {
    static void get(backend::source_generator &src,
            const inline_spmv<M, V>&,
            const backend::command_queue &queue, const std::string &prm_name,
            detail::kernel_generator_state_ptr state)
    {
        M::inline_expression(src, queue, prm_name, state);
    }
};

template <class M, class V>
struct kernel_arg_setter< inline_spmv<M, V> > {
    static void set(const inline_spmv<M, V> &term,
            backend::kernel &kernel, unsigned part, size_t index_offset,
            detail::kernel_generator_state_ptr state)
    {
        M::inline_arguments(
                kernel, part, index_offset, term.A, term.x, state
                );
    }
};

template <class M, class V>
struct expression_properties< inline_spmv<M, V> > {
    static void get(const inline_spmv<M, V> &term,
            std::vector<backend::command_queue> &queue_list,
            std::vector<size_t> &partition,
            size_t &size
            )
    {
        extract_expression_properties(term.x, queue_list, partition, size);
    }
};

} // namespace traits

#ifdef VEXCL_MULTIVECTOR_HPP
template <size_t I, class M, class V>
inline_spmv<M, vector<typename V::sub_value_type> >
get(const mv_inline_spmv<M, V> &t) {
    return make_inline(t.A * t.x(I));
}
#endif

} // namespace vex

#endif
