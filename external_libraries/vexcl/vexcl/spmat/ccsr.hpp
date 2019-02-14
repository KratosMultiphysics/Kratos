#ifndef VEXCL_SPMAT_CCSR_HPP
#define VEXCL_SPMAT_CCSR_HPP

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
 * \file   vexcl/spmat/ccsr.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Sparse matrix in Compressed CSR format.
 */

namespace vex {

/// Sparse matrix in CCSR format.
/**
 * Compressed CSR format. row, col, and val arrays contain unique rows of the
 * matrix. Column numbers in col array are relative to diagonal. idx array
 * contains index into row vector, corresponding to each row of the matrix. So
 * that matrix-vector multiplication may be performed as follows:
 \code
 for(unsigned i = 0; i < n; i++) {
     val_t sum = 0;
     for(unsigned j = row[idx[i]]; j < row[idx[i] + 1]; j++)
         sum += val[j] * x[i + col[j]];
     y[i] = sum;
 }
 \endcode
 * This format does not support multi-device computation, so it accepts single
 * queue at initialization. Vectors x and y should also be single-queued and
 * reside on the same device with matrix.
 */
template <typename val_t, typename col_t = ptrdiff_t, typename idx_t = size_t>
struct SpMatCCSR {
    static_assert(std::is_signed<col_t>::value,
            "Column type for CCSR format has to be signed.");

    /// Constructor for CCSR format.
    /**
     * Constructs GPU representation of the CCSR matrix.
     * \param queue single queue.
     * \param n     number of rows in the matrix.
     * \param m     number of unique rows in the matrix.
     * \param idx   index into row vector.
     * \param row   row index into col and val vectors.
     * \param col   column positions of nonzero elements wrt to diagonal.
     * \param val   values of nonzero elements of the matrix.
     */
    SpMatCCSR(const backend::command_queue &queue, size_t n, size_t m,
            const idx_t *idx, const idx_t *row, const col_t *col, const val_t *val
            )
        : queue(queue), n(n),
          idx(queue, n,      idx, backend::MEM_READ_ONLY),
          row(queue, (m+1),  row, backend::MEM_READ_ONLY),
          col(queue, row[m], col, backend::MEM_READ_ONLY),
          val(queue, row[m], val, backend::MEM_READ_ONLY)
    { }

    backend::command_queue queue;
    size_t n;
    backend::device_vector<idx_t> idx;
    backend::device_vector<idx_t> row;
    backend::device_vector<col_t> col;
    backend::device_vector<val_t> val;
};

struct ccsr_product_terminal {};

typedef vector_expression<
    typename boost::proto::terminal< ccsr_product_terminal >::type
    > ccsr_product_terminal_expression;

template <typename val_t, typename col_t, typename idx_t, typename T>
struct ccsr_product : public ccsr_product_terminal_expression
{
    typedef val_t value_type;
    typedef SpMatCCSR<val_t, col_t, idx_t> matrix;

    const matrix    &A;
    const vector<T> &x;

    ccsr_product(const matrix &A, const vector<T> &x) : A(A), x(x) {}
};

template <typename val_t, typename col_t, typename idx_t, typename T>
ccsr_product<val_t, col_t, idx_t, T> operator*(
        const SpMatCCSR<val_t, col_t, idx_t> &A,
        const vector<T> &x)
{
    return ccsr_product<val_t, col_t, idx_t, T>(A, x);
}


#ifdef VEXCL_MULTIVECTOR_HPP
struct mv_ccsr_product_terminal {};

typedef multivector_expression<
    typename boost::proto::terminal< mv_ccsr_product_terminal >::type
    > mv_ccsr_product_terminal_expression;

template <typename val_t, typename col_t, typename idx_t, class MV>
struct mv_ccsr_product : public mv_ccsr_product_terminal_expression
{
    typedef SpMatCCSR<val_t, col_t, idx_t> matrix;

    const matrix &A;
    const MV     &x;

    mv_ccsr_product(const matrix &A, const MV &x) : A(A), x(x) {}
};

template <typename val_t, typename col_t, typename idx_t, class MV>
typename std::enable_if<
    std::is_base_of<multivector_terminal_expression, MV>::value,
    mv_ccsr_product<val_t, col_t, idx_t, MV>
>::type
operator*(
        const SpMatCCSR<val_t, col_t, idx_t> &A,
        const MV &x)
{
    return mv_ccsr_product<val_t, col_t, idx_t, MV>(A, x);
}
#endif


// Allow ccsr_product to participate in vector expressions:
namespace traits {

template <>
struct is_vector_expr_terminal< ccsr_product_terminal > : std::true_type {};

template <>
struct proto_terminal_is_value< ccsr_product_terminal > : std::true_type {};

#ifdef VEXCL_MULTIVECTOR_HPP
template <>
struct proto_terminal_is_value< mv_ccsr_product_terminal >
    : std::true_type
{ };

template <>
struct is_multivector_expr_terminal< mv_ccsr_product_terminal >
    : std::true_type
{ };

template <size_t I, typename val_t, typename col_t, typename idx_t, typename MV>
struct component< I, mv_ccsr_product<val_t, col_t, idx_t, MV> > {
    typedef
        ccsr_product<val_t, col_t, idx_t, typename MV::sub_value_type>
        type;
};
#endif

template <typename val_t, typename col_t, typename idx_t, typename T>
struct terminal_preamble< ccsr_product<val_t, col_t, idx_t, T> > {
    static void get(backend::source_generator &src,
            const ccsr_product<val_t, col_t, idx_t, T>&,
            const backend::command_queue&, const std::string &prm_name,
            detail::kernel_generator_state_ptr)
    {
        typedef decltype(val_t() * T()) res_t;

        src.begin_function<res_t>(prm_name + "_spmv");
        src.begin_function_parameters();
        src.template parameter< global_ptr<const idx_t> >("idx");
        src.template parameter< global_ptr<const idx_t> >("row");
        src.template parameter< global_ptr<const col_t> >("col");
        src.template parameter< global_ptr<const val_t> >("val");
        src.template parameter< global_ptr<const T>     >("vec");
        src.template parameter< size_t >("i");
        src.end_function_parameters();

        src.new_line() << type_name<res_t>() << " sum = 0;";
        src.new_line() << "for(size_t pos = idx[i], j = row[pos], end = row[pos+1]; j < end; ++j)";
        src.open("{");
        src.new_line() << "sum += val[j] * vec[i + col[j]];";
        src.close("}");
        src.new_line() << "return sum;";
        src.end_function();
    }
};

template <typename val_t, typename col_t, typename idx_t, typename T>
struct kernel_param_declaration< ccsr_product<val_t, col_t, idx_t, T> > {
    static void get(backend::source_generator &src,
            const ccsr_product<val_t, col_t, idx_t, T>&,
            const backend::command_queue&, const std::string &prm_name,
            detail::kernel_generator_state_ptr)
    {
        src.template parameter< global_ptr<const idx_t> >(prm_name + "_idx");
        src.template parameter< global_ptr<const idx_t> >(prm_name + "_row");
        src.template parameter< global_ptr<const col_t> >(prm_name + "_col");
        src.template parameter< global_ptr<const val_t> >(prm_name + "_val");
        src.template parameter< global_ptr<const T    > >(prm_name + "_vec");
    }
};

template <typename val_t, typename col_t, typename idx_t, typename T>
struct partial_vector_expr< ccsr_product<val_t, col_t, idx_t, T> > {
    static void get(backend::source_generator &src,
            const ccsr_product<val_t, col_t, idx_t, T>&,
            const backend::command_queue&, const std::string &prm_name,
            detail::kernel_generator_state_ptr)
    {
        src << prm_name << "_spmv" << "("
            << prm_name << "_idx, "
            << prm_name << "_row, "
            << prm_name << "_col, "
            << prm_name << "_val, "
            << prm_name << "_vec, idx)";
    }
};

template <typename val_t, typename col_t, typename idx_t, typename T>
struct kernel_arg_setter< ccsr_product<val_t, col_t, idx_t, T> > {
    static void set(const ccsr_product<val_t, col_t, idx_t, T> &term,
            backend::kernel &kernel, unsigned part, size_t/*index_offset*/,
            detail::kernel_generator_state_ptr)
    {
        assert(part == 0);

        kernel.push_arg(term.A.idx);
        kernel.push_arg(term.A.row);
        kernel.push_arg(term.A.col);
        kernel.push_arg(term.A.val);
        kernel.push_arg(term.x(part));
    }
};

template <typename val_t, typename col_t, typename idx_t, typename T>
struct expression_properties< ccsr_product<val_t, col_t, idx_t, T> > {
    static void get(const ccsr_product<val_t, col_t, idx_t, T> &term,
            std::vector<backend::command_queue> &queue_list,
            std::vector<size_t> &partition,
            size_t &size
            )
    {
        queue_list = term.x.queue_list();
        partition  = term.x.partition();
        size       = term.x.size();

        assert(partition.size() == 2);
        partition.back() = size;
    }
};

} // namespace traits

#ifdef VEXCL_MULTIVECTOR_HPP
template <size_t I, typename val_t, typename col_t, typename idx_t, typename MV>
ccsr_product<val_t, col_t, idx_t, typename MV::sub_value_type>
get(const mv_ccsr_product<val_t, col_t, idx_t, MV> &t) {
    return t.A * t.x(I);
}
#endif

} // namespace vex

#endif
