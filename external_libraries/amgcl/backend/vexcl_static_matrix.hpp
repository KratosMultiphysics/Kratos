#ifndef AMGCL_BACKEND_VEXCL_STATIC_MATRIX_HPP
#define AMGCL_BACKEND_VEXCL_STATIC_MATRIX_HPP

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
 * \file   amgcl/backend/vexcl_static_matrix.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Static matrix support for the VexCL backend.
 */

#include <amgcl/backend/vexcl.hpp>
#include <amgcl/value_type/static_matrix.hpp>

namespace vex {

template <typename T, int N, int M>
struct is_cl_native< amgcl::static_matrix<T, N, M> > : std::true_type {};

template <typename T, int N, int M>
struct type_name_impl< amgcl::static_matrix<T, N, M> >
{
    static std::string get() {
        std::ostringstream s;
        s << "amgcl_matrix_" << type_name<T>() << "_" << N << "x" << M;
        return s.str();
    }
};

template <typename T, int N, int M>
struct cl_scalar_of< amgcl::static_matrix<T, N, M> > {
    typedef T type;
};

namespace sparse {

template <typename T, int N>
struct rhs_of< amgcl::static_matrix<T, N, N> > {
    typedef amgcl::static_matrix<T, N, 1> type;
};

template <typename T, int N>
struct spmv_ops_impl<amgcl::static_matrix<T,N,N>, amgcl::static_matrix<T,N,1>> {
    typedef amgcl::static_matrix<T,N,N> matrix_value;
    typedef amgcl::static_matrix<T,N,1> vector_value;

    static void decl_accum_var(backend::source_generator &src, const std::string &name)
    {
        src.new_line() << type_name<vector_value>() << " " << name << ";";
        for(int i = 0; i < N; ++i) {
            src.new_line() << name << ".data[" << i << "][0] = 0;";
        }
    }

    static void append(backend::source_generator &src,
            const std::string &sum, const std::string &val)
    {
        for(int i = 0; i < N; ++i)
            src.new_line() << sum << ".data[" << i << "][0] += " << val << ".data[" << i << "][0];";
    }

    static void append_product(backend::source_generator &src,
            const std::string &sum, const std::string &mat_val, const std::string &vec_val)
    {
        src.open("{");
        src.new_line() << type_name<vector_value>() << " v = " << vec_val << ";";
        for(int i = 0; i < N; ++i) {
            src.new_line() << sum << ".data[" << i << "][0] += ";
            for(int j = 0; j < N; ++j) {
                if (j) src << " + ";
                src << mat_val << ".data[" << i << "][" << j << "] * v.data[" << j << "][0]";
            }
            src << ";";
        }
        src.close("}");
    }
};

template <typename T, int N, typename Col, typename Ptr>
class ell<amgcl::static_matrix<T, N, N>, Col, Ptr> {
    public:
        typedef amgcl::static_matrix<T, N, N> Val;
        typedef Val value_type;

        typedef Val val_type;
        typedef Col col_type;
        typedef Ptr ptr_type;

        template <class PtrRange, class ColRange, class ValRange>
        ell(
                const std::vector<backend::command_queue> &q,
                size_t nrows, size_t ncols,
                const PtrRange &ptr,
                const ColRange &col,
                const ValRange &val,
                bool fast = true
           ) : q(q[0]), n(nrows), m(ncols),
               nnz(std::distance(std::begin(val), std::end(val))),
               ell_pitch(alignup(nrows, 16U)), csr_nnz(0)
        {
            precondition(q.size() == 1,
                    "sparse::ell is only supported for single-device contexts");

            if (fast) {
                convert(ptr, col, val);
                return;
            }

            /* 1. Get optimal ELL widths for local and remote parts. */
            // Speed of ELL relative to CSR:
            const double ell_vs_csr = 3.0;

            // Find maximum widths for local and remote parts:
            size_t max_width = 0;
            for(size_t i = 0; i < n; ++i)
                max_width = std::max<size_t>(max_width, ptr[i+1] - ptr[i]);

            // Build width distribution histogram.
            std::vector<Ptr> hist(max_width + 1, 0);
            for(size_t i = 0; i < n; ++i)
                ++hist[ptr[i+1] - ptr[i]];

            // Estimate optimal width for ELL part of the matrix.
            ell_width = max_width;
            for(size_t i = 0, rows = n; i < max_width; ++i) {
                rows -= hist[i]; // Number of rows wider than i.
                if (ell_vs_csr * rows < n) {
                    ell_width = i;
                    break;
                }
            }

            if (ell_width == 0) {
                assert(csr_nnz == nnz);

                csr_ptr = backend::device_vector<Col>(q[0], n + 1,   &ptr[0]);
                csr_col = backend::device_vector<Col>(q[0], csr_nnz, &col[0]);
                csr_val = create_device_vector       (q[0], csr_nnz, &val[0], false);

                return;
            }

            size_t ell_nnz = ell_pitch * ell_width;

            // Count nonzeros in CSR part of the matrix.
            for(size_t i = ell_width + 1; i <= max_width; ++i)
                csr_nnz += hist[i] * (i - ell_width);

            /* 3. Split the input matrix into ELL and CSR submatrices. */
            std::vector<Col> _ell_col(ell_nnz, static_cast<Col>(-1));
            std::vector<T>   _ell_val(ell_nnz * N * N);
            std::vector<Ptr> _csr_ptr;
            std::vector<Col> _csr_col;
            std::vector<T>   _csr_val;

            if (csr_nnz) {
                _csr_ptr.resize(n + 1);
                _csr_col.resize(csr_nnz);
                _csr_val.resize(csr_nnz * N * N);

                _csr_ptr[0] = 0;
                for(size_t i = 0; i < n; ++i) {
                    size_t w = ptr[i+1] - ptr[i];
                    _csr_ptr[i+1] = _csr_ptr[i] + (w > ell_width ? w - ell_width : 0);
                }
            }


            for(size_t i = 0; i < n; ++i) {
                size_t w = 0;
                Ptr csr_head = csr_nnz ? _csr_ptr[i] : 0;
                for(Ptr j = ptr[i], e = ptr[i+1]; j < e; ++j, ++w) {
                    Col c = col[j];
                    Val v = val[j];

                    if (w < ell_width) {
                        _ell_col[i + w * ell_pitch] = c;
                        for(int k = 0, ii = 0; ii < N; ++ii)
                            for(int jj = 0; jj < N; ++jj, ++k)
                                _ell_val[k * ell_nnz + w * ell_pitch + i] = v(ii,jj);
                    } else {
                        _csr_col[csr_head] = c;
                        for(int k = 0, ii = 0; ii < N; ++ii)
                            for(int jj = 0; jj < N; ++jj, ++k)
                                _csr_val[k * csr_nnz + csr_head] = v(ii,jj);
                        ++csr_head;
                    }
                }
            }

            {
                size_t ell_size = ell_pitch * ell_width;
                ell_col = backend::device_vector<Col>(q[0], ell_size, _ell_col.data());
                ell_val = backend::device_vector<T>  (q[0], ell_nnz * N * N, _ell_val.data());
            }

            if (csr_nnz) {
                csr_ptr = backend::device_vector<Ptr>(q[0], n + 1,   _csr_ptr.data());
                csr_col = backend::device_vector<Col>(q[0], csr_nnz, _csr_col.data());
                csr_val = backend::device_vector<T>  (q[0], csr_nnz * N * N, _csr_val.data());
            }
        }

        // Dummy matrix; used internally to pass empty parameters to kernels.
        ell(const backend::command_queue &q)
            : q(q), n(0), m(0), nnz(0), ell_width(0), ell_pitch(0), csr_nnz(0)
        {}

        template <class Expr>
        friend
        typename std::enable_if<
            boost::proto::matches<
                typename boost::proto::result_of::as_expr<Expr>::type,
                vector_expr_grammar
            >::value,
            matrix_vector_product<ell, Expr>
        >::type
        operator*(const ell &A, const Expr &x) {
            return matrix_vector_product<ell, Expr>(A, x);
        }

        template <class Vector>
        static void terminal_preamble(const Vector &x, backend::source_generator &src,
            const backend::command_queue &q, const std::string &prm_name,
            detail::kernel_generator_state_ptr state)
        {
            detail::output_terminal_preamble tp(src, q, prm_name + "_x", state);
            boost::proto::eval(boost::proto::as_child(x), tp);
        }

        template <class Vector>
        static void local_terminal_init(const Vector &x, backend::source_generator &src,
            const backend::command_queue &q, const std::string &prm_name,
            detail::kernel_generator_state_ptr state)
        {
            typedef typename detail::return_type<Vector>::type x_type;
            typedef amgcl::static_matrix<T,N,1> vector_value;

            src.new_line() << type_name<vector_value>() << " " << prm_name << "_sum;";
            for(int i = 0; i < N; ++i) {
                src.new_line() << prm_name << "_sum.data[" << i << "][0] = 0;";
            }
            src.open("{");

            // ELL part
            src.new_line() << type_name<size_t>() << " " << prm_name << "_ell_size = " << prm_name << "_ell_width * " << prm_name << "_ell_pitch;";
            src.new_line() << "for(size_t j = 0; j < " << prm_name << "_ell_width; ++j)";
            src.open("{");
            src.new_line() << type_name<Col>() << " nnz_idx = idx + j * " << prm_name << "_ell_pitch;";
            src.new_line() << type_name<Col>() << " c = " << prm_name << "_ell_col[nnz_idx];";
            src.new_line() << "if (c != (" << type_name<Col>() << ")(-1))";
            src.open("{");

            src.new_line() << type_name<Col>() << " idx = c;";

            {
                detail::output_local_preamble init_x(src, q, prm_name + "_x", state);
                boost::proto::eval(boost::proto::as_child(x), init_x);

                backend::source_generator vec_value;
                detail::vector_expr_context expr_x(vec_value, q, prm_name + "_x", state);
                boost::proto::eval(boost::proto::as_child(x), expr_x);

                src.open("{");
                src.new_line() << type_name<x_type>() << " v = " << vec_value.str() << ";";
                for(int k = 0, j = 0; j < N; ++j) {
                    src.new_line() << prm_name << "_sum.data[" << j << "][0] += ";
                    for(int i = 0; i < N; ++i, ++k) {
                        if (i) src << " + ";
                        src << prm_name << "_ell_val[" << k << " * " << prm_name << "_ell_size + nnz_idx] * v.data[" << i << "][0]";
                    }
                    src << ";";
                }
                src.close("}");
            }

            src.close("} else break;");
            src.close("}");

            // CSR part
            src.new_line() << "if (" << prm_name << "_csr_ptr)";
            src.open("{");
            src.new_line() << type_name<size_t>() << " " << prm_name << "_csr_size = " << prm_name << "_csr_ptr[n];";
            src.new_line() << type_name<Ptr>() << " csr_beg = " << prm_name << "_csr_ptr[idx];";
            src.new_line() << type_name<Ptr>() << " csr_end = " << prm_name << "_csr_ptr[idx+1];";
            src.new_line() << "for(" << type_name<Ptr>() << " j = csr_beg; j < csr_end; ++j)";
            src.open("{");

            src.new_line() << type_name<Col>() << " idx = " << prm_name << "_csr_col[j];";

            {
                detail::output_local_preamble init_x(src, q, prm_name + "_x", state);
                boost::proto::eval(boost::proto::as_child(x), init_x);

                backend::source_generator vec_value;
                detail::vector_expr_context expr_x(vec_value, q, prm_name + "_x", state);
                boost::proto::eval(boost::proto::as_child(x), expr_x);

                src.open("{");
                src.new_line() << type_name<x_type>() << " v = " << vec_value.str() << ";";
                for(int k = 0, j = 0; j < N; ++j) {
                    src.new_line() << prm_name << "_sum.data[" << j << "][0] += ";
                    for(int i = 0; i < N; ++i, ++k) {
                        if (i) src << " + ";
                        src << prm_name << "_csr_val[" << k << " * " << prm_name << "_csr_size + j] * v.data[" << i << "][0]";
                    }
                    src << ";";
                }
                src.close("}");
            }

            src.close("}");
            src.close("}");
            src.close("}");
        }

        template <class Vector>
        static void kernel_param_declaration(const Vector &x, backend::source_generator &src,
            const backend::command_queue &q, const std::string &prm_name,
            detail::kernel_generator_state_ptr state)
        {
            src.parameter< size_t >(prm_name + "_ell_width");
            src.parameter< size_t >(prm_name + "_ell_pitch");

            src.parameter< global_ptr<Col> >(prm_name + "_ell_col");
            src.parameter< global_ptr<T  > >(prm_name + "_ell_val");
            src.parameter< global_ptr<Ptr> >(prm_name + "_csr_ptr");
            src.parameter< global_ptr<Col> >(prm_name + "_csr_col");
            src.parameter< global_ptr<T  > >(prm_name + "_csr_val");

            detail::declare_expression_parameter decl_x(src, q, prm_name + "_x", state);
            detail::extract_terminals()(boost::proto::as_child(x), decl_x);
        }

        template <class Vector>
        static void partial_vector_expr(const Vector&, backend::source_generator &src,
            const backend::command_queue&, const std::string &prm_name,
            detail::kernel_generator_state_ptr)
        {
            src << prm_name << "_sum";
        }

        template <class Vector>
        void kernel_arg_setter(const Vector &x,
            backend::kernel &kernel, unsigned part, size_t index_offset,
            detail::kernel_generator_state_ptr state) const
        {
            kernel.push_arg(ell_width);
            kernel.push_arg(ell_pitch);
            if (ell_width) {
                kernel.push_arg(ell_col);
                kernel.push_arg(ell_val);
            } else {
                kernel.push_arg(static_cast<size_t>(0));
                kernel.push_arg(static_cast<size_t>(0));
            }
            if (csr_nnz) {
                kernel.push_arg(csr_ptr);
                kernel.push_arg(csr_col);
                kernel.push_arg(csr_val);
            } else {
                kernel.push_arg(static_cast<size_t>(0));
                kernel.push_arg(static_cast<size_t>(0));
                kernel.push_arg(static_cast<size_t>(0));
            }

            detail::set_expression_argument x_args(kernel, part, index_offset, state);
            detail::extract_terminals()( boost::proto::as_child(x), x_args);
        }

        template <class Vector>
        void expression_properties(const Vector&,
            std::vector<backend::command_queue> &queue_list,
            std::vector<size_t> &partition,
            size_t &size) const
        {
            queue_list = std::vector<backend::command_queue>(1, q);
            partition  = std::vector<size_t>(2, 0);
            partition.back() = size = n;
        }

        size_t rows()     const { return n; }
        size_t cols()     const { return m; }
        size_t nonzeros() const { return nnz; }
    private:
        backend::command_queue q;

        size_t n, m, nnz, ell_width, ell_pitch, csr_nnz;

        backend::device_vector<Col> ell_col;
        backend::device_vector<T>   ell_val;

        backend::device_vector<Ptr> csr_ptr;
        backend::device_vector<Col> csr_col;
        backend::device_vector<T>   csr_val;

        backend::device_vector<T> create_device_vector(const backend::command_queue &q,
                size_t nnz, const Val *host_data, bool fast = true)
        {
            backend::device_vector<T> val(q, nnz * N * N);

            if (fast) {
                backend::device_vector<T> tmp(q, nnz * N * N, reinterpret_cast<const T*>(host_data));

                VEX_FUNCTION(T, transpose, (int,k)(int,m)(int,nnz)(T*, v),
                        int i = k / nnz;
                        int j = k % nnz;
                        return v[j * m + i];
                        );

                vex::vector<T>(q,val) = transpose(vex::element_index(), N*N, nnz, raw_pointer(vex::vector<T>(q, tmp)));
            } else {
                auto v = val.map(q);

                for(int k = 0, i = 0; i < N; ++i)
                    for(int j = 0; j < N; ++j, ++k)
                        for(size_t m = 0; m < nnz; ++m)
                            v[k * nnz + m] = host_data[m](i,j);
            }

            return val;
        }

        backend::kernel& csr2ell_kernel() const {
            using namespace vex::detail;
            static kernel_cache cache;

            auto kernel = cache.find(q);
            if (kernel == cache.end()) {
                backend::source_generator src(q);

                src.begin_kernel("convert_csr2ell");
                src.begin_kernel_parameters();
                src.template parameter<size_t>("n");
                src.template parameter<size_t>("ell_width");
                src.template parameter<size_t>("ell_pitch");
                src.template parameter< global_ptr<const ptr_type> >("ptr");
                src.template parameter< global_ptr<const col_type> >("col");
                src.template parameter< global_ptr<const T> >("val");
                src.template parameter< global_ptr<col_type> >("ell_col");
                src.template parameter< global_ptr<T> >("ell_val");
                src.template parameter< global_ptr<const ptr_type> >("csr_ptr");
                src.template parameter< global_ptr<col_type> >("csr_col");
                src.template parameter< global_ptr<T> >("csr_val");
                src.end_kernel_parameters();
                src.new_line() << type_name<size_t>() << " nnz = ptr[n];";
                src.new_line() << type_name<size_t>() << " ell_nnz = ell_width * ell_pitch;";
                src.new_line() << type_name<size_t>() << " csr_nnz = csr_ptr ? csr_ptr[n] : 0;";
                src.grid_stride_loop().open("{");

                src.new_line() << type_name<int>() << " w = 0;";
                src.new_line() << type_name<ptr_type>() << " csr_head = 0;";
                src.new_line() << "if (csr_ptr) csr_head = csr_ptr[idx];";
                src.new_line() << "for(" << type_name<ptr_type>() << " j = ptr[idx], e = ptr[idx+1]; j < e; ++j, ++w)";
                src.open("{");
                src.new_line() << type_name<col_type>() << " c = col[j];";
                src.new_line() << "if (w < ell_width) {";
                src.new_line() << "  ell_col[idx + w * ell_pitch] = c;";
                for(int i = 0; i < N * N; ++i)
                src.new_line() << "  ell_val[" << i << " * ell_nnz + w * ell_pitch + idx] = val[" << i << " * nnz + j];";
                src.new_line() << "} else {";
                src.new_line() << "  csr_col[csr_head] = c;";
                for(int i = 0; i < N * N; ++i)
                src.new_line() << "  csr_val[" << i << " * csr_nnz + csr_head] = val[" << i << " * nnz + j];";
                src.new_line() << "  ++csr_head;";
                src.new_line() << "}";
                src.close("}");
                src.close("}");
                src.end_kernel();

                kernel = cache.insert(q, backend::kernel(q, src.str(), "convert_csr2ell"));
            }

            return kernel->second;
        }

        template <class PtrRange, class ColRange, class ValRange>
        void convert(
                const PtrRange &host_ptr,
                const ColRange &host_col,
                const ValRange &host_val
                )
        {
            size_t nnz = host_ptr[n];

            backend::device_vector<Ptr> Aptr(q, n + 1, &host_ptr[0]);
            backend::device_vector<Col> Acol(q, nnz, &host_col[0]);
            backend::device_vector<T>   Aval = create_device_vector(q, nnz, &host_val[0]);

            /* 1. Get optimal ELL widths for local and remote parts. */
            // Speed of ELL relative to CSR:
            const double ell_vs_csr = 3.0;

            // Find maximum widths for local and remote parts:
            std::vector<backend::command_queue> ctx(1, q);
            Reductor<int, MAX> max(ctx);

            vex::vector<Ptr> ptr(q, Aptr);

            VEX_FUNCTION(Ptr, row_width, (size_t, i)(const Ptr*, ptr),
                    return ptr[i+1] - ptr[i];
                    );

            int max_width = max(row_width(element_index(0, n), raw_pointer(ptr)));

            // Build width distribution histogram.
            vex::vector<int> hist(ctx, max_width + 1);
            hist = 0;
            eval(atomic_add(&permutation(row_width(element_index(0, n), raw_pointer(ptr)))(hist), 1));

            // Estimate optimal width for ELL part of the matrix,
            // count nonzeros in CSR part of the matrix
            ell_width = max_width;
            {
                auto h = hist.map(0);

                for(int i = 0, rows = n; i < max_width; ++i) {
                    rows -= h[i]; // Number of rows wider than i.
                    if (ell_vs_csr * rows < n) {
                        ell_width = i;
                        break;
                    }
                }

                for(int i = ell_width + 1; i <= max_width; ++i)
                    csr_nnz += h[i] * (i - ell_width);
            }

            if (ell_width == 0) {
                assert(csr_nnz == nnz);

                csr_ptr = Aptr;
                csr_col = Acol;
                csr_val = Aval;

                return;
            }

            if (csr_nnz) {
                VEX_FUNCTION(int, csr_width, (size_t, ell_width)(size_t, i)(const Ptr*, ptr),
                        if (i == 0) return 0;
                        int w = ptr[i] - ptr[i-1];
                        return (w > ell_width) ? (w - ell_width) : 0;
                        );

                vex::vector<ptr_type> csr_w(ctx, n+1);

                csr_ptr = backend::device_vector<Ptr>(q, n + 1);
                csr_col = backend::device_vector<Col>(q, csr_nnz);
                csr_val = backend::device_vector<T>  (q, csr_nnz * N * N);

                csr_w = csr_width(ell_width, element_index(), raw_pointer(ptr));
                vector<ptr_type> csr_p(q, csr_ptr);
                inclusive_scan(csr_w, csr_p);
            }


            /* 3. Split the input matrix into ELL and CSR submatrices. */
            ell_col = backend::device_vector<Col>(q, ell_pitch * ell_width);
            ell_val = backend::device_vector<T>  (q, ell_pitch * ell_width * N * N);

            vex::vector<Col>(q, ell_col) = -1;

            auto &convert = csr2ell_kernel();

            convert.push_arg(n);
            convert.push_arg(ell_width);
            convert.push_arg(ell_pitch);
            convert.push_arg(Aptr);
            convert.push_arg(Acol);
            convert.push_arg(Aval);
            convert.push_arg(ell_col);
            convert.push_arg(ell_val);
            if (csr_nnz) {
                convert.push_arg(csr_ptr);
                convert.push_arg(csr_col);
                convert.push_arg(csr_val);
            } else {
                convert.push_arg(static_cast<size_t>(0));
                convert.push_arg(static_cast<size_t>(0));
                convert.push_arg(static_cast<size_t>(0));
            }
            convert(q);
        }
};

} // namespace sparse
} // namespace vex

namespace amgcl {
namespace backend {

template <typename T, int N>
std::string vexcl_static_matrix_declaration() {
    std::ostringstream s;
    s << "typedef struct { " << vex::type_name<T>() << " data[" << N << "][" << N << "]; } "
         "amgcl_matrix_" << vex::type_name<T>() << "_" << N << "x" << N << ";\n";
    if (N != 1)
    s << "typedef struct { " << vex::type_name<T>() << " data[" << N << "][" << 1 << "]; } "
         "amgcl_matrix_" << vex::type_name<T>() << "_" << N << "x" << 1 << ";\n";
    return s.str();
}

template <typename T, int N>
struct vex_scale {
    typedef static_matrix<T,N,1> vector;

    struct apply_type : vex::UserFunction<apply_type, vector(T, vector)> {
        apply_type() {}

        static std::string name() {
            return "scale_" + vex::type_name<vector>();
        }

        static void define(vex::backend::source_generator &src, const std::string &fname = name()) {
            src.begin_function<vector>(fname);
            src.begin_function_parameters();
            src.parameter<T>("a");
            src.parameter<vector>("m");
            src.end_function_parameters();
            for(int i = 0; i < N; ++i)
                src.new_line() << "m.data[" << i << "][0] *= a;";
            src.new_line() << "return m;";
            src.end_function();
        }
    } const apply;
};

template <typename T, int N>
struct vex_add {
    typedef static_matrix<T,N,1> vector;

    struct apply_type : vex::UserFunction<apply_type, vector(vector, vector)> {
        apply_type() {}

        static std::string name() {
            return "add_" + vex::type_name<vector>();
        }

        static void define(vex::backend::source_generator &src, const std::string &fname = name()) {
            src.begin_function<vector>(fname);
            src.begin_function_parameters();
            src.parameter<vector>("a");
            src.parameter<vector>("b");
            src.end_function_parameters();
            for(int i = 0; i < N; ++i)
                src.new_line() << "a.data[" << i << "][0] += "
                               << "b.data[" << i << "][0];";
            src.new_line() << "return a;";
            src.end_function();
        }
    } const apply;
};

template <typename T, int N>
struct vex_sub {
    typedef static_matrix<T,N,1> vector;

    struct apply_type : vex::UserFunction<apply_type, vector(vector, vector)> {
        apply_type() {}

        static std::string name() {
            return "sub_" + vex::type_name<vector>();
        }

        static void define(vex::backend::source_generator &src, const std::string &fname = name()) {
            src.begin_function<vector>(fname);
            src.begin_function_parameters();
            src.parameter<vector>("a");
            src.parameter<vector>("b");
            src.end_function_parameters();
            for(int i = 0; i < N; ++i)
                src.new_line() << "a.data[" << i << "][0] -= "
                               << "b.data[" << i << "][0];";
            src.new_line() << "return a;";
            src.end_function();
        }
    } const apply;
};

template <typename T, int N>
struct vex_mul {
    typedef static_matrix<T,N,N> matrix;
    typedef static_matrix<T,N,1> vector;

    struct apply_type : vex::UserFunction<apply_type, vector(matrix, vector)> {
        apply_type() {}

        static std::string name() {
            return "mul_" + vex::type_name<matrix>();
        }

        static void define(vex::backend::source_generator &src, const std::string &fname = name()) {
            src.begin_function<vector>(fname);
            src.begin_function_parameters();
            src.parameter<matrix>("a");
            src.parameter<vector>("b");
            src.end_function_parameters();
            src.new_line() << vex::type_name<vector>() << " c;";
            for(int i = 0; i < N; ++i) {
                src.new_line() << "c.data[" << i << "][0] = ";
                for(int j = 0; j < N; ++j) {
                    if (j) src << " + ";
                    src << "a.data[" << i << "][" << j << "] * b.data[" << j << "][0]";
                }
                src << ";";
            }
            src.new_line() << "return c;";
            src.end_function();
        }
    } const apply;
};

template <typename Alpha, typename Beta, typename T, int B>
struct spmv_impl<Alpha,
    vex::sparse::distributed<vex::sparse::matrix<static_matrix<T,B,B>, ptrdiff_t, ptrdiff_t>>,
    vex::vector<static_matrix<T,B,1>>, Beta, vex::vector<static_matrix<T,B,1>>>
{
    typedef vex::sparse::distributed<vex::sparse::matrix<static_matrix<T,B,B>, ptrdiff_t, ptrdiff_t>> matrix;
    typedef vex::vector<static_matrix<T,B,1>> vector;

    static void apply(Alpha alpha, const matrix &A, const vector &x, Beta beta, vector &y)
    {
        if (beta)
            y = vex_add<T,B>().apply(vex_scale<T,B>().apply(alpha, A * x), vex_scale<T,B>().apply(beta, y));
        else
            y = vex_scale<T,B>().apply(alpha, A * x);
    }
};

template <typename T, int B>
struct residual_impl<
    vex::sparse::distributed<vex::sparse::matrix<static_matrix<T,B,B>, ptrdiff_t, ptrdiff_t>>,
    vex::vector<static_matrix<T,B,1>>,
    vex::vector<static_matrix<T,B,1>>,
    vex::vector<static_matrix<T,B,1>>
    >
{
    typedef vex::sparse::distributed<vex::sparse::matrix<static_matrix<T,B,B>, ptrdiff_t, ptrdiff_t>> matrix;
    typedef vex::vector<static_matrix<T,B,1>> vector;

    static void apply(const vector &rhs, const matrix &A, const vector &x, vector &r)
    {
        r = vex_sub<T,B>().apply(rhs, A * x);
    }
};

template < typename Alpha, typename Beta, typename T, int B >
struct vmul_impl<
    Alpha, vex::vector< static_matrix<T,B,B> >,
    vex::vector< static_matrix<T,B,1> >,
    Beta, vex::vector< static_matrix<T,B,1> >
    >
{
    typedef vex::vector< static_matrix<T,B,B> > matrix;
    typedef vex::vector< static_matrix<T,B,1> > vector;

    static void apply(Alpha a, const matrix &x, const vector &y, Beta b, vector &z)
    {
        if (b)
            z = vex_add<T,B>().apply(vex_scale<T,B>().apply(a, vex_mul<T,B>().apply(x, y)), vex_scale<T,B>().apply(b, z));
        else
            z = vex_scale<T,B>().apply(a, vex_mul<T,B>().apply(x, y));
    }
};

template < typename T, int B >
struct clear_impl< vex::vector< static_matrix<T,B,1> > >
{
    typedef static_matrix<T,B,1> vector_value;
    typedef vex::vector<vector_value> vector;

    static void apply(vector &x) {
        x.template reinterpret<T>() = 0;
    }
};

template < typename T, int B >
struct copy_impl<
    vex::vector< static_matrix<T,B,1> >,
    vex::vector< static_matrix<T,B,1> >
    >
{
    typedef vex::vector< static_matrix<T,B,1> > vector;

    static void apply(const vector &x, vector &y) {
        auto X = x.template reinterpret<T>();
        auto Y = y.template reinterpret<T>();
        Y = X;
    }
};

template < typename A, typename B, typename T, int N >
struct axpby_impl<
    A, vex::vector< static_matrix<T, N, 1> >,
    B, vex::vector< static_matrix<T, N, 1> >
    >
{
    typedef vex::vector< static_matrix<T,N,1> > vector;

    static void apply(A a, const vector &x, B b, vector &y) {
        if (b)
            y.template reinterpret<T>() =
                a * x.template reinterpret<T>() +
                b * y.template reinterpret<T>();
        else
            y.template reinterpret<T>() =
                a * x.template reinterpret<T>();
    }
};

template < typename A, typename B, typename C, typename T, int N >
struct axpbypcz_impl<
    A, vex::vector< static_matrix<T, N, 1> >,
    B, vex::vector< static_matrix<T, N, 1> >,
    C, vex::vector< static_matrix<T, N, 1> >
    >
{
    typedef vex::vector< static_matrix<T,N,1> > vector;

    static void apply(A a, const vector &x, B b, const vector &y, C c, vector &z) {
        if (c)
            z.template reinterpret<T>() =
                a * x.template reinterpret<T>() +
                b * y.template reinterpret<T>() +
                c * z.template reinterpret<T>();
        else
            z.template reinterpret<T>() =
                a * x.template reinterpret<T>() +
                b * y.template reinterpret<T>();
    }
};

template < typename T, int B >
struct inner_product_impl<
    vex::vector< static_matrix<T,B,1> >,
    vex::vector< static_matrix<T,B,1> >
    >
{
    typedef T return_type;
    typedef static_matrix<T,B,1> vector_value;
    typedef vex::vector<vector_value> vector;

    static return_type get(const vector &x, const vector &y) {
        vex::Reductor<T, vex::SUM_Kahan> sum( x.queue_list() );
        return sum( x.template reinterpret<T>() * y.template reinterpret<T>() );
    }
};

} // namespace backend
} // namespace amgcl

#endif
