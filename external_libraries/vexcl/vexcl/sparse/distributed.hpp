#ifndef VEXCL_SPARSE_DISTRIBUTED_HPP
#define VEXCL_SPARSE_DISTRIBUTED_HPP

#include <vector>
#include <unordered_map>
#include <algorithm>

#include <vexcl/util.hpp>
#include <vexcl/backend.hpp>
#include <vexcl/operations.hpp>
#include <vexcl/sparse/product.hpp>
#include <vexcl/sparse/spmv_ops.hpp>

namespace vex {
namespace sparse {

/// Value-type of a vector corresponding to the matrix with value-type T.
template <class T, class Enable = void>
struct rhs_of {
    typedef T type;
};

template <class Matrix, typename rhs_type = typename rhs_of<typename Matrix::val_type>::type>
class distributed {
    public:
        typedef typename Matrix::value_type value_type;

        template <class PtrRange, class ColRange, class ValRange>
        distributed(
                const std::vector<backend::command_queue> &q,
                size_t nrows, size_t ncols,
                const PtrRange &ptr,
                const ColRange &col,
                const ValRange &val,
                bool fast_setup = true
           )
            : q(q), n(nrows), m(ncols), nnz(boost::size(val)),
              row_part(partition(n, q)), col_part(partition(m, q)),
              A_loc(q.size()), A_rem(q.size())
        {
            if (q.size() == 1) {
                A_loc[0] = std::make_shared<Matrix>(q, nrows, ncols, ptr, col, val, fast_setup);
                return;
            }

            std::vector<std::vector<col_type>> rcols(q.size());

#ifdef _OPENMP
#  pragma omp parallel for schedule(static,1)
#endif
            for(int d = 0; d < static_cast<int>(q.size()); ++d) {
                size_t loc_rows = row_part[d+1] - row_part[d];

                if (!loc_rows) continue;

                col_type col_beg = static_cast<col_type>(col_part[d]);
                col_type col_end = static_cast<col_type>(col_part[d+1]);

                std::vector<ptr_type> loc_ptr(loc_rows + 1); loc_ptr[0] = 0;
                std::vector<ptr_type> rem_ptr(loc_rows + 1); rem_ptr[0] = 0;

                // Count nonzeros in local and remote parts of the matrix.
                for(size_t i = row_part[d], ii = 0; i < row_part[d+1]; ++i, ++ii) {
                    loc_ptr[ii+1] = loc_ptr[ii];
                    rem_ptr[ii+1] = rem_ptr[ii];

                    for(ptr_type j = ptr[i]; j < ptr[i+1]; ++j) {
                        col_type c = col[j];

                        if (col_beg <= c && c < col_end) {
                            ++loc_ptr[ii+1];
                        } else {
                            ++rem_ptr[ii+1];
                        }
                    }
                }

                // Fill local and remote parts of the matrix.
                size_t loc_nnz = loc_ptr.back();
                size_t rem_nnz = rem_ptr.back();

                std::vector<col_type> loc_col; loc_col.reserve(loc_nnz);
                std::vector<val_type> loc_val; loc_val.reserve(loc_nnz);

                std::vector<col_type> rem_col; rem_col.reserve(rem_nnz);
                std::vector<val_type> rem_val; rem_val.reserve(rem_nnz);

                for(size_t i = row_part[d], ii = 0; i < row_part[d+1]; ++i, ++ii) {
                    for(ptr_type j = ptr[i]; j < ptr[i+1]; ++j) {
                        col_type c = col[j];
                        val_type v = val[j];

                        if (col_beg <= c && c < col_end) {
                            loc_col.push_back(c - col_beg);
                            loc_val.push_back(v);
                        } else {
                            rem_col.push_back(c);
                            rem_val.push_back(v);
                        }
                    }
                }

                // Get list of unique remote columns.
                rcols[d] = rem_col;
                std::sort(rcols[d].begin(), rcols[d].end());
                rcols[d].erase(std::unique(rcols[d].begin(), rcols[d].end()), rcols[d].end());

                // Renumber remote columns.
                size_t nrcols = rcols[d].size();
                std::unordered_map<col_type, col_type> idx(2 * nrcols);

                for(size_t i = 0; i < nrcols; ++i) {
                    idx.insert(idx.end(), std::make_pair(rcols[d][i], static_cast<col_type>(i)));
                }

                for(size_t i = 0; i < rem_nnz; ++i) {
                    rem_col[i] = idx[rem_col[i]];
                }

                // Create local and remote parts of the matrix on the current
                // device.
                std::vector<backend::command_queue> qd = {q[d]};
                if (loc_nnz)
                    A_loc[d] = std::make_shared<Matrix>(
                            qd, loc_rows, col_end - col_beg,
                            loc_ptr, loc_col, loc_val, fast_setup);

                if (nrcols)
                    A_rem[d] = std::make_shared<Matrix>(
                            qd, loc_rows, nrcols, rem_ptr, rem_col, rem_val,
                            fast_setup);
            }


            // Setup exchange.
            // 1. Build the combined vector of ghost points across all GPUs
            size_t nrcols = 0;
            for(size_t d = 0; d < q.size(); ++d)
                nrcols += rcols[d].size();

            // K-way merge of rcols[...] into rem_cols:
            std::vector<col_type> rem_cols; rem_cols.reserve(nrcols);
            {
                std::vector<typename std::vector<col_type>::const_iterator> rc(q.size());
                for(size_t d = 0; d < q.size(); ++d) rc[d] = rcols[d].begin();

                while(true) {
                    bool found = false;
                    size_t winner;

                    for(size_t d = 0; d < q.size(); ++d) {
                        if (rc[d] == rcols[d].end()) continue;

                        if (!found) {
                            found  = true;
                            winner = d;
                        } else {
                            if (*rc[d] < *rc[winner]) winner = d;
                        }
                    }
                    if (!found) break;

                    col_type c = *rc[winner]++;

                    if (rem_cols.empty() || rem_cols.back() != c)
                        rem_cols.push_back(c);
                }
            }

            // 2. Each GPU has to update the global ghost vector with its data,
            // and each GPU has to get the points it needs from the vector.
            ex.resize(q.size());

            // See what elements of rem_vals each GPU needs to receive:
            for(size_t d = 0; d < q.size(); ++d) {
                size_t nrecv = rcols[d].size();
                if (!nrecv) continue;

                ex[d].vals_to_recv.resize(nrecv);
                ex[d].rem_x = backend::device_vector<rhs_type>(q[d], nrecv);

                // See where in rem_vals each element of rcols[d] is placed.
                ex[d].cols_to_recv.reserve(nrecv);
                auto rc = std::lower_bound(rem_cols.begin(), rem_cols.end(), rcols[d].front());
                for(auto c = rcols[d].begin(), e = rcols[d].end(); c != e; ++c) {
                    while(*rc < *c) ++rc;
                    assert(*rc == *c);
                    ex[d].cols_to_recv.push_back(static_cast<col_type>(std::distance(rem_cols.begin(), rc)));
                }
            }

            // See what elements of rem_vals each GPU needs to send:
            rval_ptr.resize(q.size() + 1);
            rval_ptr[0] = 0;
            for(size_t d = 0; d < q.size(); ++d) {
                col_type col_beg = static_cast<col_type>(col_part[d]);
                col_type col_end = static_cast<col_type>(col_part[d+1]);

                rval_ptr[d+1] = std::distance(rem_cols.begin(),
                        std::lower_bound(
                            rem_cols.begin() + rval_ptr[d], rem_cols.end(), col_end)
                        );
                for(size_t i = rval_ptr[d]; i < rval_ptr[d+1]; ++i)
                    rem_cols[i] -= col_beg;

                size_t nsend = rval_ptr[d+1] - rval_ptr[d];

                if (nsend) {
                    ex[d].vals_to_send = backend::device_vector<rhs_type>(q[d], nsend);
                    ex[d].cols_to_send = backend::device_vector<col_type>(q[d], nsend,
                            &rem_cols[rval_ptr[d]]);
                }
            }

            rem_vals.resize(rem_cols.size());
        }

        template <class Expr>
        friend
        typename std::enable_if<
            boost::proto::matches<
                typename boost::proto::result_of::as_expr<Expr>::type,
                vector_expr_grammar
            >::value,
            matrix_vector_product<distributed, Expr>
        >::type
        operator*(const distributed &A, const Expr &x) {
            A.exchange(x);
            return matrix_vector_product<distributed, Expr>(A, x);
        }

        template <class Vector>
        static void terminal_preamble(const Vector &x, backend::source_generator &src,
            const backend::command_queue &q, const std::string &prm_name,
            detail::kernel_generator_state_ptr state)
        {
            vex::vector<rhs_type> dummy;

            Matrix::terminal_preamble(x,     src, q, prm_name + "_loc", state);
            Matrix::terminal_preamble(dummy, src, q, prm_name + "_rem", state);
        }

        template <class Vector>
        static void local_terminal_init(const Vector &x, backend::source_generator &src,
            const backend::command_queue &q, const std::string &prm_name,
            detail::kernel_generator_state_ptr state)
        {
            typedef typename detail::return_type<Vector>::type x_type;
            typedef spmv_ops_impl<value_type, x_type> spmv_ops;

            vex::vector<rhs_type> dummy;

            Matrix::local_terminal_init(x,     src, q, prm_name + "_loc", state);
            Matrix::local_terminal_init(dummy, src, q, prm_name + "_rem", state);

            src.new_line() << type_name<rhs_type>() << " " << prm_name << "_sum = ";
            Matrix::partial_vector_expr(x, src, q, prm_name + "_loc", state);
            src << ";";

            backend::source_generator rem_value;
            Matrix::partial_vector_expr(dummy, rem_value, q, prm_name + "_rem", state);
            spmv_ops::append(src, prm_name + "_sum", rem_value.str());
        }

        template <class Vector>
        static void kernel_param_declaration(const Vector &x, backend::source_generator &src,
            const backend::command_queue &q, const std::string &prm_name,
            detail::kernel_generator_state_ptr state)
        {
            vex::vector<rhs_type> dummy;

            Matrix::kernel_param_declaration(x,     src, q, prm_name + "_loc", state);
            Matrix::kernel_param_declaration(dummy, src, q, prm_name + "_rem", state);
        }

        template <class Vector>
        static void partial_vector_expr(const Vector &, backend::source_generator &src,
            const backend::command_queue &, const std::string &prm_name,
            detail::kernel_generator_state_ptr)
        {
            src << prm_name << "_sum";
        }

        template <class Vector>
        void kernel_arg_setter(const Vector &x,
            backend::kernel &kernel, unsigned part, size_t index_offset,
            detail::kernel_generator_state_ptr state) const
        {
            if (A_loc[part]) {
                A_loc[part]->kernel_arg_setter(x, kernel, part, index_offset, state);
            } else {
                Matrix dummy_A(q[part]);
                backend::device_vector<rhs_type> dummy_x;
                dummy_A.kernel_arg_setter(dummy_x, kernel, part, index_offset, state);
            }

            if (A_rem[part]) {
                A_rem[part]->kernel_arg_setter(ex[part].rem_x, kernel, part, index_offset, state);
            } else {
                Matrix dummy_A(q[part]);
                backend::device_vector<rhs_type> dummy_x;
                dummy_A.kernel_arg_setter(dummy_x, kernel, part, index_offset, state);
            }
        }

        template <class Vector>
        void expression_properties(const Vector &x,
            std::vector<backend::command_queue> &queue_list,
            std::vector<size_t> &partition,
            size_t &size) const
        {
            queue_list = q;
            partition  = row_part;
            size       = n;
        }

        size_t rows()     const { return n;   }
        size_t cols()     const { return m;   }
        size_t nonzeros() const { return nnz; }
    private:
        typedef typename Matrix::ptr_type       ptr_type;
        typedef typename Matrix::col_type       col_type;
        typedef typename Matrix::val_type       val_type;

        mutable std::vector<backend::command_queue> q;

        size_t n, m, nnz;
        std::vector<size_t> row_part, col_part;
        std::vector<std::shared_ptr<Matrix>> A_loc, A_rem;

        mutable std::vector<rhs_type> rem_vals;
        std::vector<size_t>   rval_ptr;

        struct exdata {
            std::vector<col_type> cols_to_recv;
            mutable std::vector<rhs_type> vals_to_recv;

            backend::device_vector<col_type> cols_to_send;
            mutable backend::device_vector<rhs_type> vals_to_send;

            mutable backend::device_vector<rhs_type> rem_x;
        };

        std::vector<exdata> ex;

        template <class Expr>
        void exchange(const Expr &expr) const {
            using namespace vex::detail;
            static kernel_cache cache;

            if (q.size() == 1) return;

            // Gather values to send on the GPUs:
            for(unsigned d = 0; d < q.size(); ++d) {
                size_t nsend = rval_ptr[d+1] - rval_ptr[d];
                if (nsend == 0) continue;

                auto K = cache.find(q[d]);
                backend::select_context(q[d]);

                if (K == cache.end()) {
                    backend::source_generator src(q[d]);

                    output_terminal_preamble otp(src, q[d], "prm", empty_state());
                    boost::proto::eval(boost::proto::as_child(expr), otp);

                    src.begin_kernel("vexcl_sparse_gather");
                    src.begin_kernel_parameters();
                    src.template parameter<size_t>("n");
                    src.template parameter< global_ptr<const col_type> >("cols_to_send");
                    src.template parameter< global_ptr<rhs_type> >("vals_to_send");

                    extract_terminals()(boost::proto::as_child(expr),
                            declare_expression_parameter(src, q[d], "prm", empty_state()));

                    src.end_kernel_parameters();
                    src.grid_stride_loop().open("{");

                    src.new_line() << type_name<rhs_type>() << " cur_val;";
                    src.open("{");

                    src.new_line() << type_name<size_t>() << " pos = cols_to_send[idx];";
                    src.new_line() << type_name<size_t>() << " idx = pos;";

                    output_local_preamble olp(src, q[d], "prm", empty_state());
                    boost::proto::eval(boost::proto::as_child(expr), olp);

                    src.new_line() << "cur_val = ";
                    vector_expr_context vec(src, q[d], "prm", empty_state());
                    boost::proto::eval(boost::proto::as_child(expr), vec);
                    src << ";";

                    src.close("}");
                    src.new_line() << "vals_to_send[idx] = cur_val;";

                    src.close("}");
                    src.end_kernel();

                    K = cache.insert(q[d], backend::kernel(
                                q[d], src.str(), "vexcl_sparse_gather"));
                }

                auto &krn = K->second;

                krn.push_arg(nsend);
                krn.push_arg(ex[d].cols_to_send);
                krn.push_arg(ex[d].vals_to_send);

                extract_terminals()(boost::proto::as_child(expr),
                        set_expression_argument(krn, d, col_part[d], empty_state()));

                krn(q[d]);

                ex[d].vals_to_send.read(q[d], 0, nsend, &rem_vals[rval_ptr[d]]);
            }

            for(unsigned d = 0; d < q.size(); ++d)
                if (rval_ptr[d+1] > rval_ptr[d]) q[d].finish();

            for(unsigned d = 0; d < q.size(); ++d) {
                for(size_t i = 0; i < ex[d].cols_to_recv.size(); ++i)
                    ex[d].vals_to_recv[i] = rem_vals[ex[d].cols_to_recv[i]];

                ex[d].rem_x.write(q[d], 0, ex[d].cols_to_recv.size(), ex[d].vals_to_recv.data());
            }
        }
};

} // namespace sparse
} // namespace vex

#endif
