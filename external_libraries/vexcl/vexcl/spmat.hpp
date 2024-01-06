#ifndef VEXCL_SPMAT_HPP
#define VEXCL_SPMAT_HPP

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
 * \file   vexcl/spmat.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  OpenCL sparse matrix.
 */

#include <vector>
#include <set>
#include <unordered_map>
#include <string>
#include <memory>
#include <algorithm>
#include <iostream>
#include <type_traits>

#include <vexcl/vector.hpp>
#include <vexcl/vector_view.hpp>

#if defined(VEXCL_BACKEND_CUDA)
#  include <cuda.h>
#  if (CUDA_VERSION > 10010) && defined(VEXCL_USE_CUSPARSE)
#    undef VEXCL_USE_CUSPARSE
#  endif
#endif

namespace vex {

/// Sparse matrix in hybrid ELL-CSR format.
template <typename val_t, typename col_t = size_t, typename idx_t = size_t>
class SpMat {
    public:
        typedef val_t value_type;
        typedef typename cl_scalar_of<val_t>::type scalar_type;

        /// Empty constructor.
        SpMat() : nrows(0), ncols(0), nnz(0) {}

        /// Constructor.
        /**
         * Constructs GPU representation of the \f$n \times m\f$matrix. Input
         * matrix is in CSR format. GPU matrix utilizes ELL format and is split
         * equally across all compute devices.
         */
        SpMat(const std::vector<backend::command_queue> &queue,
              size_t n, size_t m, const idx_t *row, const col_t *col, const val_t *val
              )
            : queue(queue), part(partition(n, queue)),
              mtx(queue.size()), exc(queue.size()),
              nrows(n), ncols(m), nnz(row[n])
        {
            auto col_part = partition(m, queue);

            // Create secondary queues.
            for(auto q = queue.begin(); q != queue.end(); q++)
                squeue.push_back(backend::duplicate_queue(*q));

            std::vector<std::set<col_t>> ghost_cols = setup_exchange(col_part, row, col);

            // Each device get it's own strip of the matrix.
#ifdef _OPENMP
#  pragma omp parallel for schedule(static,1)
#endif
            for(int d = 0; d < static_cast<int>(queue.size()); d++) {
                if (part[d + 1] > part[d]) {
                    if ( backend::is_cpu(queue[d]) )
                        mtx[d].reset(
                                new SpMatCSR(queue[d],
                                    row + part[d], row + part[d+1], col, val,
                                    static_cast<col_t>(col_part[d]), static_cast<col_t>(col_part[d+1]), ghost_cols[d])
                                );
                    else
                        mtx[d].reset(
                                new SpMatHELL(queue[d],
                                    row + part[d], row + part[d + 1], col, val,
                                    static_cast<col_t>(col_part[d]), static_cast<col_t>(col_part[d+1]), ghost_cols[d])
                                );
                }
            }
        }


        // Matrix-vector multiplication.
        /*
         * Matrix vector multiplication (\f$y = \alpha Ax\f$ or \f$y += \alpha
         * Ax\f$) is performed in parallel on all registered compute devices.
         * Ghost values of x are transfered across GPU boundaries as needed.
         * \param x      input vector.
         * \param y      output vector.
         * \param alpha  coefficient in front of matrix-vector product
         * \param append if set, matrix-vector product is appended to y.
         *               Otherwise, y is replaced with matrix-vector product.
         */
        void apply(const vex::vector<val_t> &x, vex::vector<val_t> &y,
                 scalar_type alpha = 1, bool append = false) const
        {
            using namespace detail;

            if (rx.size()) {
                // Gather values to send to neighbors.
                for(unsigned d = 0; d < queue.size(); d++) {
                    if (cidx[d + 1] > cidx[d]) {
                        vex::vector<col_t> cols(queue[d], exc[d].cols_to_send);
                        vex::vector<val_t> vals(queue[d], exc[d].vals_to_send);
                        vex::vector<val_t> xloc(queue[d], x(d));

                        vals = permutation(cols)(xloc);
                    }
                }

                for(unsigned d = 0; d < queue.size(); d++)
                    if (cidx[d + 1] > cidx[d]) queue[d].finish();
            }

            // Start computing contribution from local part of the matrix.
            for(unsigned d = 0; d < queue.size(); d++)
                if (mtx[d]) {
                    backend::select_context(queue[d]);
                    mtx[d]->mul_local(x(d), y(d), alpha, append);
                }


            if (rx.size()) {
                // Meanwhile, get gathered values to host, ...
                for(unsigned d = 0; d < queue.size(); d++) {
                    if (cidx[d + 1] > cidx[d]) {
                        backend::select_context(squeue[d]);
                        vex::vector<val_t> vals(squeue[d], exc[d].vals_to_send);
                        vex::copy(vals.begin(), vals.end(), &rx[cidx[d]], /*blocking=*/false);
                    }
                }

                for(unsigned d = 0; d < queue.size(); d++)
                    if (cidx[d + 1] > cidx[d]) squeue[d].finish();

                // ... send ghost points from our neighbors to device, ...
                for(unsigned d = 0; d < queue.size(); d++) {
                    if (exc[d].cols_to_recv.size()) {
                        for(size_t i = 0; i < exc[d].cols_to_recv.size(); i++)
                            exc[d].vals_to_recv[i] = rx[exc[d].cols_to_recv[i]];

                        exc[d].rx.write(squeue[d], 0, exc[d].vals_to_recv.size(),
                                exc[d].vals_to_recv.data()
                                );
                    }
                }

                for(unsigned d = 0; d < queue.size(); d++)
                    if (exc[d].cols_to_recv.size()) squeue[d].finish();

                // Compute contribution from remote part of the matrix.
                for(unsigned d = 0; d < queue.size(); d++) {
                    if (exc[d].cols_to_recv.size()) {
                        backend::select_context(queue[d]);
                        mtx[d]->mul_remote(exc[d].rx, y(d), alpha);
                    }
                }
            }
        }

        /// Number of rows.
        size_t rows() const { return nrows; }
        /// Number of columns.
        size_t cols() const { return ncols; }
        /// Number of non-zero entries.
        size_t nonzeros() const { return nnz;   }

#if !defined(VEXCL_BACKEND_CUDA) || !defined(VEXCL_USE_CUSPARSE)
        static void inline_preamble(backend::source_generator &src,
                const backend::command_queue &queue, const std::string &prm_name,
                detail::kernel_generator_state_ptr)
        {
            if (backend::is_cpu(queue))
                SpMatCSR::inline_preamble(src, prm_name);
            else
                SpMatHELL::inline_preamble(src, prm_name);
        }

        static void inline_expression(backend::source_generator &src,
                const backend::command_queue &queue, const std::string &prm_name,
                detail::kernel_generator_state_ptr)
        {
            if (backend::is_cpu(queue))
                SpMatCSR::inline_expression(src, prm_name);
            else
                SpMatHELL::inline_expression(src, prm_name);
        }

        static void inline_parameters(backend::source_generator &src,
                const backend::command_queue &queue, const std::string &prm_name,
                detail::kernel_generator_state_ptr)
        {
            if (backend::is_cpu(queue))
                SpMatCSR::inline_parameters(src, prm_name);
            else
                SpMatHELL::inline_parameters(src, prm_name);
        }

        static void inline_arguments(backend::kernel &kernel, unsigned part,
                size_t /*index_offset*/, const SpMat &A, const vector<val_t> &x,
                detail::kernel_generator_state_ptr)
        {
            A.mtx[part]->setArgs(kernel, part, x);
        }
#endif
    private:
        template <typename T>
        static inline size_t bytes(const std::vector<T> &v) {
            return v.size() * sizeof(T);
        }

        struct sparse_matrix {
            virtual void mul_local(
                    const backend::device_vector<val_t> &x,
                    backend::device_vector<val_t> &y,
                    scalar_type alpha, bool append
                    ) const = 0;

            virtual void mul_remote(
                    const backend::device_vector<val_t> &x,
                    backend::device_vector<val_t> &y,
                    scalar_type alpha
                    ) const = 0;

#if !defined(VEXCL_BACKEND_CUDA) || !defined(VEXCL_USE_CUSPARSE)
            virtual void setArgs(backend::kernel &kernel, unsigned part, const vector<val_t> &x) const = 0;
#endif

            virtual ~sparse_matrix() {}
        };


#if !defined(VEXCL_BACKEND_CUDA) || !defined(VEXCL_USE_CUSPARSE)
#  include <vexcl/spmat/hybrid_ell.inl>
#  include <vexcl/spmat/csr.inl>
#else
#  include <vexcl/backend/cuda/cusparse.hpp>
#  include <vexcl/backend/cuda/hybrid_ell.inl>
#  include <vexcl/backend/cuda/csr.inl>
#endif

        struct exdata {
            std::vector<col_t> cols_to_recv;
            mutable std::vector<val_t> vals_to_recv;

            backend::device_vector<col_t> cols_to_send;
            backend::device_vector<val_t> vals_to_send;
            backend::device_vector<val_t> rx;
        };

        mutable std::vector<backend::command_queue> queue;
        mutable std::vector<backend::command_queue> squeue;
        const std::vector<size_t>           part;

        std::vector< std::unique_ptr<sparse_matrix> > mtx;

        mutable std::vector<exdata> exc;
        std::vector<size_t> cidx;
        mutable std::vector<val_t> rx;

        size_t nrows;
        size_t ncols;
        size_t nnz;

        std::vector<std::set<col_t>> setup_exchange(
                const std::vector<size_t> &col_part,
                const idx_t *row, const col_t *col
                )
        {
            auto is_local = [col_part](size_t c, int part) {
                return c >= col_part[part] && c < col_part[part + 1];
            };

            std::vector<std::set<col_t>> ghost_cols(queue.size());

            if (queue.size() <= 1) return ghost_cols;

            // Build sets of ghost points.
#ifdef _OPENMP
#  pragma omp parallel for schedule(static,1)
#endif
            for(int d = 0; d < static_cast<int>(queue.size()); d++) {
                for(size_t i = part[d]; i < part[d + 1]; i++) {
                    for(idx_t j = row[i]; j < row[i + 1]; j++) {
                        if (!is_local(col[j], d)) {
                            ghost_cols[d].insert(col[j]);
                        }
                    }
                }
            }

            // Complete set of points to be exchanged between devices.
            std::vector<col_t> cols_to_send;
            {
                std::set<col_t> cols_to_send_s;
                for(unsigned d = 0; d < queue.size(); d++)
                    cols_to_send_s.insert(ghost_cols[d].begin(), ghost_cols[d].end());

                cols_to_send.insert(cols_to_send.begin(), cols_to_send_s.begin(), cols_to_send_s.end());
            }

            // Build local structures to facilitate exchange.
            if (cols_to_send.size()) {
#ifdef _OPENMP
#  pragma omp parallel for schedule(static,1)
#endif
                for(int d = 0; d < static_cast<int>(queue.size()); d++) {
                    if (size_t rcols = ghost_cols[d].size()) {
                        exc[d].cols_to_recv.resize(rcols);
                        exc[d].vals_to_recv.resize(rcols);

                        exc[d].rx = backend::device_vector<val_t>(queue[d], rcols,
                                static_cast<const val_t*>(0), backend::MEM_READ_ONLY);

                        for(size_t i = 0, j = 0; i < cols_to_send.size(); i++)
                            if (ghost_cols[d].count(cols_to_send[i]))
                                exc[d].cols_to_recv[j++] = static_cast<col_t>(i);
                    }
                }

                rx.resize(cols_to_send.size());
                cidx.resize(queue.size() + 1);

                {
                    auto beg = cols_to_send.begin();
                    auto end = cols_to_send.end();
                    for(unsigned d = 0; d <= queue.size(); d++) {
                        cidx[d] = std::lower_bound(beg, end, static_cast<col_t>(col_part[d]))
                                - cols_to_send.begin();
                        beg = cols_to_send.begin() + cidx[d];
                    }
                }

                for(unsigned d = 0; d < queue.size(); d++) {
                    if (size_t ncols = cidx[d + 1] - cidx[d]) {
                        exc[d].vals_to_send = backend::device_vector<val_t>(
                                queue[d], ncols);

                        for(size_t i = cidx[d]; i < cidx[d + 1]; i++)
                            cols_to_send[i] -= static_cast<col_t>(col_part[d]);

                        exc[d].cols_to_send = backend::device_vector<col_t>(
                                queue[d], ncols, &cols_to_send[cidx[d]], backend::MEM_READ_ONLY);
                    }
                }

                for(unsigned d = 0; d < queue.size(); d++)
                    if (cidx[d + 1] > cidx[d]) queue[d].finish();
            }

            return ghost_cols;
        }
};

template <typename val_t, typename col_t, typename idx_t>
additive_operator< SpMat<val_t, col_t, idx_t >, vector<val_t> >
operator*(const SpMat<val_t, col_t, idx_t> &A, const vector<val_t> &x)
{
    return additive_operator< SpMat<val_t, col_t, idx_t >, vector<val_t> >(A, x);
}

#ifdef VEXCL_MULTIVECTOR_HPP
template <typename val_t, typename col_t, typename idx_t, class V>
typename std::enable_if<
    std::is_base_of<multivector_terminal_expression, V>::value &&
    std::is_same<val_t, typename V::sub_value_type>::value,
    multiadditive_operator< SpMat<val_t, col_t, idx_t>, V >
>::type
operator*(const SpMat<val_t, col_t, idx_t> &A, const V &x) {
    return multiadditive_operator< SpMat<val_t, col_t, idx_t>, V >(A, x);
}
#endif

/// Weights device wrt to additive_operator performance.
/**
 * Launches the following kernel on each device:
 \code
 y = A * x;
 \endcode
 * where x and y are vectors, and A is matrix for 3D Poisson problem in square
 * domain. Each device gets portion of the vector proportional to the
 * performance of this operation.
 */
inline double device_spmv_perf(const backend::command_queue &q) {
    static const size_t test_size = 64U;

    std::vector<backend::command_queue> queue(1, q);

    // Construct matrix for 3D Poisson problem in cubic domain.
    const size_t n   = test_size;
    const float  h2i = (n - 1.0f) * (n - 1.0f);

    std::vector<size_t> row;
    std::vector<size_t> col;
    std::vector<float>  val;

    row.reserve(n * n * n + 1);
    col.reserve(6 * (n - 2) * (n - 2) * (n - 2) + n * n * n);
    val.reserve(6 * (n - 2) * (n - 2) * (n - 2) + n * n * n);

    row.push_back(0);
    for(size_t k = 0, idx = 0; k < n; k++) {
        for(size_t j = 0; j < n; j++) {
            for(size_t i = 0; i < n; i++, idx++) {
                if (
                        i == 0 || i == (n - 1) ||
                        j == 0 || j == (n - 1) ||
                        k == 0 || k == (n - 1)
                   )
                {
                    col.push_back(idx);
                    val.push_back(1);
                    row.push_back(row.back() + 1);
                } else {
                    col.push_back(idx - n * n);
                    val.push_back(-h2i);

                    col.push_back(idx - n);
                    val.push_back(-h2i);

                    col.push_back(idx - 1);
                    val.push_back(-h2i);

                    col.push_back(idx);
                    val.push_back(6 * h2i);

                    col.push_back(idx + 1);
                    val.push_back(-h2i);

                    col.push_back(idx + n);
                    val.push_back(-h2i);

                    col.push_back(idx + n * n);
                    val.push_back(-h2i);

                    row.push_back(row.back() + 7);
                }
            }
        }
    }

    // Create device vectors and copy of the matrix.
    size_t n3 = n * n * n;
    vex::SpMat<float>  A(queue, n3, n3, row.data(), col.data(), val.data());
    vex::vector<float> x(queue, n3);
    vex::vector<float> y(queue, n3);

    // Warming run.
    x = 1;
    A.apply(x, y);

    // Measure performance.
    profiler<> prof(queue);
    prof.tic_cl("");
    A.apply(x, y);
    double time = prof.toc("");
    return 1.0 / time;
}

} // namespace vex

#include <vexcl/spmat/ccsr.hpp>
#include <vexcl/spmat/inline_spmv.hpp>

#endif
