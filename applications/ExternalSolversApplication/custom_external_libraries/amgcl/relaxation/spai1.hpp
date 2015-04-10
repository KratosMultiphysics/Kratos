#ifndef AMGCL_RELAXATION_SPAI1_HPP
#define AMGCL_RELAXATION_SPAI1_HPP

/*
The MIT License

Copyright (c) 2012-2015 Denis Demidov <dennis.demidov@gmail.com>

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
 * \file   amgcl/relaxation/spai1.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Sparse approximate inverse relaxation scheme.
 */

#include <vector>
#ifdef _OPENMP
#  include <omp.h>
#endif

#include <boost/shared_ptr.hpp>
#include <amgcl/backend/interface.hpp>
#include <amgcl/util.hpp>
#include <amgcl/detail/qr.hpp>

namespace amgcl {
namespace relaxation {

/// Sparse approximate interface smoother.
/**
 * Sparsity pattern of the approximate inverse matrix coincides with that of A.
 *
 * \tparam Backend Backend for temporary structures allocation.
 * \ingroup relaxation
 * \sa \cite Broker2002
 */
template <class Backend>
struct spai1 {
    typedef typename Backend::value_type value_type;
    typedef typename Backend::vector     vector;

    /// Relaxation parameters.
    struct params {
        params() {}
        params(const boost::property_tree::ptree&) {}
        void get(boost::property_tree::ptree&, const std::string&) const {}
    };

    /// \copydoc amgcl::relaxation::damped_jacobi::damped_jacobi
    template <class Matrix>
    spai1( const Matrix &A, const params &, const typename Backend::params &backend_prm)
    {
        typedef typename backend::value_type<Matrix>::type   value_type;
        typedef typename backend::row_iterator<Matrix>::type row_iterator;

        const size_t n = rows(A);
        const size_t m = cols(A);

        boost::shared_ptr<Matrix> Ainv = boost::make_shared<Matrix>();
        Ainv->nrows = n;
        Ainv->ncols = m;

        Ainv->ptr = A.ptr;
        Ainv->col = A.col;
        Ainv->val.assign(A.ptr.back(), 0);

#pragma omp parallel
        {
#ifdef _OPENMP
            int nt  = omp_get_num_threads();
            int tid = omp_get_thread_num();

            size_t chunk_size  = (n + nt - 1) / nt;
            size_t chunk_start = tid * chunk_size;
            size_t chunk_end   = std::min(n, chunk_start + chunk_size);
#else
            size_t chunk_start = 0;
            size_t chunk_end   = n;
#endif

            std::vector<ptrdiff_t> marker(m, -1);
            std::vector<ptrdiff_t> I, J;
            std::vector<value_type> B, ek;
            amgcl::detail::QR<value_type> qr;

            for(size_t i = chunk_start; i < chunk_end; ++i) {
                ptrdiff_t row_beg = A.ptr[i];
                ptrdiff_t row_end = A.ptr[i + 1];

                I.assign(A.col.begin() + row_beg, A.col.begin() + row_end);

                J.clear();
                for(ptrdiff_t j = row_beg; j < row_end; ++j) {
                    ptrdiff_t c = A.col[j];

                    for(ptrdiff_t jj = A.ptr[c], ee = A.ptr[c + 1]; jj < ee; ++jj) {
                        ptrdiff_t cc = A.col[jj];
                        if (marker[cc] < 0) {
                            marker[cc] = 1;
                            J.push_back(cc);
                        }
                    }
                }
                std::sort(J.begin(), J.end());
                B.assign(I.size() * J.size(), 0);
                ek.assign(J.size(), 0);
                for(size_t j = 0; j < J.size(); ++j) {
                    marker[J[j]] = j;
                    if (J[j] == static_cast<ptrdiff_t>(i)) ek[j] = 1;
                }

                for(ptrdiff_t j = row_beg; j < row_end; ++j) {
                    ptrdiff_t c = A.col[j];

                    for(row_iterator a = row_begin(A, c); a; ++a)
                        B[marker[a.col()] * I.size() + j - row_beg] = a.value();
                }

                qr.compute(J.size(), I.size(), B.data(), /*need Q: */false);
                qr.solve(ek.data(), &Ainv->val[row_beg]);

                for(size_t j = 0; j < J.size(); ++j)
                    marker[J[j]] = -1;
            }
        }

        M = Backend::copy_matrix(Ainv, backend_prm);
    }

    /// \copydoc amgcl::relaxation::damped_jacobi::apply_pre
    template <class Matrix, class VectorRHS, class VectorX, class VectorTMP>
    void apply_pre(
            const Matrix &A, const VectorRHS &rhs, VectorX &x, VectorTMP &tmp,
            const params&
            ) const
    {
        apply(A, rhs, x, tmp);
    }

    /// \copydoc amgcl::relaxation::damped_jacobi::apply_post
    template <class Matrix, class VectorRHS, class VectorX, class VectorTMP>
    void apply_post(
            const Matrix &A, const VectorRHS &rhs, VectorX &x, VectorTMP &tmp,
            const params&
            ) const
    {
        apply(A, rhs, x, tmp);
    }

    private:
        boost::shared_ptr<typename Backend::matrix> M;

        template <class Matrix, class VectorRHS, class VectorX, class VectorTMP>
        void apply(
                const Matrix &A, const VectorRHS &rhs, VectorX &x, VectorTMP &tmp
                ) const
        {
            backend::residual(rhs, A, x, tmp);
            backend::spmv(1, *M, tmp, 1, x);
        }

};

} // namespace relaxation
} // namespace amgcl


#endif
