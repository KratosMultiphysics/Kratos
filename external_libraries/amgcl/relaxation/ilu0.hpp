#ifndef AMGCL_RELAXATION_ILU0_HPP
#define AMGCL_RELAXATION_ILU0_HPP

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
 * \file   amgcl/relaxation/ilu0.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Incomplete LU with zero fill-in relaxation scheme.
 */

#include <boost/typeof/typeof.hpp>

#include <amgcl/backend/builtin.hpp>
#include <amgcl/util.hpp>
#include <amgcl/relaxation/detail/ilu_solve.hpp>

namespace amgcl {
namespace relaxation {

/// ILU(0) smoother.
/**
 * \note ILU(0) is a serial algorithm and is only applicable to backends that
 * support matrix row iteration (e.g. amgcl::backend::builtin or
 * amgcl::backend::eigen).
 *
 * \param Backend Backend for temporary structures allocation.
 * \ingroup relaxation
 */
template <class Backend>
struct ilu0 {
    typedef typename Backend::value_type      value_type;
    typedef typename Backend::vector          vector;
    typedef typename Backend::matrix          matrix;
    typedef typename Backend::matrix_diagonal matrix_diagonal;

    typedef typename math::scalar_of<value_type>::type scalar_type;

    /// Relaxation parameters.
    struct params {
        /// Damping factor.
        scalar_type damping;

        /// Number of Jacobi iterations.
        /** \note Used for approximate solution of triangular systems on parallel backends */
        unsigned jacobi_iters;

        params() : damping(1), jacobi_iters(2) {}

        params(const boost::property_tree::ptree &p)
            : AMGCL_PARAMS_IMPORT_VALUE(p, damping)
            , AMGCL_PARAMS_IMPORT_VALUE(p, jacobi_iters)
        {
            AMGCL_PARAMS_CHECK(p, (damping)(jacobi_iters));
        }

        void get(boost::property_tree::ptree &p, const std::string &path) const {
            AMGCL_PARAMS_EXPORT_VALUE(p, path, damping);
            AMGCL_PARAMS_EXPORT_VALUE(p, path, jacobi_iters);
        }
    };

    /// \copydoc amgcl::relaxation::damped_jacobi::damped_jacobi
    template <class Matrix>
    ilu0( const Matrix &A, const params &, const typename Backend::params &bprm)
    {
        typedef typename backend::builtin<value_type>::matrix build_matrix;
        const size_t n = backend::rows(A);
        BOOST_AUTO(Aptr, A.ptr_data());
        BOOST_AUTO(Acol, A.col_data());
        BOOST_AUTO(Aval, A.val_data());

        size_t Lnz = 0, Unz = 0;

        for(ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(n); ++i) {
            ptrdiff_t row_beg = Aptr[i];
            ptrdiff_t row_end = Aptr[i + 1];

            for(ptrdiff_t j = row_beg; j < row_end; ++j) {
                ptrdiff_t c = Acol[j];
                if (c < i)
                    ++Lnz;
                else if (c > i)
                    ++Unz;
            }
        }

        boost::shared_ptr<build_matrix> L = boost::make_shared<build_matrix>();
        boost::shared_ptr<build_matrix> U = boost::make_shared<build_matrix>();

        L->nrows = L->ncols = n;
        L->ptr.reserve(n+1); L->ptr.push_back(0);
        L->col.reserve(Lnz);
        L->val.reserve(Lnz);

        U->nrows = U->ncols = n;
        U->ptr.reserve(n+1); U->ptr.push_back(0);
        U->col.reserve(Unz);
        U->val.reserve(Unz);

        std::vector<value_type> D;
        D.reserve(n);

        std::vector<value_type*> work(n, NULL);

        for(ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(n); ++i) {
            ptrdiff_t row_beg = Aptr[i];
            ptrdiff_t row_end = Aptr[i + 1];

            for(ptrdiff_t j = row_beg; j < row_end; ++j) {
                ptrdiff_t  c = Acol[j];
                value_type v = Aval[j];

                if (c < i) {
                    L->col.push_back(c);
                    L->val.push_back(v);
                    work[c] = &L->val.back();
                } else if (c == i) {
                    D.push_back(v);
                    work[c] = &D.back();
                } else {
                    U->col.push_back(c);
                    U->val.push_back(v);
                    work[c] = &U->val.back();
                }
            }

            L->ptr.push_back(L->val.size());
            U->ptr.push_back(U->val.size());

            for(ptrdiff_t j = row_beg; j < row_end; ++j) {
                ptrdiff_t c = Acol[j];

                // Exit if diagonal is reached
                if (c >= i) {
                    precondition(c == i, "No diagonal value in system matrix");
                    precondition(!math::is_zero(D[i]), "Zero pivot in ILU");

                    D[i] = math::inverse(D[i]);
                    break;
                }

                // Compute the multiplier for jrow
                value_type tl = (*work[c]) * D[c];
                *work[c] = tl;

                // Perform linear combination
                for(ptrdiff_t k = U->ptr[c]; k < U->ptr[c+1]; ++k) {
                    value_type *w = work[U->col[k]];
                    if (w) *w -= tl * U->val[k];
                }
            }

            // Refresh work
            for(ptrdiff_t j = row_beg; j < row_end; ++j)
                work[Acol[j]] = NULL;
        }

        this->D = Backend::copy_vector(D, bprm);
        this->L = Backend::copy_matrix(L, bprm);
        this->U = Backend::copy_matrix(U, bprm);

        if (!serial_backend::value) {
            t1 = Backend::create_vector(n, bprm);
            t2 = Backend::create_vector(n, bprm);
        }
    }

    /// \copydoc amgcl::relaxation::damped_jacobi::apply_pre
    template <class Matrix, class VectorRHS, class VectorX, class VectorTMP>
    void apply_pre(
            const Matrix &A, const VectorRHS &rhs, VectorX &x, VectorTMP &tmp,
            const params &prm
            ) const
    {
        backend::residual(rhs, A, x, tmp);
        solve(tmp, prm, serial_backend());
        backend::axpby(prm.damping, tmp, math::identity<scalar_type>(), x);
    }

    /// \copydoc amgcl::relaxation::damped_jacobi::apply_post
    template <class Matrix, class VectorRHS, class VectorX, class VectorTMP>
    void apply_post(
            const Matrix &A, const VectorRHS &rhs, VectorX &x, VectorTMP &tmp,
            const params &prm
            ) const
    {
        backend::residual(rhs, A, x, tmp);
        solve(tmp, prm, serial_backend());
        backend::axpby(prm.damping, tmp, math::identity<scalar_type>(), x);
    }

    /// \copydoc amgcl::relaxation::damped_jacobi::apply_post
    template <class Matrix, class VectorRHS, class VectorX>
    void apply(const Matrix &A, const VectorRHS &rhs, VectorX &x, const params &prm) const
    {
        backend::copy(rhs, x);
        solve(x, prm, serial_backend());
    }

    private:
        typedef typename boost::is_same<
                Backend, backend::builtin<value_type>
            >::type serial_backend;

        boost::shared_ptr<matrix> L, U;
        boost::shared_ptr<matrix_diagonal> D;
        boost::shared_ptr<vector> t1, t2;

        template <class VectorX>
        void solve(VectorX &x, const params &prm, boost::true_type) const
        {
            relaxation::detail::serial_ilu_solve(*L, *U, *D, x);
        }

        template <class VectorX>
        void solve(VectorX &x, const params &prm, boost::false_type) const
        {
            relaxation::detail::parallel_ilu_solve(
                    *L, *U, *D, x, *t1, *t2, prm.jacobi_iters);
        }

};

} // namespace relaxation
} // namespace amgcl



#endif
