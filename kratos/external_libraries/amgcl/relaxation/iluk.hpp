#ifndef AMGCL_RELAXATION_ILUK_HPP
#define AMGCL_RELAXATION_ILUK_HPP

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
 * \file   amgcl/relaxation/iluk.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Incomplete LU with fill-in level.
 */

#include <vector>
#include <deque>
#include <queue>
#include <cmath>

#include <boost/typeof/typeof.hpp>
#include <boost/foreach.hpp>

#include <amgcl/backend/builtin.hpp>
#include <amgcl/util.hpp>
#include <amgcl/relaxation/detail/ilu_solve.hpp>

namespace amgcl {
namespace relaxation {

/// ILU(k) smoother.
template <class Backend>
struct iluk {
    typedef typename Backend::value_type      value_type;
    typedef typename Backend::matrix          matrix;
    typedef typename Backend::matrix_diagonal matrix_diagonal;
    typedef typename Backend::vector          vector;

    typedef typename math::scalar_of<value_type>::type scalar_type;

    /// Relaxation parameters.
    struct params {
        /// Level of fill-in.
        int k;

        /// Damping factor.
        scalar_type damping;

        /// Number of Jacobi iterations.
        /** \note Used for approximate solution of triangular systems on parallel backends */
        unsigned jacobi_iters;

        params() : k(1), damping(1), jacobi_iters(2) {}

        params(const boost::property_tree::ptree &p)
            : AMGCL_PARAMS_IMPORT_VALUE(p, k)
            , AMGCL_PARAMS_IMPORT_VALUE(p, damping)
            , AMGCL_PARAMS_IMPORT_VALUE(p, jacobi_iters)
        {
            AMGCL_PARAMS_CHECK(p, (k)(damping)(jacobi_iters));
        }

        void get(boost::property_tree::ptree &p, const std::string &path) const {
            AMGCL_PARAMS_EXPORT_VALUE(p, path, k);
            AMGCL_PARAMS_EXPORT_VALUE(p, path, damping);
            AMGCL_PARAMS_EXPORT_VALUE(p, path, jacobi_iters);
        }
    };

    /// \copydoc amgcl::relaxation::damped_jacobi::damped_jacobi
    template <class Matrix>
    iluk( const Matrix &A, const params &prm, const typename Backend::params &bprm)
    {
        typedef typename backend::row_iterator<Matrix>::type row_iterator;
        const size_t n = backend::rows(A);

        boost::shared_ptr<build_matrix> L = boost::make_shared<build_matrix>();
        boost::shared_ptr<build_matrix> U = boost::make_shared<build_matrix>();

        L->nrows = L->ncols = n;
        L->ptr.reserve(n+1); L->ptr.push_back(0);

        L->col.reserve(backend::nonzeros(A) / 3);
        L->val.reserve(backend::nonzeros(A) / 3);

        U->nrows = U->ncols = n;
        U->ptr.reserve(n+1); U->ptr.push_back(0);

        U->col.reserve(backend::nonzeros(A) / 3);
        U->val.reserve(backend::nonzeros(A) / 3);

        std::vector<int> Ulev; Ulev.reserve(backend::nonzeros(A) / 3);

        std::vector<value_type> D;
        D.reserve(n);

        sparse_vector w(n, prm.k);

        for(ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(n); ++i) {
            w.reset(i);

            for(row_iterator a = backend::row_begin(A, i); a; ++a) {
                w.add(a.col(), a.value(), 0);
            }

            while(!w.q.empty()) {
                nonzero &a = w.next_nonzero();
                a.val = a.val * D[a.col];

                for(ptrdiff_t j = U->ptr[a.col], e = U->ptr[a.col+1]; j < e; ++j) {
                    int lev = std::max(a.lev, Ulev[j]) + 1;
                    w.add(U->col[j], -a.val * U->val[j], lev);
                }
            }

            w.sort();

            BOOST_FOREACH(const nonzero &e, w.nz) {
                if (e.col < i) {
                    L->col.push_back(e.col);
                    L->val.push_back(e.val);
                } else if (e.col == i) {
                    D.push_back(math::inverse(e.val));
                } else {
                    U->col.push_back(e.col);
                    U->val.push_back(e.val);
                    Ulev.push_back(e.lev);
                }
            }

            L->ptr.push_back(L->col.size());
            U->ptr.push_back(U->col.size());
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

        typedef typename backend::builtin<value_type>::matrix build_matrix;

        boost::shared_ptr<matrix> L, U;
        boost::shared_ptr<matrix_diagonal> D;
        boost::shared_ptr<vector> t1, t2;

        struct nonzero {
            ptrdiff_t  col;
            value_type val;
            int        lev;

            nonzero() : col(-1) {}

            nonzero(ptrdiff_t col, value_type val, int lev)
                : col(col), val(val), lev(lev) {}

            friend bool operator<(const nonzero &a, const nonzero &b) {
                return a.col < b.col;
            }
        };

        struct sparse_vector {
            struct comp_indices {
                const std::vector<nonzero> &nz;

                comp_indices(const std::vector<nonzero> &nz) : nz(nz) {}

                bool operator()(int a, int b) const {
                    return nz[a].col > nz[b].col;
                }
            };

            typedef
                std::priority_queue<int, std::vector<int>, comp_indices>
                priority_queue;

            int lfil;

            std::vector<nonzero>   nz;
            std::vector<ptrdiff_t> idx;
            priority_queue q;

            ptrdiff_t dia;

            sparse_vector(size_t n, int lfil)
                : lfil(lfil), idx(n, -1), q(comp_indices(nz)), dia(0)
            {
                nz.reserve(16);
            }

            void add(ptrdiff_t col, value_type val, int lev) {
                if (idx[col] < 0) {
                    if (lev <= lfil) {
                        int p = nz.size();
                        idx[col] = p;
                        nz.push_back(nonzero(col, val, lev));
                        if (col < dia) q.push(p);
                    }
                } else {
                    nonzero &a = nz[idx[col]];
                    a.val += val;
                    a.lev = std::min(a.lev, lev);
                }
            }

            typename std::vector<nonzero>::iterator begin() {
                return nz.begin();
            }

            typename std::vector<nonzero>::iterator end() {
                return nz.end();
            }

            nonzero& next_nonzero() {
                int p = q.top();
                q.pop();
                return nz[p];
            }

            void sort() {
                std::sort(nz.begin(), nz.end());
            }

            void reset(ptrdiff_t d) {
                BOOST_FOREACH(const nonzero &e, nz) idx[e.col] = -1;
                nz.clear();
                dia = d;
            }
        };

        template <class VectorX>
        void solve(VectorX &x, const params &prm, boost::true_type) const
        {
            relaxation::detail::serial_ilu_solve(*L, *U, *D, x);
        }

        template <class VectorX>
        void solve(VectorX &x, const params &prm, boost::false_type) const
        {
            relaxation::detail::parallel_ilu_solve(
                    *L, *U, *D, x, *t1, *t2, prm.jacobi_iters
                    );
        }
};

} // namespace relaxation
} // namespace amgcl

#endif
