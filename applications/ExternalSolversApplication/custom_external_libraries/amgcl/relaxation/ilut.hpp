#ifndef AMGCL_RELAXATION_ILUT_HPP
#define AMGCL_RELAXATION_ILUT_HPP

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
 * \file   amgcl/relaxation/ilut.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Incomplete LU with thresholding relaxation scheme.
 */

#include <vector>
#include <queue>
#include <cmath>

#include <boost/typeof/typeof.hpp>
#include <boost/foreach.hpp>

#include <amgcl/backend/builtin.hpp>
#include <amgcl/util.hpp>
#include <amgcl/relaxation/detail/ilu_solve.hpp>

namespace amgcl {
namespace relaxation {

/// ILUT(p, tau) smoother.
/**
 * \note ILUT is a serial algorithm and is only applicable to backends that
 * support matrix row iteration (e.g. amgcl::backend::builtin or
 * amgcl::backend::eigen).
 *
 * \param Backend Backend for temporary structures allocation.
 * \ingroup relaxation
 */
template <class Backend>
struct ilut {
    typedef typename Backend::value_type      value_type;
    typedef typename Backend::matrix          matrix;
    typedef typename Backend::matrix_diagonal matrix_diagonal;
    typedef typename Backend::vector          vector;

    typedef typename math::scalar_of<value_type>::type scalar_type;

    /// Relaxation parameters.
    struct params {
        /// Maximum fill-in.
        int p;

        /// Minimum magnitude of non-zero elements relative to the current row norm.
        scalar_type tau;

        /// Damping factor.
        scalar_type damping;

        /// Number of Jacobi iterations.
        /** \note Used for approximate solution of triangular systems on parallel backends */
        unsigned jacobi_iters;

        params() : p(2), tau(1e-2f), damping(1), jacobi_iters(2) {}

        params(const boost::property_tree::ptree &p)
            : AMGCL_PARAMS_IMPORT_VALUE(p, p)
            , AMGCL_PARAMS_IMPORT_VALUE(p, tau)
            , AMGCL_PARAMS_IMPORT_VALUE(p, damping)
            , AMGCL_PARAMS_IMPORT_VALUE(p, jacobi_iters)
        {}

        void get(boost::property_tree::ptree &p, const std::string &path) const {
            AMGCL_PARAMS_EXPORT_VALUE(p, path, p);
            AMGCL_PARAMS_EXPORT_VALUE(p, path, tau);
            AMGCL_PARAMS_EXPORT_VALUE(p, path, damping);
            AMGCL_PARAMS_EXPORT_VALUE(p, path, jacobi_iters);
        }
    };

    /// \copydoc amgcl::relaxation::damped_jacobi::damped_jacobi
    template <class Matrix>
    ilut( const Matrix &A, const params &prm, const typename Backend::params &bprm)
    {
        typedef typename backend::row_iterator<Matrix>::type row_iterator;
        const size_t n = backend::rows(A);

        BOOST_AUTO(Aptr, A.ptr_data());
        BOOST_AUTO(Acol, A.col_data());

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
        L->col.reserve(Lnz + n * prm.p);
        L->val.reserve(Lnz + n * prm.p);

        U->nrows = U->ncols = n;
        U->ptr.reserve(n+1); U->ptr.push_back(0);
        U->col.reserve(Unz + n * prm.p);
        U->val.reserve(Unz + n * prm.p);

        std::vector<value_type> D;
        D.reserve(n);

        sparse_vector w(n);

        for(ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(n); ++i) {
            w.dia = i;

            int lenL = 0;
            int lenU = 0;

            scalar_type tol = math::zero<scalar_type>();

            for(row_iterator a = backend::row_begin(A, i); a; ++a) {
                w[a.col()] = a.value();
                tol += math::norm(a.value());

                if (a.col() <  i) ++lenL;
                if (a.col() >= i) ++lenU;
            }
            tol = prm.tau / (lenL + lenU);

            while(!w.q.empty()) {
                ptrdiff_t k = w.next_nonzero();
                w[k] = D[k] * w[k];
                value_type wk = w[k];

                if (math::norm(wk) > tol) {
                    for(ptrdiff_t j = U->ptr[k]; j < U->ptr[k+1]; ++j)
                        w[U->col[j]] -= wk * U->val[j];
                }
            }

            w.move_to(lenL + prm.p, lenU + prm.p, tol, *L, *U, D);
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
        apply_dispatch(A, rhs, x, tmp, prm, serial_backend());
    }

    /// \copydoc amgcl::relaxation::damped_jacobi::apply_post
    template <class Matrix, class VectorRHS, class VectorX, class VectorTMP>
    void apply_post(
            const Matrix &A, const VectorRHS &rhs, VectorX &x, VectorTMP &tmp,
            const params &prm
            ) const
    {
        apply_dispatch(A, rhs, x, tmp, prm, serial_backend());
    }

    private:
        typedef typename boost::is_same<
                Backend, backend::builtin<value_type>
            >::type serial_backend;

        typedef typename backend::builtin<value_type>::matrix build_matrix;

        boost::shared_ptr<matrix> L, U;
        boost::shared_ptr<matrix_diagonal> D;
        boost::shared_ptr<vector> t1, t2;

        struct sparse_vector {
            struct nonzero {
                ptrdiff_t  col;
                value_type val;

                nonzero() : col(-1) {}

                nonzero(ptrdiff_t col, value_type val = value_type())
                    : col(col), val(val) {}
            };

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

            std::vector<nonzero>   nz;
            std::vector<ptrdiff_t> idx;
            priority_queue q;

            ptrdiff_t dia;

            sparse_vector(size_t n) : idx(n, -1), q(comp_indices(nz)), dia(0) {
                nz.reserve(16);
            }

            value_type operator[](ptrdiff_t i) const {
                if (idx[i] >= 0) return nz[idx[i]].val;
                return value_type();
            }

            value_type& operator[](ptrdiff_t i) {
                if (idx[i] == -1) {
                    int p = nz.size();
                    idx[i] = p;
                    nz.push_back(nonzero(i));
                    if (i < dia) q.push(p);
                }
                return nz[idx[i]].val;
            }

            typename std::vector<nonzero>::iterator begin() {
                return nz.begin();
            }

            typename std::vector<nonzero>::iterator end() {
                return nz.end();
            }

            ptrdiff_t next_nonzero() {
                int p = q.top();
                q.pop();
                return nz[p].col;
            }

            struct higher_than {
                scalar_type tol;
                ptrdiff_t   dia;

                higher_than(scalar_type tol, ptrdiff_t dia)
                    : tol(tol), dia(dia) {}

                bool operator()(const nonzero &v) const {
                    return v.col == dia || math::norm(v.val) > tol;
                }
            };

            struct L_first {
                ptrdiff_t dia;

                L_first(ptrdiff_t dia) : dia(dia) {}

                bool operator()(const nonzero &v) const {
                    return v.col < dia;
                }
            };

            struct by_abs_val {
                ptrdiff_t dia;

                by_abs_val(ptrdiff_t dia) : dia(dia) {}

                bool operator()(const nonzero &a, const nonzero &b) const {
                    if (a.col == dia) return true;
                    if (b.col == dia) return false;

                    return math::norm(a.val) > math::norm(b.val);
                }
            };

            struct by_col {
                bool operator()(const nonzero &a, const nonzero &b) const {
                    return a.col < b.col;
                }
            };

            void move_to(
                    int lp, int up, scalar_type tol,
                    build_matrix &L, build_matrix &U, std::vector<value_type> &D
                    )
            {
                typedef typename std::vector<nonzero>::iterator ptr;

                ptr b = nz.begin();
                ptr e = nz.end();

                // Move zeros to back:
                e = std::partition(b, e, higher_than(tol, dia));

                // Split L and U:
                ptr m = std::partition(b, e, L_first(dia));

                // Get largest p elements in L and U.
                ptr lend = std::min(b + lp, m);
                ptr uend = std::min(m + up, e);

                if (lend != m) std::nth_element(b, lend, m, by_abs_val(dia));
                if (uend != e) std::nth_element(m, uend, e, by_abs_val(dia));

                // Sort entries by column number
                std::sort(b, lend, by_col());
                std::sort(m, uend, by_col());

                // copy L to the output matrix.
                for(ptr a = b; a != lend; ++a) {
                    L.col.push_back(a->col);
                    L.val.push_back(a->val);
                }
                L.ptr.push_back(L.val.size());

                // Store inverted diagonal.
                D.push_back(math::inverse(m->val));
                ++m;

                // copy U to the output matrix.
                for(ptr a = m; a != uend; ++a) {
                    U.col.push_back(a->col);
                    U.val.push_back(a->val);
                }
                U.ptr.push_back(U.val.size());

                BOOST_FOREACH(const nonzero &e, nz) idx[e.col] = -1;
                nz.clear();
            }
        };

        template <class Matrix, class VectorRHS, class VectorX, class VectorTMP>
        void apply_dispatch(
                const Matrix &A, const VectorRHS &rhs, VectorX &x, VectorTMP &tmp,
                const params &prm, boost::true_type
                ) const
        {
            backend::residual(rhs, A, x, tmp);
            relaxation::detail::serial_ilu_solve(*L, *U, *D, tmp);
            backend::axpby(prm.damping, tmp, math::identity<scalar_type>(), x);
        }

        template <class Matrix, class VectorRHS, class VectorX, class VectorTMP>
        void apply_dispatch(
                const Matrix &A, const VectorRHS &rhs, VectorX &x, VectorTMP &tmp,
                const params &prm, boost::false_type
                ) const
        {
            backend::residual(rhs, A, x, tmp);
            relaxation::detail::parallel_ilu_solve(
                    *L, *U, *D, tmp, *t1, *t2, prm.jacobi_iters
                    );
            backend::axpby(prm.damping, tmp, math::identity<scalar_type>(), x);
        }
};

} // namespace relaxation
} // namespace amgcl

#endif
