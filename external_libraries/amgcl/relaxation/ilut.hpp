#ifndef AMGCL_RELAXATION_ILUT_HPP
#define AMGCL_RELAXATION_ILUT_HPP

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
 * \file   amgcl/relaxation/ilut.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Incomplete LU with thresholding relaxation scheme.
 */

#include <vector>
#include <queue>
#include <cmath>


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

    typedef detail::ilu_solve<Backend> ilu_solve;

    /// Relaxation parameters.
    struct params {
        /// Fill factor.
        scalar_type p;

        /// Minimum magnitude of non-zero elements relative to the current row norm.
        scalar_type tau;

        /// Damping factor.
        scalar_type damping;

        /// Parameters for sparse triangular system solver
        typename ilu_solve::params solve;

        params() : p(2), tau(1e-2f), damping(1) {}

#ifndef AMGCL_NO_BOOST
        params(const boost::property_tree::ptree &p)
            : AMGCL_PARAMS_IMPORT_VALUE(p, p)
            , AMGCL_PARAMS_IMPORT_VALUE(p, tau)
            , AMGCL_PARAMS_IMPORT_VALUE(p, damping)
            , AMGCL_PARAMS_IMPORT_CHILD(p, solve)
        {
            check_params(p, {"p", "tau", "damping", "solve"});
        }

        void get(boost::property_tree::ptree &p, const std::string &path) const {
            AMGCL_PARAMS_EXPORT_VALUE(p, path, p);
            AMGCL_PARAMS_EXPORT_VALUE(p, path, tau);
            AMGCL_PARAMS_EXPORT_VALUE(p, path, damping);
            AMGCL_PARAMS_EXPORT_CHILD(p, path, solve);
        }
#endif
    } prm;

    /// \copydoc amgcl::relaxation::damped_jacobi::damped_jacobi
    template <class Matrix>
    ilut( const Matrix &A, const params &prm, const typename Backend::params &bprm)
      : prm(prm)
    {
        const size_t n = backend::rows(A);

        size_t Lnz = 0, Unz = 0;

        for(ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(n); ++i) {
            ptrdiff_t row_beg = A.ptr[i];
            ptrdiff_t row_end = A.ptr[i + 1];

            int lenL = 0, lenU = 0;
            for(ptrdiff_t j = row_beg; j < row_end; ++j) {
                ptrdiff_t c = A.col[j];
                if (c < i)
                    ++lenL;
                else if (c > i)
                    ++lenU;
            }

            Lnz += static_cast<size_t>(lenL * prm.p);
            Unz += static_cast<size_t>(lenU * prm.p);
        }

        auto L = std::make_shared<build_matrix>();
        auto U = std::make_shared<build_matrix>();

        L->set_size(n, n); L->set_nonzeros(Lnz); L->ptr[0] = 0;
        U->set_size(n, n); U->set_nonzeros(Unz); U->ptr[0] = 0;

        auto D = std::make_shared<backend::numa_vector<value_type> >(n, false);

        sparse_vector w(n);

        for(ptrdiff_t i = 0, Lhead = 0, Uhead = 0; i < static_cast<ptrdiff_t>(n); ++i) {
            w.dia = i;

            int lenL = 0;
            int lenU = 0;

            scalar_type tol = math::zero<scalar_type>();

            for(auto a = backend::row_begin(A, i); a; ++a) {
                w[a.col()] = a.value();
                tol += math::norm(a.value());

                if (a.col() < i) ++lenL;
                if (a.col() > i) ++lenU;
            }
            tol *= prm.tau / (lenL + lenU);

            while(!w.q.empty()) {
                ptrdiff_t k = w.next_nonzero();
                w[k] = w[k] * (*D)[k];
                value_type wk = w[k];

                if (math::norm(wk) > tol) {
                    for(ptrdiff_t j = U->ptr[k]; j < U->ptr[k+1]; ++j)
                        w[U->col[j]] -= wk * U->val[j];
                }
            }

            w.move_to(
                    static_cast<int>(lenL * prm.p),
                    static_cast<int>(lenU * prm.p),
                    tol, Lhead, *L, Uhead, *U, *D
                    );

            L->ptr[i+1] = Lhead;
            U->ptr[i+1] = Uhead;
        }

        L->nnz = L->ptr[n];
        U->nnz = U->ptr[n];

        ilu = std::make_shared<ilu_solve>(L, U, D, prm.solve, bprm);
    }

    /// \copydoc amgcl::relaxation::damped_jacobi::apply_pre
    template <class Matrix, class VectorRHS, class VectorX, class VectorTMP>
    void apply_pre(
            const Matrix &A, const VectorRHS &rhs, VectorX &x, VectorTMP &tmp
            ) const
    {
        backend::residual(rhs, A, x, tmp);
        ilu->solve(tmp);
        backend::axpby(prm.damping, tmp, math::identity<scalar_type>(), x);
    }

    /// \copydoc amgcl::relaxation::damped_jacobi::apply_post
    template <class Matrix, class VectorRHS, class VectorX, class VectorTMP>
    void apply_post(
            const Matrix &A, const VectorRHS &rhs, VectorX &x, VectorTMP &tmp
            ) const
    {
        backend::residual(rhs, A, x, tmp);
        ilu->solve(tmp);
        backend::axpby(prm.damping, tmp, math::identity<scalar_type>(), x);
    }

    template <class Matrix, class VectorRHS, class VectorX>
    void apply(const Matrix&, const VectorRHS &rhs, VectorX &x) const
    {
        backend::copy(rhs, x);
        ilu->solve(x);
    }

    size_t bytes() const {
        return ilu->bytes();
    }

    private:
        typedef typename backend::builtin<value_type>::matrix build_matrix;
        std::shared_ptr<ilu_solve> ilu;

        struct sparse_vector {
            struct nonzero {
                ptrdiff_t  col;
                value_type val;

                nonzero() : col(-1) {}

                nonzero(ptrdiff_t col, const value_type &val = value_type())
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
                    ptrdiff_t &Lhead, build_matrix &L,
                    ptrdiff_t &Uhead, build_matrix &U,
                    backend::numa_vector<value_type> &D
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
                    L.col[Lhead] = a->col;
                    L.val[Lhead] = a->val;

                    ++Lhead;
                }

                // Store inverted diagonal.
                D[dia] = math::inverse(m->val);

                if (m != uend) {
                    ++m;

                    // copy U to the output matrix.
                    for(ptr a = m; a != uend; ++a) {
                        U.col[Uhead] = a->col;
                        U.val[Uhead] = a->val;

                        ++Uhead;
                    }
                }

                for(const nonzero &e : nz) idx[e.col] = -1;
                nz.clear();
            }
        };
};

} // namespace relaxation
} // namespace amgcl

#endif
