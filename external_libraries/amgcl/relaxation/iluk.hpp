#ifndef AMGCL_RELAXATION_ILUK_HPP
#define AMGCL_RELAXATION_ILUK_HPP

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
 * \file   amgcl/relaxation/iluk.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Incomplete LU with fill-in level.
 */

#include <vector>
#include <deque>
#include <queue>
#include <cmath>


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

    typedef detail::ilu_solve<Backend> ilu_solve;

    /// Relaxation parameters.
    struct params {
        /// Level of fill-in.
        int k;

        /// Damping factor.
        scalar_type damping;

        /// Parameters for sparse triangular system solver
        typename ilu_solve::params solve;

        params() : k(1), damping(1) {}

#ifndef AMGCL_NO_BOOST
        params(const boost::property_tree::ptree &p)
            : AMGCL_PARAMS_IMPORT_VALUE(p, k)
            , AMGCL_PARAMS_IMPORT_VALUE(p, damping)
            , AMGCL_PARAMS_IMPORT_CHILD(p, solve)
        {
            check_params(p, {"k", "damping", "solve"});
        }

        void get(boost::property_tree::ptree &p, const std::string &path) const {
            AMGCL_PARAMS_EXPORT_VALUE(p, path, k);
            AMGCL_PARAMS_EXPORT_VALUE(p, path, damping);
            AMGCL_PARAMS_EXPORT_CHILD(p, path, solve);
        }
#endif
    } prm;

    /// \copydoc amgcl::relaxation::damped_jacobi::damped_jacobi
    template <class Matrix>
    iluk( const Matrix &A, const params &prm, const typename Backend::params &bprm)
      : prm(prm)
    {
        typedef typename backend::builtin<value_type>::matrix build_matrix;

        const size_t n = backend::rows(A);

        size_t Anz = backend::nonzeros(A);

        std::vector<ptrdiff_t>  Lptr; Lptr.reserve(n+1); Lptr.push_back(0);
        std::vector<ptrdiff_t>  Lcol; Lcol.reserve(Anz / 3);
        std::vector<value_type> Lval; Lval.reserve(Anz / 3);

        std::vector<ptrdiff_t>  Uptr; Uptr.reserve(n+1); Uptr.push_back(0);
        std::vector<ptrdiff_t>  Ucol; Ucol.reserve(Anz / 3);
        std::vector<value_type> Uval; Uval.reserve(Anz / 3);

        std::vector<int> Ulev; Ulev.reserve(Anz / 3);

        auto D = std::make_shared<backend::numa_vector<value_type> >(n, false);

        sparse_vector w(n, prm.k);

        for(ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(n); ++i) {
            w.reset(i);

            for(auto a = backend::row_begin(A, i); a; ++a) {
                w.add(a.col(), a.value(), 0);
            }

            while(!w.q.empty()) {
                nonzero &a = w.next_nonzero();
                a.val = a.val * (*D)[a.col];

                for(ptrdiff_t j = Uptr[a.col], e = Uptr[a.col+1]; j < e; ++j) {
                    int lev = std::max(a.lev, Ulev[j]) + 1;
                    w.add(Ucol[j], -a.val * Uval[j], lev);
                }
            }

            w.sort();

            for(const nonzero &e : w.nz) {
                if (e.col < i) {
                    Lcol.push_back(e.col);
                    Lval.push_back(e.val);
                } else if (e.col == i) {
                    (*D)[i] = math::inverse(e.val);
                } else {
                    Ucol.push_back(e.col);
                    Uval.push_back(e.val);
                    Ulev.push_back(e.lev);
                }
            }

            Lptr.push_back(Lcol.size());
            Uptr.push_back(Ucol.size());
        }

        ilu = std::make_shared<ilu_solve>(
                std::make_shared<build_matrix>(n, n, Lptr, Lcol, Lval),
                std::make_shared<build_matrix>(n, n, Uptr, Ucol, Uval),
                D, prm.solve, bprm);
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
        std::shared_ptr<ilu_solve> ilu;

        struct nonzero {
            ptrdiff_t  col;
            value_type val;
            int        lev;

            nonzero() : col(-1) {}

            nonzero(ptrdiff_t col, const value_type &val, int lev)
                : col(col), val(val), lev(lev) {}

            friend bool operator<(const nonzero &a, const nonzero &b) {
                return a.col < b.col;
            }
        };

        struct sparse_vector {
            struct comp_indices {
                const std::deque<nonzero> &nz;

                comp_indices(const std::deque<nonzero> &nz) : nz(nz) {}

                bool operator()(int a, int b) const {
                    return nz[a].col > nz[b].col;
                }
            };

            typedef
                std::priority_queue<int, std::vector<int>, comp_indices>
                priority_queue;

            int lfil;

            std::deque<nonzero>    nz;
            std::vector<ptrdiff_t> idx;
            priority_queue q;

            ptrdiff_t dia;

            sparse_vector(size_t n, int lfil)
                : lfil(lfil), idx(n, -1), q(comp_indices(nz)), dia(0)
            {}

            void add(ptrdiff_t col, const value_type &val, int lev) {
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
                for(const nonzero &e : nz) idx[e.col] = -1;
                nz.clear();
                dia = d;
            }
        };
};

} // namespace relaxation

namespace backend {

template <class Backend>
struct bytes_impl< relaxation::iluk<Backend> > {
    static size_t get(const relaxation::iluk<Backend> &R) {
        return R.bytes();
    }
};

} // namespace backend
} // namespace amgcl

#endif
