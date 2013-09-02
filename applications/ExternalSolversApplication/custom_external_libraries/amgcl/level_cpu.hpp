#ifndef AMGCL_LEVEL_CPU_HPP
#define AMGCL_LEVEL_CPU_HPP

/*
The MIT License

Copyright (c) 2012-2013 Denis Demidov <ddemidov@ksu.ru>

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
 * \file   level_cpu.hpp
 * \author Denis Demidov <ddemidov@ksu.ru>
 * \brief  Level of an AMG hierarchy for use with arrays located in main (CPU) memory.
 */

#include <vector>

#include <boost/array.hpp>

#include <amgcl/common.hpp>
#include <amgcl/level_params.hpp>
#include <amgcl/spmat.hpp>
#include <amgcl/spai.hpp>
#include <amgcl/tictoc.hpp>

namespace amgcl {

/// Storage schemes for levels of AMG hierarchy.
namespace level {

/**
 * \defgroup levels Level storage backends
 * \brief Provided storage and acceleration backends.
 */

struct cpu_damped_jacobi {
    struct params {
        float damping;
        params(float w = 0.72) : damping(w) {}
    };

    template <typename value_t, typename index_t>
    struct instance {
        instance() {}

        template <class spmat>
        instance(const spmat&) {}

        template <class spmat, class vector1, class vector2, class vector3>
        void apply_pre(const spmat &A, const vector1 &rhs, vector2 &x, vector3 &tmp, const params &prm) const {
            const index_t n = sparse::matrix_rows(A);

            const index_t *Arow = sparse::matrix_outer_index(A);
            const index_t *Acol = sparse::matrix_inner_index(A);
            const value_t *Aval = sparse::matrix_values(A);

#pragma omp parallel for schedule(dynamic, 1024)
            for(index_t i = 0; i < n; ++i) {
                value_t temp = rhs[i];
                value_t diag = 1;

                for(index_t j = Arow[i], e = Arow[i + 1]; j < e; ++j) {
                    index_t c = Acol[j];
                    value_t v = Aval[j];

                    temp -= v * x[c];

                    if (c == i) diag = v;
                }

                tmp[i] = x[i] + prm.damping * (temp / diag);
            }

            vector_copy(tmp, x);
        }

        template <class spmat, class vector1, class vector2, class vector3>
        void apply_post(const spmat &A, const vector1 &rhs, vector2 &x, vector3 &tmp, const params &prm) const {
            apply_pre(A, rhs, x, tmp, prm);
        }

        template <class U>
        inline static void vector_copy(U &u, U &v) {
            using namespace std;
            swap(u, v);
        }

        template <class U, class V>
        inline static void vector_copy(U &u, V &v) {
            std::copy(u.begin(), u.end(), &v[0]);
        }
    };
};

struct cpu_gauss_seidel {
    struct params { };

    template <typename value_t, typename index_t>
    struct instance {
        instance() {}

        template <class spmat>
        instance(const spmat&) {}

#define GS_INNER_LOOP \
                value_t temp = rhs[i]; \
                value_t diag = 1; \
                for(index_t j = Arow[i], e = Arow[i + 1]; j < e; ++j) { \
                    index_t c = Acol[j]; \
                    value_t v = Aval[j]; \
                    if (c == i) \
                        diag = v; \
                    else \
                        temp -= v * x[c]; \
                } \
                x[i] = temp / diag

        template <class spmat, class vector1, class vector2, class vector3>
        void apply_pre(const spmat &A, const vector1 &rhs, vector2 &x, vector3& /*tmp*/, const params&) const {
            const index_t n = sparse::matrix_rows(A);

            const index_t *Arow = sparse::matrix_outer_index(A);
            const index_t *Acol = sparse::matrix_inner_index(A);
            const value_t *Aval = sparse::matrix_values(A);

            for(index_t i = 0; i < n; ++i) {
                GS_INNER_LOOP;
            }
        }

        template <class spmat, class vector1, class vector2, class vector3>
        void apply_post(const spmat &A, const vector1 &rhs, vector2 &x, vector3& /*tmp*/, const params&) const {
            const index_t n = A.rows;

            const index_t *Arow = sparse::matrix_outer_index(A);
            const index_t *Acol = sparse::matrix_inner_index(A);
            const value_t *Aval = sparse::matrix_values(A);

            for(index_t i = n - 1; i >= 0; --i) {
                GS_INNER_LOOP;
            }
        }

#undef GS_INNER_LOOP
    };
};

struct cpu_ilu0 {
    struct params {
        float damping;
        params(float w = 0.72) : damping(w) {}
    };

    template <typename value_t, typename index_t>
    struct instance {
        instance() {}

        template <class spmat>
        instance(spmat &A) : diag(sparse::matrix_rows(A)) {
            sparse::sort_rows(A);

            const index_t n = sparse::matrix_rows(A);

            const index_t *Arow = sparse::matrix_outer_index(A);
            const index_t *Acol = sparse::matrix_inner_index(A);
            const value_t *Aval = sparse::matrix_values(A);

            luval.assign(Aval, Aval + Arow[n]);

            std::vector<index_t> work(n, -1);

            for(index_t i = 0; i < n; ++i) {
                index_t row_beg = Arow[i];
                index_t row_end = Arow[i + 1];

                for(index_t j = row_beg; j < row_end; ++j)
                    work[Acol[j]] = j;

                for(index_t j = row_beg; j < row_end; ++j) {
                    index_t c = Acol[j];

                    // Exit if diagonal is reached
                    if (c >= i) {
                        if (c != i)
                            throw std::runtime_error("No diagonal?");
                        if (fabs(luval[j]) < 1e-32)
                            throw std::runtime_error("Zero pivot in ILU");

                        diag[i]  = j;
                        luval[j] = 1 / luval[j];
                        break;
                    }

                    // Compute the multiplier for jrow
                    value_t tl = luval[j] * luval[diag[c]];
                    luval[j] = tl;

                    // Perform linear combination
                    for(index_t k = diag[c] + 1; k < Arow[c + 1]; ++k) {
                        index_t w = work[Acol[k]];
                        if (w >= 0) luval[w] -= tl * luval[k];
                    }
                }

                // Refresh work
                for(index_t j = row_beg; j < row_end; ++j)
                    work[Acol[j]] = -1;
            }
        }

        template <class spmat, class vector1, class vector2, class vector3>
        void apply_pre(const spmat &A, const vector1 &rhs, vector2 &x, vector3 &tmp, const params &prm) const {
            const index_t n = sparse::matrix_rows(A);

            const index_t *Arow = sparse::matrix_outer_index(A);
            const index_t *Acol = sparse::matrix_inner_index(A);
            const value_t *Aval = sparse::matrix_values(A);

#pragma omp parallel for schedule(dynamic, 1024)
            for(index_t i = 0; i < n; i++) {
                value_t buf = rhs[i];
                for(index_t j = Arow[i], e = Arow[i + 1]; j < e; ++j)
                    buf -= Aval[j] * x[Acol[j]];
                tmp[i] = buf;
            }

            for(index_t i = 0; i < n; i++) {
                for(index_t j = Arow[i], e = diag[i]; j < e; ++j)
                    tmp[i] -= luval[j] * tmp[Acol[j]];
            }

            for(index_t i = n - 1; i >= 0; --i) {
                for(index_t j = diag[i] + 1, e = Arow[i + 1]; j < e; ++j)
                    tmp[i] -= luval[j] * tmp[Acol[j]];
                tmp[i] *= luval[diag[i]];
            }

#pragma omp parallel for schedule(dynamic, 1024)
            for(index_t i = 0; i < n; i++) x[i] += prm.damping * tmp[i];
        }

        template <class spmat, class vector1, class vector2, class vector3>
        void apply_post(const spmat &A, const vector1 &rhs, vector2 &x, vector3 &tmp, const params &prm) const {
            apply_pre(A, rhs, x, tmp, prm);
        }

        std::vector<value_t> luval;
        std::vector<index_t> diag;
    };
};

struct cpu_spai0 {
    struct params { };

    template <typename value_t, typename index_t>
    struct instance {
        instance() {}

        template <class spmat>
        instance(const spmat &A) : m(spai::level0(A)) { }

        template <class spmat, class vector1, class vector2, class vector3>
        void apply_pre(const spmat &A, const vector1 &rhs, vector2 &x, vector3 &tmp, const params&) const {
            const index_t n = sparse::matrix_rows(A);

            const index_t *Arow = sparse::matrix_outer_index(A);
            const index_t *Acol = sparse::matrix_inner_index(A);
            const value_t *Aval = sparse::matrix_values(A);

#pragma omp parallel for schedule(dynamic, 1024)
            for(index_t i = 0; i < n; i++) {
                value_t buf = rhs[i];
                for(index_t j = Arow[i], e = Arow[i + 1]; j < e; ++j)
                    buf -= Aval[j] * x[Acol[j]];
                tmp[i] = buf;
            }

#pragma omp parallel for schedule(dynamic, 1024)
            for(index_t i = 0; i < n; ++i) {
                x[i] += m[i] * tmp[i];
            }
        }

        template <class spmat, class vector1, class vector2, class vector3>
        void apply_post(const spmat &A, const vector1 &rhs, vector2 &x, vector3 &tmp, const params &prm) const {
            apply_pre(A, rhs, x, tmp, prm);
        }


        std::vector<value_t> m;
    };
};

template <relax::scheme Relaxation>
struct cpu_relax_scheme;

AMGCL_REGISTER_RELAX_SCHEME(cpu, damped_jacobi);
AMGCL_REGISTER_RELAX_SCHEME(cpu, gauss_seidel);
AMGCL_REGISTER_RELAX_SCHEME(cpu, ilu0);
AMGCL_REGISTER_RELAX_SCHEME(cpu, spai0);

/// CPU-based AMG hierarchy.
/**
 * Level of an AMG hierarchy for use with arrays located in main (CPU) memory.
 * Employs OpenMP parallelization.
 * \ingroup levels
 *
 * \param Relaxation Relaxation \ref relax::scheme "scheme" (smoother) to use
 *                   inside V-cycles.
 */
template <relax::scheme Relaxation = relax::damped_jacobi>
struct cpu {

/// Parameters for CPU-based level storage scheme.
struct params
    : public amgcl::level::params
{
    typename cpu_relax_scheme<Relaxation>::type::params relax;
};

template <typename value_t, typename index_t>
class instance {
    public:
        typedef sparse::matrix<value_t, index_t> matrix;

        // Construct complete multigrid level from system matrix (a),
        // prolongation (p) and restriction (r) operators.
        // The matrices are moved into the local members.
        instance(matrix &a, matrix &p, matrix &r, const params &prm, unsigned nlevel)
            : t(a.rows), relax(a)
        {
            A.swap(a);
            P.swap(p);
            R.swap(r);

            if (nlevel) {
                u.resize(A.rows);
                f.resize(A.rows);

                if (prm.kcycle && nlevel % prm.kcycle == 0)
                    for(size_t i = 0; i < cg.size(); ++i) cg[i].resize(A.rows);
            }

            t.resize(A.rows);
        }

        // Construct the coarsest hierarchy level from system matrix (a) and
        // its inverse (ai).
        instance(matrix &a, matrix &ai, const params&, unsigned /*nlevel*/)
            : u(a.rows), f(a.rows), t(a.rows)
        {
            A.swap(a);
            Ai.swap(ai);
        }

        // Returns reference to the system matrix
        const matrix& get_matrix() const {
            return A;
        }

        // Compute residual value.
        template <class vector1, class vector2>
        value_t resid(const vector1 &rhs, vector2 &x) const {
            TIC("residual");
            const index_t n = A.rows;
            value_t norm = 0;

#pragma omp parallel for reduction(+:norm) schedule(dynamic, 1024)
            for(index_t i = 0; i < n; ++i) {
                value_t temp = rhs[i];

                for(index_t j = A.row[i], e = A.row[i + 1]; j < e; ++j)
                    temp -= A.val[j] * x[A.col[j]];

                norm += temp * temp;
            }

            TOC("residual");
            return sqrt(norm);
        }

        // Perform one V-cycle. Coarser levels are cycled recursively. The
        // coarsest level is solved directly.
        template <class Iterator, class vector1, class vector2>
        static void cycle(Iterator plvl, Iterator end, const params &prm,
                const vector1 &rhs, vector2 &x)
        {
            Iterator pnxt = plvl; ++pnxt;

            instance *lvl = plvl->get();
            instance *nxt = pnxt->get();

            const index_t n = lvl->A.rows;

            if (pnxt != end) {
                const index_t nc = nxt->A.rows;

                for(unsigned j = 0; j < prm.ncycle; ++j) {
                    TIC("relax");
                    for(unsigned i = 0; i < prm.npre; ++i)
                        lvl->relax.apply_pre(lvl->A, rhs, x, lvl->t, prm.relax);
                    TOC("relax");

                    //lvl->t = rhs - lvl->A * x;
                    TIC("residual");
#pragma omp parallel for schedule(dynamic, 1024)
                    for(index_t i = 0; i < n; ++i) {
                        value_t temp = rhs[i];

                        for(index_t j = lvl->A.row[i], e = lvl->A.row[i + 1]; j < e; ++j)
                            temp -= lvl->A.val[j] * x[lvl->A.col[j]];

                        lvl->t[i] = temp;
                    }
                    TOC("residual");

                    //nxt->f = lvl->R * lvl->t;
                    TIC("restrict");
#pragma omp parallel for schedule(dynamic, 1024)
                    for(index_t i = 0; i < nc; ++i) {
                        value_t temp = 0;

                        for(index_t j = lvl->R.row[i], e = lvl->R.row[i + 1]; j < e; ++j)
                            temp += lvl->R.val[j] * lvl->t[lvl->R.col[j]];

                        nxt->f[i] = temp;
                    }
                    TOC("restrict");

                    std::fill(nxt->u.begin(), nxt->u.end(), static_cast<value_t>(0));

                    if (nxt->cg[0].empty())
                        cycle(pnxt, end, prm, nxt->f, nxt->u);
                    else
                        kcycle(pnxt, end, prm, nxt->f, nxt->u);

                    //x += lvl->P * nxt->u;
                    TIC("prolongate");
#pragma omp parallel for schedule(dynamic, 1024)
                    for(index_t i = 0; i < n; ++i) {
                        value_t temp = 0;

                        for(index_t j = lvl->P.row[i], e = lvl->P.row[i + 1]; j < e; ++j)
                            temp += lvl->P.val[j] * nxt->u[lvl->P.col[j]];

                        x[i] += temp;
                    }
                    TOC("prolongate");

                    TIC("relax");
                    for(unsigned i = 0; i < prm.npost; ++i)
                        lvl->relax.apply_post(lvl->A, rhs, x, lvl->t, prm.relax);
                    TOC("relax");
                }
            } else {
                TIC("coarse");
#pragma omp parallel for
                for(index_t i = 0; i < n; ++i) {
                    value_t temp = 0;
                    for(index_t j = lvl->Ai.row[i], e = lvl->Ai.row[i + 1]; j < e; ++j)
                        temp += lvl->Ai.val[j] * rhs[lvl->Ai.col[j]];
                    x[i] = temp;
                }
                TOC("coarse");
            }
        }

        template <class Iterator, class vector1, class vector2>
        static void kcycle(Iterator plvl, Iterator end, const params &prm,
                const vector1 &rhs, vector2 &x)
        {
            Iterator pnxt = plvl; ++pnxt;

            instance *lvl = plvl->get();

            const index_t n = lvl->A.rows;

            if (pnxt != end) {
                std::vector<value_t> &r = lvl->cg[0];
                std::vector<value_t> &s = lvl->cg[1];
                std::vector<value_t> &p = lvl->cg[2];
                std::vector<value_t> &q = lvl->cg[3];

                std::copy(&rhs[0], &rhs[0] + n, &r[0]);

                value_t rho1 = 0, rho2 = 0;

                for(int iter = 0; iter < 2; ++iter) {
                    std::fill(&s[0], &s[0] + n, static_cast<value_t>(0));
                    cycle(plvl, end, prm, r, s);

                    TIC("kcycle");
                    rho2 = rho1;
                    rho1 = lvl->inner_prod(r, s);

                    if (iter) {
                        value_t beta = rho1 / rho2;
#pragma omp parallel for schedule(dynamic, 1024)
                        for(index_t i = 0; i < n; ++i) {
                            p[i] = s[i] + beta * p[i];
                        }
                    } else {
                        std::copy(&s[0], &s[0] + n, &p[0]);
                    }

#pragma omp parallel for schedule(dynamic, 1024)
                    for(index_t i = 0; i < n; ++i) {
                        value_t temp = 0;

                        for(index_t j = lvl->A.row[i], e = lvl->A.row[i + 1]; j < e; ++j)
                            temp += lvl->A.val[j] * p[lvl->A.col[j]];

                        q[i] = temp;
                    }

                    value_t alpha = rho1 / lvl->inner_prod(q, p);

#pragma omp parallel for schedule(dynamic, 1024)
                    for(index_t i = 0; i < n; ++i) {
                        x[i] += alpha * p[i];
                        r[i] -= alpha * q[i];
                    }
                    TOC("kcycle");
                }
            } else {
                TIC("coarse");
#pragma omp parallel for
                for(index_t i = 0; i < n; ++i) {
                    value_t temp = 0;
                    for(index_t j = lvl->Ai.row[i], e = lvl->Ai.row[i + 1]; j < e; ++j)
                        temp += lvl->Ai.val[j] * rhs[lvl->Ai.col[j]];
                    x[i] = temp;
                }
                TOC("coarse");
            }
        }

        index_t size() const {
            return t.size();
        }

        index_t nonzeros() const {
            return sparse::matrix_nonzeros(A);
        }
    private:
        matrix A;
        matrix P;
        matrix R;

        matrix Ai;

        mutable std::vector<value_t> u;
        mutable std::vector<value_t> f;
        mutable std::vector<value_t> t;

        typename cpu_relax_scheme<Relaxation>::type::template instance<value_t, index_t> relax;

        mutable boost::array<std::vector<value_t>, 4> cg;

        template <class vector1, class vector2>
        value_t inner_prod(const vector1 &x, const vector2 &y) const {
            const index_t n = A.rows;

            value_t sum = 0;

#pragma omp parallel for reduction(+:sum) schedule(dynamic, 1024)
            for(index_t i = 0; i < n; ++i)
                sum += x[i] * y[i];

            return sum;
        }
};

};

} // namespace level
} // namespace amgcl

#endif
