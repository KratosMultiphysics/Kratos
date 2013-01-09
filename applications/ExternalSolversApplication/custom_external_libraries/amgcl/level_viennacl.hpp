#ifndef AMGCL_LEVEL_VIENNACL_HPP
#define AMGCL_LEVEL_VIENNACL_HPP

/*
The MIT License

Copyright (c) 2012 Denis Demidov <ddemidov@ksu.ru>

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
 * \file   level_vexcl.hpp
 * \author Denis Demidov <ddemidov@ksu.ru>
 * \brief  Level of an AMG hierarchy for use with VexCL vectors.
 */


#include <boost/array.hpp>
#include <boost/typeof/typeof.hpp>

#include <amgcl/common.hpp>
#include <amgcl/level_params.hpp>
#include <amgcl/spmat.hpp>
#include <amgcl/spai.hpp>
#include <amgcl/operations_viennacl.hpp>

#include <viennacl/vector.hpp>
#include <viennacl/compressed_matrix.hpp>
#include <viennacl/ell_matrix.hpp>
#include <viennacl/hyb_matrix.hpp>
#include <viennacl/linalg/inner_prod.hpp>
#include <viennacl/linalg/prod.hpp>

namespace amgcl {
namespace level {

template <gpu_matrix_format Format, typename value_type>
struct matrix_format;

template <typename value_type>
struct matrix_format<GPU_MATRIX_CRS, value_type> {
    typedef viennacl::compressed_matrix<value_type> type;
};

template <typename value_type>
struct matrix_format<GPU_MATRIX_ELL, value_type> {
    typedef viennacl::ell_matrix<value_type> type;
};

template <typename value_type>
struct matrix_format<GPU_MATRIX_HYB, value_type> {
    typedef viennacl::hyb_matrix<value_type> type;
};

struct viennacl_damped_jacobi {
    struct params {
        float damping;
        params(float w = 0.72) : damping(w) {}
    };

    template <typename value_t, typename index_t>
    struct instance {
        instance() {}

        template <class spmat>
        instance(const spmat &A) : dia(sparse::matrix_rows(A)) {
            ::viennacl::fast_copy(sparse::diagonal(A), dia);
        }

        template <class spmat, class vector>
        void apply(const spmat &A, const vector &rhs, vector &x, vector &tmp, const params &prm) const {
            tmp = ::viennacl::linalg::prod(A, x);
            tmp = rhs - tmp;
            x += prm.damping * (::viennacl::linalg::element_div(tmp, dia));
        }

        ::viennacl::vector<value_t> dia;
    };
};

struct viennacl_spai0 {
    struct params { };

    template <typename value_t, typename index_t>
    struct instance {
        instance() {}

        template <class spmat>
        instance(const spmat &A) : M(sparse::matrix_rows(A)) {
            ::viennacl::fast_copy(spai::level0(A), M);
        }

        template <class spmat, class vector>
        void apply(const spmat &A, const vector &rhs, vector &x, vector &tmp, const params &prm) const {
            tmp = ::viennacl::linalg::prod(A, x);
            tmp = rhs - tmp;
            x += ::viennacl::linalg::element_prod(M, tmp);
        }

        ::viennacl::vector<value_t> M;
    };
};

template <relax::scheme Relaxation>
struct viennacl_relax_scheme;

AMGCL_REGISTER_RELAX_SCHEME(viennacl, damped_jacobi);
AMGCL_REGISTER_RELAX_SCHEME(viennacl, spai0);

/// ViennaCL-based AMG hierarchy.
/**
 * Level of an AMG hierarchy for use with ViennaCL vectors. ViennaCL provides
 * several backends (OpenCL, CUDA, OpenMP) and is thus able to run on various
 * hardware.
 * \ingroup levels
 *
 * \param Format Matrix storage \ref gpu_matrix_format "format" to use on each
 *               level.
 * \param Relaxation Relaxation \ref relax::scheme "scheme" (smoother) to use
 *               inside V-cycles.
 */
template <
    gpu_matrix_format Format = GPU_MATRIX_HYB,
    relax::scheme Relaxation = relax::spai0
    >
struct viennacl {

/// Parameters for CPU-based level storage scheme.
struct params : public amgcl::level::params
{
    typename viennacl_relax_scheme<Relaxation>::type::params relax;
};

template <typename value_t, typename index_t = long long>
class instance {
    public:
        typedef sparse::matrix<value_t, index_t>              cpu_matrix;
        typedef typename matrix_format<Format, value_t>::type matrix;
        typedef ::viennacl::vector<value_t>                   vector;

        // Construct complete multigrid level from system matrix (a),
        // prolongation (p) and restriction (r) operators.
        // The matrices are moved into the local members.
        instance(cpu_matrix &a, cpu_matrix &p, cpu_matrix &r, const params &prm, unsigned nlevel)
            : t(a.rows), nnz(sparse::matrix_nonzeros(a)), relax(a)
        {
            ::viennacl::copy(sparse::viennacl_map(a), A);
            ::viennacl::copy(sparse::viennacl_map(p), P);
            ::viennacl::copy(sparse::viennacl_map(r), R);

            if (nlevel) {
                u.resize(a.rows);
                f.resize(a.rows);

                if (prm.kcycle && nlevel % prm.kcycle == 0)
                    for(BOOST_AUTO(v, cg.begin()); v != cg.end(); v++)
                        v->resize(a.rows);
            }

            a.clear();
            p.clear();
            r.clear();
        }

        // Construct the coarsest hierarchy level from system matrix (a) and
        // its inverse (ai).
        instance(cpu_matrix &a, cpu_matrix &ai, const params &prm, unsigned nlevel)
            : u(a.rows), f(a.rows), t(a.rows),
              nnz(sparse::matrix_nonzeros(a))
        {
            ::viennacl::copy(sparse::viennacl_map(a),  A);
            ::viennacl::copy(sparse::viennacl_map(ai), Ainv);

            a.clear();
            ai.clear();
        }

        // Returns reference to the system matrix
        const matrix& get_matrix() const {
            return A;
        }

        // Compute residual value.
        value_t resid(const vector &rhs, vector &x) const {
            t = ::viennacl::linalg::prod(A, x);
            t = rhs - t;

            return sqrt(::viennacl::linalg::inner_prod(t, t));
        }

        // Perform one V-cycle. Coarser levels are cycled recursively. The
        // coarsest level is solved directly.
        template <class Iterator>
        static void cycle(Iterator plvl, Iterator end, const params &prm,
                const vector &rhs, vector &x)
        {
            Iterator pnxt = plvl; ++pnxt;

            instance *lvl = plvl->get();
            instance *nxt = pnxt->get();

            if (pnxt != end) {
                for(unsigned j = 0; j < prm.ncycle; ++j) {
                    for(unsigned i = 0; i < prm.npre; ++i)
                        lvl->relax.apply(lvl->A, rhs, x, lvl->t, prm.relax);

                    lvl->t = ::viennacl::linalg::prod(lvl->A, x);
                    lvl->t = rhs - lvl->t;
                    nxt->f = ::viennacl::linalg::prod(lvl->R, lvl->t);
                    nxt->u.clear();

                    if (nxt->cg[0].size())
                        kcycle(pnxt, end, prm, nxt->f, nxt->u);
                    else
                        cycle(pnxt, end, prm, nxt->f, nxt->u);

                    lvl->t = ::viennacl::linalg::prod(lvl->P, nxt->u);
                    x += lvl->t;

                    for(unsigned i = 0; i < prm.npost; ++i)
                        lvl->relax.apply(lvl->A, rhs, x, lvl->t, prm.relax);
                }
            } else {
                x = ::viennacl::linalg::prod(lvl->Ainv, rhs);
            }
        }

        template <class Iterator>
        static void kcycle(Iterator plvl, Iterator end, const params &prm,
                const vector &rhs, vector &x)
        {
            Iterator pnxt = plvl; ++pnxt;

            instance *lvl = plvl->get();
            instance *nxt = pnxt->get();

            if (pnxt != end) {
                vector &r = lvl->cg[0];
                vector &s = lvl->cg[1];
                vector &p = lvl->cg[2];
                vector &q = lvl->cg[3];

                r = rhs;

                value_t rho1 = 0, rho2 = 0;

                for(int iter = 0; iter < 2; ++iter) {
                    s.clear();
                    cycle(plvl, end, prm, r, s);

                    rho2 = rho1;
                    rho1 = ::viennacl::linalg::inner_prod(r, s);

                    if (iter)
                        p = s + (rho1 / rho2) * p;
                    else
                        p = s;

                    q = ::viennacl::linalg::prod(lvl->A, p);

                    value_t alpha = rho1 / ::viennacl::linalg::inner_prod(q, p);

                    x += alpha * p;
                    r -= alpha * q;
                }
            } else {
                x = ::viennacl::linalg::prod(lvl->Ainv, rhs);
            }
        }

        index_t size() const {
            return A.size1();
        }

        index_t nonzeros() const {
            return nnz;
        }
    private:
        matrix A;
        matrix P;
        matrix R;
        matrix Ainv;

        mutable vector u;
        mutable vector f;
        mutable vector t;

        typename viennacl_relax_scheme<Relaxation>::type::template instance<value_t, index_t> relax;

        mutable boost::array<vector, 4> cg;

        index_t nnz;
};

};

} // namespace level

} // namespace amgcl

#endif
