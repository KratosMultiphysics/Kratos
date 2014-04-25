#ifndef AMGCL_LEVEL_VIENNACL_HPP
#define AMGCL_LEVEL_VIENNACL_HPP

/*
The MIT License

Copyright (c) 2012-2014 Denis Demidov <dennis.demidov@gmail.com>

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
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Level of an AMG hierarchy for use with VexCL vectors.
 */


#include <boost/smart_ptr/scoped_ptr.hpp>

#include <amgcl/common.hpp>
#include <amgcl/level_params.hpp>
#include <amgcl/spmat.hpp>
#include <amgcl/spai.hpp>
#include <amgcl/operations_viennacl.hpp>
#include <amgcl/gmres.hpp>

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
        instance(const spmat &A, const params&) : dia(sparse::matrix_rows(A)) {
            ::viennacl::fast_copy(sparse::diagonal(A), dia);
        }

        template <class spmat, class vector>
        void apply(const spmat &A, const vector &rhs, vector &x, vector &tmp, const params &prm) const {
            tmp = ::viennacl::linalg::prod(A, x);
            x += static_cast<typename vector::value_type>(prm.damping) *
                 ::viennacl::linalg::element_div(rhs - tmp, dia);
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
        instance(const spmat &A, const params&) : M(sparse::matrix_rows(A)) {
            ::viennacl::fast_copy(spai::level0(A), M);
        }

        template <class spmat, class vector>
        void apply(const spmat &A, const vector &rhs, vector &x, vector &tmp, const params&) const {
            tmp = ::viennacl::linalg::prod(A, x);
            tmp = rhs - tmp;
            x += ::viennacl::linalg::element_prod(M, tmp);
        }

        ::viennacl::vector<value_t> M;
    };
};

struct viennacl_chebyshev {
    struct params {
        unsigned degree;
        float    lower;

        params(unsigned degree = 5, float lower = 1.0f / 30.0f)
            : degree(degree), lower(lower)
        {}
    };

    template <typename value_t, typename index_t>
    struct instance {
        instance() {}

        template <class spmat>
        instance(const spmat &A, const params &prm)
            : p(sparse::matrix_rows(A)), q(sparse::matrix_rows(A))
        {
            value_t r = spectral_radius(A);
            C = chebyshev_coefficients(prm.degree, r * prm.lower, r);
        }

        template <class spmat, class vector>
        void apply(const spmat &A, const vector &rhs, vector &x, vector &res, const params&) const {
            res = ::viennacl::linalg::prod(A, x);
            res = rhs - res;

            p = C[0] * res;

            typedef typename std::vector<value_t>::const_iterator ci;
            for(ci c = C.begin() + 1; c != C.end(); ++c) {
                q = ::viennacl::linalg::prod(A, p);
                p = (*c) * res + q;
            }

            x += p;
        }

        std::vector<value_t> C;
        mutable ::viennacl::vector<value_t> p, q;
    };
};

template <relax::scheme Relaxation>
struct viennacl_relax_scheme;

AMGCL_REGISTER_RELAX_SCHEME(viennacl, damped_jacobi);
AMGCL_REGISTER_RELAX_SCHEME(viennacl, spai0);
AMGCL_REGISTER_RELAX_SCHEME(viennacl, chebyshev);

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
            : t(a.rows), relax(a, prm.relax), nnz(sparse::matrix_nonzeros(a))
        {
            ::viennacl::copy(sparse::viennacl_map(a), A);
            ::viennacl::copy(sparse::viennacl_map(p), P);
            ::viennacl::copy(sparse::viennacl_map(r), R);

            if (nlevel) {
                u.resize(a.rows);
                f.resize(a.rows);

                if (prm.kcycle && nlevel % prm.kcycle == 0)
                    gmres.reset(new gmres_data<vector>(prm.kcycle_iterations, a.rows));
            }

            a.clear();
            p.clear();
            r.clear();
        }

        // Construct the coarsest hierarchy level from system matrix (a) and
        // its inverse (ai).
        instance(cpu_matrix &a, cpu_matrix &ai, const params&, unsigned /*nlevel*/)
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

                    if (nxt->gmres)
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

            if (pnxt != end) {
                cycle_precond<Iterator> p(plvl, end, prm);

                lvl->gmres->restart(lvl->A, rhs, p, x);

                for(int i = 0; i < lvl->gmres->M; ++i)
                    lvl->gmres->iteration(lvl->A, p, i);

                lvl->gmres->update(x, lvl->gmres->M - 1);
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

        mutable boost::scoped_ptr< gmres_data<vector> > gmres;

        index_t nnz;

        template <class Iterator>
        struct cycle_precond {
            cycle_precond(Iterator lvl, Iterator end, const params &prm)
                : lvl(lvl), end(end), prm(prm) {}

            void apply(const vector &r, vector &x) const {
                cycle(lvl, end, prm, r, x);
            }

            Iterator lvl, end;
            const params &prm;
        };
};

};

} // namespace level

} // namespace amgcl

#endif
