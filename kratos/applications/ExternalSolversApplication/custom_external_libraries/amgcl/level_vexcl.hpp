#ifndef AMGCL_LEVEL_VEXCL_HPP
#define AMGCL_LEVEL_VEXCL_HPP

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
 * \file   level_vexcl.hpp
 * \author Denis Demidov <ddemidov@ksu.ru>
 * \brief  Level of an AMG hierarchy for use with VexCL vectors.
 */

#include <array>
#include <memory>

#include <amgcl/level_params.hpp>
#include <amgcl/spmat.hpp>
#include <amgcl/spai.hpp>
#include <amgcl/operations_vexcl.hpp>
#include <amgcl/gmres.hpp>

#include <vexcl/vexcl.hpp>

namespace amgcl {
namespace level {

struct vexcl_damped_jacobi {
    struct params {
        float damping;
        params(float w = 0.72) : damping(w) {}
    };

    template <typename value_t, typename index_t>
    struct instance {
        instance() {}

        template <class spmat>
        instance(const std::vector<cl::CommandQueue> &queue, const spmat &A)
            : dia(queue, sparse::diagonal(A))
        {}

        template <class spmat, class vector>
        void apply(const spmat &A, const vector &rhs, vector &x, vector &tmp, const params &prm) const {
            tmp = rhs - A * x;
            x += prm.damping * tmp / dia;
        }

        vex::vector<value_t> dia;
    };
};

struct vexcl_spai0 {
    struct params { };

    template <typename value_t, typename index_t>
    struct instance {
        instance() {}

        template <class spmat>
        instance(const std::vector<cl::CommandQueue> &queue, const spmat &A)
            : M(queue, spai::level0(A))
        { }

        template <class spmat, class vector>
        void apply(const spmat &A, const vector &rhs, vector &x, vector &tmp, const params&) const {
            tmp = rhs - A * x;
            x += M * tmp;
        }

        vex::vector<value_t> M;
    };
};

template <relax::scheme Relaxation>
struct vexcl_relax_scheme;

AMGCL_REGISTER_RELAX_SCHEME(vexcl, damped_jacobi);
AMGCL_REGISTER_RELAX_SCHEME(vexcl, spai0);

/// VexCL-based AMG hierarchy.
/**
 * Level of an AMG hierarchy for use with VexCL vectors. Uses OpenCL technology
 * for acceleration and is able to run on OpenCL-compatible devices (GPUs or
 * CPUs).
 * \ingroup levels
 *
 * \param Relaxation Relaxation \ref relax::scheme "scheme" (smoother) to use
 *                   inside V-cycles.
 */
template <relax::scheme Relaxation = relax::damped_jacobi>
struct vexcl {

/// Parameters for VexCL-based level storage scheme.
struct params : public amgcl::level::params
{
    vex::Context *ctx;  ///< VexCL Context for VexCL objects creation.
    typename vexcl_relax_scheme<Relaxation>::type::params relax;

    params() : ctx(0) { }
};

template <typename value_t, typename index_t = long long>
class instance {
    public:
        typedef sparse::matrix<value_t, index_t>      cpu_matrix;
        typedef vex::SpMat<value_t, index_t, index_t> matrix;
        typedef vex::vector<value_t>                  vector;

        // Construct complete multigrid level from system matrix (a),
        // prolongation (p) and restriction (r) operators.
        // The matrices are moved into the local members.
        instance(cpu_matrix &a, cpu_matrix &p, cpu_matrix &r, const params &prm, unsigned nlevel)
            : A(prm.ctx ? prm.ctx->queue() : vex::StaticContext<>::get().queue(),
                    a.rows, a.cols, a.row.data(), a.col.data(), a.val.data()),
              P(prm.ctx ? prm.ctx->queue() : vex::StaticContext<>::get().queue(),
                      p.rows, p.cols, p.row.data(), p.col.data(), p.val.data()),
              R(prm.ctx ? prm.ctx->queue() : vex::StaticContext<>::get().queue(),
                          r.rows, r.cols, r.row.data(), r.col.data(), r.val.data()),
              t(prm.ctx ? prm.ctx->queue() : vex::StaticContext<>::get().queue(), a.rows),
              sum(prm.ctx ? prm.ctx->queue() : vex::StaticContext<>::get().queue()),
              relax(prm.ctx ? prm.ctx->queue() : vex::StaticContext<>::get().queue(), a)
        {
            if (nlevel) {
                u.resize(prm.ctx ? prm.ctx->queue() : vex::StaticContext<>::get().queue(), a.rows);
                f.resize(prm.ctx ? prm.ctx->queue() : vex::StaticContext<>::get().queue(), a.rows);

                if (prm.kcycle && nlevel % prm.kcycle == 0)
                    gmres.reset(new gmres_data<vector>(prm.kcycle_iterations, a.rows));
            }

            a.clear();
            p.clear();
            r.clear();
        }

        // Construct the coarsest hierarchy level from system matrix (a) and
        // its inverse (ai).
        instance(cpu_matrix &a, cpu_matrix &ai, const params &prm, unsigned /*nlevel*/)
            : A(prm.ctx ? prm.ctx->queue() : vex::StaticContext<>::get().queue(),
                        a.rows, a.cols, a.row.data(), a.col.data(), a.val.data()),
              Ainv(prm.ctx ? prm.ctx->queue() : vex::StaticContext<>::get().queue(),
                          ai.rows, ai.cols, ai.row.data(), ai.col.data(), ai.val.data()),
              u(prm.ctx ? prm.ctx->queue() : vex::StaticContext<>::get().queue(), a.rows),
              f(prm.ctx ? prm.ctx->queue() : vex::StaticContext<>::get().queue(), a.rows),
              t(prm.ctx ? prm.ctx->queue() : vex::StaticContext<>::get().queue(), a.rows),
              sum(prm.ctx ? prm.ctx->queue() : vex::StaticContext<>::get().queue())
        {
            a.clear();
            ai.clear();
        }

        // Returns reference to the system matrix
        const matrix& get_matrix() const {
            return A;
        }

        // Compute residual value.
        value_t resid(const vector &rhs, vector &x) const {
            t = rhs - A * x;

            return sqrt(sum(t * t));
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

                    lvl->t = rhs - lvl->A * x;
                    nxt->f = lvl->R * lvl->t;
                    nxt->u = 0;

                    if (nxt->gmres)
                        kcycle(pnxt, end, prm, nxt->f, nxt->u);
                    else
                        cycle(pnxt, end, prm, nxt->f, nxt->u);

                    x += lvl->P * nxt->u;

                    for(unsigned i = 0; i < prm.npost; ++i)
                        lvl->relax.apply(lvl->A, rhs, x, lvl->t, prm.relax);
                }
            } else {
                x = lvl->Ainv * rhs;
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
                x = lvl->Ainv * rhs;
            }
        }

        index_t size() const {
            return A.rows();
        }

        index_t nonzeros() const {
            return A.nonzeros();
        }
    private:
        matrix A;
        matrix P;
        matrix R;
        matrix Ainv;

        mutable vector u;
        mutable vector f;
        mutable vector t;

        vex::Reductor<value_t, vex::SUM> sum;

        typename vexcl_relax_scheme<Relaxation>::type::template instance<value_t, index_t> relax;

        mutable std::unique_ptr< gmres_data<vector> > gmres;

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
