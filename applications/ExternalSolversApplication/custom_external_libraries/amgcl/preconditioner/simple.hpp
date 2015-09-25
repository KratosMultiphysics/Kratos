#ifndef AMGCL_PRECONDITIONER_SIMPLE_HPP
#define AMGCL_PRECONDITIONER_SIMPLE_HPP

/*
The MIT License

Copyright (c) 2012-2015 Denis Demidov <dennis.demidov@gmail.com>
Copyright (c) 2015, Riccardo Rossi, CIMNE (International Center for Numerical Methods in Engineering)

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
 * \file   amgcl/preconditioner/cpr.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  CPR (constrained pressure residual) preconditioner implementation.
 */

#include <vector>
#include <cassert>

#include <boost/shared_ptr.hpp>
#include <amgcl/backend/builtin.hpp>

#include <boost/shared_ptr.hpp>
#include <boost/function.hpp>
#include <boost/range/irange.hpp>
#include <boost/range/algorithm.hpp>
#include <amgcl/amgcl.hpp>
#include <amgcl/relaxation/ilu0.hpp>
#include <amgcl/util.hpp>

namespace amgcl {
namespace preconditioner {

template <
    class Backend,
    class Coarsening,
    template <class> class Relax
    >
class simple {
    public:
        typedef amg<Backend, Coarsening, Relax> AMG;

        typedef typename Backend::value_type value_type;
        typedef typename Backend::matrix     matrix;
        typedef typename Backend::vector     vector;

        typedef typename AMG::params amg_params;

        template <class Matrix>
        simple(
                const Matrix &M,
                boost::function<bool(size_t)> pmask,
                const amg_params &prm = amg_params()
           ) : aprm(prm), pmask(pmask), n(backend::rows(M))
        {
            typedef typename backend::row_iterator<Matrix>::type row_iterator;

            // Extract system matrix subblocks:
            //   - App,
            //   - Asp,
            //   - Dss = Dia(Ass),
            //   - Aps * Dss^-1
            const size_t np = boost::count_if(boost::irange<size_t>(0, n), pmask);
            const size_t ns = n - np;

            boost::shared_ptr<build_matrix> App = boost::make_shared<build_matrix>();
            boost::shared_ptr<build_matrix> Aps = boost::make_shared<build_matrix>();
            boost::shared_ptr<build_matrix> Asp = boost::make_shared<build_matrix>();
            boost::shared_ptr<build_matrix> Ass = boost::make_shared<build_matrix>();

            std::vector<value_type> Dss(ns);

            App->nrows = np;
            App->ncols = np;
            App->ptr.resize(np + 1, 0);

            Aps->nrows = np;
            Aps->ncols = ns;
            Aps->ptr.resize(np + 1, 0);

            Asp->nrows = ns;
            Asp->ncols = np;
            Asp->ptr.resize(ns + 1, 0);

            Ass->nrows = ns;
            Ass->ncols = ns;
            Ass->ptr.resize(ns + 1, 0);

            idx.resize(n);

            for(size_t i = 0, ip = 0, is = 0; i < n; ++i) {
                if (pmask(i)) {
                    idx[i] = ip++;
                } else {
                    idx[i] = is++;
                }

                for(row_iterator a = backend::row_begin(M, i); a; ++a) {
                    size_t j = a.col();

                    if (pmask(i)) {
                        if (pmask(j)) {
                            ++App->ptr[ip];
                        } else {
                            ++Aps->ptr[ip];
                        }
                    } else {
                        if (pmask(j)) {
                            ++Asp->ptr[is];
                        } else {
                            ++Ass->ptr[is];
                            if (j == i) {
                                Dss[idx[i]] = a.value();
                            }
                        }
                    }
                }
            }

            boost::partial_sum(App->ptr, App->ptr.begin());
            App->col.resize(App->ptr.back());
            App->val.resize(App->ptr.back());

            boost::partial_sum(Aps->ptr, Aps->ptr.begin());
            Aps->col.resize(Aps->ptr.back());
            Aps->val.resize(Aps->ptr.back());

            boost::partial_sum(Asp->ptr, Asp->ptr.begin());
            Asp->col.resize(Asp->ptr.back());
            Asp->val.resize(Asp->ptr.back());

            boost::partial_sum(Ass->ptr, Ass->ptr.begin());
            Ass->col.resize(Ass->ptr.back());
            Ass->val.resize(Ass->ptr.back());

            for(size_t i = 0; i < n; ++i) {
                size_t ii = idx[i];

                if (pmask(i)) {
                    size_t p_head = App->ptr[ii];
                    size_t s_head = Aps->ptr[ii];

                    for(row_iterator a = backend::row_begin(M, i); a; ++a) {
                        size_t j  = a.col();
                        size_t jj = idx[j];

                        if (pmask(j)) {
                            App->col[p_head] = jj;
                            App->val[p_head] = a.value();
                            ++p_head;
                        } else {
                            Aps->col[s_head] = jj;
                            Aps->val[s_head] = a.value();
                            ++s_head;
                        }
                    }
                } else {
                    size_t p_head = Asp->ptr[ii];
                    size_t s_head = Ass->ptr[ii];

                    for(row_iterator a = backend::row_begin(M, i); a; ++a) {
                        size_t j  = a.col();
                        size_t jj = idx[j];

                        if (pmask(j)) {
                            Asp->col[p_head] = jj;
                            Asp->val[p_head] = a.value() / Dss[ii];
                            ++p_head;
                        } else {
                            Ass->col[s_head] = jj;
                            Ass->val[s_head] = a.value();
                            ++s_head;
                        }
                    }
                }
            }

            // Now we need to compute the matrix Ap that will be used for
            // construction of amg preconditioner. In the original CPR method
            // the matrix is written as
            //   Ap = App - Aps * Dss^-1 * Asp,
            // but we use the simplified form of
            //   Ap = App - Dia(Aps * Dss*-1 * Asp)
            // This allows us to modify App in place.
            for(size_t i = 0; i < np; ++i) {
                typedef typename backend::row_iterator<build_matrix>::type row_iterator;

                // Find the diagonal of App (we need to update it):
                size_t pp_dia = App->ptr[i];
                for(; pp_dia < App->ptr[i+1]; ++pp_dia)
                    if (App->col[pp_dia] == i) break;
                assert(App->col[pp_dia] == i && "No diagonal in App?");

                // Compute i-th diagonal of Aps * Asp:
                value_type apsp = 0;
                for(row_iterator a = backend::row_begin(*Aps, i); a; ++a)
                    for(row_iterator b = backend::row_begin(*Asp, a.col()); b; ++b)
                        if (b.col() == i) apsp += a.value() * b.value();

                App->val[pp_dia] -= apsp;
            }

            // We are ready to create AMG and ILU preconditioners.
            K   = Backend::copy_matrix(Ass, aprm.backend);
            D   = Backend::copy_matrix(Aps, aprm.backend);
            G   = Backend::copy_matrix(Asp, aprm.backend);
            bp  = Backend::create_vector(np, aprm.backend);
            xp  = Backend::create_vector(np, aprm.backend);
            bs  = Backend::create_vector(ns, aprm.backend);
            xs  = Backend::create_vector(ns, aprm.backend);
            tmp = Backend::create_vector(ns, aprm.backend);

            A = Backend::copy_matrix(boost::make_shared<build_matrix>(M), aprm.backend);
            P = boost::make_shared<AMG>(*App, aprm);

            iprm.damping = 1;

            I = boost::make_shared<ILU>(*Ass, iprm, aprm.backend);
        }

        template <class Vec1, class Vec2>
        void apply(
                const Vec1 &rhs,
#ifdef BOOST_NO_CXX11_RVALUE_REFERENCES
                Vec2       &x
#else
                Vec2       &&x
#endif
                ) const
        {
            typedef typename backend::row_iterator<matrix>::type row_iterator;

            // Split RHS into p and s parts:
            for(size_t i = 0; i < n; ++i) {
                if (pmask(i))
                    (*bp)[idx[i]] = rhs[i];
                else
                    (*bs)[idx[i]] = rhs[i];
            }

            // Solve for s part:
            backend::clear(*xs);
            I->apply_pre(*K, *bs, *xs, *tmp, iprm);

            // Compute RHS for the reduced pressure problem:
            backend::spmv(-1, *D, *xs, 1, *bp);

            // Do single V-cycle of amg on the computed RHS:
            P->apply(*bp, *xp);

            // Update xs:
            backend::spmv(-1, *G, *xp, 1, *xs);

            // Expand partial vectors onto complete solution:
            // TODO: this only works for host-addressable backends now.
            for(size_t i = 0; i < n; ++i) {
                if (pmask(i))
                    x[i] = (*xp)[idx[i]];
                else
                    x[i] = (*xs)[idx[i]];
            }
        }

        const matrix& system_matrix() const {
            return *A;
        }
    private:
        typedef typename backend::builtin<value_type>::matrix build_matrix;
        typedef relaxation::ilu0<Backend> ILU;
        typedef typename ILU::params ilu_params;

        amg_params aprm;
        ilu_params iprm;

        boost::function<bool(size_t)> pmask;

        size_t n;
        std::vector<size_t> idx;
        boost::shared_ptr<matrix> A, K, D, G;
        boost::shared_ptr<vector> xp, bp, xs, bs, tmp;
        boost::shared_ptr<AMG>    P;
        boost::shared_ptr<ILU>    I;
};

} // namespace preconditioner
} // namespace amgcl

#endif
