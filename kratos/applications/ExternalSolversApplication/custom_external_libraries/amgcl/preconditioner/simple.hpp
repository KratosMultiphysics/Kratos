#ifndef AMGCL_PRECONDITIONER_SIMPLE_HPP
#define AMGCL_PRECONDITIONER_SIMPLE_HPP

/*
The MIT License

Copyright (c) 2012-2016 Denis Demidov <dennis.demidov@gmail.com>
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
 * \file   amgcl/preconditioner/simple.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  SIMPLE preconditioner (Semi-Implicit Method for Pressure-Linked Equations).
 */

#include <vector>
#include <cassert>

#include <boost/shared_ptr.hpp>
#include <amgcl/backend/builtin.hpp>

#include <boost/shared_ptr.hpp>
#include <boost/function.hpp>
#include <boost/range/irange.hpp>
#include <boost/range/algorithm.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/util.hpp>

namespace amgcl {
namespace preconditioner {


/// SIMPLE preconditioner (Semi-Implicit Method for Pressure-Linked Equations).
/** \cite caretto1973two */
template <
    class PressurePrecond,
    class FlowPrecond
    >
class simple {
    public:
        BOOST_STATIC_ASSERT_MSG(
                (
                 boost::is_same<
                     typename PressurePrecond::backend_type,
                     typename FlowPrecond::backend_type
                     >::value
                ),
                "Backends for pressure and global preconditioners should coinside!"
                );

        typedef typename PressurePrecond::backend_type backend_type;

        typedef typename backend_type::value_type value_type;
        typedef typename backend_type::matrix     matrix;
        typedef typename backend_type::vector     vector;
        typedef typename backend_type::params     backend_params;

        struct params {
            typedef typename PressurePrecond::params pressure_params;
            typedef typename FlowPrecond::params     flow_params;

            pressure_params pressure;
            flow_params     flow;

            std::vector<char> pmask;

            params() {}

            params(const boost::property_tree::ptree &p)
                : AMGCL_PARAMS_IMPORT_CHILD(p, pressure),
                  AMGCL_PARAMS_IMPORT_CHILD(p, flow)
            {
                void *pm = 0;
                size_t n = 0;

                pm = p.get("pmask",     pm);
                n  = p.get("pmask_size", n);

                precondition(
                        pm,
                        "Error in SIMPLE parameters: pmask is not set"
                        );

                precondition(
                        n > 0,
                        "Error in SIMPLE parameters: "
                        "pmask is set, but pmask_size is not"
                        );

                pmask.assign(static_cast<char*>(pm), static_cast<char*>(pm) + n);
            }

            void get(
                    boost::property_tree::ptree &p,
                    const std::string &path = ""
                    ) const
            {
                AMGCL_PARAMS_EXPORT_CHILD(p, path, pressure);
                AMGCL_PARAMS_EXPORT_CHILD(p, path, flow);
            }
        } prm;

        template <class Matrix>
        simple(
                const Matrix &M,
                const params &prm = params(),
                const backend_params &bprm = backend_params()
           ) : prm(prm), n(backend::rows(M))
        {
            typedef typename backend::row_iterator<Matrix>::type row_iterator;

            // Extract system matrix subblocks:
            //   - App,
            //   - Asp,
            //   - Dss = Dia(Ass),
            //   - Aps * Dss^-1
            const size_t ns = boost::count(prm.pmask, false);
            const size_t np = n - ns;

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

            std::vector<size_t> idx(n);
            std::vector<ptrdiff_t> pix(np);
            std::vector<ptrdiff_t> six(ns);

            for(size_t i = 0, ip = 0, is = 0; i < n; ++i) {
                if (prm.pmask[i]) {
                    pix[ip] = i;
                    idx[i] = ip;
                    ++ip;
                } else {
                    six[is] = i;
                    idx[i] = is;
                    ++is;
                }
            }

#pragma omp parallel for
            for(ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(n); ++i) {
                size_t ii = idx[i], in = ii + 1;

                for(row_iterator a = backend::row_begin(M, i); a; ++a) {
                    ptrdiff_t j = a.col();

                    if (prm.pmask[i]) {
                        if (prm.pmask[j]) {
                            ++App->ptr[in];
                        } else {
                            ++Aps->ptr[in];
                        }
                    } else {
                        if (prm.pmask[j]) {
                            ++Asp->ptr[in];
                        } else {
                            ++Ass->ptr[in];
                            if (j == i) {
                                Dss[ii] = a.value();
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

#pragma omp parallel for
            for(ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(n); ++i) {
                size_t ii = idx[i];

                if (prm.pmask[i]) {
                    size_t p_head = App->ptr[ii];
                    size_t s_head = Aps->ptr[ii];

                    for(row_iterator a = backend::row_begin(M, i); a; ++a) {
                        size_t j  = a.col();
                        size_t jj = idx[j];

                        if (prm.pmask[j]) {
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

                        if (prm.pmask[j]) {
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
#pragma omp parallel for
            for(ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(np); ++i) {
                typedef typename backend::row_iterator<build_matrix>::type row_iterator;

                // Find the diagonal of App (we need to update it):
                ptrdiff_t pp_dia = App->ptr[i];
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
            K   = backend_type::copy_matrix(Ass,  bprm);
            D   = backend_type::copy_matrix(Aps,  bprm);
            G   = backend_type::copy_matrix(Asp,  bprm);
            bp  = backend_type::create_vector(np, bprm);
            xp  = backend_type::create_vector(np, bprm);
            bs  = backend_type::create_vector(ns, bprm);
            xs  = backend_type::create_vector(ns, bprm);

            A = backend_type::copy_matrix(boost::make_shared<build_matrix>(M), bprm);
            P = boost::make_shared<PressurePrecond>(*App, prm.pressure, bprm);
            I = boost::make_shared<FlowPrecond>(*Ass, prm.flow, bprm);

            p_gather  = boost::make_shared<typename backend_type::gather >(n, pix, bprm);
            p_scatter = boost::make_shared<typename backend_type::scatter>(n, pix, bprm);
            s_gather  = boost::make_shared<typename backend_type::gather >(n, six, bprm);
            s_scatter = boost::make_shared<typename backend_type::scatter>(n, six, bprm);
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
            // Split RHS into p and s parts:
            (*p_gather)(rhs, *bp);
            (*s_gather)(rhs, *bs);

            // Solve for s part:
            I->apply(*bs, *xs);

            // Compute RHS for the reduced pressure problem:
            backend::spmv(-1, *D, *xs, 1, *bp);

            // Do single V-cycle of amg on the computed RHS:
            P->apply(*bp, *xp);

            // Update xs:
            backend::spmv(-1, *G, *xp, 1, *xs);

            // Expand partial vectors onto complete solution:
            (*p_scatter)(*xp, x);
            (*s_scatter)(*xs, x);
        }

        const matrix& system_matrix() const {
            return *A;
        }
    private:
        typedef typename backend::builtin<value_type>::matrix build_matrix;

        size_t n;

        boost::shared_ptr<typename backend_type::gather>  p_gather;
        boost::shared_ptr<typename backend_type::gather>  s_gather;
        boost::shared_ptr<typename backend_type::scatter> p_scatter;
        boost::shared_ptr<typename backend_type::scatter> s_scatter;

        boost::shared_ptr<matrix> A, K, D, G;
        boost::shared_ptr<vector> xp, bp, xs, bs;
        boost::shared_ptr<PressurePrecond> P;
        boost::shared_ptr<FlowPrecond>     I;
};

} // namespace preconditioner
} // namespace amgcl

#endif
