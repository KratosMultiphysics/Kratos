#ifndef AMGCL_PRECONDITIONER_CPR_HPP
#define AMGCL_PRECONDITIONER_CPR_HPP

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

#include <amgcl/util.hpp>

namespace amgcl {
namespace preconditioner {

/// CPR (constrained pressure residual) preconditioner.
/** \cite stueben2007algebraic */
template <
    class PressurePrecond,
    class FlowPrecond
    >
class cpr {
    public:
        BOOST_STATIC_ASSERT_MSG(
                (
                 boost::is_same<
                     typename PressurePrecond::backend_type,
                     typename FlowPrecond::backend_type
                     >::value
                ),
                "Backends for pressure and flow preconditioners should coinside!"
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
                        "Error in CPR parameters: pmask is not set"
                        );

                precondition(
                        n > 0,
                        "Error in CPR parameters: "
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
        cpr(
                const Matrix &M,
                const params &prm = params(),
                const backend_params &bprm = backend_params()
           ) : prm(prm)
        {
            typedef typename backend::row_iterator<Matrix>::type row_iterator;

            const size_t n = backend::rows(M);

            precondition(n == prm.pmask.size(), "Incorrect pmask size in CPR");

            // Extract system matrix subblocks:
            //   - App,
            //   - Asp,
            //   - Dss = Dia(Ass),
            //   - Aps * Dss^-1
            size_t ns = boost::count(prm.pmask, false);
            np = n - ns;

            build_matrix App, Aps, Asp;
            boost::shared_ptr<build_matrix> B = boost::make_shared<build_matrix>();
            std::vector<value_type> Dss(ns);

            App.nrows = np;
            App.ncols = np;
            App.ptr.resize(np + 1, 0);

            Aps.nrows = np;
            Aps.ncols = ns;
            Aps.ptr.resize(np + 1, 0);

            Asp.nrows = ns;
            Asp.ncols = np;
            Asp.ptr.resize(ns + 1, 0);

            B->nrows = np;
            B->ncols = n;
            B->ptr.resize(np + 1, 1); B->ptr[0] = 0;

            std::vector<ptrdiff_t> idx(n);
            std::vector<ptrdiff_t> pix(np);

            for(size_t i = 0, ip = 0, is = 0; i < n; ++i) {
                if (prm.pmask[i]) {
                    pix[ip] = i;
                    idx[i] = ip++;
                } else {
                    idx[i] = is++;
                }
            }

#pragma omp parallel for
            for(ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(n); ++i) {
                size_t ii = idx[i], in = ii + 1;

                for(row_iterator a = backend::row_begin(M, i); a; ++a) {
                    ptrdiff_t j = a.col();

                    if (prm.pmask[i]) {
                        if (prm.pmask[j]) {
                            ++App.ptr[in];
                        } else {
                            ++Aps.ptr[in];
                            ++B->ptr[in];
                        }
                    } else {
                        if (j == i) {
                            Dss[ii] = a.value();
                        } else if (prm.pmask[j]) {
                            ++Asp.ptr[in];
                        }
                    }
                }
            }

            boost::partial_sum(App.ptr, App.ptr.begin());
            App.col.resize(App.ptr.back());
            App.val.resize(App.ptr.back());

            boost::partial_sum(Aps.ptr, Aps.ptr.begin());
            Aps.col.resize(Aps.ptr.back());
            Aps.val.resize(Aps.ptr.back());

            boost::partial_sum(Asp.ptr, Asp.ptr.begin());
            Asp.col.resize(Asp.ptr.back());
            Asp.val.resize(Asp.ptr.back());

            boost::partial_sum(B->ptr, B->ptr.begin());
            B->col.resize(B->ptr.back());
            B->val.resize(B->ptr.back());

#pragma omp parallel for
            for(ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(n); ++i) {
                size_t ii = idx[i];

                if (prm.pmask[i]) {
                    size_t p_head = App.ptr[ii];
                    size_t s_head = Aps.ptr[ii];
                    size_t b_head = B->ptr[ii];

                    B->col[b_head] = i;
                    B->val[b_head] = 1;
                    ++b_head;

                    for(row_iterator a = backend::row_begin(M, i); a; ++a) {
                        size_t j  = a.col();
                        size_t jj = idx[j];

                        if (prm.pmask[j]) {
                            App.col[p_head] = jj;
                            App.val[p_head] = a.value();
                            ++p_head;
                        } else {
                            value_type aps = a.value() / Dss[jj];

                            Aps.col[s_head] = jj;
                            Aps.val[s_head] = aps;
                            ++s_head;

                            B->col[b_head] = j;
                            B->val[b_head] = -aps;
                            ++b_head;
                        }
                    }
                } else {
                    size_t p_head = Asp.ptr[ii];

                    for(row_iterator a = backend::row_begin(M, i); a; ++a) {
                        size_t j  = a.col();

                        if (prm.pmask[j]) {
                            Asp.col[p_head] = idx[j];
                            Asp.val[p_head] = a.value();
                            ++p_head;
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
                ptrdiff_t pp_dia = App.ptr[i];
                for(; pp_dia < App.ptr[i+1]; ++pp_dia)
                    if (App.col[pp_dia] == i) break;
                assert(App.col[pp_dia] == i && "No diagonal in App?");

                // Compute i-th diagonal of Aps * Asp:
                value_type apsp = 0;
                for(row_iterator a = backend::row_begin(Aps, i); a; ++a)
                    for(row_iterator b = backend::row_begin(Asp, a.col()); b; ++b)
                        if (b.col() == i) apsp += a.value() * b.value();

                App.val[pp_dia] -= apsp;
            }

            // We are ready to create AMG and ILU preconditioners, and also
            // transfer the Br matrix to the backend.
            Br  = backend_type::copy_matrix(B, bprm);
            bp  = backend_type::create_vector(np, bprm);
            xp  = backend_type::create_vector(np, bprm);
            bu  = backend_type::create_vector(n, bprm);
            xu  = backend_type::create_vector(n, bprm);

            P = boost::make_shared<PressurePrecond>(App, prm.pressure, bprm);
            I = boost::make_shared<FlowPrecond>(M, prm.flow, bprm);

            scatter = boost::make_shared<typename backend_type::scatter>(n, pix, bprm);
            tmp.resize(np);
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
            backend::clear(x);

            // Compute RHS for the reduced pressure problem:
            backend::spmv(1, *Br, rhs, 0, *bp);

            // Do single V-cycle of amg on the computed RHS:
            P->apply(*bp, *xp);

            // Expand pressure vector onto complete solution:
            (*scatter)(*xp, x);

            // Postprocess the complete solution with an ILU sweep:
            backend::residual(rhs, system_matrix(), x, *bu);
            I->apply(*bu, *xu);
            backend::axpby(1, *xu, 1, x);
        }

        const matrix& system_matrix() const {
            return I->system_matrix();
        }
    private:
        typedef typename backend::builtin<value_type>::matrix build_matrix;

        size_t np;
        std::vector<value_type> tmp;
        boost::shared_ptr<matrix> Br;
        boost::shared_ptr<vector> xp, bp, xu, bu;
        boost::shared_ptr<PressurePrecond> P;
        boost::shared_ptr<FlowPrecond>     I;

        boost::shared_ptr<typename backend_type::scatter> scatter;
};

} // namespace preconditioner
} // namespace amgcl

#endif
