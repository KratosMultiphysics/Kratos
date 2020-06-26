#ifndef AMGCL_SOLVER_RUNTIME_HPP
#define AMGCL_SOLVER_RUNTIME_HPP

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
 * \file   amgcl/solver/runtime.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Runtime-configurable wrappers around amgcl iterative solvers.
 */

#include <iostream>
#include <stdexcept>
#include <type_traits>

#ifdef AMGCL_NO_BOOST
#  error Runtime interface relies on Boost.PropertyTree!
#endif

#include <boost/property_tree/ptree.hpp>

#include <amgcl/util.hpp>
#include <amgcl/solver/cg.hpp>
#include <amgcl/solver/bicgstab.hpp>
#include <amgcl/solver/bicgstabl.hpp>
#include <amgcl/solver/gmres.hpp>
#include <amgcl/solver/lgmres.hpp>
#include <amgcl/solver/fgmres.hpp>
#include <amgcl/solver/idrs.hpp>
#include <amgcl/solver/detail/default_inner_product.hpp>

namespace amgcl {
namespace runtime {
namespace solver {

enum type {
    cg,         ///< Conjugate gradients method
    bicgstab,   ///< BiConjugate Gradient Stabilized
    bicgstabl,  ///< BiCGStab(ell)
    gmres,      ///< GMRES
    lgmres,     ///< LGMRES
    fgmres,     ///< FGMRES
    idrs        ///< IDR(s)
};

inline std::ostream& operator<<(std::ostream &os, type s)
{
    switch (s) {
        case cg:
            return os << "cg";
        case bicgstab:
            return os << "bicgstab";
        case bicgstabl:
            return os << "bicgstabl";
        case gmres:
            return os << "gmres";
        case lgmres:
            return os << "lgmres";
        case fgmres:
            return os << "fgmres";
        case idrs:
            return os << "idrs";
        default:
            return os << "???";
    }
}

inline std::istream& operator>>(std::istream &in, type &s)
{
    std::string val;
    in >> val;

    if (val == "cg")
        s = cg;
    else if (val == "bicgstab")
        s = bicgstab;
    else if (val == "bicgstabl")
        s = bicgstabl;
    else if (val == "gmres")
        s = gmres;
    else if (val == "lgmres")
        s = lgmres;
    else if (val == "fgmres")
        s = fgmres;
    else if (val == "idrs")
        s = idrs;
    else
        throw std::invalid_argument("Invalid solver value. Valid choices are: "
                "cg, bicgstab, bicgstabl, gmres, lgmres, fgmres, idrs.");

    return in;
}

template <
    class Backend,
    class InnerProduct = amgcl::solver::detail::default_inner_product
    >
struct wrapper {
    typedef boost::property_tree::ptree                params;
    typedef typename Backend::params                   backend_params;
    typedef typename Backend::value_type               value_type;
    typedef typename math::scalar_of<value_type>::type scalar_type;
    typedef Backend                                    backend_type;

    type s;
    void *handle;

    wrapper(size_t n, params prm = params(),
            const backend_params &bprm = backend_params(),
            const InnerProduct &inner_product = InnerProduct()
            )
        : s(prm.get("type", runtime::solver::bicgstab)), handle(0)
    {
        if (!prm.erase("type")) AMGCL_PARAM_MISSING("type");

        switch(s) {

#define AMGCL_RUNTIME_SOLVER(type) \
            case type: \
                handle = static_cast<void*>(new amgcl::solver::type<Backend, InnerProduct>(n, prm, bprm, inner_product)); \
                break

            AMGCL_RUNTIME_SOLVER(cg);
            AMGCL_RUNTIME_SOLVER(bicgstab);
            AMGCL_RUNTIME_SOLVER(bicgstabl);
            AMGCL_RUNTIME_SOLVER(gmres);
            AMGCL_RUNTIME_SOLVER(lgmres);
            AMGCL_RUNTIME_SOLVER(fgmres);
            AMGCL_RUNTIME_SOLVER(idrs);

#undef AMGCL_RUNTIME_SOLVER

            default:
                throw std::invalid_argument("Unsupported solver type");
        }
    }

    ~wrapper() {
        switch(s) {

#define AMGCL_RUNTIME_SOLVER(type) \
            case type: \
                delete static_cast<amgcl::solver::type<Backend, InnerProduct>*>(handle); \
                break

            AMGCL_RUNTIME_SOLVER(cg);
            AMGCL_RUNTIME_SOLVER(bicgstab);
            AMGCL_RUNTIME_SOLVER(bicgstabl);
            AMGCL_RUNTIME_SOLVER(gmres);
            AMGCL_RUNTIME_SOLVER(lgmres);
            AMGCL_RUNTIME_SOLVER(fgmres);
            AMGCL_RUNTIME_SOLVER(idrs);

#undef AMGCL_RUNTIME_SOLVER
        }
    }

    template <class Matrix, class Precond, class Vec1, class Vec2>
    std::tuple<size_t, scalar_type> operator()(
            const Matrix &A, const Precond &P, const Vec1 &rhs, Vec2 &&x) const
    {
        switch(s) {

#define AMGCL_RUNTIME_SOLVER(type) \
            case type: \
                return static_cast<amgcl::solver::type<Backend, InnerProduct>*>(handle)->operator()(A, P, rhs, x)

            AMGCL_RUNTIME_SOLVER(cg);
            AMGCL_RUNTIME_SOLVER(bicgstab);
            AMGCL_RUNTIME_SOLVER(bicgstabl);
            AMGCL_RUNTIME_SOLVER(gmres);
            AMGCL_RUNTIME_SOLVER(lgmres);
            AMGCL_RUNTIME_SOLVER(fgmres);
            AMGCL_RUNTIME_SOLVER(idrs);

#undef AMGCL_RUNTIME_SOLVER

            default:
                throw std::invalid_argument("Unsupported solver type");
        }
    }

    template <class Precond, class Vec1, class Vec2>
    std::tuple<size_t, scalar_type> operator()(
            const Precond &P, const Vec1 &rhs, Vec2 &&x) const
    {
        return (*this)(P.system_matrix(), P, rhs, x);
    }

    friend std::ostream& operator<<(std::ostream &os, const wrapper &w) {
        switch(w.s) {

#define AMGCL_RUNTIME_SOLVER(type) \
            case type: \
                return os << *static_cast<amgcl::solver::type<Backend, InnerProduct>*>(w.handle)

            AMGCL_RUNTIME_SOLVER(cg);
            AMGCL_RUNTIME_SOLVER(bicgstab);
            AMGCL_RUNTIME_SOLVER(bicgstabl);
            AMGCL_RUNTIME_SOLVER(gmres);
            AMGCL_RUNTIME_SOLVER(lgmres);
            AMGCL_RUNTIME_SOLVER(fgmres);
            AMGCL_RUNTIME_SOLVER(idrs);

#undef AMGCL_RUNTIME_SOLVER

            default:
                throw std::invalid_argument("Unsupported solver type");
        }
    }

    size_t bytes() const {
        switch(s) {

#define AMGCL_RUNTIME_SOLVER(type) \
            case type: \
                return backend::bytes(*static_cast<amgcl::solver::type<Backend, InnerProduct>*>(handle))

            AMGCL_RUNTIME_SOLVER(cg);
            AMGCL_RUNTIME_SOLVER(bicgstab);
            AMGCL_RUNTIME_SOLVER(bicgstabl);
            AMGCL_RUNTIME_SOLVER(gmres);
            AMGCL_RUNTIME_SOLVER(lgmres);
            AMGCL_RUNTIME_SOLVER(fgmres);
            AMGCL_RUNTIME_SOLVER(idrs);

#undef AMGCL_RUNTIME_SOLVER

            default:
                throw std::invalid_argument("Unsupported solver type");
        }
    }
};

} // namespace solver
} // namespace runtime
} // namespace amgcl
#endif
