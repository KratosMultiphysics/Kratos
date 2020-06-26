#ifndef AMGCL_RELAXATION_RUNTIME_HPP
#define AMGCL_RELAXATION_RUNTIME_HPP

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
 * \file   amgcl/relaxation/runtime.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Runtime-configurable smoother as standalone preconditioner.
 */

#include <type_traits>

#ifdef AMGCL_NO_BOOST
#  error Runtime interface relies on Boost.PropertyTree!
#endif

#include <boost/property_tree/ptree.hpp>

#include <amgcl/util.hpp>
#include <amgcl/relaxation/gauss_seidel.hpp>
#include <amgcl/relaxation/ilu0.hpp>
#include <amgcl/relaxation/iluk.hpp>
#include <amgcl/relaxation/ilut.hpp>
#include <amgcl/relaxation/damped_jacobi.hpp>
#include <amgcl/relaxation/spai0.hpp>
#include <amgcl/relaxation/spai1.hpp>
#include <amgcl/relaxation/chebyshev.hpp>

namespace amgcl {
namespace runtime {
namespace relaxation {

/// Relaxation schemes.
enum type {
    gauss_seidel,               ///< Gauss-Seidel smoothing
    ilu0,                       ///< Incomplete LU with zero fill-in
    iluk,                       ///< Level-based incomplete LU
    ilut,                       ///< Incomplete LU with thresholding
    damped_jacobi,              ///< Damped Jacobi
    spai0,                      ///< Sparse approximate inverse of 0th order
    spai1,                      ///< Sparse approximate inverse of 1st order
    chebyshev                   ///< Chebyshev relaxation
};

inline std::ostream& operator<<(std::ostream &os, type r)
{
    switch (r) {
        case gauss_seidel:
            return os << "gauss_seidel";
        case ilu0:
            return os << "ilu0";
        case iluk:
            return os << "iluk";
        case ilut:
            return os << "ilut";
        case damped_jacobi:
            return os << "damped_jacobi";
        case spai0:
            return os << "spai0";
        case spai1:
            return os << "spai1";
        case chebyshev:
            return os << "chebyshev";
        default:
            return os << "???";
    }
}

inline std::istream& operator>>(std::istream &in, type &r)
{
    std::string val;
    in >> val;

    if (val == "gauss_seidel")
        r = gauss_seidel;
    else if (val == "ilu0")
        r = ilu0;
    else if (val == "iluk")
        r = iluk;
    else if (val == "ilut")
        r = ilut;
    else if (val == "damped_jacobi")
        r = damped_jacobi;
    else if (val == "spai0")
        r = spai0;
    else if (val == "spai1")
        r = spai1;
    else if (val == "chebyshev")
        r = chebyshev;
    else
        throw std::invalid_argument("Invalid relaxation value. Valid choices are:"
                "gauss_seidel, ilu0, iluk, ilut, damped_jacobi, spai0, spai1, chebyshev.");

    return in;
}

template <class Backend>
struct wrapper {
    typedef boost::property_tree::ptree params;
    typedef typename Backend::params    backend_params;
    type r;
    void *handle;

    template <class Matrix>
    wrapper(const Matrix &A, params prm = params(),
            const backend_params &bprm = backend_params()
            )
      : r(prm.get("type", runtime::relaxation::spai0)), handle(0)
    {
        if (!prm.erase("type")) AMGCL_PARAM_MISSING("type");

        switch(r) {

#define AMGCL_RUNTIME_RELAXATION(type) \
            case type: \
                handle = call_constructor<amgcl::relaxation::type>(A, prm, bprm); \
                break

            AMGCL_RUNTIME_RELAXATION(gauss_seidel);
            AMGCL_RUNTIME_RELAXATION(ilu0);
            AMGCL_RUNTIME_RELAXATION(iluk);
            AMGCL_RUNTIME_RELAXATION(ilut);
            AMGCL_RUNTIME_RELAXATION(damped_jacobi);
            AMGCL_RUNTIME_RELAXATION(spai0);
            AMGCL_RUNTIME_RELAXATION(spai1);
            AMGCL_RUNTIME_RELAXATION(chebyshev);

#undef AMGCL_RUNTIME_RELAXATION

            default:
                throw std::invalid_argument("Unsupported relaxation type");
        }
    }

    ~wrapper() {
        switch(r) {

#define AMGCL_RUNTIME_RELAXATION(type) \
            case type: \
                delete static_cast<amgcl::relaxation::type<Backend>*>(handle); \
                break

            AMGCL_RUNTIME_RELAXATION(gauss_seidel);
            AMGCL_RUNTIME_RELAXATION(ilu0);
            AMGCL_RUNTIME_RELAXATION(iluk);
            AMGCL_RUNTIME_RELAXATION(ilut);
            AMGCL_RUNTIME_RELAXATION(damped_jacobi);
            AMGCL_RUNTIME_RELAXATION(spai0);
            AMGCL_RUNTIME_RELAXATION(spai1);
            AMGCL_RUNTIME_RELAXATION(chebyshev);

#undef AMGCL_RUNTIME_RELAXATION
        }
    }

    template <class Matrix, class VectorRHS, class VectorX, class VectorTMP>
    void apply_pre(
            const Matrix &A, const VectorRHS &rhs, VectorX &x, VectorTMP &tmp
            ) const
    {
        switch(r) {

#define AMGCL_RUNTIME_RELAXATION(type) \
            case type: \
                call_apply_pre<amgcl::relaxation::type>(A, rhs, x, tmp); \
                break

            AMGCL_RUNTIME_RELAXATION(gauss_seidel);
            AMGCL_RUNTIME_RELAXATION(ilu0);
            AMGCL_RUNTIME_RELAXATION(iluk);
            AMGCL_RUNTIME_RELAXATION(ilut);
            AMGCL_RUNTIME_RELAXATION(damped_jacobi);
            AMGCL_RUNTIME_RELAXATION(spai0);
            AMGCL_RUNTIME_RELAXATION(spai1);
            AMGCL_RUNTIME_RELAXATION(chebyshev);

#undef AMGCL_RUNTIME_RELAXATION

            default:
                throw std::invalid_argument("Unsupported relaxation type");
        }
    }

    template <class Matrix, class VectorRHS, class VectorX, class VectorTMP>
    void apply_post(
            const Matrix &A, const VectorRHS &rhs, VectorX &x, VectorTMP &tmp
            ) const
    {
        switch(r) {

#define AMGCL_RUNTIME_RELAXATION(type) \
            case type: \
                call_apply_post<amgcl::relaxation::type>(A, rhs, x, tmp); \
                break

            AMGCL_RUNTIME_RELAXATION(gauss_seidel);
            AMGCL_RUNTIME_RELAXATION(ilu0);
            AMGCL_RUNTIME_RELAXATION(iluk);
            AMGCL_RUNTIME_RELAXATION(ilut);
            AMGCL_RUNTIME_RELAXATION(damped_jacobi);
            AMGCL_RUNTIME_RELAXATION(spai0);
            AMGCL_RUNTIME_RELAXATION(spai1);
            AMGCL_RUNTIME_RELAXATION(chebyshev);

#undef AMGCL_RUNTIME_RELAXATION

            default:
                throw std::invalid_argument("Unsupported relaxation type");
        }
    }

    template <class Matrix, class VectorRHS, class VectorX>
    void apply( const Matrix &A, const VectorRHS &rhs, VectorX &x) const
    {
        switch(r) {

#define AMGCL_RUNTIME_RELAXATION(type) \
            case type: \
                call_apply<amgcl::relaxation::type>(A, rhs, x); \
                break

            AMGCL_RUNTIME_RELAXATION(gauss_seidel);
            AMGCL_RUNTIME_RELAXATION(ilu0);
            AMGCL_RUNTIME_RELAXATION(iluk);
            AMGCL_RUNTIME_RELAXATION(ilut);
            AMGCL_RUNTIME_RELAXATION(damped_jacobi);
            AMGCL_RUNTIME_RELAXATION(spai0);
            AMGCL_RUNTIME_RELAXATION(spai1);
            AMGCL_RUNTIME_RELAXATION(chebyshev);

#undef AMGCL_RUNTIME_RELAXATION

            default:
                throw std::invalid_argument("Unsupported relaxation type");
        }
    }

    size_t bytes() const {
        switch(r) {

#define AMGCL_RUNTIME_RELAXATION(type) \
            case type: \
                return backend::bytes(*static_cast<amgcl::relaxation::type<Backend>*>(handle))

            AMGCL_RUNTIME_RELAXATION(gauss_seidel);
            AMGCL_RUNTIME_RELAXATION(ilu0);
            AMGCL_RUNTIME_RELAXATION(iluk);
            AMGCL_RUNTIME_RELAXATION(ilut);
            AMGCL_RUNTIME_RELAXATION(damped_jacobi);
            AMGCL_RUNTIME_RELAXATION(spai0);
            AMGCL_RUNTIME_RELAXATION(spai1);
            AMGCL_RUNTIME_RELAXATION(chebyshev);

#undef AMGCL_RUNTIME_RELAXATION

            default:
                throw std::invalid_argument("Unsupported relaxation type");
        }
    }

    template <template <class> class Relaxation, class Matrix>
    typename std::enable_if<
        backend::relaxation_is_supported<Backend, Relaxation>::value,
        void*
    >::type
    call_constructor(
            const Matrix &A, const params &prm, const backend_params &bprm)
    {
        return static_cast<void*>(new Relaxation<Backend>(A, prm, bprm));
    }

    template <template <class> class Relaxation, class Matrix>
    typename std::enable_if<
        !backend::relaxation_is_supported<Backend, Relaxation>::value,
        void*
    >::type
    call_constructor(const Matrix&, const params&, const backend_params&)
    {
        throw std::logic_error("The relaxation is not supported by the backend");
    }

    template <template <class> class Relaxation, class Matrix, class VectorRHS, class VectorX, class VectorTMP>
    typename std::enable_if<
        backend::relaxation_is_supported<Backend, Relaxation>::value,
        void
    >::type
    call_apply_pre(
            const Matrix &A, const VectorRHS &rhs, VectorX &x, VectorTMP &tmp) const
    {
        static_cast<Relaxation<Backend>*>(handle)->apply_pre(A, rhs, x, tmp);
    }

    template <template <class> class Relaxation, class Matrix, class VectorRHS, class VectorX, class VectorTMP>
    typename std::enable_if<
        !backend::relaxation_is_supported<Backend, Relaxation>::value,
        void
    >::type
    call_apply_pre(const Matrix&, const VectorRHS&, VectorX&, VectorTMP&) const {
        throw std::logic_error("The relaxation is not supported by the backend");
    }

    template <template <class> class Relaxation, class Matrix, class VectorRHS, class VectorX, class VectorTMP>
    typename std::enable_if<
        backend::relaxation_is_supported<Backend, Relaxation>::value,
        void
    >::type
    call_apply_post(
            const Matrix &A, const VectorRHS &rhs, VectorX &x, VectorTMP &tmp) const
    {
        static_cast<Relaxation<Backend>*>(handle)->apply_post(A, rhs, x, tmp);
    }

    template <template <class> class Relaxation, class Matrix, class VectorRHS, class VectorX, class VectorTMP>
    typename std::enable_if<
        !backend::relaxation_is_supported<Backend, Relaxation>::value,
        void
    >::type
    call_apply_post(const Matrix&, const VectorRHS&, VectorX&, VectorTMP&) const {
        throw std::logic_error("The relaxation is not supported by the backend");
    }

    template <template <class> class Relaxation, class Matrix, class VectorRHS, class VectorX>
    typename std::enable_if<
        backend::relaxation_is_supported<Backend, Relaxation>::value,
        void
    >::type
    call_apply(
            const Matrix &A, const VectorRHS &rhs, VectorX &x) const
    {
        static_cast<Relaxation<Backend>*>(handle)->apply(A, rhs, x);
    }

    template <template <class> class Relaxation, class Matrix, class VectorRHS, class VectorX>
    typename std::enable_if<
        !backend::relaxation_is_supported<Backend, Relaxation>::value,
        void
    >::type
    call_apply(const Matrix&, const VectorRHS&, VectorX&) const {
        throw std::logic_error("The relaxation is not supported by the backend");
    }

};

} // namespace relaxation
} // namespace runtime
} // namespace amgcl

#endif
