#ifndef AMGCL_RELAXATION_RUNTIME_HPP
#define AMGCL_RELAXATION_RUNTIME_HPP

/*
The MIT License

Copyright (c) 2012-2017 Denis Demidov <dennis.demidov@gmail.com>

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

#include <amgcl/relaxation/gauss_seidel.hpp>
#include <amgcl/relaxation/multicolor_gauss_seidel.hpp>
#include <amgcl/relaxation/ilu0.hpp>
#include <amgcl/relaxation/parallel_ilu0.hpp>
#include <amgcl/relaxation/iluk.hpp>
#include <amgcl/relaxation/ilut.hpp>
#include <amgcl/relaxation/damped_jacobi.hpp>
#include <amgcl/relaxation/spai0.hpp>
#include <amgcl/relaxation/spai1.hpp>
#include <amgcl/relaxation/chebyshev.hpp>
#include <amgcl/relaxation/as_preconditioner.hpp>

namespace amgcl {
namespace runtime {
namespace relaxation {

/// Relaxation schemes.
enum type {
    gauss_seidel,               ///< Gauss-Seidel smoothing
    multicolor_gauss_seidel,    ///< Multicolor Gauss-seidel
    ilu0,                       ///< Incomplete LU with zero fill-in
    parallel_ilu0,              ///< Parallel version of ILU(0)
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
        case multicolor_gauss_seidel:
            return os << "multicolor_gauss_seidel";
        case ilu0:
            return os << "ilu0";
        case parallel_ilu0:
            return os << "parallel_ilu0";
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

    if (val == "multicolor_gauss_seidel")
        r = multicolor_gauss_seidel;
    else if (val == "gauss_seidel")
        r = gauss_seidel;
    else if (val == "ilu0")
        r = ilu0;
    else if (val == "parallel_ilu0")
        r = parallel_ilu0;
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
        throw std::invalid_argument("Invalid relaxation value");

    return in;
}

namespace detail {

template <
    class Backend,
    template <class> class Relaxation,
    class Func
    >
inline
typename boost::enable_if<
    typename backend::relaxation_is_supported<Backend, Relaxation>::type,
    void
>::type
process_rap(const Func &func) {
    typedef amgcl::relaxation::as_preconditioner<Backend, Relaxation> RAP;
    func.template process<RAP>();
}

template <
    class Backend,
    template <class> class Relaxation,
    class Func
    >
inline
typename boost::disable_if<
    typename backend::relaxation_is_supported<Backend, Relaxation>::type,
    void
>::type
process_rap(const Func &func) {
    throw std::logic_error("The relaxation scheme is not supported by the backend");
}

template <class Backend, class Func>
void process_rap(runtime::relaxation::type relaxation, const Func &func) {
    switch (relaxation) {
        case runtime::relaxation::gauss_seidel:
            process_rap<Backend, amgcl::relaxation::gauss_seidel>(func);
            break;
        case runtime::relaxation::multicolor_gauss_seidel:
            process_rap<Backend, amgcl::relaxation::multicolor_gauss_seidel>(func);
            break;
        case runtime::relaxation::ilu0:
            process_rap<Backend, amgcl::relaxation::ilu0>(func);
            break;
        case runtime::relaxation::parallel_ilu0:
            process_rap<Backend, amgcl::relaxation::parallel_ilu0>(func);
            break;
        case runtime::relaxation::iluk:
            process_rap<Backend, amgcl::relaxation::iluk>(func);
            break;
        case runtime::relaxation::ilut:
            process_rap<Backend, amgcl::relaxation::ilut>(func);
            break;
        case runtime::relaxation::damped_jacobi:
            process_rap<Backend, amgcl::relaxation::damped_jacobi>(func);
            break;
        case runtime::relaxation::spai0:
            process_rap<Backend, amgcl::relaxation::spai0>(func);
            break;
        case runtime::relaxation::spai1:
            process_rap<Backend, amgcl::relaxation::spai1>(func);
            break;
        case runtime::relaxation::chebyshev:
            process_rap<Backend, amgcl::relaxation::chebyshev>(func);
            break;
    }
}

template <class Backend, class Matrix>
struct rap_create {
    typedef boost::property_tree::ptree params;
    typedef typename Backend::params backend_params;

    void * &handle;
    const Matrix &A;
    const params &p;
    const backend_params &bp;

    rap_create(void* &handle, const Matrix &A, const params &p, const backend_params &bp)
        : handle(handle), A(A), p(p), bp(bp) {}

    template <class RAP>
    void process() const {
        handle = static_cast<void*>(new RAP(A, p, bp));
    }
};

struct rap_destroy {
    void * handle;

    rap_destroy(void *handle) : handle(handle) {}

    template <class RAP>
    void process() const {
        delete static_cast<RAP*>(handle);
    }
};

template <class Vec1, class Vec2>
struct rap_apply {
    void * handle;
    Vec1 const &rhs;
    Vec2 &x;

    rap_apply(void *handle, const Vec1 &rhs, Vec2 &x)
        : handle(handle), rhs(rhs), x(x) {}

    template <class RAP>
    void process() const {
        static_cast<RAP*>(handle)->apply(rhs, x);
    }
};

template <class Matrix>
struct rap_matrix {
    void * handle;
    const Matrix * &A;

    rap_matrix(void *handle, const Matrix * &A) : handle(handle), A(A) {}

    template <class RAP>
    void process() const {
        A = &(static_cast<RAP*>(handle)->system_matrix());
    }
};

}

/// Use one of AMGCL smoothers as standalone preconditioner.
/**
 * The exact smoother is selected at runtime through prm.type parameter
 */
template <class Backend>
class as_preconditioner {
    public:
        typedef Backend backend_type;

        typedef typename Backend::matrix  matrix;

        typedef boost::property_tree::ptree params;
        typedef typename Backend::params  backend_params;


        template <class Matrix>
        as_preconditioner(
                const Matrix &A,
                params prm = params(),
                const backend_params &bprm = backend_params()
                )
            : relaxation(prm.get("type", runtime::relaxation::spai0)),
              handle(0)
        {
            if (!prm.erase("type")) AMGCL_PARAM_MISSING("type");

            runtime::relaxation::detail::process_rap<Backend>(
                    relaxation,
                    runtime::relaxation::detail::rap_create<Backend, Matrix>(
                        handle, A, prm, bprm
                        )
                    );
        }

        ~as_preconditioner() {
            runtime::relaxation::detail::process_rap<Backend>(
                    relaxation,
                    runtime::relaxation::detail::rap_destroy(handle)
                    );
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
            runtime::relaxation::detail::process_rap<Backend>(
                    relaxation,
                    runtime::relaxation::detail::rap_apply<Vec1, Vec2>(
                        handle, rhs, x
                        )
                    );
        }

        const matrix& system_matrix() const {
            const matrix *A = 0;

            runtime::relaxation::detail::process_rap<Backend>(
                    relaxation,
                    runtime::relaxation::detail::rap_matrix<matrix>(handle, A)
                    );

            return *A;
        }
    private:
        runtime::relaxation::type relaxation;
        void * handle;

        friend std::ostream& operator<<(std::ostream &os, const as_preconditioner &p)
        {
            os << "Using " << p.relaxation << " as preconditioner" << std::endl;
            os << "  unknowns: " << backend::rows(p.system_matrix()) << std::endl;
            os << "  nonzeros: " << backend::nonzeros(p.system_matrix()) << std::endl;
            return os;
        }
};

}
}
}

#endif
