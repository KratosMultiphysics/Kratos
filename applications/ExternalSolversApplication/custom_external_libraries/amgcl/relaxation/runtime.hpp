#ifndef AMGCL_RELAXATION_RUNTIME_HPP
#define AMGCL_RELAXATION_RUNTIME_HPP

/*
The MIT License

Copyright (c) 2012-2016 Denis Demidov <dennis.demidov@gmail.com>

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

#include <amgcl/runtime.hpp>
#include <amgcl/relaxation/as_preconditioner.hpp>

namespace amgcl {
namespace runtime {
namespace relaxation {

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
                const params &prm = params(),
                const backend_params &bprm = backend_params()
                )
            : relaxation(prm.get("type", runtime::relaxation::spai0)),
              handle(0)
        {
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
};

}
}
}

#endif
