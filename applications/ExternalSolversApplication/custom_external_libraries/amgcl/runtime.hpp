#ifndef AMGCL_RUNTIME_HPP
#define AMGCL_RUNTIME_HPP

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
 * \file   amgcl/runtime.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Runtime-configurable wrappers around amgcl classes.
 */

#include <iostream>
#include <stdexcept>

#include <boost/type_traits.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/noncopyable.hpp>

#include <amgcl/amg.hpp>

#include <amgcl/backend/interface.hpp>
#include <amgcl/coarsening/ruge_stuben.hpp>
#include <amgcl/coarsening/aggregation.hpp>
#include <amgcl/coarsening/smoothed_aggregation.hpp>
#include <amgcl/coarsening/smoothed_aggr_emin.hpp>

#include <amgcl/relaxation/gauss_seidel.hpp>
#include <amgcl/relaxation/multicolor_gauss_seidel.hpp>
#include <amgcl/relaxation/ilu0.hpp>
#include <amgcl/relaxation/parallel_ilu0.hpp>
#include <amgcl/relaxation/ilut.hpp>
#include <amgcl/relaxation/damped_jacobi.hpp>
#include <amgcl/relaxation/spai0.hpp>
#include <amgcl/relaxation/spai1.hpp>
#include <amgcl/relaxation/chebyshev.hpp>

#include <amgcl/solver/cg.hpp>
#include <amgcl/solver/bicgstab.hpp>
#include <amgcl/solver/bicgstabl.hpp>
#include <amgcl/solver/gmres.hpp>
#include <amgcl/solver/detail/default_inner_product.hpp>


namespace amgcl {

/// Runtime-configurable interface to AMGCL.
namespace runtime {

/// Coarsening kinds.
namespace coarsening {
enum type {
    ruge_stuben,            ///< Ruge-Stueben coarsening
    aggregation,            ///< Aggregation
    smoothed_aggregation,   ///< Smoothed aggregation
    smoothed_aggr_emin      ///< Smoothed aggregation with energy minimization
};

inline std::ostream& operator<<(std::ostream &os, type c) {
    switch (c) {
        case ruge_stuben:
            return os << "ruge_stuben";
        case aggregation:
            return os << "aggregation";
        case smoothed_aggregation:
            return os << "smoothed_aggregation";
        case smoothed_aggr_emin:
            return os << "smoothed_aggr_emin";
        default:
            return os << "???";
    }
}

inline std::istream& operator>>(std::istream &in, type &c)
{
    std::string val;
    in >> val;

    if (val == "ruge_stuben")
        c = ruge_stuben;
    else if (val == "aggregation")
        c = aggregation;
    else if (val == "smoothed_aggregation")
        c = smoothed_aggregation;
    else if (val == "smoothed_aggr_emin")
        c = smoothed_aggr_emin;
    else
        throw std::invalid_argument("Invalid coarsening value");

    return in;
}

} // namespace coarsening

/// Relaxation schemes.
namespace relaxation {
enum type {
    gauss_seidel,               ///< Gauss-Seidel smoothing
    multicolor_gauss_seidel,    ///< Multicolor Gauss-seidel
    ilu0,                       ///< Incoplete LU with zero fill-in
    parallel_ilu0,              ///< Parallel version of ILU(0)
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

} // namespace relaxation

namespace detail {

//---------------------------------------------------------------------------
template <
    class Backend,
    class Coarsening,
    template <class> class Relaxation,
    class Func
    >
inline
typename boost::enable_if<
    typename backend::relaxation_is_supported<Backend, Relaxation>::type,
    void
>::type
process_amg(const Func &func) {
    typedef amgcl::amg<Backend, Coarsening, Relaxation> AMG;
    func.template process<AMG>();
}

template <
    class Backend,
    class Coarsening,
    template <class> class Relaxation,
    class Func
    >
inline
typename boost::disable_if<
    typename backend::relaxation_is_supported<Backend, Relaxation>::type,
    void
>::type
process_amg(const Func&) {
    throw std::logic_error("The relaxation scheme is not supported by the backend");
}

//---------------------------------------------------------------------------
template <
    class Backend,
    class Coarsening,
    class Func
    >
inline
typename boost::disable_if<
    typename backend::coarsening_is_supported<Backend, Coarsening>::type,
    void
    >::type
process_amg(
        runtime::relaxation::type relaxation,
        const Func &func
        )
{
    throw std::logic_error("The coarsening is not supported by the backend");
}

template <
    class Backend,
    class Coarsening,
    class Func
    >
inline
typename boost::enable_if<
    typename backend::coarsening_is_supported<Backend, Coarsening>::type,
    void
    >::type
process_amg(
        runtime::relaxation::type relaxation,
        const Func &func
        )
{
    switch (relaxation) {
        case runtime::relaxation::gauss_seidel:
            process_amg<
                Backend,
                Coarsening,
                amgcl::relaxation::gauss_seidel
                >(func);
            break;
        case runtime::relaxation::multicolor_gauss_seidel:
            process_amg<
                Backend,
                Coarsening,
                amgcl::relaxation::multicolor_gauss_seidel
                >(func);
            break;
        case runtime::relaxation::ilu0:
            process_amg<
                Backend,
                Coarsening,
                amgcl::relaxation::ilu0
                >(func);
            break;
        case runtime::relaxation::parallel_ilu0:
            process_amg<
                Backend,
                Coarsening,
                amgcl::relaxation::parallel_ilu0
                >(func);
            break;
        case runtime::relaxation::ilut:
            process_amg<
                Backend,
                Coarsening,
                amgcl::relaxation::ilut
                >(func);
            break;
        case runtime::relaxation::damped_jacobi:
            process_amg<
                Backend,
                Coarsening,
                amgcl::relaxation::damped_jacobi
                >(func);
            break;
        case runtime::relaxation::spai0:
            process_amg<
                Backend,
                Coarsening,
                amgcl::relaxation::spai0
                >(func);
            break;
        case runtime::relaxation::spai1:
            process_amg<
                Backend,
                Coarsening,
                amgcl::relaxation::spai1
                >(func);
            break;
        case runtime::relaxation::chebyshev:
            process_amg<
                Backend,
                Coarsening,
                amgcl::relaxation::chebyshev
                >(func);
            break;
    }
}

//---------------------------------------------------------------------------
template <
    class Backend,
    class Func
    >
inline void process_amg(
        runtime::coarsening::type coarsening,
        runtime::relaxation::type relaxation,
        const Func &func
        )
{
    switch (coarsening) {
        case runtime::coarsening::ruge_stuben:
            process_amg<
                Backend,
                amgcl::coarsening::ruge_stuben
                >(relaxation, func);
            break;
        case runtime::coarsening::aggregation:
            process_amg<
                Backend,
                amgcl::coarsening::aggregation
                >(relaxation, func);
            break;
        case runtime::coarsening::smoothed_aggregation:
            process_amg<
                Backend,
                amgcl::coarsening::smoothed_aggregation
                >(relaxation, func);
            break;
        case runtime::coarsening::smoothed_aggr_emin:
            process_amg<
                Backend,
                amgcl::coarsening::smoothed_aggr_emin
                >(relaxation, func);
            break;
    }
}

template <class Backend, class Matrix>
struct amg_create {
    typedef boost::property_tree::ptree params;
    typedef typename Backend::params backend_params;

    void* &handle;

    const Matrix &A;
    const params &p;
    const backend_params &bp;

    amg_create(void* &handle, const Matrix &A, const params &p, const backend_params &bp)
        : handle(handle), A(A), p(p), bp(bp) {}

    template <class AMG>
    void process() const {
        handle = static_cast<void*>( new AMG(A, p, bp) );
    }
};

struct amg_destroy {
    void *handle;

    amg_destroy(void *handle) : handle(handle) {}

    template <class AMG>
    void process() const {
        delete static_cast<AMG*>(handle);
    }
};

template <class Vec1, class Vec2>
struct amg_cycle {
    void       *handle;

    Vec1 const &rhs;
    Vec2       &x;

    amg_cycle(void *handle, const Vec1 &rhs, Vec2 &x)
        : handle(handle), rhs(rhs), x(x) {}

    template <class AMG>
    void process() const {
        static_cast<AMG*>(handle)->cycle(rhs, x);
    }
};

template <class Vec1, class Vec2>
struct amg_apply {
    void       *handle;

    Vec1 const &rhs;
    Vec2       &x;

    amg_apply(void *handle, const Vec1 &rhs, Vec2 &x)
        : handle(handle), rhs(rhs), x(x) {}

    template <class AMG>
    void process() const {
        static_cast<AMG*>(handle)->apply(rhs, x);
    }
};

template <class Matrix>
struct amg_system_matrix {
    void * handle;
    mutable const Matrix * matrix;

    amg_system_matrix(void * handle) : handle(handle), matrix(0) {}

    template <class AMG>
    void process() const {
        matrix = &(static_cast<AMG*>(handle)->system_matrix());
    }
};

struct amg_print {
    void * handle;
    std::ostream &os;

    amg_print(void * handle, std::ostream &os) : handle(handle), os(os) {}

    template <class AMG>
    void process() const {
        os << *static_cast<AMG*>(handle);
    }
};

struct amg_get_params {
    void * handle;
    boost::property_tree::ptree &p;

    amg_get_params(void * handle, boost::property_tree::ptree &p)
        : handle(handle), p(p) {}

    template <class AMG>
    void process() const {
        static_cast<AMG*>(handle)->prm.get(p, "precond.");
    }
};

} // namespace detail

/// Runtime-configurable AMG preconditioner.
template <class Backend>
class amg : boost::noncopyable {
    public:
        typedef Backend backend_type;

        typedef typename Backend::value_type value_type;
        typedef typename Backend::matrix     matrix;
        typedef typename Backend::vector     vector;
        typedef typename Backend::params     backend_params;

        typedef boost::property_tree::ptree params;

        /** Constructs the AMG hierarchy for the system matrix \p A.
         * \rst
         * ``prm`` is an instance of ``boost::property_tree::ptree`` class.
         * The property tree may contain parameters "coarsening.type" and
         * "relax.type". Default values are
         * ``amgcl::runtime::coarsening::smoothed_aggregation`` and
         * ``runtime::relaxation::spai0``.
         * The rest of the property tree should copy the structure of
         * the corresponding ``amgcl::amg::params`` struct. For example, when
         * smoothed aggregation is selected for coarsening, one could:
         *
         * .. code-block:: cpp
         *
         *   prm.put("coarsening.aggr.eps_strong", 1e-2);
         *
         * .. note::
         *
         *   Any parameters that are not relevant to the selected AMG
         *   components are silently ignored.
         * \endrst
         */
        template <class Matrix>
        amg(
                const Matrix &A,
                const params &prm = params(),
                const backend_params &backend_prm = backend_params()
           )
          : coarsening(prm.get("coarsening.type", runtime::coarsening::smoothed_aggregation)),
            relaxation(prm.get("relax.type", runtime::relaxation::spai0)),
            handle(0)
        {
            runtime::detail::process_amg<Backend>(
                    coarsening, relaxation,
                    runtime::detail::amg_create<Backend, Matrix>(handle, A, prm, backend_prm)
                    );
        }

        // Destructor.
        ~amg() {
            runtime::detail::process_amg<Backend>(
                    coarsening, relaxation,
                    runtime::detail::amg_destroy(handle)
                    );
        }

        // Fills the property tree with the actual parameters used.
        void get_params(boost::property_tree::ptree &p) const {
            runtime::detail::process_amg<Backend>(
                    coarsening, relaxation,
                    runtime::detail::amg_get_params(handle, p)
                    );
        }

        /** Performs single V-cycle for the given right-hand side \p rhs and
         * solution \p x.
         */
        template <class Vec1, class Vec2>
        void cycle(const Vec1 &rhs, Vec2 &x) const {
            runtime::detail::process_amg<Backend>(
                    coarsening, relaxation,
                    runtime::detail::amg_cycle<Vec1, Vec2>(handle, rhs, x)
                    );
        }

        /** Performs single V-cycle for the given right-hand side \p rhs after
         * clearing \p x.  This is intended for use as a preconditioning
         * procedure.
         */
        template <class Vec1, class Vec2>
        void apply(const Vec1 &rhs, Vec2 &x) const {
            runtime::detail::process_amg<Backend>(
                    coarsening, relaxation,
                    runtime::detail::amg_apply<Vec1, Vec2>(handle, rhs, x)
                    );
        }

        /** Returns the system matrix in the backend format */
        const matrix& system_matrix() const {
            runtime::detail::amg_system_matrix<matrix> top(handle);
            runtime::detail::process_amg<Backend>(
                    coarsening, relaxation, top
                    );
            return *top.matrix;
        }

        /** Returns the problem size at the finest level. */
        size_t size() const {
            return backend::rows( system_matrix() );
        }

        /// Prints some info about the AMG hierarchy to the output stream.
        friend std::ostream& operator<<(std::ostream &os, const amg &a)
        {
            os << "coarsening:          " << a.coarsening << std::endl;
            os << "relaxation:          " << a.relaxation << std::endl;

            runtime::detail::process_amg<Backend>(
                    a.coarsening, a.relaxation,
                    runtime::detail::amg_print(a.handle, os)
                    );
            return os;
        }

    private:
        const runtime::coarsening::type coarsening;
        const runtime::relaxation::type relaxation;

        void *handle;
};

/// Iterative solvers.
namespace solver {
enum type {
    cg,         ///< Conjugate gradients method
    bicgstab,   ///< BiConjugate Gradient Stabilized
    bicgstabl,  ///< BiCGStab(ell)
    gmres       ///< GMRES
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
    else
        throw std::invalid_argument("Invalid solver value");

    return in;
}

} // namespace solver

namespace detail {

template <
    class Backend,
    class InnerProduct,
    class Func
    >
inline void process_solver(
        runtime::solver::type solver,
        const Func &func
        )
{
    switch (solver) {
        case runtime::solver::cg:
            {
                typedef amgcl::solver::cg<Backend, InnerProduct> Solver;
                func.template process<Solver>();
            }
            break;
        case runtime::solver::bicgstab:
            {
                typedef amgcl::solver::bicgstab<Backend, InnerProduct> Solver;
                func.template process<Solver>();
            }
            break;
        case runtime::solver::bicgstabl:
            {
                typedef amgcl::solver::bicgstabl<Backend, InnerProduct> Solver;
                func.template process<Solver>();
            }
            break;
        case runtime::solver::gmres:
            {
                typedef amgcl::solver::gmres<Backend, InnerProduct> Solver;
                func.template process<Solver>();
            }
            break;
    }
}

template <class Backend>
struct solver_create {
    typedef boost::property_tree::ptree params;
    typedef typename Backend::params backend_params;

    void * &handle;
    size_t n;
    const params &sprm;
    const backend_params &bprm;

    solver_create(void * &handle, size_t n, const params &sprm, const backend_params &bprm)
        : handle(handle), n(n), sprm(sprm), bprm(bprm) {}

    template <class Solver>
    void process() const {
        handle = static_cast<void*>(new Solver(n, sprm, bprm));
    }
};

struct solver_destroy {
    void * handle;

    solver_destroy(void * handle) : handle(handle) {}

    template <class Solver>
    void process() const {
        delete static_cast<Solver*>(handle);
    }
};

struct solver_get_params {
    void * handle;
    boost::property_tree::ptree &p;

    solver_get_params(void * handle, boost::property_tree::ptree &p)
        : handle(handle), p(p) {}

    template <class Solver>
    void process() const {
        static_cast<Solver*>(handle)->prm.get(p, "");
    }
};

template <
    class Matrix,
    class Precond,
    class Vec1,
    class Vec2
    >
struct solver_solve {
    typedef typename Precond::backend_type::value_type value_type;
    typedef typename math::scalar_of<value_type>::type scalar_type;

    void * handle;

    Matrix  const &A;
    Precond const &P;
    Vec1    const &rhs;
    Vec2          &x;

    size_t      &iters;
    scalar_type &resid;

    solver_solve(void * handle, const Matrix &A, const Precond &P,
            const Vec1 &rhs, Vec2 &x, size_t &iters, scalar_type &resid
            )
        : handle(handle), A(A), P(P), rhs(rhs), x(x),
          iters(iters), resid(resid)
    {}

    template <class Solver>
    void process() const {
        boost::tie(iters, resid) = static_cast<Solver*>(handle)->operator()(
                A, P, rhs, x);
    }
};

} // namespace detail

/** This is runtime wrapper around AMGCL iterative solver types. Allows to
 * select the actual solver at runtime.
 */
template <
    class Backend,
    class InnerProduct = amgcl::solver::detail::default_inner_product
    >
class iterative_solver {
    public:
        typedef Backend backend_type;
        typedef boost::property_tree::ptree   params;
        typedef typename Backend::params backend_params;

        typedef typename Backend::value_type value_type;
        typedef typename Backend::vector     vector;

        typedef typename math::scalar_of<value_type>::type scalar_type;

        /** Constructs the iterative solver for the problem size \p n.
         * \rst
         * The property tree ``solver_prm`` may contain "type" entry that would
         * determine the actual type of the iterative solver. Default value:
         * ``amgcl::runtime::solver::bicgstab``.
         *
         * .. note::
         *
         *   Any parameters that are not relevant to the selected solver are
         *   silently ignored.
         * \endrst
         */
        iterative_solver(
                size_t n,
                const params &solver_prm = params(),
                const backend_params &backend_prm = backend_params()
                )
            : solver(solver_prm.get("type", runtime::solver::bicgstab)),
              handle(0)
        {
            runtime::detail::process_solver<Backend, InnerProduct>(
                    solver,
                    runtime::detail::solver_create<Backend>(
                        handle, n, solver_prm, backend_prm
                        )
                    );
        }

        // Destructor.
        ~iterative_solver() {
            runtime::detail::process_solver<Backend, InnerProduct>(
                    solver,
                    runtime::detail::solver_destroy(handle)
                    );
        }

        /** Computes the solution for the given system matrix \p A and the
         * right-hand side \p rhs.  Returns the number of iterations made and
         * the achieved residual as a ``boost::tuple``. The solution vector
         * \p x provides initial approximation in input and holds the computed
         * solution on output.
         *
         * \rst
         * The system matrix may differ from the matrix used during
         * initialization. This may be used for the solution of non-stationary
         * problems with slowly changing coefficients. There is a strong chance
         * that a preconditioner built for a time step will act as a reasonably
         * good preconditioner for several subsequent time steps [DeSh12]_.
         * \endrst
         */
        template <class Matrix, class Precond, class Vec1, class Vec2>
        boost::tuple<size_t, scalar_type> operator()(
                Matrix  const &A,
                Precond const &P,
                Vec1    const &rhs,
#ifdef BOOST_NO_CXX11_RVALUE_REFERENCES
                Vec2          &x
#else
                Vec2          &&x
#endif
                ) const
        {
            size_t      iters = 0;
            scalar_type resid = 0;

            runtime::detail::process_solver<Backend, InnerProduct>(
                    solver,
                    runtime::detail::solver_solve<Matrix, Precond, Vec1, Vec2>(
                        handle, A, P, rhs, x, iters, resid)
                    );

            return boost::make_tuple(iters, resid);
        }

        /** Computes the solution for the given right-hand side \p rhs. The
         * system matrix is the same that was used for the setup of the
         * preconditioner \p P.  Returns the number of iterations made and the
         * achieved residual as a ``boost::tuple``. The solution vector \p x
         * provides initial approximation in input and holds the computed
         * solution on output.
         */
        template <class Precond, class Vec1, class Vec2>
        boost::tuple<size_t, scalar_type> operator()(
                Precond const &P,
                Vec1    const &rhs,
#ifdef BOOST_NO_CXX11_RVALUE_REFERENCES
                Vec2          &x
#else
                Vec2          &&x
#endif
                ) const
        {
            return (*this)(P.system_matrix(), P, rhs, x);
        }

    private:
        const runtime::solver::type  solver;
        void                        *handle;
};

} // namespace runtime
} // namespace amgcl


#endif
