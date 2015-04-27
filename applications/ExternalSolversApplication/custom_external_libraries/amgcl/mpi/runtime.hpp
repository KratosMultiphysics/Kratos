#ifndef AMGCL_MPI_RUNTIME_HPP
#define AMGCL_MPI_RUNTIME_HPP

/*
The MIT License

Copyright (c) 2012-2015 Denis Demidov <dennis.demidov@gmail.com>

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
 * \file   amgcl/mpi/runtime.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Runtime-configurable wrappers around amgcl mpi classes.
 */

#include <amgcl/mpi/subdomain_deflation.hpp>
#include <amgcl/runtime.hpp>

#ifdef AMGCL_HAVE_PASTIX
#  include <amgcl/mpi/pastix.hpp>
#endif

namespace amgcl {

namespace runtime {

/// Direct solvers.
namespace direct_solver {
enum type {
    skyline_lu
#ifdef AMGCL_HAVE_PASTIX
  , pastix
#endif
};

std::ostream& operator<<(std::ostream &os, type s)
{
    switch (s) {
        case skyline_lu:
            return os << "skyline_lu";
#ifdef AMGCL_HAVE_PASTIX
        case pastix:
            return os << "pastix";
#endif
        default:
            return os << "???";
    }
}

std::istream& operator>>(std::istream &in, type &s)
{
    std::string val;
    in >> val;

    if (val == "skyline_lu")
        s = skyline_lu;
#ifdef AMGCL_HAVE_PASTIX
    else if (val == "pastix")
        s = pastix;
#endif
    else
        throw std::invalid_argument("Invalid direct solver value");

    return in;
}

} // namespace direct_solver

/// Distributed algorithms and structures.
namespace mpi {

namespace detail {

template <
    class Backend,
    class Coarsening,
    template <class> class Relaxation,
    template <class, class> class IterativeSolver,
    class Func
    >
inline
typename boost::enable_if<
    typename backend::relaxation_is_supported<Backend, Relaxation>::type,
    void
>::type
process_sdd(
        runtime::direct_solver::type direct_solver,
        const Func &func
        )
{
    switch (direct_solver) {
        case amgcl::runtime::direct_solver::skyline_lu:
            {
                typedef amgcl::mpi::subdomain_deflation<
                    Backend,
                    Coarsening,
                    Relaxation,
                    IterativeSolver,
                    amgcl::mpi::skyline_lu<typename Backend::value_type>
                    > SDD;
                func.template process<SDD>();
            }
            break;
#ifdef AMGCL_HAVE_PASTIX
        case amgcl::runtime::direct_solver::pastix:
            {
                typedef amgcl::mpi::subdomain_deflation<
                    Backend,
                    Coarsening,
                    Relaxation,
                    IterativeSolver,
                    amgcl::mpi::PaStiX<typename Backend::value_type>
                    > SDD;
                func.template process<SDD>();
            }
            break;
#endif
    }
}

template <
    class Backend,
    class Coarsening,
    template <class> class Relaxation,
    template <class, class> class IterativeSolver,
    class Func
    >
inline
typename boost::disable_if<
    typename backend::relaxation_is_supported<Backend, Relaxation>::type,
    void
>::type
process_sdd(
        runtime::direct_solver::type,
        const Func&
        )
{
    throw std::logic_error("The relaxation scheme is not supported by the backend");
}

template <
    class Backend,
    class Coarsening,
    template <class> class Relaxation,
    class Func
    >
inline void process_sdd(
        runtime::solver::type        iterative_solver,
        runtime::direct_solver::type direct_solver,
        const Func &func
        )
{
    switch (iterative_solver) {
        case runtime::solver::cg:
            process_sdd<
                Backend,
                Coarsening,
                Relaxation,
                amgcl::solver::cg
                >(direct_solver, func);
            break;
        case runtime::solver::bicgstab:
            process_sdd<
                Backend,
                Coarsening,
                Relaxation,
                amgcl::solver::bicgstab
                >(direct_solver, func);
            break;
        case runtime::solver::bicgstabl:
            process_sdd<
                Backend,
                Coarsening,
                Relaxation,
                amgcl::solver::bicgstabl
                >(direct_solver, func);
            break;
        case runtime::solver::gmres:
            process_sdd<
                Backend,
                Coarsening,
                Relaxation,
                amgcl::solver::gmres
                >(direct_solver, func);
            break;
    }
}

template <
    class Backend,
    class Coarsening,
    class Func
    >
inline void process_sdd(
        runtime::relaxation::type    relaxation,
        runtime::solver::type        iterative_solver,
        runtime::direct_solver::type direct_solver,
        const Func &func
        )
{
    switch (relaxation) {
        case runtime::relaxation::gauss_seidel:
            process_sdd<
                Backend,
                Coarsening,
                amgcl::relaxation::gauss_seidel
                >(iterative_solver, direct_solver, func);
            break;
        case runtime::relaxation::multicolor_gauss_seidel:
            process_sdd<
                Backend,
                Coarsening,
                amgcl::relaxation::gauss_seidel
                >(iterative_solver, direct_solver, func);
            break;
        case runtime::relaxation::ilu0:
            process_sdd<
                Backend,
                Coarsening,
                amgcl::relaxation::ilu0
                >(iterative_solver, direct_solver, func);
            break;
        case runtime::relaxation::damped_jacobi:
            process_sdd<
                Backend,
                Coarsening,
                amgcl::relaxation::damped_jacobi
                >(iterative_solver, direct_solver, func);
            break;
        case runtime::relaxation::spai0:
            process_sdd<
                Backend,
                Coarsening,
                amgcl::relaxation::spai0
                >(iterative_solver, direct_solver, func);
            break;
        case runtime::relaxation::spai1:
            process_sdd<
                Backend,
                Coarsening,
                amgcl::relaxation::spai1
                >(iterative_solver, direct_solver, func);
            break;
        case runtime::relaxation::chebyshev:
            process_sdd<
                Backend,
                Coarsening,
                amgcl::relaxation::chebyshev
                >(iterative_solver, direct_solver, func);
            break;
    }
}

template <
    class Backend,
    class Func
    >
inline void process_sdd(
        runtime::coarsening::type    coarsening,
        runtime::relaxation::type    relaxation,
        runtime::solver::type        iterative_solver,
        runtime::direct_solver::type direct_solver,
        const Func &func
        )
{
    switch (coarsening) {
        case runtime::coarsening::ruge_stuben:
            process_sdd<
                Backend,
                amgcl::coarsening::ruge_stuben
                >(relaxation, iterative_solver, direct_solver, func);
            break;
        case runtime::coarsening::aggregation:
            process_sdd<
                Backend,
                amgcl::coarsening::aggregation
                >(relaxation, iterative_solver, direct_solver, func);
            break;
        case runtime::coarsening::smoothed_aggregation:
            process_sdd<
                Backend,
                amgcl::coarsening::smoothed_aggregation
                >(relaxation, iterative_solver, direct_solver, func);
            break;
        case runtime::coarsening::smoothed_aggr_emin:
            process_sdd<
                Backend,
                amgcl::coarsening::smoothed_aggr_emin
                >(relaxation, iterative_solver, direct_solver, func);
            break;
    }
}

template <class Matrix, class DefVec>
struct sdd_create {
    typedef boost::property_tree::ptree params;

    void * &handle;

    MPI_Comm comm;

    const Matrix &A;
    const DefVec &def_vec;
    const params &prm;

    sdd_create(
            void * &handle, MPI_Comm comm, const Matrix &A,
            const DefVec &def_vec, const params &prm
            )
        : handle(handle), comm(comm), A(A), def_vec(def_vec), prm(prm)
    {}

    template <class SDD>
    void process() const {
        handle = static_cast<void*>( new SDD(comm, A, def_vec, prm) );
    }
};

struct sdd_destroy {
    void * handle;

    sdd_destroy(void * handle) : handle(handle) {}

    template <class SDD>
    void process() const {
        delete static_cast<SDD*>(handle);
    }
};

struct sdd_get_params {
    void * handle;
    boost::property_tree::ptree &p;

    sdd_get_params(void * handle, boost::property_tree::ptree &p)
        : handle(handle), p(p) {}

    template <class SDD>
    void process() const {
        static_cast<SDD*>(handle)->get_params(p);
    }
};

template <class Backend, class Vec1, class Vec2>
struct sdd_solve {
    typedef typename Backend::value_type value_type;

    void * handle;

    Vec1 const &rhs;
    Vec2       &x;

    size_t     &iters;
    value_type &resid;

    sdd_solve(void * handle, const Vec1 &rhs, Vec2 &x,
            size_t &iters, value_type &resid)
        : handle(handle), rhs(rhs), x(x), iters(iters), resid(resid)
    {}

    template <class SDD>
    void process() const {
        boost::tie(iters, resid) = static_cast<SDD*>(handle)->operator()(rhs, x);
    }
};

} // namespace detail

/// Runtime-configurable distributed solver based on subdomain deflation.
/**
 * \sa \cite Frank2001
 */
template <class Backend>
class subdomain_deflation : boost::noncopyable {
    public:
        typedef typename Backend::value_type value_type;
        typedef boost::property_tree::ptree params;

        template <class Matrix, class DefVec>
        subdomain_deflation(
                runtime::coarsening::type    coarsening,
                runtime::relaxation::type    relaxation,
                runtime::solver::type        iterative_solver,
                runtime::direct_solver::type direct_solver,
                MPI_Comm comm,
                const Matrix &A,
                const DefVec &def_vec,
                const params &prm
                )
            : coarsening(coarsening),
              relaxation(relaxation),
              iterative_solver(iterative_solver),
              direct_solver(direct_solver),
              n( backend::rows(A) ),
              handle(0)
        {
            runtime::mpi::detail::process_sdd<Backend>(
                    coarsening,
                    relaxation,
                    iterative_solver,
                    direct_solver,
                    runtime::mpi::detail::sdd_create<Matrix, DefVec>(
                        handle, comm, A, def_vec, prm
                        )
                    );
        }

        template <class Matrix, class DefVec>
        subdomain_deflation(
                MPI_Comm comm,
                const Matrix &A,
                const DefVec &def_vec,
                const params &prm
                )
            : coarsening(amgcl::runtime::coarsening::smoothed_aggregation),
              relaxation(amgcl::runtime::relaxation::spai0),
              iterative_solver(amgcl::runtime::solver::bicgstabl),
              direct_solver(amgcl::runtime::direct_solver::skyline_lu),
              n( backend::rows(A) ), handle(0)
        {
            runtime::mpi::detail::process_sdd<Backend>(
                    coarsening,
                    relaxation,
                    iterative_solver,
                    direct_solver,
                    runtime::mpi::detail::sdd_create<Matrix, DefVec>(
                        handle, comm, A, def_vec, prm
                        )
                    );
        }

        ~subdomain_deflation() {
            runtime::mpi::detail::process_sdd<Backend>(
                    coarsening,
                    relaxation,
                    iterative_solver,
                    direct_solver,
                    runtime::mpi::detail::sdd_destroy(handle)
                    );
        }

        void get_params(boost::property_tree::ptree &p) const {
            runtime::mpi::detail::process_sdd<Backend>(
                    coarsening,
                    relaxation,
                    iterative_solver,
                    direct_solver,
                    runtime::mpi::detail::sdd_get_params(handle, p)
                    );

        }

        template <class Vec1, class Vec2>
        boost::tuple<size_t, value_type>
        operator()(const Vec1 &rhs, Vec2 &x) const {
            size_t     iters;
            value_type resid;

            runtime::mpi::detail::process_sdd<Backend>(
                    coarsening,
                    relaxation,
                    iterative_solver,
                    direct_solver,
                    runtime::mpi::detail::sdd_solve<Backend, Vec1, Vec2>(
                        handle, rhs, x, iters, resid
                        )
                    );

            return boost::make_tuple(iters, resid);
        }

        size_t local_size() const {
            return n;
        }
    private:
        const runtime::coarsening::type    coarsening;
        const runtime::relaxation::type    relaxation;
        const runtime::solver::type        iterative_solver;
        const runtime::direct_solver::type direct_solver;

        size_t n;
        void * handle;
};

} // namespace runtime
} // namespace mpi
} // namespace amgcl

#endif
