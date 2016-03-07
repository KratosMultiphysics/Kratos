#ifndef AMGCL_MPI_RUNTIME_HPP
#define AMGCL_MPI_RUNTIME_HPP

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
    class LocalPrecond,
    template <class, class> class IterativeSolver,
    class Func
    >
inline void process_sdd(
        runtime::direct_solver::type direct_solver,
        const Func &func
        )
{
    typedef typename LocalPrecond::backend_type::value_type value_type;

    switch (direct_solver) {
        case amgcl::runtime::direct_solver::skyline_lu:
            {
                typedef amgcl::mpi::subdomain_deflation<
                    LocalPrecond,
                    IterativeSolver,
                    amgcl::mpi::skyline_lu<value_type>
                    > SDD;
                func.template process<SDD>();
            }
            break;
#ifdef AMGCL_HAVE_PASTIX
        case amgcl::runtime::direct_solver::pastix:
            {
                typedef amgcl::mpi::subdomain_deflation<
                    LocalPrecond,
                    IterativeSolver,
                    amgcl::mpi::PaStiX<value_type>
                    > SDD;
                func.template process<SDD>();
            }
            break;
#endif
    }
}

template <
    class LocalPrecond,
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
                LocalPrecond,
                amgcl::solver::cg
                >(direct_solver, func);
            break;
        case runtime::solver::bicgstab:
            process_sdd<
                LocalPrecond,
                amgcl::solver::bicgstab
                >(direct_solver, func);
            break;
        case runtime::solver::bicgstabl:
            process_sdd<
                LocalPrecond,
                amgcl::solver::bicgstabl
                >(direct_solver, func);
            break;
        case runtime::solver::gmres:
            process_sdd<
                LocalPrecond,
                amgcl::solver::gmres
                >(direct_solver, func);
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

template <typename value_type, class Vec1, class Vec2>
struct sdd_solve {
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
template <class LocalPrecond>
class subdomain_deflation : boost::noncopyable {
    public:
        typedef typename LocalPrecond::backend_type Backend;
        typedef typename Backend::value_type value_type;
        typedef boost::property_tree::ptree params;

        template <class Matrix, class DefVec>
        subdomain_deflation(MPI_Comm comm, const Matrix &A, const DefVec &def_vec, const params &prm)
            : iterative_solver(prm.get("solver.type", amgcl::runtime::solver::bicgstabl)),
              direct_solver(prm.get("direct_solver.type",
#ifdef AMGCL_HAVE_PASTIX
                          amgcl::runtime::direct_solver::pastix
#else
                          amgcl::runtime::direct_solver::skyline_lu
#endif
                          )),
              n( backend::rows(A) ), handle(0)
        {
            runtime::mpi::detail::process_sdd<LocalPrecond>(
                    iterative_solver,
                    direct_solver,
                    runtime::mpi::detail::sdd_create<Matrix, DefVec>(
                        handle, comm, A, def_vec, prm
                        )
                    );
        }

        ~subdomain_deflation() {
            runtime::mpi::detail::process_sdd<LocalPrecond>(
                    iterative_solver,
                    direct_solver,
                    runtime::mpi::detail::sdd_destroy(handle)
                    );
        }

        void get_params(boost::property_tree::ptree &p) const {
            runtime::mpi::detail::process_sdd<LocalPrecond>(
                    iterative_solver,
                    direct_solver,
                    runtime::mpi::detail::sdd_get_params(handle, p)
                    );

        }

        template <class Vec1, class Vec2>
        boost::tuple<size_t, value_type>
        operator()(const Vec1 &rhs, Vec2 &x) const {
            size_t     iters = 0;
            value_type resid = 0;

            runtime::mpi::detail::process_sdd<LocalPrecond>(
                    iterative_solver,
                    direct_solver,
                    runtime::mpi::detail::sdd_solve<value_type, Vec1, Vec2>(
                        handle, rhs, x, iters, resid
                        )
                    );

            return boost::make_tuple(iters, resid);
        }

        size_t local_size() const {
            return n;
        }
    private:
        const runtime::solver::type        iterative_solver;
        const runtime::direct_solver::type direct_solver;

        size_t n;
        void * handle;
};

} // namespace runtime
} // namespace mpi
} // namespace amgcl

#endif
