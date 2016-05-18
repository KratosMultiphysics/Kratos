#ifndef AMGCL_MPI_DIRECT_SOLVER_HPP
#define AMGCL_MPI_DIRECT_SOLVER_HPP

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
 * \file   amgcl/mpi/direct_solver.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Runtime wrapper for distributed direct solvers.
 */

#include <amgcl/util.hpp>
#include <amgcl/mpi/skyline_lu.hpp>
#ifdef AMGCL_HAVE_PASTIX
#  include <amgcl/mpi/pastix.hpp>
#endif

namespace amgcl {
namespace runtime {
namespace mpi {

/// Direct solvers.
namespace dsolver {
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
} // namespace dsolver

template <class value_type>
class direct_solver {
    public:
        typedef boost::property_tree::ptree params;

        static int comm_size(int n_global_rows, params prm = params()) {
            dsolver::type s = prm.get("type", dsolver::skyline_lu);
            if (!prm.erase("type")) AMGCL_PARAM_MISSING("type");

            switch (s) {
                case dsolver::skyline_lu:
                    {
                        typedef amgcl::mpi::skyline_lu<value_type> S;
                        return S::comm_size(n_global_rows, prm);
                    }
#ifdef AMGCL_HAVE_PASTIX
                case dsolver::pastix:
                    {
                        typedef amgcl::mpi::PaStiX<value_type> S;
                        return S::comm_size(n_global_rows, prm);
                    }
#endif
                default:
                    amgcl::precondition(false, "Unsupported direct solver type");
            }

            return 1;
        }

        template <class PRng, class CRng, class VRng>
        direct_solver(
                MPI_Comm mpi_comm,
                int n_local_rows,
                const PRng &p_ptr,
                const CRng &p_col,
                const VRng &p_val,
                params prm = params()
                )
            : solver(prm.get("type", dsolver::skyline_lu))
        {
            if (!prm.erase("type")) AMGCL_PARAM_MISSING("type");

            switch (solver) {
                case dsolver::skyline_lu:
                    {
                        typedef amgcl::mpi::skyline_lu<value_type> S;
                        handle = static_cast<void*>(
                                new S(mpi_comm, n_local_rows, p_ptr, p_col, p_val, prm));
                    }
                    break;
#ifdef AMGCL_HAVE_PASTIX
                case dsolver::pastix:
                    {
                        typedef amgcl::mpi::PaStiX<value_type> S;
                        handle = static_cast<void*>(
                                new S(mpi_comm, n_local_rows, p_ptr, p_col, p_val, prm));
                    }
                    break;
#endif
                default:
                    amgcl::precondition(false, "Unsupported direct solver type");
            }
        }

        template <class Vec1, class Vec2>
        void operator()(const Vec1 &rhs, Vec2 &x) const {
            switch (solver) {
                case dsolver::skyline_lu:
                    {
                        typedef amgcl::mpi::skyline_lu<value_type> S;
                        static_cast<const S*>(handle)->operator()(rhs, x);
                    }
                    break;
#ifdef AMGCL_HAVE_PASTIX
                case dsolver::pastix:
                    {
                        typedef amgcl::mpi::PaStiX<value_type> S;
                        static_cast<const S*>(handle)->operator()(rhs, x);
                    }
                    break;
#endif
                default:
                    amgcl::precondition(false, "Unsupported direct solver type");
            }
        }

        ~direct_solver() {
            switch (solver) {
                case dsolver::skyline_lu:
                    {
                        typedef amgcl::mpi::skyline_lu<value_type> S;
                        delete static_cast<S*>(handle);
                    }
                    break;
#ifdef AMGCL_HAVE_PASTIX
                case dsolver::pastix:
                    {
                        typedef amgcl::mpi::PaStiX<value_type> S;
                        delete static_cast<S*>(handle);
                    }
                    break;
#endif
                default:
                    amgcl::precondition(false, "Unsupported direct solver type");
            }
        }
    private:
        dsolver::type solver;
        void *handle;
};

} // namespace mpi
} // namespace runtime
} // namespace amgcl

#endif
