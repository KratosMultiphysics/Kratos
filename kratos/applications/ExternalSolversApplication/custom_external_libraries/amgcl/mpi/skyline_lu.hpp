#ifndef AMGCL_MPI_SKYLINE_LU_HPP
#define AMGCL_MPI_SKYLINE_LU_HPP

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
\file   amgcl/mpi/skyline_lu.hpp
\author Denis Demidov <dennis.demidov@gmail.com>
\brief  MPI wrapper for Skyline LU factorization solver.

This is a wrapper around Skyline LU factorization solver that provides a
distributed direct solver interface but always works sequentially.
*/

#include <amgcl/adapter/crs_tuple.hpp>
#include <amgcl/solver/skyline_lu.hpp>

namespace amgcl {
namespace mpi {

/// Provides distributed direct solver interface for Skyline LU solver.
template <typename value_type>
class skyline_lu {
    public:
        typedef typename solver::skyline_lu<value_type>::params params;

        /// The number of processes optimal for the given problem size.
        static int comm_size(int /*n_global_rows*/) {
            return 1;
        }

        /// Constructor.
        /**
         * \param comm MPI communicator containing processes to participate in
         *        solution of the problem. The number of processes in
         *        communicator should be (but not necessarily) equal to the
         *        result of comm_size().
         * \param n_local_rows Number of matrix rows belonging to the calling
         *        process.
         * \param ptr Start of each row in col and val arrays.
         * \param col Column numbers of nonzero elements.
         * \param val Values of nonzero elements.
         * \param prm Solver parameters.
         */
        template <class PRng, class CRng, class VRng>
        skyline_lu(
                MPI_Comm,
                int n_local_rows,
                const PRng &ptr,
                const CRng &col,
                const VRng &val,
                const params &prm = params()
                ) : S( boost::tie(n_local_rows, ptr, col, val), prm )
        {}

        /// Solves the problem for the given right-hand side.
        /**
         * \param rhs The right-hand side.
         * \param x   The solution.
         */
        template <class Vec1, class Vec2>
        void operator()(const Vec1 &rhs, Vec2 &x) const {
            S(rhs, x);
        }
    private:
        solver::skyline_lu<value_type> S;
};

}
}

#endif
