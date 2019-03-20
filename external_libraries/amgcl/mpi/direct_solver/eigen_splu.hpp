#ifndef AMGCL_MPI_DIRECT_SOLVER_EIGEN_SPLU_HPP
#define AMGCL_MPI_DIRECT_SOLVER_EIGEN_SPLU_HPP

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
\file   amgcl/mpi/direct_solver/eigen_splu.hpp
\author Denis Demidov <dennis.demidov@gmail.com>
\brief  MPI wrapper for Eigen::SparseLU solver.

This is a wrapper around Eigen SparseLU solver that provides a
distributed direct solver interface but always works sequentially.
*/

#include <mpi.h>

#include <memory>

#include <Eigen/SparseLU>
#include <amgcl/backend/builtin.hpp>
#include <amgcl/adapter/crs_tuple.hpp>
#include <amgcl/solver/eigen.hpp>
#include <amgcl/mpi/util.hpp>
#include <amgcl/mpi/direct_solver/solver_base.hpp>

namespace amgcl {
namespace mpi {
namespace direct {

/// Provides distributed direct solver interface for Skyline LU solver.
template <typename value_type>
class eigen_splu : public solver_base< value_type, eigen_splu<value_type> > {
    public:
        typedef
            amgcl::solver::EigenSolver<
                Eigen::SparseLU<
                    Eigen::SparseMatrix<value_type, Eigen::ColMajor, int>
                    >
                >
            Solver;
        typedef typename Solver::params params;
        typedef backend::crs<value_type> build_matrix;

        /// Constructor.
        template <class Matrix>
        eigen_splu(communicator comm, const Matrix &A,
                const params &prm = params()) : prm(prm)
        {
            static_cast<Base*>(this)->init(comm, A);
        }

        static size_t coarse_enough() {
            return Base::coarse_enough();
        }

        int comm_size(int /*n*/) const {
            return 1;
        }

        void init(communicator, const build_matrix &A) {
            S = std::make_shared<Solver>(A, prm);
        }

        /// Solves the problem for the given right-hand side.
        /**
         * \param rhs The right-hand side.
         * \param x   The solution.
         */
        template <class Vec1, class Vec2>
        void solve(const Vec1 &rhs, Vec2 &x) const {
            (*S)(rhs, x);
        }
    private:
        typedef solver_base< value_type, eigen_splu<value_type> > Base;
        params prm;
        std::shared_ptr<Solver> S;
};

} // namespace direct
} // namespace mpi
} // namespace amgcl

#endif
