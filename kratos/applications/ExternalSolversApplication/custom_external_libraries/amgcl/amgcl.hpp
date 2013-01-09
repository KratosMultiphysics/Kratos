#ifndef AMGCL_AMGCL_HPP
#define AMGCL_AMGCL_HPP

/*
The MIT License

Copyright (c) 2012 Denis Demidov <ddemidov@ksu.ru>

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
 * \file   amgcl.hpp
 * \author Denis Demidov <ddemidov@ksu.ru>
 * \brief  Generic algebraic multigrid framework.
 */

/**
\mainpage amgcl Generic algebraic multigrid framework.

amgcl is a simple and generic algebraic
<a href="http://en.wikipedia.org/wiki/Multigrid_method">multigrid</a> (AMG)
hierarchy builder (and a work in progress).  The constructed hierarchy may be
used as a standalone solver or as a preconditioner with some iterative solver.
Several \ref iterative "iterative solvers" are provided, and it is also
possible to use generic solvers from other libraries, e.g.
<a href="http://viennacl.sourceforge.net">ViennaCL</a>.

The setup phase is completely CPU-based. The constructed levels of AMG
hierarchy may be stored and used through several \ref levels "backends". This
allows for transparent acceleration of the solution phase with help of OpenCL,
CUDA, or OpenMP technologies.  See
<a href="https://github.com/ddemidov/amgcl/blob/master/examples/vexcl.cpp">examples/vexcl.cpp</a>,
<a href="https://github.com/ddemidov/amgcl/blob/master/examples/viennacl.cpp">examples/viennacl.cpp</a> and
<a href="https://github.com/ddemidov/amgcl/blob/master/examples/eigen.cpp">examples/eigen.cpp</a>
for examples of
using amgcl with
<a href="https://github.com/ddemidov/vexcl">VexCL</a>,
<a href="http://viennacl.sourceforge.net">ViennaCL</a>, or
CPU backends.

\section overview Overview

You can use amgcl to solve large sparse system of linear equations in three
simple steps: first, you have to select method components (this is a compile
time decision); second, the AMG hierarchy has to be constructed from a system
matrix; and third, the hierarchy is used to solve the equation system for a
given right-hand side.

The list of interpolation schemes and available backends may be found in
\ref interpolation "Interpolation" and \ref levels "Level Storage Backends" 
modules.  The aggregation and smoothed-aggregation interpolation schemes use
less memory and are set up faster than classic interpolation, but their
convergence rate is slower. They are well suited for GPU-accelerated backends,
where the cost of the setup phase is much more important.

\code
// First, we need to include relevant headers. Each header basically
// corresponds to an AMG component. Let's say we want to use conjugate gradient
// method preconditioned with smoothed aggregation AMG with VexCL backend:

// This is generic hierarchy builder.
#include <amgcl/amgcl.hpp>
// It will use the following components:

// Interpolation scheme based on smoothed aggregation.
#include <amgcl/interp_smoothed_aggr.hpp>
// Aggregates will be constructed with plain aggregation:
#include <amgcl/aggr_plain.hpp>
// VexCL will be used as a backend:
#include <amgcl/level_vexcl.hpp>
// The definition of conjugate gradient method:
#include <amgcl/cg.hpp>

int main() {
    // VexCL context initialization (let's use all GPUs that support double precision):
    vex::Context ctx( vex::Filter::Type(CL_DEVICE_TYPE_GPU) && vex::Filter::DoublePrecision );

    // Here, the system matrix and right-hand side are somehow constructed. The
    // system matrix data is stored in compressed row storage format in vectors
    // row, col, and val.
    int size;
    std::vector<int>    row, col;
    std::vector<double> val, rhs;

    // We wrap the matrix data into amgcl-compatible type.
    // No data is copied here:
    auto A = amgcl::sparse::map(size, size, row.data(), col.data(), val.data());

    // The AMG builder type. Note the use of damped Jacobi relaxation (smoothing) on each level.
    typedef amgcl::solver<
        double, int,
        amgcl::interp::smoothed_aggregation<amgcl::aggr::plain>,
        amgcl::level::vexcl<amgcl::relax::damped_jacobi>
    > AMG;

    // The parameters. Most of the parameters have some reasonable defaults.
    // VexCL backend needs to know what context to use:
    AMG::params prm;
    prm.level.ctx = &ctx;

    // Here we construct the hierarchy:
    AMG amg(A, prm);

    // Now let's solve the system of equations. We need to transfer matrix,
    // right-hand side, and initial approximation to GPUs. The matrix part may
    // be omitted though, since AMG already has it as part of the hierarchy:
    std::vector<double> x(size, 0.0);

    vex::vector<double> f(ctx.queue(), rhs);
    vex::vector<double> u(ctx.queue(), x);

    // Call AMG-preconditioned CG method:
    auto cnv = amgcl::solve(amg.top_matrix(), f, amg, u, amgcl::cg_tag());

    std::cout << "Iterations: " << std::get<0>(cnv) << std::endl
              << "Error:      " << std::get<1>(cnv) << std::endl;

    // Copy the solution back to host:
    vex::copy(u, x);
}
\endcode

The following command line would compile the example:
\verbatim
g++ -o example -std=c++0x -O3 -fopenmp example.cpp -I<path/to/vexcl> -I<path/to/amgcl> -lOpenCL -lboost_chrono
\endverbatim

The C++11 support is enabled here (by -std=c++0x flag) because it is required
by VexCL library. amgcl relies on Boost instead. Also note the use of
`-fopenmp` switch. It enables an OpenMP-based parallelization of the setup
stage.


\section install Installation

The library is header-only, so there is nothing to compile or link to. You just
need to copy amgcl folder somewhere and tell your compiler to scan it for
include files.

\section references References
 -# \anchor Trottenberg_2001 <em>U. Trottenberg, C. Oosterlee, A. Shuller,</em>
    Multigrid, Academic Press, London, 2001.
 -# \anchor Stuben_1999 <em>K. Stuben,</em> Algebraic multigrid (AMG): an
    introduction with applications, Journal of Computational and Applied
     Mathematics,  2001, Vol. 128, Pp. 281-309.
 -# \anchor Vanek_1996 <em>P. Vanek, J. Mandel, M. Brezina,</em> Algebraic multigrid
    by smoothed aggregation for second and fourth order elliptic problems,
    Computing 56, 1996, Pp. 179-196.
 -# \anchor Notay_2008 <em>Y. Notay, P. Vassilevski,</em> Recursive
    Krylov-based multigrid cycles, Numer. Linear Algebra Appl. 2008; 15:473-487.
 -# \anchor Templates_1994 <em>R. Barrett, M. Berry,
    T. F. Chan et al.</em> Templates for the Solution of Linear Systems:
    Building Blocks for Iterative Methods, 2nd Edition, SIAM, Philadelphia, PA,
    1994.
 -# \anchor spai_2002 <em>O. Broeker, M. Grote,</em> Sparse approximate inverse
    smoothers for geometric and algebraic multigrid, Applied Numerical
    Mathematics, Volume 41, Issue 1, April 2002, Pages 61–80.
 -# \anchor Sala_2008 <em>M. Sala, R. Tuminaro,</em> A new Petrov-Galerkin
    smoothed aggregation preconditioner for nonsymmetric linear systems.
    SIAM J. Sci. Comput. 2008, Vol. 31, No.1, pp. 143-166.
*/

#include <iostream>
#include <iomanip>
#include <utility>
#include <list>

#include <boost/static_assert.hpp>
#include <boost/smart_ptr/shared_ptr.hpp>
#include <boost/smart_ptr/make_shared.hpp>
#include <boost/typeof/typeof.hpp>
#include <boost/type_traits/is_signed.hpp>

#include <amgcl/spmat.hpp>
#include <amgcl/tictoc.hpp>

/// Primary namespace for the library.
namespace amgcl {

/// Interpolation-related types and functions.
namespace interp {

/// Galerkin operator.
struct galerkin_operator {
    template <class spmat, class Params>
    static spmat apply(const spmat &R, const spmat &A, const spmat &P,
            const Params &prm)
    {
        return sparse::prod(sparse::prod(R, A), P);
    }
};

/// Returns coarse level construction scheme for a given interpolation scheme.
/**
 * By default, Galerkin operator is used to construct coarse level from system
 * matrix, restriction and prolongation operators:
 * \f[A^H = R A^h P.\f] Usually, \f$R = P^T.\f$
 *
 * \param Interpolation interpolation scheme.
 */
template <class Interpolation>
struct coarse_operator {
    typedef galerkin_operator type;
};

} // namespace interp

/// Possible relaxation (smoothing) schemes.
namespace relax {

/// Possible relaxation (smoothing) schemes.
/**
 * Each backend may support only a limited subset of these.
 * \sa relax_vs_backend
 */
enum scheme {
    damped_jacobi, ///< Damped Jacobi.
    spai0,         ///< SPAI-0 algorithm from \ref spai_2002 "Broeker (2002)".
    gauss_seidel,  ///< Gauss-Seidel.
    ilu0           ///< Incomplete LU decomposition with zero fill-in.
};

} // namespace relax.

/// Algebraic multigrid method.
/**
 * \param value_t  Type for matrix entries (double/float).
 * \param index_t  Type for matrix indices. Should be signed integral type.
 * \param interp_t \ref interpolation "Interpolation scheme".
 * \param level_t  Hierarchy level \ref levels "storage backend".
 */
template <
    typename value_t, typename index_t, typename interp_t, typename level_t
    >
class solver {
    private:
        typedef sparse::matrix<value_t, index_t> matrix;
        typedef typename level_t::template instance<value_t, index_t> level_type;

    public:
        typedef value_t value_type;
        typedef index_t index_type;

        /// Parameters for AMG components.
        struct params {
            /// When level is coarse enough to be solved directly?
            /**
             * If number of variables at a next level in hierarchy becomes
             * lower than this threshold, then the hierarchy construction is
             * stopped and the linear system is solved explicitly at this
             * level.
             */
            unsigned coarse_enough;

            typename interp_t::params interp; ///< Interpolation parameters.
            typename level_t::params  level;  ///< Level/Solution parameters.

            params() : coarse_enough(300) { }
        };

        /// Constructs the AMG hierarchy from the system matrix.
        /** 
         * The input matrix is copied here and may be freed afterwards.
         *
         * \param A   The system matrix. Should be convertible to
         *            amgcl::sparse::matrix<>.
         * \param prm Parameters controlling the setup and solution phases.
         *
         * \sa amgcl::sparse::map()
         */
        template <typename spmat>
        solver(const spmat &A, const params &prm = params()) : prm(prm)
        {
            BOOST_STATIC_ASSERT_MSG(boost::is_signed<index_t>::value,
                    "Matrix index type should be signed");

            matrix copy(A);
            build_level(copy, prm);
        }

        /// The AMG hierarchy is used as a standalone solver.
        /** 
         * The vector types should be compatible with level_t:
         *
         * -# Any type with operator[] should work on a CPU.
         * -# vex::vector<value_t> should be used with level::vexcl.
         * -# viennacl::vector<value_t> should be used with level::ViennaCL.
         *
         * \param rhs Right-hand side.
         * \param x   Solution. Contains an initial approximation on input, and
         *            the approximated solution on output.
         */
        template <class vector1, class vector2>
        std::pair< int, value_t > solve(const vector1 &rhs, vector2 &x) const {
            int     iter = 0;
            value_t res  = 2 * prm.level.tol;

            for(; res > prm.level.tol && iter < prm.level.maxiter; ++iter) {
                apply(rhs, x);
                res = hier.front()->resid(rhs, x);
            }

            return std::make_pair(iter, res);
        }

        /// Performs single multigrid cycle.
        /**
         * Is intended to be used as a preconditioner with iterative methods.
         *
         * The vector types should be compatible with level_t:
         *
         * -# Any type with operator[] should work on a CPU.
         * -# vex::vector<value_t> should be used with level::vexcl.
         * -# viennacl::vector<value_t> should be used with level::ViennaCL.
         *
         * \param rhs Right-hand side.
         * \param x   Solution. Contains an initial approximation on input, and
         *            the approximated solution on output.
         */
        template <class vector1, class vector2>
        void apply(const vector1 &rhs, vector2 &x) const {
            level_type::cycle(hier.begin(), hier.end(), prm.level, rhs, x);
        }

        /// Output some general information about the AMG hierarchy.
        std::ostream& print(std::ostream &os) const {
            BOOST_AUTO(ff, os.flags());
            BOOST_AUTO(pp, os.precision());

            index_t sum_dof = 0;
            index_t sum_nnz = 0;
            for(BOOST_AUTO(lvl, hier.begin()); lvl != hier.end(); ++lvl) {
                sum_dof += (*lvl)->size();
                sum_nnz += (*lvl)->nonzeros();
            }

            os << "Number of levels:    "   << hier.size()
               << "\nOperator complexity: " << std::fixed << std::setprecision(2)
                                            << 1.0 * sum_nnz / hier.front()->nonzeros()
               << "\nGrid complexity:     " << std::fixed << std::setprecision(2)
                                            << 1.0 * sum_dof / hier.front()->size()
               << "\n\nlevel     unknowns       nonzeros\n"
               << "---------------------------------\n";

            index_t depth = 0;
            for(BOOST_AUTO(lvl, hier.begin()); lvl != hier.end(); ++lvl, ++depth)
                os << std::setw(5)  << depth
                   << std::setw(13) << (*lvl)->size()
                   << std::setw(15) << (*lvl)->nonzeros() << " ("
                   << std::setw(5) << std::fixed << std::setprecision(2)
                   << 100.0 * (*lvl)->nonzeros() / sum_nnz
                   << "%)" << std::endl;

            os.flags(ff);
            os.precision(pp);
            return os;
        }

        /// Number of unknowns at the finest level.
        index_t size() const {
            return hier.front()->size();
        }

        /// Returns reference to the system matrix at the finest level.
        const typename level_type::matrix& top_matrix() const {
            return hier.front()->get_matrix();
        }
    private:
        void build_level(matrix &A, const params &prm, unsigned nlevel = 0)
        {
            if (A.rows <= prm.coarse_enough) {
                TIC("coarsest level");
                matrix Ai = sparse::inverse(A);
                hier.push_back( boost::shared_ptr<level_type>(new level_type(A, Ai, prm.level, nlevel) ) );
                TOC("coarsest level");
            } else {
                TIC("construct level");

                TIC("interp");
                BOOST_AUTO(PR, interp_t::interp(A, prm.interp));
                matrix &P = PR.first;
                matrix &R = PR.second;
                TOC("interp");

                TIC("coarse operator");
                matrix a = interp::coarse_operator<interp_t>::type::apply(
                        R, A, P, prm.interp);
                TOC("coarse operator");

                TOC("construct level");

                TIC("transfer level data");
                hier.push_back( boost::shared_ptr<level_type>(new level_type(A, P, R, prm.level, nlevel) ) );
                TOC("transfer level data");

                build_level(a, prm, nlevel + 1);
            }
        }

        params prm;
        std::list< boost::shared_ptr<level_type> > hier;
};

} // namespace amgcl

/// Output some general information about the AMG hierarchy.
template <
    typename value_t,
    typename index_t,
    typename interp_t,
    typename level_t
    >
std::ostream& operator<<(std::ostream &os, const amgcl::solver<value_t, index_t, interp_t, level_t> &amg) {
    return amg.print(os);
}

#endif
