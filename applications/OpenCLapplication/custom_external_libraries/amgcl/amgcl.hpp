#ifndef AMGCL_AMGCL_HPP
#define AMGCL_AMGCL_HPP

/*
The MIT License

Copyright (c) 2012-2014 Denis Demidov <dennis.demidov@gmail.com>

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
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Generic algebraic multigrid framework.
 */

#include <iostream>
#include <iomanip>
#include <utility>
#include <list>

#include <boost/static_assert.hpp>
#include <boost/smart_ptr/shared_ptr.hpp>
#include <boost/smart_ptr/make_shared.hpp>
#include <boost/type_traits/is_signed.hpp>
#include <boost/io/ios_state.hpp>

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
            const Params&)
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
    ilu0,          ///< Incomplete LU decomposition with zero fill-in.
    chebyshev      ///< Chebyshev polynomial smoother.
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
            unsigned iter = 0;
            value_t  res  = 2 * prm.level.tol;

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
            boost::io::ios_all_saver stream_state(os);

            index_t sum_dof = 0;
            index_t sum_nnz = 0;
            for(typename std::list< boost::shared_ptr<level_type> >::const_iterator lvl = hier.begin(); lvl != hier.end(); ++lvl) {
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
            for(typename std::list< boost::shared_ptr<level_type> >::const_iterator lvl = hier.begin(); lvl != hier.end(); ++lvl, ++depth)
                os << std::setw(5)  << depth
                   << std::setw(13) << (*lvl)->size()
                   << std::setw(15) << (*lvl)->nonzeros() << " ("
                   << std::setw(5) << std::fixed << std::setprecision(2)
                   << 100.0 * (*lvl)->nonzeros() / sum_nnz
                   << "%)" << std::endl;

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
            if (static_cast<size_t>(A.rows) <= prm.coarse_enough) {
                TIC("coarsest level");
                matrix Ai = sparse::inverse(A);
                hier.push_back( boost::shared_ptr<level_type>(new level_type(A, Ai, prm.level, nlevel) ) );
                TOC("coarsest level");
            } else {
                TIC("construct level");

                TIC("interp");
                std::pair<sparse::matrix<value_t, index_t>, sparse::matrix<value_t, index_t> > PR = interp_t::interp(A, prm.interp);
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
