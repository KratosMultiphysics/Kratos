#ifndef AMGCL_MAKE_BLOCK_SOLVER_HPP
#define AMGCL_MAKE_BLOCK_SOLVER_HPP

#include <amgcl/backend/interface.hpp>
#include <amgcl/adapter/block_matrix.hpp>
#include <amgcl/value_type/static_matrix.hpp>
#include <amgcl/make_solver.hpp>
#include <amgcl/util.hpp>

namespace amgcl {

namespace backend {

} // namespace backend

/* Creates solver that operates in non-scalar domain but may take scalar inputs
 * for the system matrix and the rhs/solution vectors.
 */
template <class Precond, class IterativeSolver>
class make_block_solver {
    public:
        typedef typename Precond::backend_type             backend_type;
        typedef typename backend_type::value_type          value_type;
        typedef typename backend_type::params              backend_params;
        typedef typename backend_type::vector              vector;
        typedef typename math::scalar_of<value_type>::type scalar_type;

        typedef typename make_solver<Precond, IterativeSolver>::params params;

        template <class Matrix>
        make_block_solver(
                const Matrix &A,
                const params &prm = params(),
                const backend_params &bprm = backend_params()
                )
        {
            S = std::make_shared<Solver>(adapter::block_matrix<value_type>(A), prm, bprm);
        }

        template <class Matrix, class Vec1, class Vec2>
        std::tuple<size_t, scalar_type> operator()(
                const Matrix &A, const Vec1 &rhs, Vec2 &&x) const
        {
            auto F = backend::reinterpret_as_rhs<value_type>(rhs);
            auto X = backend::reinterpret_as_rhs<value_type>(x);

            return (*S)(A, F, X);
        }

        template <class Vec1, class Vec2>
        std::tuple<size_t, scalar_type>
        operator()(const Vec1 &rhs, Vec2 &&x) const {
            auto F = backend::reinterpret_as_rhs<value_type>(rhs);
            auto X = backend::reinterpret_as_rhs<value_type>(x);

            return (*S)(F, X);
        }

        std::shared_ptr<typename Precond::matrix> system_matrix_ptr() const {
            return S->system_matrix_ptr();
        }

        typename Precond::matrix const& system_matrix() const {
            return S->system_matrix();
        }

        friend std::ostream& operator<<(std::ostream &os, const make_block_solver &p) {
            return os << *p.S << std::endl;
        }

        size_t bytes() const {
            return backend::bytes(*S);
        }
    private:
        typedef make_solver<Precond, IterativeSolver> Solver;
        std::shared_ptr<Solver> S;
};

} // namespace amgcl

#endif
