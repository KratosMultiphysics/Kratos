#ifndef AMGCL_MAKE_BLOCK_SOLVER_HPP
#define AMGCL_MAKE_BLOCK_SOLVER_HPP

#include <amgcl/backend/interface.hpp>
#include <amgcl/adapter/block_matrix.hpp>
#include <amgcl/value_type/static_matrix.hpp>
#include <amgcl/make_solver.hpp>

namespace amgcl {

namespace backend {

/* Allows to do matrix-vector products with mixed scalar/nonscalar types.
 * Reinterprets pointers to the vectors data into appropriate types.
 */
template <class Alpha, class Matrix, class Vector1, class Beta, class Vector2>
struct spmv_impl<
    Alpha, Matrix, Vector1, Beta, Vector2,
    typename std::enable_if<
            detail::use_builtin_matrix_ops<Matrix>::value && (
            math::static_rows<typename value_type<Matrix>::type>::value != math::static_rows<typename value_type<Vector1>::type>::value ||
            math::static_rows<typename value_type<Matrix>::type>::value != math::static_rows<typename value_type<Vector2>::type>::value)
        >::type
    >
{
    static void apply(
            Alpha alpha, const Matrix &A, const Vector1 &x, Beta beta, Vector2 &y
            )
    {
        typedef typename value_type<Matrix>::type     val_type;
        typedef typename math::rhs_of<val_type>::type rhs_type;

        const size_t n = backend::rows(A);
        const size_t m = backend::cols(A);

        rhs_type const * xptr = reinterpret_cast<rhs_type const *>(&x[0]);
        rhs_type       * yptr = reinterpret_cast<rhs_type       *>(&y[0]);

        boost::iterator_range<rhs_type const *> xrng(xptr, xptr + m);
        boost::iterator_range<rhs_type       *> yrng(yptr, yptr + n);

        spmv(alpha, A, xrng, beta, yrng);
    }
};

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
        typedef typename math::rhs_of<value_type>::type    rhs_type;

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
            const size_t n = backend::rows(system_matrix());

            rhs_type const * fptr = reinterpret_cast<rhs_type const *>(&rhs[0]);
            rhs_type       * xptr = reinterpret_cast<rhs_type       *>(&x[0]);

            boost::iterator_range<rhs_type const *> frng(fptr, fptr + n);
            boost::iterator_range<rhs_type       *> xrng(xptr, xptr + n);

            return (*S)(A, frng, xrng);
        }

        template <class Vec1, class Vec2>
        std::tuple<size_t, scalar_type>
        operator()(const Vec1 &rhs, Vec2 &&x) const {
            const size_t n = backend::rows(system_matrix());

            rhs_type const * fptr = reinterpret_cast<rhs_type const *>(&rhs[0]);
            rhs_type       * xptr = reinterpret_cast<rhs_type       *>(&x[0]);

            boost::iterator_range<rhs_type const *> frng(fptr, fptr + n);
            boost::iterator_range<rhs_type       *> xrng(xptr, xptr + n);

            return (*S)(frng, xrng);
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
    private:
        typedef make_solver<Precond, IterativeSolver> Solver;
        std::shared_ptr<Solver> S;
};

} // namespace amgcl

#endif
