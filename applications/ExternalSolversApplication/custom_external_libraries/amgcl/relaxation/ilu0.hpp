#ifndef AMGCL_RELAXATION_ILU0_HPP
#define AMGCL_RELAXATION_ILU0_HPP

/**
 * \file   amgcl/relaxation/ilu0.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Incomplete LU relaxation scheme.
 */

#include <boost/shared_ptr.hpp>
#include <amgcl/backend/interface.hpp>
#include <amgcl/util.hpp>

namespace amgcl {
namespace relaxation {

/// ILU(0) smoother.
/**
 * \note ILU(0) is a serial algorithm and is only applicable to backends that
 * support matrix row iteration (e.g. amgcl::backend::builtin or
 * amgcl::backend::eigen).
 *
 * \param Backend Backend for temporary structures allocation.
 * \ingroup relaxation
 */
template <class Backend>
struct ilu0 {
    typedef typename Backend::value_type value_type;
    typedef typename Backend::vector     vector;

    /// Relaxation parameters.
    struct params {
        /// Damping factor.
        float damping;

        params(float damping = 0.72) : damping(damping) {}

        params(const boost::property_tree::ptree &p)
            : AMGCL_PARAMS_IMPORT_VALUE(p, damping)
        {}

        void get(boost::property_tree::ptree &p, const std::string &path) const {
            AMGCL_PARAMS_EXPORT_VALUE(p, path, damping);
        }
    };

    /// \copydoc amgcl::relaxation::damped_jacobi::damped_jacobi
    template <class Matrix>
    ilu0( const Matrix &A, const params &, const typename Backend::params&)
        : luval( A.val ),
          dia  ( backend::rows(A) )
    {
        const size_t n = backend::rows(A);
        const value_type eps = amgcl::detail::eps<value_type>(1);

        std::vector<ptrdiff_t> work(n, -1);

        for(size_t i = 0; i < n; ++i) {
            ptrdiff_t row_beg = A.ptr[i];
            ptrdiff_t row_end = A.ptr[i + 1];

            for(ptrdiff_t j = row_beg; j < row_end; ++j)
                work[ A.col[j] ] = j;

            for(ptrdiff_t j = row_beg; j < row_end; ++j) {
                ptrdiff_t c = A.col[j];

                // Exit if diagonal is reached
                if (static_cast<size_t>(c) >= i) {
                    precondition(
                            static_cast<size_t>(c) == i,
                            "No diagonal value in system matrix"
                            );
                    precondition(
                            fabs(luval[j]) > eps,
                            "Zero pivot in ILU"
                            );

                    dia[i]   = j;
                    luval[j] = 1 / luval[j];
                    break;
                }

                // Compute the multiplier for jrow
                value_type tl = luval[j] * luval[dia[c]];
                luval[j] = tl;

                // Perform linear combination
                for(ptrdiff_t k = dia[c] + 1; k < A.ptr[c + 1]; ++k) {
                    ptrdiff_t w = work[A.col[k]];
                    if (w >= 0) luval[w] -= tl * luval[k];
                }
            }

            // Refresh work
            for(ptrdiff_t j = row_beg; j < row_end; ++j)
                work[A.col[j]] = -1;
        }
    }

    /// \copydoc amgcl::relaxation::damped_jacobi::apply_pre
    template <class Matrix, class VectorRHS, class VectorX, class VectorTMP>
    void apply_pre(
            const Matrix &A, const VectorRHS &rhs, VectorX &x, VectorTMP &tmp,
            const params &prm
            ) const
    {
        apply(A, rhs, x, tmp, prm);
    }

    /// \copydoc amgcl::relaxation::damped_jacobi::apply_post
    template <class Matrix, class VectorRHS, class VectorX, class VectorTMP>
    void apply_post(
            const Matrix &A, const VectorRHS &rhs, VectorX &x, VectorTMP &tmp,
            const params &prm
            ) const
    {
        apply(A, rhs, x, tmp, prm);
    }

    private:
        std::vector<value_type> luval;
        std::vector<ptrdiff_t>  dia;

        template <class Matrix, class VectorRHS, class VectorX, class VectorTMP>
        void apply(
                const Matrix &A, const VectorRHS &rhs, VectorX &x, VectorTMP &tmp,
                const params &prm
                ) const
        {
            const size_t n = backend::rows(A);

#pragma omp parallel for
            for(ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(n); ++i) {
                value_type buf = rhs[i];
                for(ptrdiff_t j = A.ptr[i], e = A.ptr[i + 1]; j < e; ++j)
                    buf -= A.val[j] * x[A.col[j]];
                tmp[i] = buf;
            }

            for(size_t i = 0; i < n; i++) {
                for(ptrdiff_t j = A.ptr[i], e = dia[i]; j < e; ++j)
                    tmp[i] -= luval[j] * tmp[A.col[j]];
            }

            for(size_t i = n; i-- > 0;) {
                for(ptrdiff_t j = dia[i] + 1, e = A.ptr[i + 1]; j < e; ++j)
                    tmp[i] -= luval[j] * tmp[A.col[j]];
                tmp[i] *= luval[dia[i]];
            }

#pragma omp parallel for
            for(ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(n); ++i)
                x[i] += prm.damping * tmp[i];
        }

};

} // namespace relaxation

namespace backend {

template <class Backend>
struct relaxation_is_supported<
    Backend,
    relaxation::ilu0,
    typename boost::disable_if<
            typename boost::is_same<
                Backend,
                builtin<typename Backend::value_type>
            >::type
        >::type
    > : boost::false_type
{};

} // namespace backend
} // namespace amgcl



#endif
