#ifndef AMGCL_RELAXATION_DETAIL_ILU_SOLVE_HPP
#define AMGCL_RELAXATION_DETAIL_ILU_SOLVE_HPP

#include <amgcl/backend/interface.hpp>

namespace amgcl {
namespace relaxation {
namespace detail {

template <class Matrix, class VecD, class VecX>
void serial_ilu_solve(
        const Matrix &L, const Matrix &U, const VecD &D, VecX &x
        )
{
    const size_t n = backend::rows(L);

    for(size_t i = 0; i < n; i++) {
        for(ptrdiff_t j = L.ptr[i], e = L.ptr[i+1]; j < e; ++j)
            x[i] -= L.val[j] * x[L.col[j]];
    }

    for(size_t i = n; i-- > 0;) {
        for(ptrdiff_t j = U.ptr[i], e = U.ptr[i+1]; j < e; ++j)
            x[i] -= U.val[j] * x[U.col[j]];
        x[i] = D[i] * x[i];
    }
}

template <class Matrix, class VecD, class VecX, class VecT>
void parallel_ilu_solve(
        const Matrix &L, const Matrix &U, const VecD &D,
        VecX &x, VecT &t1, VecT &t2, unsigned jacobi_iters
        )
{
    typedef typename backend::value_type<Matrix>::type value_type;
    typedef typename math::scalar_of<value_type>::type scalar_type;

    VecT *y0 = &t1;
    VecT *y1 = &t2;

    backend::copy(x, *y0);
    for(unsigned i = 0; i < jacobi_iters; ++i) {
        backend::residual(x, L, *y0, *y1);
        std::swap(y0, y1);
    }

    backend::copy(*y0, x);
    for(unsigned i = 0; i < jacobi_iters; ++i) {
        backend::residual(*y0, U, x, *y1);
        backend::vmul(math::identity<scalar_type>(), D, *y1, math::zero<scalar_type>(), x);
    }
}

} // namespace detail
} // namespace relaxation
} // namespace amgcl

#endif
