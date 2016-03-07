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

template <class Matrix, class VecD, class VecX>
void parallel_ilu_solve(
        const Matrix &L, const Matrix &U, const VecD &D,
        VecX &x, VecX &t1, VecX &t2, unsigned jacobi_iters
        )
{
    typedef typename backend::value_type<Matrix>::type value_type;
    typedef typename math::scalar_of<value_type>::type scalar_type;

    VecX *b  = &x;
    VecX *y0 = &t1;
    VecX *y1 = &t2;

    backend::copy(*b, *y0);
    for(unsigned i = 0; i < jacobi_iters; ++i) {
        backend::residual(*b, L, *y0, *y1);
        std::swap(y0, y1);
    }

    b  = y0;
    y0 = &x;
    backend::copy(*b, *y0);
    for(unsigned i = 0; i < jacobi_iters; ++i) {
        backend::residual(*b, U, *y0, *y1);
        backend::vmul(math::identity<scalar_type>(), D, *y1, math::zero<scalar_type>(), *y0);
    }
}

} // namespace detail
} // namespace relaxation
} // namespace amgcl

#endif
