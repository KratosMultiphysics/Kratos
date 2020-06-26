#ifndef AMGCL_MPI_RELAXATION_SPAI0_HPP
#define AMGCL_MPI_RELAXATION_SPAI0_HPP

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
 * \file   amgcl/mpi/relaxation/spai0.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Distributed memory sparse approximate inverse relaxation scheme.
 */

#include <memory>
#include <amgcl/backend/interface.hpp>
#include <amgcl/value_type/interface.hpp>
#include <amgcl/mpi/util.hpp>
#include <amgcl/mpi/distributed_matrix.hpp>

namespace amgcl {
namespace mpi {
namespace relaxation {

template <class Backend>
struct spai0 {
    typedef Backend                                    backend_type;
    typedef typename Backend::value_type               value_type;
    typedef typename Backend::matrix_diagonal          matrix_diagonal;
    typedef typename math::scalar_of<value_type>::type scalar_type;
    typedef amgcl::detail::empty_params                params;
    typedef typename Backend::params                   backend_params;

    spai0(
            const distributed_matrix<Backend> &A,
            const params &, const backend_params &bprm = backend_params()
         )
    {
        typedef backend::crs<value_type> build_matrix;

        const ptrdiff_t n = A.loc_rows();
        const build_matrix &A_loc = *A.local();
        const build_matrix &A_rem = *A.remote();

        auto m = std::make_shared< backend::numa_vector<value_type> >(n, false);
        typedef backend::crs<value_type> build_matrix;

#pragma omp parallel for
        for(ptrdiff_t i = 0; i < n; ++i) {
            value_type  num = math::zero<value_type>();
            scalar_type den = math::zero<scalar_type>();

            for(ptrdiff_t j = A_loc.ptr[i], e = A_loc.ptr[i+1]; j < e; ++j) {
                value_type v = A_loc.val[j];
                scalar_type norm_v = math::norm(v);
                den += norm_v * norm_v;
                if (A_loc.col[j] == i) num += v;
            }

            for(ptrdiff_t j = A_rem.ptr[i], e = A_rem.ptr[i+1]; j < e; ++j) {
                value_type v = A_rem.val[j];
                scalar_type norm_v = math::norm(v);
                den += norm_v * norm_v;
            }

            (*m)[i] = math::inverse(den) * num;
        }

        M = Backend::copy_vector(m, bprm);
    }

    /// \copydoc amgcl::relaxation::damped_jacobi::apply_pre
    template <class Matrix, class VectorRHS, class VectorX, class VectorTMP>
    void apply_pre(
            const Matrix &A, const VectorRHS &rhs, VectorX &x, VectorTMP &tmp
            ) const
    {
        static const scalar_type one = math::identity<scalar_type>();
        backend::residual(rhs, A, x, tmp);
        backend::vmul(one, *M, tmp, one, x);
    }

    /// \copydoc amgcl::relaxation::damped_jacobi::apply_post
    template <class Matrix, class VectorRHS, class VectorX, class VectorTMP>
    void apply_post(
            const Matrix &A, const VectorRHS &rhs, VectorX &x, VectorTMP &tmp
            ) const
    {
        static const scalar_type one = math::identity<scalar_type>();
        backend::residual(rhs, A, x, tmp);
        backend::vmul(one, *M, tmp, one, x);
    }

    template <class Matrix, class VectorRHS, class VectorX>
    void apply( const Matrix&, const VectorRHS &rhs, VectorX &x) const
    {
        backend::vmul(math::identity<scalar_type>(), *M, rhs, math::zero<scalar_type>(), x);
    }

    private:
        std::shared_ptr<matrix_diagonal> M;
};

} // namespace relaxation
} // namespace mpi
} // namespace amgcl

#endif
