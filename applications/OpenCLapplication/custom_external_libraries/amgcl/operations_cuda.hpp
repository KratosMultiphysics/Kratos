#ifndef AMGCL_OPERATIONS_CUDA_HPP
#define AMGCL_OPERATIONS_CUDA_HPP

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
 * \file   operations_cuda.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Adaptors for thrust types.
 */

#include <thrust/device_vector.h>
#include <amgcl/common.hpp>

namespace amgcl {

/// Returns value type for a thrust vector.
template <typename T>
struct value_type< thrust::device_vector<T> > {
    typedef T type;
};

/// Clears (sets elements to zero) thrust vector.
template <typename T>
void clear(thrust::device_vector<T> &x) {
    thrust::fill(x.begin(), x.end(), T());
}

/// Returns inner product of two thrust vectors.
template <typename T>
T inner_prod(const thrust::device_vector<T> &x, const thrust::device_vector<T> &y) {
    return thrust::inner_product(x.begin(), x.end(), y.begin(), T());
}

/// Returns norm of a thrust vector.
template <typename T>
T norm(const thrust::device_vector<T> &x) {
    return sqrt(thrust::inner_product(x.begin(), x.end(), x.begin(), T()));
}

/// Specialization of residual operation for thrust types.
template <class matrix, typename real>
void residual(
        const matrix &A,
        const thrust::device_vector<real> &x,
        const thrust::device_vector<real> &f,
        thrust::device_vector<real> &y
        )
{
    thrust::copy(f.begin(), f.end(), y.begin());;
    A.mul(-1, x, 1, y);
}

/// Specialization of matrix-vector product for thrust types.
template <class matrix, typename real>
void axpy(
        const matrix &A,
        const thrust::device_vector<real> &x,
        thrust::device_vector<real> &y
        )
{
    A.mul(1, x, 0, y);
}

namespace cuda_ops {

template <typename T>
struct ax {
    T a;

    ax(T a) : a(a) {}

    template <class Tuple>
    __device__ __host__ void operator()(Tuple t) const {
        thrust::get<0>(t) = a * thrust::get<1>(t);
    }
};

template <typename T>
struct axpby {
    T a;
    T b;

    axpby(T a, T b) : a(a), b(b) {}

    template <class Tuple>
    __device__ __host__ void operator()(Tuple t) const {
        thrust::get<0>(t) =
            a * thrust::get<1>(t) +
            b * thrust::get<2>(t);
    }
};

template <typename T>
struct axpbypcz {
    T a;
    T b;
    T c;

    axpbypcz(T a, T b, T c) : a(a), b(b), c(c) {}

    template <class Tuple>
    __device__ __host__ void operator()(Tuple t) const {
        thrust::get<0>(t) =
            a * thrust::get<1>(t) +
            b * thrust::get<2>(t) +
            c * thrust::get<3>(t);
    }
};

} // namespace cuda_ops

namespace sparse {
template <typename value_type, gpu_matrix_format> class cuda_matrix;
}

namespace gmres_ops {

template <class GMRES>
void update(GMRES &gmres, thrust::device_vector<typename GMRES::value_t> &x, int k) {
    std::copy(gmres.s.begin(), gmres.s.end(), gmres.y.begin());

    for (int i = k; i >= 0; --i) {
        gmres.y[i] /= gmres.H[i * gmres.M + i];
        for (int j = i - 1; j >= 0; --j)
            gmres.y[j] -= gmres.H[j * gmres.M + i] * gmres.y[i];
    }

    for (int j = 0; j <= k; j++)
        thrust::for_each(
                thrust::make_zip_iterator(thrust::make_tuple(
                        x.begin(), x.begin(), gmres.v[j].begin())),
                thrust::make_zip_iterator(thrust::make_tuple(
                        x.end(), x.end(), gmres.v[j].end())),
                cuda_ops::axpby<typename GMRES::value_t>(1, gmres.y[j])
                );
}

template <class GMRES, class matrix, class precond>
typename GMRES::value_t restart(GMRES &gmres,
        const matrix &A, const thrust::device_vector<typename GMRES::value_t> &rhs,
        const precond &P, const thrust::device_vector<typename GMRES::value_t> &x
        )
{
    typedef typename GMRES::value_t value_t;

    residual(A, x, rhs, gmres.w);
    clear(gmres.r);
    P.apply(gmres.w, gmres.r);

    gmres.s[0] = norm(gmres.r);

    thrust::for_each(
            thrust::make_zip_iterator(thrust::make_tuple(
                    gmres.v[0].begin(), gmres.r.begin())),
            thrust::make_zip_iterator(thrust::make_tuple(
                    gmres.v[0].end(), gmres.r.end())),
            cuda_ops::ax<value_t>(1 / gmres.s[0])
            );

    std::fill(gmres.s.begin() + 1, gmres.s.end(), value_t());

    return gmres.s[0];
}

template <class GMRES, typename value_t, gpu_matrix_format format, class precond>
value_t iteration(GMRES &gmres, const sparse::cuda_matrix<value_t, format> &A, const precond &P, int i)
{
    A.mul(1, gmres.v[i], 0, gmres.r);
    clear(gmres.w);
    P.apply(gmres.r, gmres.w);

    for(int k = 0; k <= i; ++k) {
        gmres.H[k * gmres.M + i] = inner_prod(gmres.w, gmres.v[k]);
        thrust::for_each(
                thrust::make_zip_iterator(thrust::make_tuple(
                    gmres.w.begin(), gmres.w.begin(), gmres.v[k].begin())),
                thrust::make_zip_iterator(thrust::make_tuple(
                    gmres.w.end(), gmres.w.end(), gmres.v[k].end())),
                cuda_ops::axpby<value_t>(1, -gmres.H[k * gmres.M + i])
                );
    }

    gmres.H[(i+1) * gmres.M + i] = norm(gmres.w);

    thrust::for_each(
            thrust::make_zip_iterator(thrust::make_tuple(
                    gmres.v[i+1].begin(), gmres.w.begin())),
            thrust::make_zip_iterator(thrust::make_tuple(
                    gmres.v[i+1].end(), gmres.w.end())),
            cuda_ops::ax<value_t>(1 / gmres.H[(i+1) * gmres.M + i])
            );

    for(int k = 0; k < i; ++k)
        apply_plane_rotation(gmres.H[k * gmres.M + i], gmres.H[(k+1) * gmres.M + i], gmres.cs[k], gmres.sn[k]);

    generate_plane_rotation(gmres.H[i * gmres.M + i], gmres.H[(i+1) * gmres.M + i], gmres.cs[i], gmres.sn[i]);
    apply_plane_rotation(gmres.H[i * gmres.M + i], gmres.H[(i+1) * gmres.M + i], gmres.cs[i], gmres.sn[i]);
    apply_plane_rotation(gmres.s[i], gmres.s[i+1], gmres.cs[i], gmres.sn[i]);

    return fabs(gmres.s[i + 1]);
}

} // namespace gmres_ops

} // namespace amgcl

#endif
