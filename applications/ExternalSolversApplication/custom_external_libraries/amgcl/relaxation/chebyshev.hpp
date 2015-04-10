#ifndef AMGCL_RELAXATION_CHEBYSHEV_HPP
#define AMGCL_RELAXATION_CHEBYSHEV_HPP

/*
The MIT License

Copyright (c) 2012-2015 Denis Demidov <dennis.demidov@gmail.com>

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
 * \file   amgcl/relaxation/chebyshev.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Chebyshev polynomial smoother.
 */

#include <vector>
#include <cmath>
#include <boost/range/iterator_range.hpp>
#include <boost/foreach.hpp>
#include <boost/multi_array.hpp>
#include <boost/math/constants/constants.hpp>

#include <amgcl/detail/inverse.hpp>
#include <amgcl/util.hpp>

namespace amgcl {
namespace relaxation {

/// Chebyshev polynomial smoother.
/**
 * \param Backend Backend for temporary structures allocation.
 * \ingroup relaxation
 */
template <class Backend>
class chebyshev {
    public:
        typedef typename Backend::value_type value_type;
        typedef typename Backend::vector     vector;

        /// Relaxation parameters.
        struct params {
            /// Chebyshev polynomial degree.
            unsigned degree;

            /// Lowest-to-highest eigen value ratio.
            float lower;

            params() : degree(5), lower(1.0f / 30) {}

            params(const boost::property_tree::ptree &p)
                : AMGCL_PARAMS_IMPORT_VALUE(p, degree),
                  AMGCL_PARAMS_IMPORT_VALUE(p, lower)
            {}

            void get(boost::property_tree::ptree &p, const std::string &path) const {
                AMGCL_PARAMS_EXPORT_VALUE(p, path, degree);
                AMGCL_PARAMS_EXPORT_VALUE(p, path, lower);
            }
        };

        /// \copydoc amgcl::relaxation::damped_jacobi::damped_jacobi
        template <class Matrix>
        chebyshev(
                const Matrix &A, const params &prm,
                const typename Backend::params &backend_prm
            ) : C( prm.degree ),
                p( Backend::create_vector(rows(A), backend_prm) ),
                q( Backend::create_vector(rows(A), backend_prm) )
        {
            typedef value_type V;

            V hi = spectral_radius(A);
            V lo = hi * prm.lower;

            // Chebyshev polynomial roots on the interval [lo, hi].
            std::vector<value_type> roots(prm.degree);
            for(unsigned i = 0; i < prm.degree; ++i) {
                using boost::math::constants::pi;
                using boost::math::constants::half;

                roots[i] = lo + half<V>() * ( hi - lo ) * (
                        1 + cos( pi<V>() * ( i + half<V>() ) / prm.degree )
                        );
            }

            // Construct linear system to determine Chebyshev coefficients.
            boost::multi_array<V, 2> S(boost::extents[prm.degree][prm.degree]);
            std::vector<V> rhs(prm.degree);
            for(unsigned i = 0; i < prm.degree; ++i) {
                V x = roots[i];
                V x_to_j = 1;
                for(unsigned j = 0; j < prm.degree; ++j) {
                    S[i][j] = x_to_j;
                    x_to_j *= x;
                }
                rhs[i] = -x_to_j;
            }

            // Invert S, compute coefficients.
            amgcl::detail::inverse(prm.degree, S.data());

            value_type const_c = 1;
            for(unsigned i = 0; i < prm.degree; ++i) {
                value_type c = 0;
                for(unsigned j = 0; j < prm.degree; ++j)
                    c += S[i][j] * rhs[j];
                if (i == 0)
                    const_c = c;
                else
                    C[prm.degree - i] = -c / const_c;
            }
            C[0] = -1 / const_c;
        }

        /// \copydoc amgcl::relaxation::damped_jacobi::apply_pre
        template <class Matrix, class VectorRHS, class VectorX, class VectorTMP>
        void apply_pre(
                const Matrix &A, const VectorRHS &rhs, VectorX &x, VectorTMP &tmp,
                const params&
                ) const
        {
            apply(A, rhs, x, tmp);
        }

        /// \copydoc amgcl::relaxation::damped_jacobi::apply_post
        template <class Matrix, class VectorRHS, class VectorX, class VectorTMP>
        void apply_post(
                const Matrix &A, const VectorRHS &rhs, VectorX &x, VectorTMP &tmp,
                const params&
                ) const
        {
            apply(A, rhs, x, tmp);
        }

    private:
        std::vector<value_type> C;
        mutable boost::shared_ptr<vector> p, q;

        template <class Matrix, class VectorRHS, class VectorX, class VectorTMP>
        void apply(
                const Matrix &A, const VectorRHS &rhs, VectorX &x, VectorTMP &res
                ) const
        {
            backend::residual(rhs, A, x, res);
            backend::axpby(C[0], res, 0, *p);

            BOOST_FOREACH(value_type c, boost::make_iterator_range(C.begin() + 1, C.end()))
            {
                backend::spmv(1, A, *p, 0, *q);
                backend::axpbypcz(c, res, 1, *q, 0, *p);
            }

            backend::axpby(1, *p, 1, x);
        }

        template <class Matrix>
        static value_type spectral_radius(const Matrix &A) {
            typedef typename backend::row_iterator<Matrix>::type row_iterator;
            const size_t n = rows(A);

            value_type emax = 0;

#pragma omp parallel
            {
#ifdef _OPENMP
                int nt  = omp_get_num_threads();
                int tid = omp_get_thread_num();

                size_t chunk_size  = (n + nt - 1) / nt;
                size_t chunk_start = tid * chunk_size;
                size_t chunk_end   = std::min(n, chunk_start + chunk_size);
#else
                size_t chunk_start = 0;
                size_t chunk_end   = n;
#endif
                value_type my_emax = 0;
                for(size_t i = chunk_start; i < chunk_end; ++i) {
                    value_type hi = 0;

                    for(row_iterator a = backend::row_begin(A, i); a; ++a)
                        hi += std::fabs( a.value() );

                    my_emax = std::max(my_emax, hi);
                }

#pragma omp critical
                emax = std::max(emax, my_emax);
            }

            return emax;
        }
};

} // namespace relaxation
} // namespace amgcl



#endif
