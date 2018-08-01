#ifndef AMGCL_RELAXATION_DAMPED_JACOBI_HPP
#define AMGCL_RELAXATION_DAMPED_JACOBI_HPP

/*
The MIT License

Copyright (c) 2012-2018 Denis Demidov <dennis.demidov@gmail.com>

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
 * \file   amgcl/relaxation/damped_jacobi.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Damped Jacobi relaxation scheme.
 */

#include <memory>
#include <amgcl/backend/interface.hpp>
#include <amgcl/util.hpp>

namespace amgcl {

/// Smoothers
namespace relaxation {

/**
 * \defgroup relaxation
 * \brief Relaxation schemes
 */

/// Damped Jacobi relaxation.
/**
 * \param Backend Backend for temporary structures allocation.
 * \ingroup relaxation
 */
template <class Backend>
struct damped_jacobi {
    typedef typename Backend::value_type               value_type;
    typedef typename math::scalar_of<value_type>::type scalar_type;

    /// Relaxation parameters.
    struct params {
        /// Damping factor.
        scalar_type damping;

        params(scalar_type damping = 0.72) : damping(damping) {}

#ifndef AMGCL_NO_BOOST
        params(const boost::property_tree::ptree &p)
            : AMGCL_PARAMS_IMPORT_VALUE(p, damping)
        {
            check_params(p, {"damping"});
        }

        void get(boost::property_tree::ptree &p, const std::string &path) const {
            AMGCL_PARAMS_EXPORT_VALUE(p, path, damping);
        }
#endif
    } prm;

    std::shared_ptr<typename Backend::matrix_diagonal> dia;

    /// Constructs smoother for the system matrix.
    /**
     * \param A           The system matrix.
     * \param prm         Relaxation parameters.
     * \param backend_prm Backend parameters.
     */
    template <class Matrix>
    damped_jacobi(
            const Matrix &A,
            const params &prm,
            const typename Backend::params &backend_prm
            )
        : prm(prm), dia( Backend::copy_vector( diagonal(A, true), backend_prm ) )
    { }

    /// Apply pre-relaxation
    /**
     * \param A   System matrix.
     * \param rhs Right-hand side.
     * \param x   Solution vector.
     * \param tmp Scratch vector.
     * \param prm Relaxation parameters.
     */
    template <class Matrix, class VectorRHS, class VectorX, class VectorTMP>
    void apply_pre(
            const Matrix &A, const VectorRHS &rhs, VectorX &x, VectorTMP &tmp
            ) const
    {
        backend::residual(rhs, A, x, tmp);
        backend::vmul(prm.damping, *dia, tmp, math::identity<scalar_type>(), x);
    }

    /// Apply post-relaxation
    /**
     * \param A   System matrix.
     * \param rhs Right-hand side.
     * \param x   Solution vector.
     * \param tmp Scratch vector.
     * \param prm Relaxation parameters.
     */
    template <class Matrix, class VectorRHS, class VectorX, class VectorTMP>
    void apply_post(
            const Matrix &A, const VectorRHS &rhs, VectorX &x, VectorTMP &tmp
            ) const
    {
        backend::residual(rhs, A, x, tmp);
        backend::vmul(prm.damping, *dia, tmp, math::identity<scalar_type>(), x);
    }

    template <class Matrix, class VectorRHS, class VectorX>
    void apply(const Matrix&, const VectorRHS &rhs, VectorX &x) const
    {
        backend::vmul(math::identity<scalar_type>(), *dia, rhs, math::zero<scalar_type>(), x);
    }
};

} // namespace relaxation

namespace backend {

template <class Backend>
struct bytes_impl< relaxation::damped_jacobi<Backend> > {
    static size_t get(const relaxation::damped_jacobi<Backend> &R) {
        return backend::bytes(*R.dia);
    }
};

} // namespace backend
} // namespace amgcl

#endif
