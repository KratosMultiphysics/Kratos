#ifndef AMGCL_MPI_COARSENING_AGGREGATION_HPP
#define AMGCL_MPI_COARSENING_AGGREGATION_HPP

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
 * \file   amgcl/mpi/coarsening/aggregation.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Distributed non-smoothed aggregation coarsening scheme.
 */

#include <tuple>
#include <memory>

#include <amgcl/backend/builtin.hpp>
#include <amgcl/util.hpp>
#include <amgcl/coarsening/detail/scaled_galerkin.hpp>
#include <amgcl/mpi/util.hpp>
#include <amgcl/mpi/distributed_matrix.hpp>
#include <amgcl/mpi/coarsening/pmis.hpp>

namespace amgcl {
namespace mpi {
namespace coarsening {

template <class Backend>
struct aggregation {
    typedef typename Backend::value_type value_type;
    typedef typename math::scalar_of<value_type>::type scalar_type;
    typedef backend::crs<value_type> build_matrix;

    struct params {
        // aggregation params
        typedef typename pmis<Backend>::params aggr_params;
        aggr_params aggr;

        /// Over-interpolation factor \f$\alpha\f$.
        /**
         * In case of aggregation coarsening, coarse-grid
         * correction of smooth error, and by this the overall convergence, can
         * often be substantially improved by using "over-interpolation", that is,
         * by multiplying the actual correction (corresponding to piecewise
         * constant interpolation) by some factor \f$\alpha > 1\f$. Equivalently,
         * this means that the coarse-level Galerkin operator is re-scaled by
         * \f$1 / \alpha\f$:
         * \f[I_h^HA_hI_H^h \to \frac{1}{\alpha}I_h^HA_hI_H^h.\f]
         *
         * \sa  \cite Stuben1999, Section 9.1 "Re-scaling of the Galerkin operator".
         */
        float over_interp;

        params() : over_interp(1.5f) { }

#ifndef AMGCL_NO_BOOST
        params(const boost::property_tree::ptree &p)
            : AMGCL_PARAMS_IMPORT_CHILD(p, aggr),
              AMGCL_PARAMS_IMPORT_VALUE(p, over_interp)
        {
            check_params(p, {"aggr", "over_interp"});
        }

        void get(boost::property_tree::ptree &p, const std::string &path) const {
            AMGCL_PARAMS_EXPORT_CHILD(p, path, aggr);
            AMGCL_PARAMS_EXPORT_VALUE(p, path, over_interp);
        }
#endif
    } prm;

    aggregation(const params &prm = params()) : prm(prm) {}

    std::tuple<
        std::shared_ptr< distributed_matrix<Backend> >,
        std::shared_ptr< distributed_matrix<Backend> >
        >
    transfer_operators(const distributed_matrix<Backend> &A) {
        pmis<Backend> aggr(A, prm.aggr);
        return std::make_tuple(aggr.p_tent, transpose(*aggr.p_tent));
    }

    std::shared_ptr< distributed_matrix<Backend> >
    coarse_operator(
            const distributed_matrix<Backend> &A,
            const distributed_matrix<Backend> &P,
            const distributed_matrix<Backend> &R
            ) const
    {
        return amgcl::coarsening::detail::scaled_galerkin(A, P, R, 1 / prm.over_interp);
    }

};

template <class Backend>
unsigned block_size(const aggregation<Backend> &c) {
    return c.prm.aggr.block_size;
}

} // namespace coarsening
} // namespace mpi
} // namespace amgcl


#endif
