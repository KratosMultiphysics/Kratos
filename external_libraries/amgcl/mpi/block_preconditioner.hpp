#ifndef AMGCL_MPI_BLOCK_PRECONDITIONER_HPP
#define AMGCL_MPI_BLOCK_PRECONDITIONER_HPP

/*
The MIT License

Copyright (c) 2012-2017 Denis Demidov <dennis.demidov@gmail.com>

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
 * \file   amgcl/mpi/block_preconditioner.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Distributed block preconditioner.
 */

#include <vector>

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/range/numeric.hpp>

#include <mpi.h>

#include <amgcl/backend/builtin.hpp>
#include <amgcl/mpi/util.hpp>
#include <amgcl/mpi/inner_product.hpp>
#include <amgcl/mpi/distributed_matrix.hpp>

namespace amgcl {
namespace mpi {

template <class Precond>
class block_preconditioner {
    public:
        typedef typename Precond::params       params;
        typedef typename Precond::backend_type backend_type;
        typedef typename backend_type::params  backend_params;

        typedef typename backend_type::value_type value_type;
        typedef typename backend_type::matrix     bmatrix;
        typedef distributed_matrix<backend_type>  matrix;

        template <class Matrix>
        block_preconditioner(
                communicator comm,
                const Matrix &Astrip,
                const params &prm = params(),
                const backend_params &bprm = backend_params()
                )
        {
            typedef backend::crs<value_type> build_matrix;
            typedef typename backend::row_iterator<Matrix>::type row_iterator;

            ptrdiff_t n = backend::rows(Astrip);

            // Get sizes of each domain in comm.
            std::vector<ptrdiff_t> domain = mpi::exclusive_sum(comm, n);
            ptrdiff_t loc_beg = domain[comm.rank];
            ptrdiff_t loc_end = domain[comm.rank + 1];

            // Split the matrix into local and remote parts.
            boost::shared_ptr<build_matrix> Aloc = boost::make_shared<build_matrix>();
            boost::shared_ptr<build_matrix> Arem = boost::make_shared<build_matrix>();

            Aloc->set_size(n, n, true);
            Arem->set_size(n, 0, true);

#pragma omp parallel for
            for(ptrdiff_t i = 0; i < n; ++i) {
                for(row_iterator a = backend::row_begin(Astrip, i); a; ++a) {
                    ptrdiff_t c = a.col();

                    if (loc_beg <= c && c < loc_end)
                        ++Aloc->ptr[i + 1];
                    else
                        ++Arem->ptr[i + 1];
                }
            }

            std::partial_sum(Aloc->ptr, Aloc->ptr + n + 1, Aloc->ptr);
            std::partial_sum(Arem->ptr, Arem->ptr + n + 1, Arem->ptr);

            Aloc->set_nonzeros(Aloc->ptr[n]);
            Arem->set_nonzeros(Arem->ptr[n]);

#pragma omp parallel for
            for(ptrdiff_t i = 0; i < n; ++i) {
                ptrdiff_t loc_head = Aloc->ptr[i];
                ptrdiff_t rem_head = Arem->ptr[i];

                for(row_iterator a = backend::row_begin(Astrip, i); a; ++a) {
                    ptrdiff_t  c = a.col();
                    value_type v = a.value();

                    if (loc_beg <= c && c < loc_end) {
                        Aloc->col[loc_head] = c - loc_beg;
                        Aloc->val[loc_head] = v;
                        ++loc_head;
                    } else {
                        Arem->col[rem_head] = c;
                        Arem->val[rem_head] = v;
                        ++rem_head;
                    }
                }
            }

            C = boost::make_shared< comm_pattern<backend_type> >(comm, n, Arem->nnz, Arem->col, bprm);
            Arem->ncols = C->renumber(Arem->nnz, Arem->col);

            P = boost::make_shared<Precond>(Aloc, prm, bprm);

            this->Arem = backend_type::copy_matrix(Arem, bprm);
            this->A = boost::make_shared<matrix>(*C, P->system_matrix(), *this->Arem);
        }

        const matrix& system_matrix() const {
            return *A;
        }

        template <class Vec1, class Vec2>
        void apply(
                const Vec1 &rhs,
#ifdef BOOST_NO_CXX11_RVALUE_REFERENCES
                Vec2       &x
#else
                Vec2       &&x
#endif
                ) const
        {
            P->apply(rhs, x);
        }
    private:
        boost::shared_ptr< comm_pattern<backend_type> > C;
        boost::shared_ptr<bmatrix>  Arem;
        boost::shared_ptr<matrix> A;
        boost::shared_ptr<Precond> P;
};

} // namespace mpi
} // namespace amgcl

#endif
