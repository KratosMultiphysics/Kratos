#ifndef AMGCL_MPI_RELAXATION_ILU0_CHOW_PATEL_HPP
#define AMGCL_MPI_RELAXATION_ILU0_CHOW_PATEL_HPP

/*
The MIT License

Copyright (c) 2026 Vicente Mataix Ferrandiz <vicente.mataix-ferrandiz@siemens.com>

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
 * \file   amgcl/mpi/relaxation/ilu0_chow_patel.hpp
 * \author Vicente Mataix Ferrandiz <vicente.mataix-ferrandiz@siemens.com>
 * \brief  Distributed memory fine-grained parallel ILU(0) (Chow-Patel).
 */

#include <amgcl/relaxation/ilu0_chow_patel.hpp>
#include <amgcl/mpi/distributed_matrix.hpp>

namespace amgcl {
namespace mpi {
namespace relaxation {

template <class Backend>
struct ilu0_chow_patel : public amgcl::relaxation::ilu0_chow_patel<Backend> {
    typedef Backend backend_type;
    typedef amgcl::relaxation::ilu0_chow_patel<Backend> Base;
    typedef typename Backend::params backend_params;
    typedef typename Base::params params;

    ilu0_chow_patel(
            const distributed_matrix<Backend> &A,
            const params &prm = params(),
            const backend_params &bprm = backend_params()
         ) : Base(*A.local(), prm, bprm)
    {}
};

} // namespace relaxation
} // namespace mpi
} // namespace amgcl

#endif
