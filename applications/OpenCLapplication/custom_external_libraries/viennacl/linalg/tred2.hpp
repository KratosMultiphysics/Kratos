#ifndef VIENNACL_LINALG_TRED2_HPP_
#define VIENNACL_LINALG_TRED2_HPP_

/* =========================================================================
   Copyright (c) 2010-2014, Institute for Microelectronics,
                            Institute for Analysis and Scientific Computing,
                            TU Wien.
   Portions of this software are copyright by UChicago Argonne, LLC.

                            -----------------
                  ViennaCL - The Vienna Computing Library
                            -----------------

   Project Head:    Karl Rupp                   rupp@iue.tuwien.ac.at

   (A list of authors and contributors can be found in the PDF manual)

   License:         MIT (X11), see file LICENSE in the base directory
============================================================================= */

/** @file viennacl/linalg/tred2.hpp
*   @brief optimized CPU tridiagonalization
*
*   Contributed by Alex Christensen.
*/

#include "viennacl/linalg/host_based/sse_kernels.hpp"
#include <boost/numeric/ublas/matrix.hpp>

namespace viennacl
{
  namespace linalg
  {
    /** @brief Inplace reduction of a hermitian (or real symmetric) to tridiagonal form using householder similarity transforms (preserving eigenvalues)
     *
     * @param A            A dense matrix to be tridiagonalized
     * @param block_size   The block size to be used
     */
    template<typename ScalarType>
    void inplace_tred2(boost::numeric::ublas::matrix<ScalarType> const & A, vcl_size_t block_size = 1)
    {

#ifdef VIENNACL_WITH_OPENMP
      vcl_size_t num_threads=(vcl_size_t)omp_get_max_threads();
      omp_set_num_threads(omp_get_max_threads());
#else
      vcl_size_t num_threads=1;
#endif

      vcl_size_t n=A.size1();
      if(n!=A.size2())
        std::cerr << "ViennaCL: Warning in inplace_tred2(): Matrix is not hermitian (or real symmetric)" << std::endl;
      block_size=std::min(n,block_size);

      //get pointers to the value arrays
      ScalarType** rows=new ScalarType*[n];
      for(vcl_size_t i=0;i<n;i++)
        rows[i]=(ScalarType*)&A(i,0);

      //call the optimized CPU code
      inplace_tred2(rows,A.size1(),block_size,num_threads);

      delete [] rows;
    }

  } //namespace linalg
} //namespace viennacl
#endif
