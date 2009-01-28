/*
 * 
 * Copyright (c) Toon Knapen, Karl Meerbergen & Kresimir Fresl 2003
 *
 * Permission to copy, modify, use and distribute this software 
 * for any non-commercial or commercial purpose is granted provided 
 * that this license appear on all copies of the software source code.
 *
 * Authors assume no responsibility whatsoever for its use and makes 
 * no guarantees about its quality, correctness or reliability.
 *
 * KF acknowledges the support of the Faculty of Civil Engineering, 
 * University of Zagreb, Croatia.
 *
 */

#ifndef BOOST_NUMERIC_BINDINGS_LAPACK_GEQRF_HPP
#define BOOST_NUMERIC_BINDINGS_LAPACK_GEQRF_HPP

#include <complex>
#include <boost/numeric/bindings/traits/traits.hpp>
#include <boost/numeric/bindings/lapack/lapack.h>
#include <boost/numeric/bindings/lapack/workspace.hpp>
#include <boost/numeric/bindings/traits/detail/array.hpp>
// #include <boost/numeric/bindings/traits/std_vector.hpp>

#ifndef BOOST_NUMERIC_BINDINGS_NO_STRUCTURE_CHECK 
#  include <boost/static_assert.hpp>
#  include <boost/type_traits.hpp>
#endif 


namespace boost { namespace numeric { namespace bindings { 

  namespace lapack {

    ///////////////////////////////////////////////////////////////////
    //
    // QR factorization of a general m x n matrix  A = Q * R
    // 
    ///////////////////////////////////////////////////////////////////

    /* 
     * geqrf() computes the QR factorization of a rectangular matrix
     * A = Q *  R,  where Q is a M x min(M,N) matrix with orthogonal
     * and normalized column (i.e. herm(Q) $ Q = I) and R is a
     * min(M,N) x N upper triangular matrix.
     *
     * On return of geqrf, the elements on and above the diagonal of
     * A contain the min(M,N) x N upper trapezoidal matrix R (R is
     * upper triangular if  m  >=  n);  the elements below the diagonal,
     * with the array TAU, represent the orthogonal matrix Q as a product
     * of min(M,N) elementary reflectors.
     */ 

    namespace detail {

      inline 
      void geqrf (int const m, int const n,
                 float* a, int const lda,
                 float* tau, float* work, int const lwork, int& info) 
      {
        LAPACK_SGEQRF (&m, &n, a, &lda, tau, work, &lwork, &info);
      }

      inline 
      void geqrf (int const m, int const n,
                 double* a, int const lda,
                 double* tau, double* work, int const lwork, int& info) 
      {
        LAPACK_DGEQRF (&m, &n, a, &lda, tau, work, &lwork, &info);
      }

      inline 
      void geqrf (int const m, int const n,
                  std::complex<float>* a, int const lda,
                  std::complex<float>* tau, std::complex<float>* work,
		  int const lwork, int& info) 
      {
        LAPACK_CGEQRF (&m, &n,
                      reinterpret_cast<fcomplex_t*> (a), &lda,
                      reinterpret_cast<fcomplex_t*> (tau),
                      reinterpret_cast<fcomplex_t*> (work), &lwork, &info );
      }
      

      inline 
      void geqrf (int const m, int const n,
                  std::complex<double>* a, int const lda,
                  std::complex<double>* tau, std::complex<double>* work,
		  int const lwork, int& info) 
      {
        LAPACK_ZGEQRF (&m, &n,
                      reinterpret_cast<dcomplex_t*> (a), &lda,
                      reinterpret_cast<dcomplex_t*> (tau),
                      reinterpret_cast<dcomplex_t*> (work), &lwork, &info );
      }
      
    } 

    template <typename A, typename Tau, typename Work>
    inline
    int geqrf (A& a, Tau& tau, Work& work) {

#ifndef BOOST_NUMERIC_BINDINGS_NO_STRUCTURE_CHECK 
      BOOST_STATIC_ASSERT((boost::is_same<
        typename traits::matrix_traits<A>::matrix_structure, 
        traits::general_t
      >::value)); 
#endif 

      int const m = traits::matrix_size1 (a);
      int const n = traits::matrix_size2 (a);
      assert (std::min(m,n) <= traits::vector_size (tau)); 
      assert (n <= traits::vector_size (work)); 

      int info; 
      detail::geqrf (m, n,
                     traits::matrix_storage (a), 
                     traits::leading_dimension (a),
                     traits::vector_storage (tau),  
                     traits::vector_storage (work),
                     traits::vector_size (work),
                     info);
      return info; 
    }

    // Computation of the QR factorization.
    // Workspace is allocated dynamically so that the optimization of
    // blas 3 calls is optimal.
    template <typename A, typename Tau>
    inline
    int geqrf (A& a, Tau& tau, optimal_workspace ) {
       typedef typename A::value_type value_type ;
       const int n = traits::matrix_size2 (a);
       traits::detail::array<value_type> work(std::max(1, n*32));
       return geqrf( a, tau, work );
    }

    // Computation of the QR factorization.
    // Workspace is allocated dynamically to its minimum size.
    // Blas 3 calls are not optimal.
    template <typename A, typename Tau>
    inline
    int geqrf (A& a, Tau& tau, minimal_workspace ) {
       typedef typename A::value_type value_type ;
       const int n = traits::matrix_size2 (a);
       traits::detail::array<value_type> work(std::max(1, n));
       return geqrf( a, tau, work );
    }

    // Computation of the QR factorization.
    // Workspace is taken from the array in workspace.
    // The calling sequence is
    // geqrf( a, tau, workspace( work ) ) where work is an array with the same value_type
    // as a.
    template <typename A, typename Tau, typename Work>
    inline
    int geqrf (A& a, Tau& tau, detail::workspace1<Work> workspace ) {
       typedef typename A::value_type value_type ;
       const int n = traits::matrix_size2 (a);
       traits::detail::array<value_type> work(std::max(1, n));
       return geqrf( a, tau, workspace.w_ );
    }

    // Function without workarray as argument
    template <typename A, typename Tau>
    inline
    int geqrf (A& a, Tau& tau) {
       return geqrf( a, tau, optimal_workspace() );
    }

  }

}}}

#endif 
