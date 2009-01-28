/*
 *
 * Copyright (c) Kresimir Fresl 2003
 * Copyright (c) Georg Baum 2004
 *
 * Permission to copy, modify, use and distribute this software
 * for any non-commercial or commercial purpose is granted provided
 * that this license appear on all copies of the software source code.
 *
 * Author assumes no responsibility whatsoever for its use and makes
 * no guarantees about its quality, correctness or reliability.
 *
 * Author acknowledges the support of the Faculty of Civil Engineering,
 * University of Zagreb, Croatia.
 *
 */
/**
 * Adapted to row_major matrices by Janosch Stascheit
 * Institute for Structural Mechanics
 * Ruhr-University Bochum
 */

#ifndef BOOST_NUMERIC_BINDINGS_SUPERLU_H
#define BOOST_NUMERIC_BINDINGS_SUPERLU_H

// #define BOOST_NUMERIC_BINDINGS_SUPERLU_PRINT

#include <boost/numeric/bindings/traits/sparse_traits.hpp>
#include <boost/numeric/bindings/superlu/superlu_overloads.hpp>

#ifndef BOOST_NUMERIC_BINDINGS_NO_STRUCTURE_CHECK
#  include <boost/static_assert.hpp>
#  include <boost/type_traits/same_traits.hpp>
#endif

#ifdef BOOST_NUMERIC_BINDINGS_SUPERLU_PRINT
#  include <iostream>
#  include <algorithm>
#  include <iterator>
#endif

namespace boost { namespace numeric { namespace bindings {

  namespace superlu {

    enum permc_spec_t {
      natural = 0,            // Pc = I
      // .. is there a bug in SuperLU ?
      // .. get_perm_c.c & mmd.c ?
      // attma_min_degree = 1,  // min. degree on the structure of A^T * A
      atpla_min_degree = 2,   // min. degree on the structure of A^T + A
      colamd = 3              // approx. min. degree for unsymmetric matrices
    };
    

    template <typename MatrA, typename MatrB>
    inline
    int gssv (MatrA const& a, MatrB& b, permc_spec_t perm_spec = colamd) {
      // adapted to row_major

      typedef traits::sparse_matrix_traits<MatrA> matraits;
      typedef traits::matrix_traits<MatrB> mbtraits;
      std::cout << "I am here" << std::endl;

// #ifndef BOOST_NUMERIC_BINDINGS_NO_STRUCTURE_CHECK
      BOOST_STATIC_ASSERT((boost::is_same<
        typename matraits::matrix_structure, traits::general_t
      >::value));
      BOOST_STATIC_ASSERT((boost::is_same<
        typename matraits::storage_format, traits::compressed_t
      >::value));
      BOOST_STATIC_ASSERT((boost::is_same<
        typename matraits::ordering_type, traits::row_major_t
      >::value));

      BOOST_STATIC_ASSERT(matraits::index_base == 0);

      BOOST_STATIC_ASSERT((boost::is_same<
        typename mbtraits::matrix_structure, traits::general_t
      >::value));
      BOOST_STATIC_ASSERT((boost::is_same<
        typename mbtraits::ordering_type, traits::row_major_t
      >::value));
// #endif

      typedef typename matraits::value_type val_t;

      MatrA& aa = const_cast<MatrA&> (a);

      int m = matraits::size1 (aa);
      assert (m == matraits::size2 (aa));
      assert (m == mbtraits::size1 (b));

      SuperMatrix A, B, L, U;
      
      /** manual index vector generation
       */
      int *index1_vector = new (std::nothrow) int[a.index1_data().size()];
      int *index2_vector = new (std::nothrow) int[a.index2_data().size()];
      for( int i = 0; i < a.index1_data().size(); i++ )
      {
          index1_vector[i] = (int)a.index1_data()[i];
//           std::cout << index1_vector[i] << " ";
      }
      std::cout << std::endl;
      for( int i = 0; i < a.index2_data().size(); i++ )
      {
          index2_vector[i] = (int)a.index2_data()[i];
//           std::cout << index2_vector[i] << " ";
      }
      std::cout << std::endl;

      
      detail::Create_CompRow_Matrix (&A, m, m, matraits::num_nonzeros (aa),
                                     matraits::value_storage (aa),
                                     index2_vector,
                                     index1_vector
//                                      (int*)matraits::index1_storage (aa),
//                                      (int*)matraits::index2_storage (aa) 
                                    );
      
      std::cout << "Matrix storage type: " << A.Mtype << std::endl;
     
      detail::Create_Dense_Matrix (&B, mbtraits::size1(b), mbtraits::size2(b),
                                   mbtraits::storage (b),
//                                    mbtraits::leading_dimension (b));
                                   mbtraits::size1(b));

#ifdef BOOST_NUMERIC_BINDINGS_SUPERLU_PRINT
      detail::Print_CompRow_Matrix (val_t(), "A", &A);
      detail::Print_Dense_Matrix (val_t(), "B", &B);
      std::cout << std::endl;
#endif

      int *perm_r = new (std::nothrow) int[m];
      if (!perm_r)
        return -100;
      int *perm_c = new (std::nothrow) int[m];
      if (!perm_c) {
        delete[] perm_r;
        return -101;
      }
#ifdef BOOST_NUMERIC_BINDINGS_SUPERLU_PRINT
      std::cout << "perm_c: ";
      std::copy (perm_c, perm_c + m,
                 std::ostream_iterator<int> (std::cout, " "));
      std::cout << std::endl << std::endl;
#endif

      superlu_options_t options;
      set_default_options(&options);
      options.RowPerm = static_cast<rowperm_t>(perm_spec);
      /**
       * testing options
       */
//       options.ColPerm = static_cast<colperm_t>(colamd);
      
      
      SuperLUStat_t stat;
      StatInit(&stat);

      int info = 0;
      detail::gssv (val_t(), &options, &A, perm_c, perm_r, &L, &U, &B, &stat, &info);

#ifdef BOOST_NUMERIC_BINDINGS_SUPERLU_PRINT
      // detail::Print_SuperNode_Matrix (val_t(), "L + U - I", &L);
      // detail::Print_CompCol_Matrix (val_t(), "U", &U);

      // std::cout << std::endl;
      // SCformat* Lstore = (SCformat*) L.Store;
      // NCformat* Ustore = (NCformat*) U.Store;
      // std::cout << "No of nonzeros in L = " << Lstore->nnz << std::endl;
      // std::cout << "No of nonzeros in U = " << Ustore->nnz << std::endl;
      // std::cout << std::endl;
      // SCformat* Lstore = (SCformat*) L.Store;
      // for (int i = 0; i < Lstore->nnz; ++i)
      //   std::cout << ((val_t*) (Lstore->nzval))[i] << " ";
      // std::cout << std::endl << std::endl;

      detail::Print_Dense_Matrix (val_t(), "X", &B);
      std::cout << std::endl;
#endif

      if (options.PrintStat) {
        StatPrint(&stat);
      }
      StatFree(&stat);
      Destroy_CompCol_Matrix (&U);
      Destroy_SuperNode_Matrix (&L);
      delete[] perm_c;
      delete[] perm_r;
      Destroy_SuperMatrix_Store (&B);
      Destroy_SuperMatrix_Store (&A);
      delete[] index1_vector;
      delete[] index2_vector;
      
      return info;
    }

    /**
     * specialisation for single equation systems Ax=b
     * note: template parameter DataType is only for distinguation
     * with AX=B specialisation
     * by Janosch Stascheit
     */
    template <typename MatrA, typename VecB>
            inline
        int gssv (MatrA const& a, VecB& b, int single, permc_spec_t perm_spec = colamd) 
    {
        // adapted to row_major
        typedef traits::sparse_matrix_traits<MatrA> matraits;
        typedef traits::vector_traits<VecB> vtraits;

#ifndef BOOST_NUMERIC_BINDINGS_NO_STRUCTURE_CHECK
      BOOST_STATIC_ASSERT((boost::is_same<
        typename matraits::matrix_structure, traits::general_t
      >::value));
      BOOST_STATIC_ASSERT((boost::is_same<
              typename matraits::storage_format, traits::compressed_t
              >::value));
      BOOST_STATIC_ASSERT((boost::is_same<
              typename matraits::ordering_type, traits::row_major_t
              >::value));

      BOOST_STATIC_ASSERT(matraits::index_base == 0);

//       BOOST_STATIC_ASSERT((boost::is_same<
//               typename mbtraits::matrix_structure, traits::general_t
//               >::value));
//       BOOST_STATIC_ASSERT((boost::is_same<
//               typename mbtraits::ordering_type, traits::row_major_t
//               >::value));
#endif

      typedef typename matraits::value_type val_t;

      MatrA& aa = const_cast<MatrA&> (a);

      int m = matraits::size1 (aa);
      assert (m == matraits::size2 (aa));
      assert (m == vtraits::size (b));

      SuperMatrix A, B, L, U;

      detail::Create_CompRow_Matrix (&A, m, m, matraits::num_nonzeros (aa),
                                      matraits::value_storage (aa),
                                      (int*)matraits::index2_storage (aa),
                                      (int*)matraits::index1_storage (aa) );

      detail::Create_Dense_Matrix (&B, m, 1, 
                                    vtraits::storage(b),
                                    b.size()
//                                     vtraits::vector_storage(b),
//                                     vtraits::leading_dimension (b)
                                  );

#ifdef BOOST_NUMERIC_BINDINGS_SUPERLU_PRINT
      detail::Print_CompCol_Matrix (val_t(), "A", &A);
      detail::Print_Dense_Matrix (val_t(), "B", &B);
      std::cout << std::endl;
#endif

      int *perm_r = new (std::nothrow) int[m];
      if (!perm_r)
          return -100;
      int *perm_c = new (std::nothrow) int[m];
      if (!perm_c) {
          delete[] perm_r;
          return -101;
      }
#ifdef BOOST_NUMERIC_BINDINGS_SUPERLU_PRINT
      std::cout << "perm_c: ";
      std::copy (perm_c, perm_c + m,
                 std::ostream_iterator<int> (std::cout, " "));
      std::cout << std::endl << std::endl;
#endif

      superlu_options_t options;
      set_default_options(&options);
      options.RowPerm = static_cast<rowperm_t>(perm_spec);
      SuperLUStat_t stat;
      StatInit(&stat);

      int info = 0;
      detail::gssv (val_t(), &options, &A, perm_c, perm_r, &L, &U, &B, &stat, &info);

#ifdef BOOST_NUMERIC_BINDINGS_SUPERLU_PRINT
      // detail::Print_SuperNode_Matrix (val_t(), "L + U - I", &L);
      // detail::Print_CompCol_Matrix (val_t(), "U", &U);

      // std::cout << std::endl;
      // SCformat* Lstore = (SCformat*) L.Store;
      // NCformat* Ustore = (NCformat*) U.Store;
      // std::cout << "No of nonzeros in L = " << Lstore->nnz << std::endl;
      // std::cout << "No of nonzeros in U = " << Ustore->nnz << std::endl;
      // std::cout << std::endl;
      // SCformat* Lstore = (SCformat*) L.Store;
      // for (int i = 0; i < Lstore->nnz; ++i)
      //   std::cout << ((val_t*) (Lstore->nzval))[i] << " ";
      // std::cout << std::endl << std::endl;

//       detail::Print_Dense_Matrix (val_t(), "X", &B);
      std::cout << std::endl;
#endif

      if (options.PrintStat) {
    StatPrint(&stat);
      }
      StatFree(&stat);
      Destroy_CompCol_Matrix (&U);
      Destroy_SuperNode_Matrix (&L);
      delete[] perm_c;
      delete[] perm_r;
      Destroy_SuperMatrix_Store (&B);
      Destroy_SuperMatrix_Store (&A);

      return info;
                    }


    template <typename MatrA, typename IVec, typename MatrL, typename MatrU>
    inline
    int gstrf (MatrA const& a, IVec& perm_c, IVec& perm_r, MatrL& l, MatrU& u, permc_spec_t perm_spec = colamd) {
      // adapted row major orderings

      typedef traits::sparse_matrix_traits<MatrA> matraits;
      typedef traits::sparse_matrix_traits<MatrL> mltraits;
      typedef traits::sparse_matrix_traits<MatrU> mutraits;
      typedef traits::vector_traits<IVec> vtraits;

#ifndef BOOST_NUMERIC_BINDINGS_NO_STRUCTURE_CHECK
      BOOST_STATIC_ASSERT((boost::is_same<
        typename matraits::matrix_structure, traits::general_t
      >::value));
      BOOST_STATIC_ASSERT((boost::is_same<
        typename matraits::storage_format, traits::compressed_t
      >::value));
      BOOST_STATIC_ASSERT((boost::is_same<
        typename matraits::ordering_type, traits::row_major_t
      >::value));

      BOOST_STATIC_ASSERT(matraits::index_base == 0);

      BOOST_STATIC_ASSERT((boost::is_same<
        typename mltraits::matrix_structure, traits::general_t
      >::value));
      BOOST_STATIC_ASSERT((boost::is_same<
        typename mltraits::storage_format, traits::compressed_t
      >::value));
      BOOST_STATIC_ASSERT((boost::is_same<
        typename mltraits::ordering_type, traits::row_major_t
      >::value));

      BOOST_STATIC_ASSERT(mltraits::index_base == 0);

      BOOST_STATIC_ASSERT((boost::is_same<
        typename mutraits::matrix_structure, traits::general_t
      >::value));
      BOOST_STATIC_ASSERT((boost::is_same<
        typename mutraits::storage_format, traits::compressed_t
      >::value));
      BOOST_STATIC_ASSERT((boost::is_same<
        typename mutraits::ordering_type, traits::row_major_t
      >::value));

      BOOST_STATIC_ASSERT(mutraits::index_base == 0);

      BOOST_STATIC_ASSERT((boost::is_same<
        typename matraits::value_type, typename mltraits::value_type
      >::value));
      BOOST_STATIC_ASSERT((boost::is_same<
        typename matraits::value_type, typename mutraits::value_type
      >::value));
#endif

      typedef typename matraits::value_type val_t;

      int m = matraits::size1 (a);
      assert (m == matraits::size2 (a));

      SuperMatrix A, L, U;

      detail::Create_CompRow_Matrix (&A, m, m, matraits::num_nonzeros (a),
                                     matraits::value_storage (a),
                                     matraits::index2_storage (a),
                                     matraits::index1_storage (a));

      detail::Create_CompRow_Matrix (&L, m, m, mltraits::num_nonzeros (l),
                                     mltraits::value_storage (l),
                                     mltraits::index2_storage (l),
                                     mltraits::index1_storage (l));

      detail::Create_CompRow_Matrix (&U, m, m, mutraits::num_nonzeros (u),
                                     mutraits::value_storage (u),
                                     mutraits::index2_storage (u),
                                     mutraits::index1_storage (u));

      int* etree = new (std::nothrow) int[m];
      if (!etree)
        return -100;
#ifdef BOOST_NUMERIC_BINDINGS_SUPERLU_PRINT
      std::cout << "perm_c: ";
      std::copy (perm_c, perm_c + m,
                 std::ostream_iterator<int> (std::cout, " "));
      std::cout << std::endl << std::endl;
#endif

      superlu_options_t options;
      set_default_options(&options);
      options.RowPerm = static_cast<rowperm_t>(perm_spec);
      SuperLUStat_t stat;
      StatInit(&stat);

      int const panel_size = sp_ienv(1);
      int const relax = sp_ienv(2);
      val_t const drop_tol = val_t();	// not used

      int info = 0;
      detail::gstrf(val_t(), options, A, drop_tol, relax, panel_size, etree, 0, 0,
                    vtraits::vector_storage(perm_c), vtraits::pointer::vector_storage(perm_r), L, U, stat, info);

      if (options.PrintStat) {
        StatPrint(&stat);
      }
      StatFree(&stat);
      Destroy_CompCol_Matrix (&U);
      Destroy_SuperNode_Matrix (&L);
      delete[] etree;
      Destroy_SuperMatrix_Store (&A);

      return info;
    }
  }

}}}

#endif
