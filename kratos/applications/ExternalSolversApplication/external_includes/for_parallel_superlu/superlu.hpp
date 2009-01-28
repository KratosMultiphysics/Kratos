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
 * Adapted to SuperLU_MT by Janosch Stascheit
 * Institute for Structural Mechanics
 * Ruhr-University Bochum
 */

#ifndef BOOST_NUMERIC_BINDINGS_SUPERLU_MT_H
#define BOOST_NUMERIC_BINDINGS_SUPERLU_MT_H

// #define BOOST_NUMERIC_BINDINGS_SUPERLU_PRINT

#include <boost/numeric/bindings/traits/sparse_traits.hpp>
#include "external_includes/superlu/superlu_mt_overloads.hpp"

#ifndef BOOST_NUMERIC_BINDINGS_NO_STRUCTURE_CHECK
#  include <boost/static_assert.hpp>
#  include <boost/type_traits/same_traits.hpp>
#endif

#ifdef BOOST_NUMERIC_BINDINGS_SUPERLU_PRINT
#  include <iostream>
#  include <algorithm>
#  include <iterator>
#endif

// #include "external_includes/superlu/slu_mt_util.h"
// #include "external_includes/superlu/pdsp_defs.h"
// #include "external_includes/superlu/slu_mt_Cnames.h"

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
    
    
    template <typename MatrA, typename MatrB> inline
            int pgssv (int nproc, MatrA const& a, MatrB& b, permc_spec_t perm_spec = colamd) 
    {
        typedef traits::sparse_matrix_traits<MatrA> matraits;
        typedef traits::matrix_traits<MatrB> mbtraits;
        std::cout << "I am in parallel pgssv, running with " << nproc << " threads" << std::endl;
        
        BOOST_STATIC_ASSERT((boost::is_same<
                typename matraits::matrix_structure, traits::general_t
                >::value));
        BOOST_STATIC_ASSERT((boost::is_same<
                typename matraits::storage_format, traits::compressed_t
                >::value));
        BOOST_STATIC_ASSERT((boost::is_same<
                typename matraits::ordering_type, traits::column_major_t
                >::value));

        BOOST_STATIC_ASSERT(matraits::index_base == 0);

        BOOST_STATIC_ASSERT((boost::is_same<
                typename mbtraits::matrix_structure, traits::general_t
                >::value));
        BOOST_STATIC_ASSERT((boost::is_same<
                typename mbtraits::ordering_type, traits::column_major_t
                >::value));
        
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
//         std::cout << std::endl;
//         std::cout << std::endl;
        for( int i = 0; i < a.index2_data().size(); i++ )
        {
            index2_vector[i] = (int)a.index2_data()[i];
//           std::cout << index2_vector[i] << " ";
        }
//         std::cout << std::endl;

        detail::Create_CompCol_Matrix (&A, m, m, matraits::num_nonzeros (aa),
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
        detail::Print_CompCol_Matrix (val_t(), "A", &A);
        detail::Print_Dense_Matrix (val_t(), "B", &B);
        std::cout << std::endl;
#endif

        int *perm_r = new (std::nothrow) int[m];
        if (!perm_r)
            return -100;
        int *perm_c = new (std::nothrow) int[m];
        if (!perm_c) {
//             std::cout << "####### deleting perm_r" << std::endl;
            delete[] perm_r;
            return -101;
        }
        
        /*
         * Get column permutation vector perm_c[], according to permc_spec:
         *   permc_spec = 0: natural ordering
         *   permc_spec = 1: minimum degree on structure of A’*A
         *   permc_spec = 2: minimum degree on structure of A’+A
         *   permc_spec = 3: approximate minimum degree for unsymmetric matrices
         */
        int permc_spec = 3;
        get_perm_c(permc_spec, &A, perm_c);
        

#ifdef BOOST_NUMERIC_BINDINGS_SUPERLU_PRINT
        std::cout << "perm_c: ";
        std::copy (perm_c, perm_c + m,
                   std::ostream_iterator<int> (std::cout, " "));
        std::cout << std::endl << std::endl;
#endif

//         superlu_mt_options_t options;
//         set_default_options(&options);
//         options.RowPerm = static_cast<rowperm_t>(perm_spec);
        /**
         * testing options
         */
//         options.ColPerm = static_cast<colperm_t>(colamd);
        
        int info = 0;
//         std::cout << "A.ncol: " << A.ncol << std::endl;
//         for( int i=0; i<4; i++ )
//         {
//             std::cout<< "####superlu.hpp#### running solver in iteration no. " << i << std::endl;
            boost::numeric::bindings::superlu::detail::pgssv( nproc, &A, perm_c, perm_r, &L, &U, &B, &info );
//         }
        
//         std::cout << "########## free memory..." << std::endl;
        //StatFree(&stat);
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
//     template <typename MatrA, typename VecB>
//             inline
//         int gssv (MatrA const& a, VecB& b, int single, permc_spec_t perm_spec = colamd) 
//     {
//         // adapted to row_major
//         typedef traits::sparse_matrix_traits<MatrA> matraits;
//         typedef traits::vector_traits<VecB> vtraits;
// 
// #ifndef BOOST_NUMERIC_BINDINGS_NO_STRUCTURE_CHECK
//       BOOST_STATIC_ASSERT((boost::is_same<
//         typename matraits::matrix_structure, traits::general_t
//       >::value));
//       BOOST_STATIC_ASSERT((boost::is_same<
//               typename matraits::storage_format, traits::compressed_t
//               >::value));
//       BOOST_STATIC_ASSERT((boost::is_same<
//               typename matraits::ordering_type, traits::row_major_t
//               >::value));
// 
//       BOOST_STATIC_ASSERT(matraits::index_base == 0);
// 
// //       BOOST_STATIC_ASSERT((boost::is_same<
// //               typename mbtraits::matrix_structure, traits::general_t
// //               >::value));
// //       BOOST_STATIC_ASSERT((boost::is_same<
// //               typename mbtraits::ordering_type, traits::row_major_t
// //               >::value));
// #endif
// 
//       typedef typename matraits::value_type val_t;
// 
//       MatrA& aa = const_cast<MatrA&> (a);
// 
//       int m = matraits::size1 (aa);
//       assert (m == matraits::size2 (aa));
//       assert (m == vtraits::size (b));
// 
//       SuperMatrix A, B, L, U;
// 
//       detail::Create_CompRow_Matrix (&A, m, m, matraits::num_nonzeros (aa),
//                                       matraits::value_storage (aa),
//                                       (int*)matraits::index2_storage (aa),
//                                       (int*)matraits::index1_storage (aa) );
// 
//       detail::Create_Dense_Matrix (&B, m, 1, 
//                                     vtraits::storage(b),
//                                     b.size()
// //                                     vtraits::vector_storage(b),
// //                                     vtraits::leading_dimension (b)
//                                   );
// 
// #ifdef BOOST_NUMERIC_BINDINGS_SUPERLU_PRINT
//       detail::Print_CompCol_Matrix (val_t(), "A", &A);
//       detail::Print_Dense_Matrix (val_t(), "B", &B);
//       std::cout << std::endl;
// #endif
// 
//       int *perm_r = new (std::nothrow) int[m];
//       if (!perm_r)
//           return -100;
//       int *perm_c = new (std::nothrow) int[m];
//       if (!perm_c) {
//           delete[] perm_r;
//           return -101;
//       }
// #ifdef BOOST_NUMERIC_BINDINGS_SUPERLU_PRINT
//       std::cout << "perm_c: ";
//       std::copy (perm_c, perm_c + m,
//                  std::ostream_iterator<int> (std::cout, " "));
//       std::cout << std::endl << std::endl;
// #endif
// 
//       superlu_options_t options;
//       set_default_options(&options);
//       options.RowPerm = static_cast<rowperm_t>(perm_spec);
//       SuperLUStat_t stat;
//       StatInit(&stat);
// 
//       int info = 0;
//       detail::gssv (val_t(), &options, &A, perm_c, perm_r, &L, &U, &B, &stat, &info);
// 
// #ifdef BOOST_NUMERIC_BINDINGS_SUPERLU_PRINT
//       // detail::Print_SuperNode_Matrix (val_t(), "L + U - I", &L);
//       // detail::Print_CompCol_Matrix (val_t(), "U", &U);
// 
//       // std::cout << std::endl;
//       // SCformat* Lstore = (SCformat*) L.Store;
//       // NCformat* Ustore = (NCformat*) U.Store;
//       // std::cout << "No of nonzeros in L = " << Lstore->nnz << std::endl;
//       // std::cout << "No of nonzeros in U = " << Ustore->nnz << std::endl;
//       // std::cout << std::endl;
//       // SCformat* Lstore = (SCformat*) L.Store;
//       // for (int i = 0; i < Lstore->nnz; ++i)
//       //   std::cout << ((val_t*) (Lstore->nzval))[i] << " ";
//       // std::cout << std::endl << std::endl;
// 
// //       detail::Print_Dense_Matrix (val_t(), "X", &B);
//       std::cout << std::endl;
// #endif
// 
//       if (options.PrintStat) {
//     StatPrint(&stat);
//       }
//       StatFree(&stat);
//       Destroy_CompCol_Matrix (&U);
//       Destroy_SuperNode_Matrix (&L);
//       delete[] perm_c;
//       delete[] perm_r;
//       Destroy_SuperMatrix_Store (&B);
//       Destroy_SuperMatrix_Store (&A);
// 
//       return info;
//                     }
// 
  }

}}}

#endif
