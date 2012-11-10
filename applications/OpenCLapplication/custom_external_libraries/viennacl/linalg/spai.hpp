#ifndef VIENNACL_LINALG_SPAI_HPP
#define VIENNACL_LINALG_SPAI_HPP

/* =========================================================================
   Copyright (c) 2010-2012, Institute for Microelectronics,
                            Institute for Analysis and Scientific Computing,
                            TU Wien.

                            -----------------
                  ViennaCL - The Vienna Computing Library
                            -----------------

   Project Head:    Karl Rupp                   rupp@iue.tuwien.ac.at
               
   (A list of authors and contributors can be found in the PDF manual)

   License:         MIT (X11), see file LICENSE in the base directory
============================================================================= */

/** @file viennacl/linalg/spai.hpp
    @brief Main include file for the sparse approximate inverse preconditioner family (SPAI and FSPAI).  Experimental in 1.2.x.
    
    Most implementation contributed by Nikolay Lukash.
*/


#include <utility>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <vector>
#include <math.h>
#include <map>

//local includes
#include "viennacl/linalg/detail/spai/spai_tag.hpp"
#include "viennacl/linalg/qr.hpp"
#include "viennacl/linalg/detail/spai/spai-dynamic.hpp"
#include "viennacl/linalg/detail/spai/spai-static.hpp"
#include "viennacl/linalg/detail/spai/sparse_vector.hpp"
#include "viennacl/linalg/detail/spai/block_matrix.hpp"
#include "viennacl/linalg/detail/spai/block_vector.hpp"
#include "viennacl/linalg/detail/spai/fspai.hpp"
#include "viennacl/linalg/detail/spai/spai.hpp"

//boost includes
#include "boost/numeric/ublas/vector.hpp"
#include "boost/numeric/ublas/matrix.hpp"
#include "boost/numeric/ublas/matrix_proxy.hpp"
#include "boost/numeric/ublas/vector_proxy.hpp"
#include "boost/numeric/ublas/storage.hpp"
#include "boost/numeric/ublas/io.hpp"
#include "boost/numeric/ublas/lu.hpp"
#include "boost/numeric/ublas/triangular.hpp"
#include "boost/numeric/ublas/matrix_expression.hpp"

// ViennaCL includes
#include "viennacl/linalg/prod.hpp"
#include "viennacl/matrix.hpp"
#include "viennacl/compressed_matrix.hpp"
#include "viennacl/linalg/compressed_matrix_operations.hpp"
#include "viennacl/linalg/matrix_operations.hpp"
#include "viennacl/scalar.hpp"
#include "viennacl/linalg/inner_prod.hpp"
#include "viennacl/linalg/ilu.hpp"
#include "viennacl/ocl/backend.hpp"
#include "viennacl/linalg/kernels/spai_source.h"
#include "viennacl/linalg/kernels/spai_kernels.h"


namespace viennacl
{
    namespace linalg
    {
        
        typedef viennacl::linalg::detail::spai::spai_tag         spai_tag;
        typedef viennacl::linalg::detail::spai::fspai_tag        fspai_tag;
        
        /** @brief Implementation of the SParse Approximate Inverse Algorithm
         * @param Matrix matrix that is used for computations
         * @param Vector vector that is used for computations
         */
        //UBLAS version
        template <typename MatrixType>
        class spai_precond
        {
        public:
            typedef typename MatrixType::value_type ScalarType;
            typedef typename boost::numeric::ublas::vector<ScalarType> VectorType;
            /** @brief Constructor
             * @param A matrix whose approximate inverse is calculated. Must be quadratic.
             * @param tag spai tag 
             */
            spai_precond(const MatrixType& A,
                         const spai_tag& tag): _tag(tag){
                
                //VCLMatrixType vcl_Ap((unsigned int)A.size2(), (unsigned int)A.size1()), vcl_A((unsigned int)A.size1(), (unsigned int)A.size2()), 
                //vcl_At((unsigned int)A.size1(), (unsigned int)A.size2());
                //UBLASDenseMatrixType dA = A;
                MatrixType pA(A.size1(), A.size2());
                MatrixType At;
                //std::cout<<A<<std::endl;
                if(!_tag.getIsRight()){
                    viennacl::linalg::detail::spai::sparse_transpose(A, At);
                }else{
                    At = A;
                }
                pA = At;
                viennacl::linalg::detail::spai::initPreconditioner(pA, _spai_m);
                viennacl::linalg::detail::spai::computeSPAI(At, _spai_m, _tag);
                //(At, pA, _tag.getIsRight(), _tag.getIsStatic(), (ScalarType)_tag.getResidualNormThreshold(), (unsigned int)_tag.getIterationLimit(),
                 //_spai_m);
                
            }
            /** @brief Application of current preconditioner, multiplication on the right-hand side vector
             * @param vec rhs vector
             */
            void apply(VectorType& vec) const {
                vec = viennacl::linalg::prod(_spai_m, vec);
            }
        private:
            // variables
            spai_tag _tag;
            // result of SPAI
            MatrixType _spai_m;
        };   
        
        //VIENNACL version
        template <typename ScalarType, unsigned int MAT_ALIGNMENT>
        class spai_precond< viennacl::compressed_matrix<ScalarType, MAT_ALIGNMENT> >
        {
            typedef viennacl::compressed_matrix<ScalarType, MAT_ALIGNMENT> MatrixType;
            typedef boost::numeric::ublas::compressed_matrix<ScalarType> UBLASSparseMatrixType;
            typedef viennacl::vector<ScalarType> VectorType;
            typedef viennacl::matrix<ScalarType> VCLDenseMatrixType;
            
            typedef boost::numeric::ublas::vector<ScalarType> UBLASVectorType;
        public:
            
            /** @brief Constructor
             * @param A matrix whose approximate inverse is calculated. Must be quadratic.
             * @param tag spai tag
             */
            spai_precond(const MatrixType& A,
                         const spai_tag& tag): _tag(tag)
            {
                viennacl::linalg::kernels::spai<ScalarType, 1>::init();
              
                MatrixType At(A.size1(), A.size2());
                UBLASSparseMatrixType ubls_A, ubls_spai_m;
                UBLASSparseMatrixType ubls_At;
                viennacl::copy(A, ubls_A);;
                if(!_tag.getIsRight()){
                    viennacl::linalg::detail::spai::sparse_transpose(ubls_A, ubls_At);
                }
                else{
                    ubls_At = ubls_A;
                }
                //current pattern is A
                //pA = ubls_At;
                //execute SPAI with ublas matrix types
                viennacl::linalg::detail::spai::initPreconditioner(ubls_At, ubls_spai_m);
                viennacl::copy(ubls_At, At);
                viennacl::linalg::detail::spai::computeSPAI(At, ubls_At, ubls_spai_m, _spai_m, _tag);
                //viennacl::copy(ubls_spai_m, _spai_m);
                
            }
            /** @brief Application of current preconditioner, multiplication on the right-hand side vector
             * @param vec rhs vector
             */
            void apply(VectorType& vec) const {
                vec = viennacl::linalg::prod(_spai_m, vec);
            }
        private:
            // variables
            spai_tag _tag;
            // result of SPAI
            MatrixType _spai_m;
        };
        
        
        //
        // FSPAI
        //
        
        /** @brief Implementation of the Factored SParse Approximate Inverse Algorithm
        * @param Matrix matrix that is used for computations
        * @param Vector vector that is used for computations
        */
        //UBLAS version
        template <typename MatrixType>
        class fspai_precond
        {
            typedef typename MatrixType::value_type ScalarType;
            typedef typename boost::numeric::ublas::vector<ScalarType> VectorType;
            typedef typename boost::numeric::ublas::matrix<ScalarType> UBLASDenseMatrixType;
            typedef typename viennacl::matrix<ScalarType> VCLMatrixType;
        public:
            
            /** @brief Constructor
            * @param A matrix whose approximate inverse is calculated. Must be quadratic.
            * @param tag SPAI configuration tag
            */
            fspai_precond(const MatrixType& A,
                        const fspai_tag& tag): tag_(tag)
            {
                MatrixType pA = A;
                viennacl::linalg::detail::spai::computeFSPAI(A, pA, L, L_trans, tag_);
            }
            
            /** @brief Application of current preconditioner, multiplication on the right-hand side vector
            * @param vec rhs vector
            */
            void apply(VectorType& vec) const 
            {
              VectorType temp = viennacl::linalg::prod(L_trans, vec);
              vec = viennacl::linalg::prod(L, temp);
            }
            
        private:
            // variables
            const fspai_tag & tag_;
            // result of SPAI
            MatrixType L;
            MatrixType L_trans;
        };   
        

        
        
        
        //
        // ViennaCL version
        //
        template <typename ScalarType, unsigned int MAT_ALIGNMENT>
        class fspai_precond< viennacl::compressed_matrix<ScalarType, MAT_ALIGNMENT> >
        {
            typedef viennacl::compressed_matrix<ScalarType, MAT_ALIGNMENT>   MatrixType;
            typedef viennacl::vector<ScalarType> VectorType;
            typedef viennacl::matrix<ScalarType> VCLDenseMatrixType;
            typedef boost::numeric::ublas::compressed_matrix<ScalarType> UBLASSparseMatrixType;
            typedef boost::numeric::ublas::vector<ScalarType> UBLASVectorType;
        public:
            
            /** @brief Constructor
            * @param A matrix whose approximate inverse is calculated. Must be quadratic.
            * @param tag SPAI configuration tag
            */
            fspai_precond(const MatrixType & A,
                        const fspai_tag & tag): tag_(tag){
                //UBLASSparseMatrixType ubls_A;
                UBLASSparseMatrixType ublas_A(A.size1(), A.size2());
                UBLASSparseMatrixType pA(A.size1(), A.size2());
                UBLASSparseMatrixType ublas_L(A.size1(), A.size2());
                UBLASSparseMatrixType ublas_L_trans(A.size1(), A.size2());
                viennacl::copy(A, ublas_A);
                //viennacl::copy(ubls_A, vcl_A);
                //vcl_At = viennacl::linalg::prod(vcl_A, vcl_A);
                //vcl_pA = viennacl::linalg::prod(vcl_A, vcl_At);
                //viennacl::copy(vcl_pA, pA);
                pA = ublas_A;
                //execute SPAI with ublas matrix types
                viennacl::linalg::detail::spai::computeFSPAI(ublas_A, pA, ublas_L, ublas_L_trans, tag_);
                //copy back to GPU
                viennacl::copy(ublas_L, L);
                viennacl::copy(ublas_L_trans, L_trans);
            }
            
            
            /** @brief Application of current preconditioner, multiplication on the right-hand side vector
            * @param vec rhs vector
            */
            void apply(VectorType& vec) const 
            {
              VectorType temp(vec.size());
              temp = viennacl::linalg::prod(L_trans, vec);
              vec = viennacl::linalg::prod(L, temp);
            }
            
        private:
            // variables
            const fspai_tag & tag_;
            MatrixType L;
            MatrixType L_trans;
        };
        
        
    }
}
#endif