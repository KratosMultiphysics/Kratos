#ifndef VIENNACL_DIRECT_SOLVE_HPP_
#define VIENNACL_DIRECT_SOLVE_HPP_

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

/** @file direct_solve.hpp
    @brief Implementations of dense direct solvers are found here.
*/

#include "viennacl/vector.hpp"
#include "viennacl/matrix.hpp"
#include "viennacl/tools/matrix_kernel_class_deducer.hpp"
#include "viennacl/tools/matrix_solve_kernel_class_deducer.hpp"
#include "viennacl/ocl/kernel.hpp"
#include "viennacl/ocl/device.hpp"
#include "viennacl/ocl/handle.hpp"


namespace viennacl
{
  namespace linalg
  {
    ////////////////// upper triangular solver (upper_tag) //////////////////////////////////////
    /** @brief Direct inplace solver for dense upper triangular systems
    *
    * @param mat    The system matrix
    * @param B      The matrix of row vectors, where the solution is directly written to
    */
    template<typename SCALARTYPE, typename F1, typename F2, unsigned int A1, unsigned int A2, typename SOLVERTAG>
    void inplace_solve(const matrix<SCALARTYPE, F1, A1> & mat,
                       matrix<SCALARTYPE, F2, A2> & B,
                       SOLVERTAG)
    {
      assert(mat.size1() == mat.size2());
      assert(mat.size2() == B.size1());
      
      typedef typename viennacl::tools::MATRIX_SOLVE_KERNEL_CLASS_DEDUCER< matrix<SCALARTYPE, F1, A1>,
                                                                           matrix<SCALARTYPE, F2, A2> >::ResultType    KernelClass;
      KernelClass::init();
      
      std::stringstream ss;
      ss << SOLVERTAG::name() << "_solve";
      viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(KernelClass::program_name(), ss.str());

      k.global_work_size(0, B.size2() * k.local_work_size());
      viennacl::ocl::enqueue(k(mat, cl_uint(mat.size1()), cl_uint(mat.size2()),
                                    cl_uint(mat.internal_size1()), cl_uint(mat.internal_size2()),
                               B,   cl_uint(B.size1()), cl_uint(B.size2()),
                                    cl_uint(B.internal_size1()), cl_uint(B.internal_size2()))
                            );        
    }
    
    /** @brief Direct inplace solver for dense upper triangular systems
    *
    * @param mat    The system matrix
    * @param B      The (transposed) matrix of row vectors, where the solution is directly written to
    */
    template<typename SCALARTYPE, typename F1, typename F2, unsigned int A1, unsigned int A2, typename SOLVERTAG>
    void inplace_solve(const matrix<SCALARTYPE, F1, A1> & mat,
                       const matrix_expression< const matrix<SCALARTYPE, F2, A2>,
                                                const matrix<SCALARTYPE, F2, A2>,
                                                op_trans> & B,
                       SOLVERTAG)
    {
      assert(mat.size1() == mat.size2());
      assert(mat.size2() == B.lhs().size2());
      
      typedef typename viennacl::tools::MATRIX_SOLVE_KERNEL_CLASS_DEDUCER< matrix<SCALARTYPE, F1, A1>,
                                                                           matrix<SCALARTYPE, F2, A2> >::ResultType    KernelClass;
      KernelClass::init();

      std::stringstream ss;
      ss << SOLVERTAG::name() << "_trans_solve";
      viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(KernelClass::program_name(), ss.str());

      k.global_work_size(0, B.lhs().size1() * k.local_work_size());
      viennacl::ocl::enqueue(k(mat, cl_uint(mat.size1()), cl_uint(mat.size2()),
                                    cl_uint(mat.internal_size1()), cl_uint(mat.internal_size2()),
                               B.lhs(), cl_uint(B.lhs().size1()), cl_uint(B.lhs().size2()),
                                        cl_uint(B.lhs().internal_size1()), cl_uint(B.lhs().internal_size2()))
                            );     
    }
    
    //upper triangular solver for transposed lower triangular matrices
    /** @brief Direct inplace solver for dense upper triangular systems that stem from transposed lower triangular systems
    *
    * @param proxy    The system matrix proxy
    * @param B        The matrix holding the load vectors, where the solution is directly written to
    */
    template<typename SCALARTYPE, typename F1, typename F2, unsigned int A1, unsigned int A2, typename SOLVERTAG>
    void inplace_solve(const matrix_expression< const matrix<SCALARTYPE, F1, A1>,
                                                const matrix<SCALARTYPE, F1, A1>,
                                                op_trans> & proxy,
                       matrix<SCALARTYPE, F2, A2> & B,
                       SOLVERTAG)
    {
      assert(proxy.lhs().size1() == proxy.lhs().size2());
      assert(proxy.lhs().size2() == B.size1());
      
      typedef typename viennacl::tools::MATRIX_SOLVE_KERNEL_CLASS_DEDUCER< matrix<SCALARTYPE, F1, A1>,
                                                                           matrix<SCALARTYPE, F2, A2> >::ResultType    KernelClass;
      KernelClass::init();

      std::stringstream ss;
      ss << "trans_" << SOLVERTAG::name() << "_solve";
      viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(KernelClass::program_name(), ss.str());

      k.global_work_size(0, B.size2() * k.local_work_size());
      viennacl::ocl::enqueue(k(proxy.lhs(), cl_uint(proxy.lhs().size1()), cl_uint(proxy.lhs().size2()),
                                            cl_uint(proxy.lhs().internal_size1()), cl_uint(proxy.lhs().internal_size2()),
                                         B, cl_uint(B.size1()), cl_uint(B.size2()),
                                            cl_uint(B.internal_size1()), cl_uint(B.internal_size2()))
                            );        
    }

    /** @brief Direct inplace solver for dense upper triangular systems that stem from transposed lower triangular systems
    *
    * @param proxy    The system matrix proxy
    * @param B        The matrix holding the load vectors, where the solution is directly written to
    */
    template<typename SCALARTYPE, typename F1, typename F2, unsigned int A1, unsigned int A2, typename SOLVERTAG>
    void inplace_solve(const matrix_expression< const matrix<SCALARTYPE, F1, A1>,
                                                const matrix<SCALARTYPE, F1, A1>,
                                                op_trans> & proxy,
                       const matrix_expression< const matrix<SCALARTYPE, F2, A2>,
                                                const matrix<SCALARTYPE, F2, A2>,
                                                op_trans> & B,
                       SOLVERTAG)
    {
      assert(proxy.lhs().size1() == proxy.lhs().size2());
      assert(proxy.lhs().size2() == B.lhs().size2());
      
      typedef typename viennacl::tools::MATRIX_SOLVE_KERNEL_CLASS_DEDUCER< matrix<SCALARTYPE, F1, A1>,
                                                                           matrix<SCALARTYPE, F2, A2> >::ResultType    KernelClass;
      KernelClass::init();

      std::stringstream ss;
      ss << "trans_" << SOLVERTAG::name() << "_trans_solve";
      viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(KernelClass::program_name(), ss.str());

      k.global_work_size(0, B.lhs().size1() * k.local_work_size());
      viennacl::ocl::enqueue(k(proxy.lhs(), cl_uint(proxy.lhs().size1()), cl_uint(proxy.lhs().size2()),
                                            cl_uint(proxy.lhs().internal_size1()), cl_uint(proxy.lhs().internal_size2()),
                               B.lhs(), cl_uint(B.lhs().size1()), cl_uint(B.lhs().size2()),
                                        cl_uint(B.lhs().internal_size1()), cl_uint(B.lhs().internal_size2()))
                            );        
    }

    template<typename SCALARTYPE, typename F, unsigned int ALIGNMENT, unsigned int VEC_ALIGNMENT, typename SOLVERTAG>
    void inplace_solve(const matrix<SCALARTYPE, F, ALIGNMENT> & mat,
                       vector<SCALARTYPE, VEC_ALIGNMENT> & vec,
                       SOLVERTAG)
    {
      assert(mat.size1() == vec.size());
      assert(mat.size2() == vec.size());
      
      typedef typename viennacl::tools::MATRIX_KERNEL_CLASS_DEDUCER< matrix<SCALARTYPE, F, ALIGNMENT> >::ResultType    KernelClass;

      std::stringstream ss;
      ss << SOLVERTAG::name() << "_triangular_substitute_inplace";
      viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(KernelClass::program_name(), ss.str());

      k.global_work_size(0, k.local_work_size());
      viennacl::ocl::enqueue(k(mat, cl_uint(mat.size1()), cl_uint(mat.size2()), 
                                    cl_uint(mat.internal_size1()), cl_uint(mat.internal_size2()), vec));        
    }

    /** @brief Direct inplace solver for dense upper triangular systems that stem from transposed lower triangular systems
    *
    * @param proxy    The system matrix proxy
    * @param vec    The load vector, where the solution is directly written to
    */
    template<typename SCALARTYPE, typename F, unsigned int ALIGNMENT, unsigned int VEC_ALIGNMENT, typename SOLVERTAG>
    void inplace_solve(const matrix_expression< const matrix<SCALARTYPE, F, ALIGNMENT>,
                                                const matrix<SCALARTYPE, F, ALIGNMENT>,
                                                op_trans> & proxy,
                       vector<SCALARTYPE, VEC_ALIGNMENT> & vec,
                       SOLVERTAG)
    {
      assert(proxy.lhs().size1() == vec.size());
      assert(proxy.lhs().size2() == vec.size());

      typedef typename viennacl::tools::MATRIX_KERNEL_CLASS_DEDUCER< matrix<SCALARTYPE, F, ALIGNMENT> >::ResultType    KernelClass;
      
      std::stringstream ss;
      ss << "trans_" << SOLVERTAG::name() << "_triangular_substitute_inplace";
      viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(KernelClass::program_name(), ss.str());
      
      k.global_work_size(0, k.local_work_size());
      viennacl::ocl::enqueue(k(proxy.lhs(), cl_uint(proxy.lhs().size1()), cl_uint(proxy.lhs().size2()),
                                            cl_uint(proxy.lhs().internal_size1()), cl_uint(proxy.lhs().internal_size2()), vec));        
    }
    
    /////////////////// general wrappers for non-inplace solution //////////////////////    

    /** @brief Convenience functions for C = solve(A, B, some_tag()); Creates a temporary result matrix and forwards the request to inplace_solve()
    *
    * @param A    The system matrix
    * @param B    The matrix of load vectors
    * @param tag    Dispatch tag
    */
    template<typename SCALARTYPE, typename F1, typename F2, unsigned int ALIGNMENT_A, unsigned int ALIGNMENT_B, typename TAG>
    matrix<SCALARTYPE, F2, ALIGNMENT_B> solve(const matrix<SCALARTYPE, F1, ALIGNMENT_A> & A,
                                        const matrix<SCALARTYPE, F2, ALIGNMENT_B> & B,
                                        TAG const & tag)
    {
      // do an inplace solve on the result vector:
      matrix<SCALARTYPE, F2, ALIGNMENT_A> result(B.size1(), B.size2());
      result = B;
    
      inplace_solve(A, result, tag);
    
      return result;
    }

    /** @brief Convenience functions for C = solve(A, B^T, some_tag()); Creates a temporary result matrix and forwards the request to inplace_solve()
    *
    * @param A    The system matrix
    * @param proxy  The transposed load vector
    * @param tag    Dispatch tag
    */
    template<typename SCALARTYPE, typename F1, typename F2, unsigned int ALIGNMENT_A, unsigned int ALIGNMENT_B, typename TAG>
    matrix<SCALARTYPE, F2, ALIGNMENT_B> solve(const matrix<SCALARTYPE, F1, ALIGNMENT_A> & A,
                                        const matrix_expression< const matrix<SCALARTYPE, F2, ALIGNMENT_B>,
                                                                     const matrix<SCALARTYPE, F2, ALIGNMENT_B>,
                                                                     op_trans> & proxy,
                                        TAG const & tag)
    {
      // do an inplace solve on the result vector:
      matrix<SCALARTYPE, F2, ALIGNMENT_B> result(proxy.lhs().size2(), proxy.lhs().size1());
      result = proxy;
    
      inplace_solve(A, result, tag);
    
      return result;
    }

    /** @brief Convenience functions for result = solve(mat, vec, some_tag()); Creates a temporary result vector and forwards the request to inplace_solve()
    *
    * @param mat    The system matrix
    * @param vec    The load vector
    * @param tag    Dispatch tag
    */
    template<typename SCALARTYPE, typename F, unsigned int ALIGNMENT, unsigned int VEC_ALIGNMENT, typename TAG>
    vector<SCALARTYPE, VEC_ALIGNMENT> solve(const matrix<SCALARTYPE, F, ALIGNMENT> & mat,
                                        const vector<SCALARTYPE, VEC_ALIGNMENT> & vec,
                                        TAG const & tag)
    {
      // do an inplace solve on the result vector:
      vector<SCALARTYPE, VEC_ALIGNMENT> result(vec.size());
      result = vec;
    
      inplace_solve(mat, result, tag);
    
      return result;
    }
    
    
    ///////////// transposed system matrix:
    /** @brief Convenience functions for result = solve(trans(mat), B, some_tag()); Creates a temporary result matrix and forwards the request to inplace_solve()
    *
    * @param proxy  The transposed system matrix proxy
    * @param B      The matrix of load vectors
    * @param tag    Dispatch tag
    */
    template<typename SCALARTYPE, typename F1, typename F2, unsigned int ALIGNMENT_A, unsigned int ALIGNMENT_B, typename TAG>
    matrix<SCALARTYPE, F2, ALIGNMENT_B> solve(const matrix_expression< const matrix<SCALARTYPE, F1, ALIGNMENT_A>,
                                                                     const matrix<SCALARTYPE, F1, ALIGNMENT_A>,
                                                                     op_trans> & proxy,
                                            const matrix<SCALARTYPE, F2, ALIGNMENT_B> & B,
                                            TAG const & tag)
    {
      // do an inplace solve on the result vector:
      matrix<SCALARTYPE, F2, ALIGNMENT_B> result(B.size1(), B.size2());
      result = B;
    
      inplace_solve(proxy, result, tag);
    
      return result;
    }
    
    
    /** @brief Convenience functions for result = solve(trans(mat), vec, some_tag()); Creates a temporary result vector and forwards the request to inplace_solve()
    *
    * @param proxy_A  The transposed system matrix proxy
    * @param proxy_B  The transposed matrix of load vectors, where the solution is directly written to
    * @param tag    Dispatch tag
    */
    template<typename SCALARTYPE, typename F1, typename F2, unsigned int ALIGNMENT_A, unsigned int ALIGNMENT_B, typename TAG>
    matrix<SCALARTYPE, F2, ALIGNMENT_B> solve(const matrix_expression< const matrix<SCALARTYPE, F1, ALIGNMENT_A>,
                                                                     const matrix<SCALARTYPE, F1, ALIGNMENT_A>,
                                                                     op_trans> & proxy_A,
                                            const matrix_expression< const matrix<SCALARTYPE, F2, ALIGNMENT_B>,
                                                                     const matrix<SCALARTYPE, F2, ALIGNMENT_B>,
                                                                     op_trans> & proxy_B,
                                            TAG const & tag)
    {
      // do an inplace solve on the result vector:
      matrix<SCALARTYPE, F2, ALIGNMENT_B> result(proxy_B.lhs().size2(), proxy_B.lhs().size1());
      result = trans(proxy_B.lhs());
    
      inplace_solve(proxy_A, result, tag);
    
      return result;
    }
    
    /** @brief Convenience functions for result = solve(trans(mat), vec, some_tag()); Creates a temporary result vector and forwards the request to inplace_solve()
    *
    * @param proxy  The transposed system matrix proxy
    * @param vec    The load vector, where the solution is directly written to
    * @param tag    Dispatch tag
    */
    template<typename SCALARTYPE, typename F, unsigned int ALIGNMENT, unsigned int VEC_ALIGNMENT, typename TAG>
    vector<SCALARTYPE, VEC_ALIGNMENT> solve(const matrix_expression< const matrix<SCALARTYPE, F, ALIGNMENT>,
                                                                     const matrix<SCALARTYPE, F, ALIGNMENT>,
                                                                     op_trans> & proxy,
                                            const vector<SCALARTYPE, VEC_ALIGNMENT> & vec,
                                            TAG const & tag)
    {
      // do an inplace solve on the result vector:
      vector<SCALARTYPE, VEC_ALIGNMENT> result(vec.size());
      result = vec;
    
      inplace_solve(proxy, result, tag);
    
      return result;
    }
    
    
    ///////////////////////////// lu factorization ///////////////////////
    /** @brief LU factorization of a dense matrix.
    *
    * @param mat    The system matrix, where the LU matrices are directly written to. The implicit unit diagonal of L is not written.
    */
    template<typename SCALARTYPE, typename F, unsigned int ALIGNMENT>
    void lu_factorize(matrix<SCALARTYPE, F, ALIGNMENT> & mat)
    {
      assert(mat.size1() == mat.size2());

      typedef typename viennacl::tools::MATRIX_KERNEL_CLASS_DEDUCER< matrix<SCALARTYPE, F, ALIGNMENT> >::ResultType    KernelClass;
      
      viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(KernelClass::program_name(), "lu_factorize");
      
      k.global_work_size(0, k.local_work_size());
      viennacl::ocl::enqueue(k(mat, cl_uint(mat.size1()), cl_uint(mat.size2()),
                                    cl_uint(mat.internal_size1()), cl_uint(mat.internal_size2())) );        
    }


    /** @brief LU substitution for the system LU = rhs.
    *
    * @param A    The system matrix, where the LU matrices are directly written to. The implicit unit diagonal of L is not written.
    * @param B    The matrix of load vectors, where the solution is directly written to
    */
    template<typename SCALARTYPE, typename F1, typename F2, unsigned int ALIGNMENT_A, unsigned int ALIGNMENT_B>
    void lu_substitute(matrix<SCALARTYPE, F1, ALIGNMENT_A> const & A,
                       matrix<SCALARTYPE, F2, ALIGNMENT_B> & B)
    {
      assert(A.size1() == A.size2());
      assert(A.size1() == A.size2());
      inplace_solve(A, B, unit_lower_tag());
      inplace_solve(A, B, upper_tag());
    }

    /** @brief LU substitution for the system LU = rhs.
    *
    * @param mat    The system matrix, where the LU matrices are directly written to. The implicit unit diagonal of L is not written.
    * @param vec    The load vector, where the solution is directly written to
    */
    template<typename SCALARTYPE, typename F, unsigned int ALIGNMENT, unsigned int VEC_ALIGNMENT>
    void lu_substitute(matrix<SCALARTYPE, F, ALIGNMENT> const & mat,
                       vector<SCALARTYPE, VEC_ALIGNMENT> & vec)
    {
      assert(mat.size1() == mat.size2());
      inplace_solve(mat, vec, unit_lower_tag());
      inplace_solve(mat, vec, upper_tag());
    }

  }
}

#endif
