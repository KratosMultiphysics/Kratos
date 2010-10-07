/* =======================================================================
   Copyright (c) 2010, Institute for Microelectronics, TU Vienna.
   http://www.iue.tuwien.ac.at
                             -----------------
                     ViennaCL - The Vienna Computing Library
                             -----------------
                            
   authors:    Karl Rupp                          rupp@iue.tuwien.ac.at
               Florian Rudolf                     flo.rudy+viennacl@gmail.com
               Josef Weinbub                      weinbub@iue.tuwien.ac.at

   license:    MIT (X11), see file LICENSE in the ViennaCL base directory

   file changelog: - May 28, 2010   New from scratch for first release
======================================================================= */

#ifndef _VIENNACL_DIRECT_SOLVE_HPP_
#define _VIENNACL_DIRECT_SOLVE_HPP_

#include "viennacl/vector.hpp"
#include "viennacl/matrix.hpp"
#include "viennacl/ocl/kernel.hpp"
#include "viennacl/ocl/device.hpp"
#include "viennacl/ocl/handle.hpp"


namespace viennacl
{
  namespace linalg
  {
    ////////////////// upper triangular solver //////////////////////////////////////
    /** @brief Direct inplace solver for dense upper triangular systems
    *
    * @param mat    The system matrix
    * @param vec    The load vector, where the solution is directly written to
    */
    template<typename SCALARTYPE, typename F, unsigned int ALIGNMENT, unsigned int VEC_ALIGNMENT>
    void inplace_solve(const matrix<SCALARTYPE, F, ALIGNMENT> & mat,
                       vector<SCALARTYPE, VEC_ALIGNMENT> & vec,
                       upper_tag)
    {
      assert(mat.size1() == vec.size());
      assert(mat.size2() == vec.size());

      unsigned int pos = 0;
      viennacl::linalg::kernels::matrix<SCALARTYPE,ALIGNMENT>::upper_triangular_substitute_inplace.setArgument(pos++, mat.handle());
      viennacl::linalg::kernels::matrix<SCALARTYPE,ALIGNMENT>::upper_triangular_substitute_inplace.setArgument(pos++, vec.handle());
      viennacl::linalg::kernels::matrix<SCALARTYPE,ALIGNMENT>::upper_triangular_substitute_inplace.setArgument(pos++, mat.internal_size2());
      viennacl::linalg::kernels::matrix<SCALARTYPE,ALIGNMENT>::upper_triangular_substitute_inplace.setArgument(pos++, vec.size());
      
      viennacl::linalg::kernels::matrix<SCALARTYPE,ALIGNMENT>::upper_triangular_substitute_inplace.start1D(
                                                viennacl::linalg::kernels::matrix<SCALARTYPE,ALIGNMENT>::upper_triangular_substitute_inplace.work_items_per_group(),
                                                viennacl::linalg::kernels::matrix<SCALARTYPE,ALIGNMENT>::upper_triangular_substitute_inplace.work_items_per_group());
    }

    /** @brief Direct inplace solver for dense upper triangular systems with unit diagonal
    *
    * @param mat    The system matrix
    * @param vec    The load vector, where the solution is directly written to
    */
    template<typename SCALARTYPE, typename F, unsigned int ALIGNMENT, unsigned int VEC_ALIGNMENT>
    void inplace_solve(const matrix<SCALARTYPE, F, ALIGNMENT> & mat,
                       vector<SCALARTYPE, VEC_ALIGNMENT> & vec,
                       unit_upper_tag)
    {
      assert(mat.size1() == vec.size());
      assert(mat.size2() == vec.size());

      unsigned int pos = 0;
      viennacl::linalg::kernels::matrix<SCALARTYPE,ALIGNMENT>::unit_upper_triangular_substitute_inplace.setArgument(pos++, mat.handle());
      viennacl::linalg::kernels::matrix<SCALARTYPE,ALIGNMENT>::unit_upper_triangular_substitute_inplace.setArgument(pos++, vec.handle());
      viennacl::linalg::kernels::matrix<SCALARTYPE,ALIGNMENT>::unit_upper_triangular_substitute_inplace.setArgument(pos++, mat.internal_size2());
      viennacl::linalg::kernels::matrix<SCALARTYPE,ALIGNMENT>::unit_upper_triangular_substitute_inplace.setArgument(pos++, vec.size());
      
      viennacl::linalg::kernels::matrix<SCALARTYPE,ALIGNMENT>::unit_upper_triangular_substitute_inplace.start1D(
                                          viennacl::linalg::kernels::matrix<SCALARTYPE,ALIGNMENT>::unit_upper_triangular_substitute_inplace.work_items_per_group(),
                                          viennacl::linalg::kernels::matrix<SCALARTYPE,ALIGNMENT>::unit_upper_triangular_substitute_inplace.work_items_per_group());
    }

    //upper triangular solver for transposed lower triangular matrices
    /** @brief Direct inplace solver for dense upper triangular systems that stem from transposed lower triangular systems
    *
    * @param proxy    The system matrix proxy
    * @param vec    The load vector, where the solution is directly written to
    */
    template<typename SCALARTYPE, typename F, unsigned int ALIGNMENT, unsigned int VEC_ALIGNMENT>
    void inplace_solve(const transposed_matrix_proxy<SCALARTYPE, F, ALIGNMENT> & proxy,
                       vector<SCALARTYPE, VEC_ALIGNMENT> & vec,
                       upper_tag)
    {
      assert(proxy.get_matrix().size1() == vec.size());
      assert(proxy.get_matrix().size2() == vec.size());

      unsigned int pos = 0;
      viennacl::linalg::kernels::matrix<SCALARTYPE,ALIGNMENT>::trans_upper_triangular_substitute_inplace.setArgument(pos++, proxy.get_matrix().handle());
      viennacl::linalg::kernels::matrix<SCALARTYPE,ALIGNMENT>::trans_upper_triangular_substitute_inplace.setArgument(pos++, vec.handle());
      viennacl::linalg::kernels::matrix<SCALARTYPE,ALIGNMENT>::trans_upper_triangular_substitute_inplace.setArgument(pos++, proxy.get_matrix().internal_size2());
      viennacl::linalg::kernels::matrix<SCALARTYPE,ALIGNMENT>::trans_upper_triangular_substitute_inplace.setArgument(pos++, vec.size());
      
      viennacl::linalg::kernels::matrix<SCALARTYPE,ALIGNMENT>::trans_upper_triangular_substitute_inplace.start1D(
                                              viennacl::linalg::kernels::matrix<SCALARTYPE,ALIGNMENT>::trans_upper_triangular_substitute_inplace.work_items_per_group(),
                                              viennacl::linalg::kernels::matrix<SCALARTYPE,ALIGNMENT>::trans_upper_triangular_substitute_inplace.work_items_per_group());
    }
    
    
    ////////////////// lower triangular solver //////////////////////////////////////
    
    /** @brief Direct inplace solver for lower triangular systems
    *
    * @param mat    The system matrix
    * @param vec    The load vector, where the solution is directly written to
    */
    template<typename SCALARTYPE, typename F, unsigned int ALIGNMENT, unsigned int VEC_ALIGNMENT>
    void inplace_solve(const matrix<SCALARTYPE, F, ALIGNMENT> & mat,
                       vector<SCALARTYPE, VEC_ALIGNMENT> & vec,
                       lower_tag)
    {
      assert(mat.size1() == vec.size());
      assert(mat.size2() == vec.size());

      unsigned int pos = 0;
      viennacl::linalg::kernels::matrix<SCALARTYPE,ALIGNMENT>::lower_triangular_substitute_inplace.setArgument(pos++, mat.handle());
      viennacl::linalg::kernels::matrix<SCALARTYPE,ALIGNMENT>::lower_triangular_substitute_inplace.setArgument(pos++, vec.handle());
      viennacl::linalg::kernels::matrix<SCALARTYPE,ALIGNMENT>::lower_triangular_substitute_inplace.setArgument(pos++, mat.internal_size2());
      viennacl::linalg::kernels::matrix<SCALARTYPE,ALIGNMENT>::lower_triangular_substitute_inplace.setArgument(pos++, vec.size());
      
      viennacl::linalg::kernels::matrix<SCALARTYPE,ALIGNMENT>::lower_triangular_substitute_inplace.start1D(
                                        viennacl::linalg::kernels::matrix<SCALARTYPE,ALIGNMENT>::lower_triangular_substitute_inplace.work_items_per_group(),
                                        viennacl::linalg::kernels::matrix<SCALARTYPE,ALIGNMENT>::lower_triangular_substitute_inplace.work_items_per_group());
    }

    /** @brief Direct inplace solver for lower triangular systems with unit diagonal
    *
    * @param mat    The system matrix
    * @param vec    The load vector, where the solution is directly written to
    */
    template<typename SCALARTYPE, typename F, unsigned int ALIGNMENT, unsigned int VEC_ALIGNMENT>
    void inplace_solve(const matrix<SCALARTYPE, F, ALIGNMENT> & mat,
                       vector<SCALARTYPE, VEC_ALIGNMENT> & vec,
                       unit_lower_tag)
    {
      assert(mat.size1() == vec.size());
      assert(mat.size2() == vec.size());

      unsigned int pos = 0;
      viennacl::linalg::kernels::matrix<SCALARTYPE,ALIGNMENT>::unit_lower_triangular_substitute_inplace.setArgument(pos++, mat.handle());
      viennacl::linalg::kernels::matrix<SCALARTYPE,ALIGNMENT>::unit_lower_triangular_substitute_inplace.setArgument(pos++, vec.handle());
      viennacl::linalg::kernels::matrix<SCALARTYPE,ALIGNMENT>::unit_lower_triangular_substitute_inplace.setArgument(pos++, mat.internal_size2());
      viennacl::linalg::kernels::matrix<SCALARTYPE,ALIGNMENT>::unit_lower_triangular_substitute_inplace.setArgument(pos++, vec.size());
      
      viennacl::linalg::kernels::matrix<SCALARTYPE,ALIGNMENT>::unit_lower_triangular_substitute_inplace.start1D(
                                                    viennacl::linalg::kernels::matrix<SCALARTYPE,ALIGNMENT>::unit_lower_triangular_substitute_inplace.work_items_per_group(),
                                                    viennacl::linalg::kernels::matrix<SCALARTYPE,ALIGNMENT>::unit_lower_triangular_substitute_inplace.work_items_per_group());
    }

    /** @brief Direct inplace solver for dense lower triangular systems that stem from transposed upper triangular systems
    *
    * @param proxy    The system matrix proxy
    * @param vec    The load vector, where the solution is directly written to
    */
    template<typename SCALARTYPE, typename F, unsigned int ALIGNMENT, unsigned int VEC_ALIGNMENT>
    void inplace_solve(const transposed_matrix_proxy<SCALARTYPE, F, ALIGNMENT> & proxy,
                       vector<SCALARTYPE, VEC_ALIGNMENT> & vec,
                       lower_tag)
    {
      assert(proxy.get_matrix().size1() == vec.size());
      assert(proxy.get_matrix().size2() == vec.size());

      unsigned int pos = 0;
      viennacl::linalg::kernels::matrix<SCALARTYPE,ALIGNMENT>::trans_lower_triangular_substitute_inplace.setArgument(pos++, proxy.get_matrix().handle());
      viennacl::linalg::kernels::matrix<SCALARTYPE,ALIGNMENT>::trans_lower_triangular_substitute_inplace.setArgument(pos++, vec.handle());
      viennacl::linalg::kernels::matrix<SCALARTYPE,ALIGNMENT>::trans_lower_triangular_substitute_inplace.setArgument(pos++, proxy.get_matrix().internal_size2());
      viennacl::linalg::kernels::matrix<SCALARTYPE,ALIGNMENT>::trans_lower_triangular_substitute_inplace.setArgument(pos++, vec.size());
      
      viennacl::linalg::kernels::matrix<SCALARTYPE,ALIGNMENT>::trans_lower_triangular_substitute_inplace.start1D(
                                            viennacl::linalg::kernels::matrix<SCALARTYPE,ALIGNMENT>::trans_lower_triangular_substitute_inplace.work_items_per_group(),
                                            viennacl::linalg::kernels::matrix<SCALARTYPE,ALIGNMENT>::trans_lower_triangular_substitute_inplace.work_items_per_group());
    }
   
   
    /////////////////// general wrappers for non-inplace solution //////////////////////    

    /** @brief Convenience functions for result = solve(mat, vec, some_tag()); Creates a temporary result vector and forwards the request to inplace_solve()
    *
    * @param mat    The system matrix
    * @param vec    The load vector, where the solution is directly written to
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

    /** @brief Convenience functions for result = solve(trans(mat), vec, some_tag()); Creates a temporary result vector and forwards the request to inplace_solve()
    *
    * @param proxy  The transposed system matrix proxy
    * @param vec    The load vector, where the solution is directly written to
    * @param tag    Dispatch tag
    */
    template<typename SCALARTYPE, typename F, unsigned int ALIGNMENT, unsigned int VEC_ALIGNMENT, typename TAG>
    vector<SCALARTYPE, VEC_ALIGNMENT> solve(const transposed_matrix_proxy<SCALARTYPE, F, ALIGNMENT> & proxy,
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

      unsigned int pos = 0;
      viennacl::linalg::kernels::matrix<SCALARTYPE,ALIGNMENT>::lu_factorize.setArgument(pos++, mat.handle());
      viennacl::linalg::kernels::matrix<SCALARTYPE,ALIGNMENT>::lu_factorize.setArgument(pos++, mat.internal_size2());
      viennacl::linalg::kernels::matrix<SCALARTYPE,ALIGNMENT>::lu_factorize.setArgument(pos++, mat.size1());
      
      viennacl::linalg::kernels::matrix<SCALARTYPE,ALIGNMENT>::lu_factorize.start1D(
                                              viennacl::linalg::kernels::matrix<SCALARTYPE,ALIGNMENT>::lu_factorize.work_items_per_group(),
                                              viennacl::linalg::kernels::matrix<SCALARTYPE,ALIGNMENT>::lu_factorize.work_items_per_group());
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
