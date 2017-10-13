#ifndef VIENNACL_ELL_MATRIX_HPP_
#define VIENNACL_ELL_MATRIX_HPP_

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

/** @file viennacl/ell_matrix.hpp
    @brief Implementation of the ell_matrix class

    Contributed by Volodymyr Kysenko.
*/


#include "viennacl/forwards.h"
#include "viennacl/vector.hpp"

#include "viennacl/tools/tools.hpp"

#include "viennacl/linalg/sparse_matrix_operations.hpp"

namespace viennacl
{
    /** @brief Sparse matrix class using the ELLPACK format for storing the nonzeros.
      *
      * This format works best for matrices where the number of nonzeros per row is mostly the same.
      * Finite element and finite difference methods on nicely shaped domains often result in such a nonzero pattern.
      * For a matrix
      *
      *   (1 2 0 0 0)
      *   (2 3 4 0 0)
      *   (0 5 6 0 7)
      *   (0 0 8 9 0)
      *
      * the entries are layed out in chunks of size 3 as
      *   (1 2 5 8; 2 3 6 9; 0 4 7 0)
      * Note that this is a 'transposed' representation in order to maximize coalesced memory access.
      */
    template<typename SCALARTYPE, unsigned int ALIGNMENT /* see forwards.h for default argument */>
    class ell_matrix
    {
      public:
        typedef viennacl::backend::mem_handle                                                              handle_type;
        typedef scalar<typename viennacl::tools::CHECK_SCALAR_TEMPLATE_ARGUMENT<SCALARTYPE>::ResultType>   value_type;
        typedef vcl_size_t                                                                                 size_type;

        ell_matrix() : rows_(0), cols_(0), maxnnz_(0) {}

        ell_matrix(viennacl::context ctx) : rows_(0), cols_(0), maxnnz_(0)
        {
            coords_.switch_active_handle_id(ctx.memory_type());
          elements_.switch_active_handle_id(ctx.memory_type());

#ifdef VIENNACL_WITH_OPENCL
          if (ctx.memory_type() == OPENCL_MEMORY)
          {
              coords_.opencl_handle().context(ctx.opencl_context());
            elements_.opencl_handle().context(ctx.opencl_context());
          }
#endif
        }

      public:
        vcl_size_t internal_size1() const { return viennacl::tools::align_to_multiple<vcl_size_t>(rows_, ALIGNMENT); }
        vcl_size_t internal_size2() const { return viennacl::tools::align_to_multiple<vcl_size_t>(cols_, ALIGNMENT); }

        vcl_size_t size1() const { return rows_; }
        vcl_size_t size2() const { return cols_; }

        vcl_size_t internal_maxnnz() const {return viennacl::tools::align_to_multiple<vcl_size_t>(maxnnz_, ALIGNMENT); }
        vcl_size_t maxnnz() const { return maxnnz_; }

        vcl_size_t nnz() const { return rows_ * maxnnz_; }
        vcl_size_t internal_nnz() const { return internal_size1() * internal_maxnnz(); }

              handle_type & handle()       { return elements_; }
        const handle_type & handle() const { return elements_; }

              handle_type & handle2()       { return coords_; }
        const handle_type & handle2() const { return coords_; }

      #if defined(_MSC_VER) && _MSC_VER < 1500          //Visual Studio 2005 needs special treatment
        template <typename CPU_MATRIX>
        friend void copy(const CPU_MATRIX & cpu_matrix, ell_matrix & gpu_matrix );
      #else
        template <typename CPU_MATRIX, typename T, unsigned int ALIGN>
        friend void copy(const CPU_MATRIX & cpu_matrix, ell_matrix<T, ALIGN> & gpu_matrix );
      #endif

      private:
        vcl_size_t rows_;
        vcl_size_t cols_;
        vcl_size_t maxnnz_;

        handle_type coords_;
        handle_type elements_;
    };

    template <typename CPU_MATRIX, typename SCALARTYPE, unsigned int ALIGNMENT>
    void copy(const CPU_MATRIX& cpu_matrix, ell_matrix<SCALARTYPE, ALIGNMENT>& gpu_matrix )
    {
      assert( (gpu_matrix.size1() == 0 || viennacl::traits::size1(cpu_matrix) == gpu_matrix.size1()) && bool("Size mismatch") );
      assert( (gpu_matrix.size2() == 0 || viennacl::traits::size2(cpu_matrix) == gpu_matrix.size2()) && bool("Size mismatch") );

      if(cpu_matrix.size1() > 0 && cpu_matrix.size2() > 0)
      {
        //determine max capacity for row
        vcl_size_t max_entries_per_row = 0;
        for (typename CPU_MATRIX::const_iterator1 row_it = cpu_matrix.begin1(); row_it != cpu_matrix.end1(); ++row_it)
        {
          vcl_size_t num_entries = 0;
          for (typename CPU_MATRIX::const_iterator2 col_it = row_it.begin(); col_it != row_it.end(); ++col_it)
          {
              ++num_entries;
          }

          max_entries_per_row = std::max(max_entries_per_row, num_entries);
        }

        //setup GPU matrix
        gpu_matrix.maxnnz_ = max_entries_per_row;
        gpu_matrix.rows_ = cpu_matrix.size1();
        gpu_matrix.cols_ = cpu_matrix.size2();

        vcl_size_t nnz = gpu_matrix.internal_nnz();

        viennacl::backend::typesafe_host_array<unsigned int> coords(gpu_matrix.handle2(), nnz);
        std::vector<SCALARTYPE> elements(nnz, 0);

        // std::cout << "ELL_MATRIX copy " << gpu_matrix.maxnnz_ << " " << gpu_matrix.rows_ << " " << gpu_matrix.cols_ << " "
        //             << gpu_matrix.internal_maxnnz() << "\n";

        for (typename CPU_MATRIX::const_iterator1 row_it = cpu_matrix.begin1(); row_it != cpu_matrix.end1(); ++row_it)
        {
          vcl_size_t data_index = 0;

          for (typename CPU_MATRIX::const_iterator2 col_it = row_it.begin(); col_it != row_it.end(); ++col_it)
          {
            coords.set(gpu_matrix.internal_size1() * data_index + col_it.index1(), col_it.index2());
            elements[gpu_matrix.internal_size1() * data_index + col_it.index1()] = *col_it;
            //std::cout << *col_it << "\n";
              data_index++;
          }
        }

        viennacl::backend::memory_create(gpu_matrix.handle2(), coords.raw_size(),                   traits::context(gpu_matrix.handle2()), coords.get());
        viennacl::backend::memory_create(gpu_matrix.handle(), sizeof(SCALARTYPE) * elements.size(), traits::context(gpu_matrix.handle()), &(elements[0]));
      }
    }

    template <typename CPU_MATRIX, typename SCALARTYPE, unsigned int ALIGNMENT>
    void copy(const ell_matrix<SCALARTYPE, ALIGNMENT>& gpu_matrix, CPU_MATRIX& cpu_matrix)
    {
      assert( (viennacl::traits::size1(cpu_matrix) == gpu_matrix.size1()) && bool("Size mismatch") );
      assert( (viennacl::traits::size2(cpu_matrix) == gpu_matrix.size2()) && bool("Size mismatch") );

      if(gpu_matrix.size1() > 0 && gpu_matrix.size2() > 0)
      {
        std::vector<SCALARTYPE> elements(gpu_matrix.internal_nnz());
        viennacl::backend::typesafe_host_array<unsigned int> coords(gpu_matrix.handle2(), gpu_matrix.internal_nnz());

        viennacl::backend::memory_read(gpu_matrix.handle(), 0, sizeof(SCALARTYPE) * elements.size(), &(elements[0]));
        viennacl::backend::memory_read(gpu_matrix.handle2(), 0, coords.raw_size(), coords.get());

        for(vcl_size_t row = 0; row < gpu_matrix.size1(); row++)
        {
          for(vcl_size_t ind = 0; ind < gpu_matrix.internal_maxnnz(); ind++)
          {
            vcl_size_t offset = gpu_matrix.internal_size1() * ind + row;

            if(elements[offset] == static_cast<SCALARTYPE>(0.0))
                continue;

            if(coords[offset] >= gpu_matrix.size2())
            {
                std::cerr << "ViennaCL encountered invalid data " << offset << " " << ind << " " << row << " " << coords[offset] << " " << gpu_matrix.size2() << std::endl;
                return;
            }

            cpu_matrix(row, coords[offset]) = elements[offset];
          }
        }
      }
    }


    //
    // Specify available operations:
    //

    /** \cond */

    namespace linalg
    {
      namespace detail
      {
        // x = A * y
        template <typename T, unsigned int A>
        struct op_executor<vector_base<T>, op_assign, vector_expression<const ell_matrix<T, A>, const vector_base<T>, op_prod> >
        {
            static void apply(vector_base<T> & lhs, vector_expression<const ell_matrix<T, A>, const vector_base<T>, op_prod> const & rhs)
            {
              // check for the special case x = A * x
              if (viennacl::traits::handle(lhs) == viennacl::traits::handle(rhs.rhs()))
              {
                viennacl::vector<T> temp(lhs);
                viennacl::linalg::prod_impl(rhs.lhs(), rhs.rhs(), temp);
                lhs = temp;
              }
              else
                viennacl::linalg::prod_impl(rhs.lhs(), rhs.rhs(), lhs);
            }
        };

        template <typename T, unsigned int A>
        struct op_executor<vector_base<T>, op_inplace_add, vector_expression<const ell_matrix<T, A>, const vector_base<T>, op_prod> >
        {
            static void apply(vector_base<T> & lhs, vector_expression<const ell_matrix<T, A>, const vector_base<T>, op_prod> const & rhs)
            {
              viennacl::vector<T> temp(lhs);
              viennacl::linalg::prod_impl(rhs.lhs(), rhs.rhs(), temp);
              lhs += temp;
            }
        };

        template <typename T, unsigned int A>
        struct op_executor<vector_base<T>, op_inplace_sub, vector_expression<const ell_matrix<T, A>, const vector_base<T>, op_prod> >
        {
            static void apply(vector_base<T> & lhs, vector_expression<const ell_matrix<T, A>, const vector_base<T>, op_prod> const & rhs)
            {
              viennacl::vector<T> temp(lhs);
              viennacl::linalg::prod_impl(rhs.lhs(), rhs.rhs(), temp);
              lhs -= temp;
            }
        };


        // x = A * vec_op
        template <typename T, unsigned int A, typename LHS, typename RHS, typename OP>
        struct op_executor<vector_base<T>, op_assign, vector_expression<const ell_matrix<T, A>, const vector_expression<const LHS, const RHS, OP>, op_prod> >
        {
            static void apply(vector_base<T> & lhs, vector_expression<const ell_matrix<T, A>, const vector_expression<const LHS, const RHS, OP>, op_prod> const & rhs)
            {
              viennacl::vector<T> temp(rhs.rhs(), viennacl::traits::context(rhs));
              viennacl::linalg::prod_impl(rhs.lhs(), temp, lhs);
            }
        };

        // x = A * vec_op
        template <typename T, unsigned int A, typename LHS, typename RHS, typename OP>
        struct op_executor<vector_base<T>, op_inplace_add, vector_expression<const ell_matrix<T, A>, const vector_expression<const LHS, const RHS, OP>, op_prod> >
        {
            static void apply(vector_base<T> & lhs, vector_expression<const ell_matrix<T, A>, const vector_expression<const LHS, const RHS, OP>, op_prod> const & rhs)
            {
              viennacl::vector<T> temp(rhs.rhs(), viennacl::traits::context(rhs));
              viennacl::vector<T> temp_result(lhs);
              viennacl::linalg::prod_impl(rhs.lhs(), temp, temp_result);
              lhs += temp_result;
            }
        };

        // x = A * vec_op
        template <typename T, unsigned int A, typename LHS, typename RHS, typename OP>
        struct op_executor<vector_base<T>, op_inplace_sub, vector_expression<const ell_matrix<T, A>, const vector_expression<const LHS, const RHS, OP>, op_prod> >
        {
            static void apply(vector_base<T> & lhs, vector_expression<const ell_matrix<T, A>, const vector_expression<const LHS, const RHS, OP>, op_prod> const & rhs)
            {
              viennacl::vector<T> temp(rhs.rhs(), viennacl::traits::context(rhs));
              viennacl::vector<T> temp_result(lhs);
              viennacl::linalg::prod_impl(rhs.lhs(), temp, temp_result);
              lhs -= temp_result;
            }
        };

     } // namespace detail
   } // namespace linalg

    /** \endcond */
}

#endif


