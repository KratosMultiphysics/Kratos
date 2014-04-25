#ifndef VIENNACL_LINALG_OPENCL_KERNELS_MATRIX_SOLVE_HPP
#define VIENNACL_LINALG_OPENCL_KERNELS_MATRIX_SOLVE_HPP

#include "viennacl/tools/tools.hpp"
#include "viennacl/ocl/kernel.hpp"
#include "viennacl/ocl/platform.hpp"
#include "viennacl/ocl/utils.hpp"

#include "viennacl/linalg/opencl/kernels/matrix.hpp"

/** @file viennacl/linalg/opencl/kernels/matrix_solve.hpp
 *  @brief OpenCL kernel file for dense matrix solves with multiple right hand side (BLAS level 3) */
namespace viennacl
{
  namespace linalg
  {
    namespace opencl
    {
      namespace kernels
      {

        template <typename StringType>
        void generate_matrix_solve_blas3(StringType & source, std::string const & numeric_string,
                                         bool row_major_A, bool row_major_B,
                                         bool transpose_A, bool transpose_B,
                                         bool upper_solve, bool unit_diagonal)
        {
          //start OpenCL code:
          source.append("__kernel void ");
          if (transpose_A)
            source.append("trans_");
          if (unit_diagonal)
            source.append("unit_");
          if (upper_solve)
            source.append("upper_");
          else
            source.append("lower_");
          if (transpose_B)
            source.append("trans_");
          source.append("solve");

          source.append("( \n");
          source.append("          __global const "); source.append(numeric_string); source.append(" * A, \n");
          source.append("          unsigned int A_start1, unsigned int A_start2, \n");
          source.append("          unsigned int A_inc1,   unsigned int A_inc2, \n");
          source.append("          unsigned int A_size1,  unsigned int A_size2, \n");
          source.append("          unsigned int A_internal_size1, unsigned int A_internal_size2, \n");
          source.append("          __global "); source.append(numeric_string); source.append(" * B, \n");
          source.append("          unsigned int B_start1, unsigned int B_start2, \n");
          source.append("          unsigned int B_inc1,   unsigned int B_inc2, \n");
          source.append("          unsigned int B_size1,  unsigned int B_size2, \n");
          source.append("          unsigned int B_internal_size1, unsigned int B_internal_size2) { \n");
          source.append("  "); source.append(numeric_string); source.append(" temp;  \n");
          if (upper_solve)
          {
            //Note: A is square, thus A_rows == A_cols and no dispatch for transposedness needed
            source.append("  for (unsigned int row_cnt = 0; row_cnt < A_size1; ++row_cnt)  \n");
            source.append("  {  \n");
            source.append("    unsigned int row = A_size1 - 1 - row_cnt; \n");
          }
          else //lower triangular solve
          {
            source.append("  for (unsigned int row = 0; row < A_size1; ++row) \n");
            source.append("  { \n");
          }

          if (!unit_diagonal)
          {
            source.append("    barrier(CLK_GLOBAL_MEM_FENCE); \n");
            source.append("    if (get_local_id(0) == 0)  \n");
            //Note: A is square, thus A_internal_rows == A_internal_cols and no dispatch for transposedness needed
            if (row_major_B && transpose_B)
              source.append("      B[(get_group_id(0) * B_inc1 + B_start1) * B_internal_size2 + (row * B_inc2 + B_start2)] /= ");
            else if (row_major_B && !transpose_B)
              source.append("      B[(row * B_inc1 + B_start1) * B_internal_size2 + (get_group_id(0) * B_inc2 + B_start2)] /= ");
            else if (!row_major_B && transpose_B)
              source.append("      B[(get_group_id(0) * B_inc1 + B_start1) + (row * B_inc2 + B_start2) * B_internal_size1] /= ");
            else if (!row_major_B && !transpose_B)
              source.append("      B[(row * B_inc1 + B_start1) + (get_group_id(0) * B_inc2 + B_start2) * B_internal_size1] /= ");

            if (row_major_A)
              source.append("A[(row * A_inc1 + A_start1) * A_internal_size2 + (row * A_inc2 + A_start2)]; \n");
            else
              source.append("A[(row * A_inc1 + A_start1) + (row * A_inc2 + A_start2)*A_internal_size1]; \n");
          }

          source.append("    barrier(CLK_GLOBAL_MEM_FENCE); \n");

          if (row_major_B && transpose_B)
            source.append("    temp = B[(get_group_id(0) * B_inc1 + B_start1) * B_internal_size2 + (row * B_inc2 + B_start2)]; \n");
          else if (row_major_B && !transpose_B)
            source.append("    temp = B[(row * B_inc1 + B_start1) * B_internal_size2 + (get_group_id(0) * B_inc2 + B_start2)]; \n");
          else if (!row_major_B && transpose_B)
            source.append("    temp = B[(get_group_id(0) * B_inc1 + B_start1) + (row * B_inc2 + B_start2) * B_internal_size1]; \n");
          else if (!row_major_B && !transpose_B)
            source.append("    temp = B[(row * B_inc1 + B_start1) + (get_group_id(0) * B_inc2 + B_start2) * B_internal_size1]; \n");

          source.append("    //eliminate column of op(A) with index 'row' in parallel: \n");
          if (upper_solve)
            source.append("    for  (unsigned int elim = get_local_id(0); elim < row; elim += get_local_size(0)) \n");
          else
            source.append("    for  (unsigned int elim = row + get_local_id(0) + 1; elim < A_size1; elim += get_local_size(0)) \n");

          if (row_major_B && transpose_B)
            source.append("      B[(get_group_id(0) * B_inc1 + B_start1) * B_internal_size2 + (elim * B_inc2 + B_start2)] -= temp * ");
          else if (row_major_B && !transpose_B)
            source.append("      B[(elim * B_inc1 + B_start1) * B_internal_size2 + (get_group_id(0) * B_inc2 + B_start2)] -= temp * ");
          else if (!row_major_B && transpose_B)
            source.append("      B[(get_group_id(0) * B_inc1 + B_start1) + (elim * B_inc2 + B_start2) * B_internal_size1] -= temp * ");
          else if (!row_major_B && !transpose_B)
            source.append("      B[(elim * B_inc1 + B_start1) + (get_group_id(0) * B_inc2 + B_start2) * B_internal_size1] -= temp * ");

          if (row_major_A && transpose_A)
            source.append("A[(row * A_inc1 + A_start1) * A_internal_size2 + (elim * A_inc2 + A_start2)]; \n");
          else if (row_major_A && !transpose_A)
            source.append("A[(elim * A_inc1 + A_start1) * A_internal_size2 + (row * A_inc2 + A_start2)]; \n");
          else if (!row_major_A && transpose_A)
            source.append("A[(row * A_inc1 + A_start1) + (elim * A_inc2 + A_start2) * A_internal_size1]; \n");
          else if (!row_major_A && !transpose_A)
            source.append("A[(elim * A_inc1 + A_start1) + (row * A_inc2 + A_start2) * A_internal_size1]; \n");

          source.append("   } \n");
          source.append("} \n");
        }


        // main kernel class
        /** @brief Main kernel class for the generation of matrix solve kernels.
          *
          * @param F1  Row/Column majority tag for the system matrix
          * @param F2  Row/Column majority tag for the right hand side matrix
          */
        template <class NumericT, typename F1, typename F2>
        struct matrix_solve
        {
          static std::string program_name()
          {
            return viennacl::ocl::type_to_string<NumericT>::apply() + "_matrix_solve_" + detail::type_to_string(F1()) + detail::type_to_string(F2());
          }

          static void init(viennacl::ocl::context & ctx)
          {
            viennacl::ocl::DOUBLE_PRECISION_CHECKER<NumericT>::apply(ctx);
            std::string numeric_string = viennacl::ocl::type_to_string<NumericT>::apply();
            bool matrix_row_major = viennacl::is_row_major<F1>::value;
            bool rhs_row_major    = viennacl::is_row_major<F2>::value;


            static std::map<cl_context, bool> init_done;
            if (!init_done[ctx.handle().get()])
            {
              std::string source;
              source.reserve(8192);

              viennacl::ocl::append_double_precision_pragma<NumericT>(ctx, source);

              // only generate for floating points (forces error for integers)
              if (numeric_string == "float" || numeric_string == "double")
              {
                generate_matrix_solve_blas3(source, numeric_string, matrix_row_major, rhs_row_major,
                                            false, false, false, false);
                generate_matrix_solve_blas3(source, numeric_string, matrix_row_major, rhs_row_major,
                                            false, false, false, true);
                generate_matrix_solve_blas3(source, numeric_string, matrix_row_major, rhs_row_major,
                                            false, false, true, false);
                generate_matrix_solve_blas3(source, numeric_string, matrix_row_major, rhs_row_major,
                                            false, false, true, true);

                generate_matrix_solve_blas3(source, numeric_string, matrix_row_major, rhs_row_major,
                                            false, true, false, false);
                generate_matrix_solve_blas3(source, numeric_string, matrix_row_major, rhs_row_major,
                                            false, true, false, true);
                generate_matrix_solve_blas3(source, numeric_string, matrix_row_major, rhs_row_major,
                                            false, true, true, false);
                generate_matrix_solve_blas3(source, numeric_string, matrix_row_major, rhs_row_major,
                                            false, true, true, true);

                generate_matrix_solve_blas3(source, numeric_string, matrix_row_major, rhs_row_major,
                                            true, false, false, false);
                generate_matrix_solve_blas3(source, numeric_string, matrix_row_major, rhs_row_major,
                                            true, false, false, true);
                generate_matrix_solve_blas3(source, numeric_string, matrix_row_major, rhs_row_major,
                                            true, false, true, false);
                generate_matrix_solve_blas3(source, numeric_string, matrix_row_major, rhs_row_major,
                                            true, false, true, true);

                generate_matrix_solve_blas3(source, numeric_string, matrix_row_major, rhs_row_major,
                                            true, true, false, false);
                generate_matrix_solve_blas3(source, numeric_string, matrix_row_major, rhs_row_major,
                                            true, true, false, true);
                generate_matrix_solve_blas3(source, numeric_string, matrix_row_major, rhs_row_major,
                                            true, true, true, false);
                generate_matrix_solve_blas3(source, numeric_string, matrix_row_major, rhs_row_major,
                                            true, true, true, true);
              }

              std::string prog_name = program_name();
              #ifdef VIENNACL_BUILD_INFO
              std::cout << "Creating program " << prog_name << std::endl;
              #endif
              ctx.add_program(source, prog_name);
              init_done[ctx.handle().get()] = true;
            } //if
          } //init
        };

      }  // namespace kernels
    }  // namespace opencl
  }  // namespace linalg
}  // namespace viennacl
#endif

