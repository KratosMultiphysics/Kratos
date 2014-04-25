#ifndef VIENNACL_LINALG_OPENCL_KERNELS_MATRIX_PROD_HPP
#define VIENNACL_LINALG_OPENCL_KERNELS_MATRIX_PROD_HPP

#include "viennacl/tools/tools.hpp"
#include "viennacl/ocl/kernel.hpp"
#include "viennacl/ocl/platform.hpp"
#include "viennacl/ocl/utils.hpp"

#include "viennacl/linalg/opencl/kernels/matrix.hpp"

/** @file viennacl/linalg/opencl/kernels/matrix_prod.hpp
 *  @brief Runtime generation of OpenCL kernels for dense matrix-matrix products */
namespace viennacl
{
  namespace linalg
  {
    namespace opencl
    {
      namespace kernels
      {

        template <typename StringType>
        void generate_matrix_prod_blas3(StringType & source, std::string const & numeric_string,
                                        bool row_major_A, bool row_major_B, bool row_major_C,
                                        bool transpose_A, bool transpose_B)
        {
          //start OpenCL code:
          source.append("__kernel void prod_");
          if (transpose_A)
            source.append("T");
          else
            source.append("A");
          if (transpose_B)
            source.append("T");
          else
            source.append("A");

          source.append("( \n");
          source.append("  "); source.append(numeric_string); source.append(" alpha, \n");
          source.append("  __global const "); source.append(numeric_string); source.append(" * A, \n");
          source.append("  unsigned int A_row_start, \n");
          source.append("  unsigned int A_col_start, \n");
          source.append("  unsigned int A_row_inc, \n");
          source.append("  unsigned int A_col_inc, \n");
          source.append("  unsigned int A_row_size, \n");   //number of elements starting from row_start!
          source.append("  unsigned int A_col_size, \n");
          source.append("  unsigned int A_internal_rows, \n");
          source.append("  unsigned int A_internal_cols, \n");

          source.append("  __global const "); source.append(numeric_string); source.append(" * B,   \n");
          source.append("  unsigned int B_row_start, \n");
          source.append("  unsigned int B_col_start, \n");
          source.append("  unsigned int B_row_inc, \n");
          source.append("  unsigned int B_col_inc, \n");
          source.append("  unsigned int B_row_size, \n");
          source.append("  unsigned int B_col_size, \n");
          source.append("  unsigned int B_internal_rows, \n");
          source.append("  unsigned int B_internal_cols, \n");

          source.append("  "); source.append(numeric_string); source.append(" beta, \n");
          source.append("  __global "); source.append(numeric_string); source.append(" * C, \n");
          source.append("  unsigned int C_row_start, \n");
          source.append("  unsigned int C_col_start, \n");
          source.append("  unsigned int C_row_inc, \n");
          source.append("  unsigned int C_col_inc, \n");
          source.append("  unsigned int C_row_size, \n");
          source.append("  unsigned int C_col_size, \n");
          source.append("  unsigned int C_internal_rows, \n");
          source.append("  unsigned int C_internal_cols)  \n");
          source.append("{  \n");

          source.append("  __local "); source.append(numeric_string); source.append(" bufA[272]; \n"); // 16 * 17
          source.append("  __local "); source.append(numeric_string); source.append(" bufB[272]; \n"); // 16 * 17

          source.append("  size_t block_size = 16; \n"); //get_local_size(0);

          source.append("  size_t row_block_id = get_group_id(0); \n");
          source.append("  size_t col_block_id = get_group_id(1); \n");
          source.append("  size_t row_thread_id = get_local_id(0); \n");
          source.append("  size_t col_thread_id = get_local_id(1); \n");

          //traverse block row of A (taking mem layout and transpose operation into account)
          if (row_major_A && transpose_A)
          {
            source.append("  size_t aBegin = (row_block_id * block_size * A_col_inc + A_col_start) + A_row_start * A_internal_cols; \n");
            source.append("  size_t aStep = block_size * A_row_inc * A_internal_cols; \n");
          }
          else if (row_major_A && !transpose_A)
          {
            source.append("  size_t aBegin = (row_block_id * block_size * A_row_inc + A_row_start) * A_internal_cols + A_col_start; \n");
            source.append("  size_t aStep = block_size * A_col_inc; \n");
          }
          else if (!row_major_A && transpose_A)
          {
            source.append("  size_t aBegin = (row_block_id * block_size * A_col_inc + A_col_start) * A_internal_rows + A_row_start; \n");
            source.append("  size_t aStep = block_size * A_row_inc; \n");
          }
          else if (!row_major_A && !transpose_A)
          {
            source.append("  size_t aBegin = (row_block_id * block_size * A_row_inc + A_row_start) + A_col_start * A_internal_rows; \n");
            source.append("  size_t aStep = block_size * A_col_inc * A_internal_rows; \n");
          }


          if (row_major_B && transpose_B)
          {
            source.append("  size_t bBegin = (col_block_id * block_size * B_row_inc + B_row_start) * B_internal_cols + B_col_start; \n");
            source.append("  size_t bStep = block_size * B_col_inc; \n");
          }
          else if (row_major_B && !transpose_B)
          {
            source.append("  size_t bBegin = (col_block_id * block_size * B_col_inc + B_col_start) + B_row_start * B_internal_cols; \n");
            source.append("  size_t bStep = block_size * B_internal_cols * B_row_inc; \n");
          }
          else if (!row_major_B && transpose_B)
          {
            source.append("  size_t bBegin = (col_block_id * block_size * B_row_inc + B_row_start) + B_col_start * B_internal_rows; \n");
            source.append("  size_t bStep = block_size * B_internal_rows * B_col_inc; \n");
          }
          else if (!row_major_B && !transpose_B)
          {
            source.append("  size_t bBegin = (col_block_id * block_size * B_col_inc + B_col_start) * B_internal_rows + B_row_start; \n");
            source.append("  size_t bStep = block_size * B_row_inc; \n");
          }


          if (transpose_A)
            source.append("  size_t block_num = (A_row_size + block_size - 1) / block_size; \n");
          else
            source.append("  size_t block_num = (A_col_size + block_size - 1) / block_size; \n");

          source.append("  "); source.append(numeric_string); source.append(" Csub = 0; \n");

          //offset of the the memory access by the thread relative to the beginning of the block:
          if (row_major_A)
            source.append("  size_t aOffset = row_thread_id * A_col_inc + col_thread_id * A_row_inc * A_internal_cols; \n");
          else
            source.append("  size_t aOffset = row_thread_id * A_row_inc + col_thread_id * A_col_inc * A_internal_rows; \n");

          if (row_major_B)
            source.append("  size_t bOffset = row_thread_id * B_col_inc + col_thread_id * B_row_inc * B_internal_cols; \n");
          else
            source.append("  size_t bOffset = row_thread_id * B_row_inc + col_thread_id * B_col_inc *  B_internal_rows; \n");

          source.append("  size_t row_thread_id_times_block_size = row_thread_id * (block_size + 1); \n");
          source.append("  size_t col_thread_id_times_block_size = col_thread_id * (block_size + 1); \n");

          source.append("  for (size_t block = 0; \n");
          source.append("           block < block_num; \n");
          source.append("           ++block) \n");
          source.append("  { \n");

          //read block from A and check for access within matrix:

          if (transpose_A && row_major_A)
            source.append("    bufA[row_thread_id_times_block_size + col_thread_id] = ((block * block_size + col_thread_id < A_row_size) && (row_block_id * block_size + row_thread_id < A_col_size)) ? A[aBegin + aOffset] : 0; \n");
          else if (transpose_A && !row_major_A)
            source.append("    bufA[col_thread_id_times_block_size + row_thread_id] = ((block * block_size + row_thread_id < A_row_size) && (row_block_id * block_size + col_thread_id < A_col_size)) ? A[aBegin + aOffset] : 0; \n");
          else if (!transpose_A && row_major_A)
            source.append("    bufA[col_thread_id_times_block_size + row_thread_id] = ((block * block_size + row_thread_id < A_col_size) && (row_block_id * block_size + col_thread_id < A_row_size)) ? A[aBegin + aOffset] : 0; \n");
          else if (!transpose_A && !row_major_A)
            source.append("    bufA[row_thread_id_times_block_size + col_thread_id] = ((block * block_size + col_thread_id < A_col_size) && (row_block_id * block_size + row_thread_id < A_row_size)) ? A[aBegin + aOffset] : 0; \n");


          if (transpose_B && row_major_B)
            source.append("    bufB[col_thread_id_times_block_size + row_thread_id] = ((block * block_size + row_thread_id < B_col_size) && (col_block_id * block_size + col_thread_id < B_row_size)) ? B[bBegin + bOffset] : 0; \n");
          else if (transpose_B && !row_major_B)
            source.append("    bufB[row_thread_id_times_block_size + col_thread_id] = ((block * block_size + col_thread_id < B_col_size) && (col_block_id * block_size + row_thread_id < B_row_size)) ? B[bBegin + bOffset] : 0; \n");
          else if (!transpose_B && row_major_B)
            source.append("    bufB[row_thread_id_times_block_size + col_thread_id] = ((block * block_size + col_thread_id < B_row_size) && (col_block_id * block_size + row_thread_id < B_col_size)) ? B[bBegin + bOffset] : 0; \n");
          else if (!transpose_B && !row_major_B)
            source.append("    bufB[col_thread_id_times_block_size + row_thread_id] = ((block * block_size + row_thread_id < B_row_size) && (col_block_id * block_size + col_thread_id < B_col_size)) ? B[bBegin + bOffset] : 0; \n");

          //computation of block-matrix-matrix product is the same for all cases:
          source.append("    barrier(CLK_LOCAL_MEM_FENCE); \n");

          //loop unrolling:
          source.append("    __local "); source.append(numeric_string); source.append(" * bufAptr = bufA + row_thread_id_times_block_size; \n");
          source.append("    __local "); source.append(numeric_string); source.append(" * bufBptr = bufB + col_thread_id_times_block_size; \n");

          for (size_t unroll = 0; unroll < 16; ++unroll) {
            source.append("      Csub += (*bufAptr) * (*bufBptr); ++bufAptr; ++bufBptr; \n");
          }

          source.append("    barrier(CLK_LOCAL_MEM_FENCE); \n");
          source.append("    aBegin += aStep; \n");
          source.append("    bBegin += bStep; \n");
          source.append("  } \n");


          if (transpose_A)
          {
            source.append("  if (get_global_id(0) < A_col_size && ");
          }
          else
          {
            source.append("  if (get_global_id(0) < A_row_size && ");
          }

          if (transpose_B)
          {
            source.append("get_global_id(1) < B_row_size) \n");
          }
          else
          {
            source.append("get_global_id(1) < B_col_size) \n");
          }

          if (row_major_C)
          {
            source.append("    C[(get_global_id(0) * C_row_inc + C_row_start) * C_internal_cols + get_global_id(1) * C_col_inc + C_col_start] = (beta == 0) ? alpha * Csub : alpha * Csub + beta * C[(get_global_id(0) * C_row_inc + C_row_start) * C_internal_cols + get_global_id(1) * C_col_inc + C_col_start]; \n");
          }
          else
          {
            source.append("    C[get_global_id(0) * C_row_inc + C_row_start + (get_global_id(1) * C_col_inc + C_col_start) * C_internal_rows] = (beta == 0) ? alpha * Csub : alpha * Csub + beta * C[get_global_id(0) * C_row_inc + C_row_start + (get_global_id(1) * C_col_inc + C_col_start) * C_internal_rows]; \n");
          }
          source.append("} \n");
        }

        template <typename StringType>
        void generate_matrix_prod16_blas3(StringType & source, std::string const & numeric_string,
                                        bool row_major_A, bool row_major_B, bool row_major_C,
                                        bool transpose_A, bool transpose_B)
        {
          //vcl_size_t vector_size =  4;
          vcl_size_t block_size  = 16;

          //start OpenCL code:
          source.append("__kernel void prod16_");
          if (transpose_A)
            source.append("T");
          else
            source.append("A");
          if (transpose_B)
            source.append("T");
          else
            source.append("A");

          source.append("( "); source.append(numeric_string); source.append(" alpha, \n");
          source.append("   __global const "); source.append(numeric_string); source.append(" * A, \n");
          source.append("   unsigned int A_row_start, \n");
          source.append("   unsigned int A_col_start, \n");
          source.append("   unsigned int A_row_inc, \n");
          source.append("   unsigned int A_col_inc, \n");
          source.append("   unsigned int A_row_size, \n");   //number of elements starting from row_start, using an increment of A_row_inc
          source.append("   unsigned int A_col_size, \n");
          source.append("   unsigned int A_internal_rows, \n");
          source.append("   unsigned int A_internal_cols, \n");
          source.append("   __global const "); source.append(numeric_string); source.append(" * B,   \n");
          source.append("   unsigned int B_row_start, \n");
          source.append("   unsigned int B_col_start, \n");
          source.append("   unsigned int B_row_inc, \n");
          source.append("   unsigned int B_col_inc, \n");
          source.append("   unsigned int B_row_size, \n");
          source.append("   unsigned int B_col_size, \n");
          source.append("   unsigned int B_internal_rows, \n");
          source.append("   unsigned int B_internal_cols, \n");
          source.append("   "); source.append(numeric_string); source.append(" beta, \n");
          source.append("   __global "); source.append(numeric_string); source.append(" * C, \n");
          source.append("   unsigned int C_row_start, \n");
          source.append("   unsigned int C_col_start, \n");
          source.append("   unsigned int C_row_inc, \n");
          source.append("   unsigned int C_col_inc, \n");
          source.append("   unsigned int C_row_size, \n");
          source.append("   unsigned int C_col_size, \n");
          source.append("   unsigned int C_internal_rows, \n");
          source.append("   unsigned int C_internal_cols)  \n");
          source.append("{  \n");
          //do not forgot to change block_size !!!
          source.append("  size_t row_block_id = get_group_id(1); \n");    //refers to the row index in op(A), op(B)
          source.append("  size_t col_block_id = get_group_id(0); \n");    //refers to the col index in op(A), op(B)
          source.append("  size_t row_thread_id = get_local_id(1); \n");
          source.append("  size_t col_thread_id = get_local_id(0); \n");

          source.append("  __local "); source.append(numeric_string); source.append(" As[256]; \n");

          source.append("  "); source.append(numeric_string); source.append(" cv[16] = {");
          for (vcl_size_t i=0; i<block_size-1; ++i)
            source.append("0,");
          source.append("0}; \n");

          //traverse block row of A (taking mem layout and transpose operation into account)
          if (row_major_A && transpose_A)
          {
            source.append("  size_t aBegin = (row_block_id * 16 * A_col_inc + A_col_start) + A_row_start * A_internal_cols; \n");
            source.append("  size_t aStep = 16 * A_internal_cols * A_row_inc; \n");
            source.append("  size_t aEnd = aBegin + A_internal_cols * A_row_inc * A_row_size; \n");
          }
          else if (row_major_A && !transpose_A)
          {
            source.append("  size_t aBegin = (row_block_id * 16 * A_row_inc + A_row_start) * A_internal_cols + A_col_start; \n");
            source.append("  size_t aStep = 16 * A_col_inc; \n");
            source.append("  size_t aEnd = aBegin + A_col_inc * A_col_size; \n");
          }
          else if (!row_major_A && transpose_A)
          {
            source.append("  size_t aBegin = (row_block_id * 16 * A_col_inc + A_col_start) * A_internal_rows + A_row_start; \n");
            source.append("  size_t aStep = 16 * A_row_inc; \n");
            source.append("  size_t aEnd = aBegin + A_row_inc * A_row_size; \n");
          }
          else if (!row_major_A && !transpose_A)
          {
            source.append("  size_t aBegin = (row_block_id * 16 * A_row_inc + A_row_start) + A_col_start * A_internal_rows; \n");
            source.append("  size_t aStep = 16 * A_internal_rows * A_col_inc; \n");
            source.append("  size_t aEnd = aBegin + A_internal_rows * A_col_inc * A_col_size; \n");
          }


          if (row_major_B && transpose_B)
          {
            source.append("  size_t bBegin = (col_block_id * 64 * B_row_inc + B_row_start) * B_internal_cols + B_col_start; \n");
            source.append("  size_t bStep = 16 * B_col_inc; \n");
          }
          else if (row_major_B && !transpose_B)
          {
            source.append("  size_t bBegin = (col_block_id * 64 * B_col_inc + B_col_start) + B_row_start * B_internal_cols; \n");
            source.append("  size_t bStep = 16 * B_row_inc * B_internal_cols; \n");
          }
          else if (!row_major_B && transpose_B)
          {
            source.append("  size_t bBegin = (col_block_id * 64 * B_row_inc + B_row_start) + B_col_start * B_internal_rows; \n");
            source.append("  size_t bStep = 16 * B_col_inc * B_internal_rows; \n");
          }
          else if (!row_major_B && !transpose_B)
          {
            source.append("  size_t bBegin = (col_block_id * 64 * B_col_inc + B_col_start) * B_internal_rows + B_row_start; \n");
            source.append("  size_t bStep = 16 * B_row_inc; \n");
          }

          source.append("  for(size_t a = aBegin, b = bBegin; a < aEnd; a += aStep, b += bStep) {  \n");

          // copy blocks of op(A) to shared memory (op(A) is column-major in shared memory then)
          source.append("    for(size_t i = 0; i < 4; i++)   \n");
          if (row_major_A && transpose_A)
            source.append("      As[ (i*4 + row_thread_id) + 16 * col_thread_id] = (A[a + A_col_inc * (i * 4 + row_thread_id) + A_internal_cols * A_row_inc * col_thread_id]);");
          else if (row_major_A && !transpose_A)
            source.append("      As[ (i*4 + row_thread_id) + 16 * col_thread_id] = (A[a + A_internal_cols * A_row_inc * (i * 4 + row_thread_id) + A_col_inc * col_thread_id]);");
          else if (!row_major_A && transpose_A)
            source.append("      As[ (i*4 + row_thread_id) + 16 * col_thread_id] = (A[a + A_internal_rows * A_col_inc * (i * 4 + row_thread_id) + A_row_inc * col_thread_id]);");
          else if (!row_major_A && !transpose_A)
            source.append("      As[ (i*4 + row_thread_id) + 16 * col_thread_id] = (A[a + A_row_inc * (i * 4 + row_thread_id) + A_internal_rows * A_col_inc * col_thread_id]);");

          source.append("    barrier(CLK_LOCAL_MEM_FENCE);  \n");

          // initialize memory pointers
          source.append("    __local  "); source.append(numeric_string); source.append(" *ap = As;  \n");
          if (row_major_B && transpose_B)
          {
            source.append("    __global const "); source.append(numeric_string); source.append(" *bp = B + (b + (16 * row_thread_id + col_thread_id) * B_row_inc * B_internal_cols);  \n");
          }
          else if (row_major_B && !transpose_B)
          {
            source.append("    __global const "); source.append(numeric_string); source.append(" *bp = B + (b + (16 * row_thread_id + col_thread_id) * B_col_inc);  \n");
          }
          else if (!row_major_B && transpose_B)
          {
            source.append("    __global const "); source.append(numeric_string); source.append(" *bp = B + (b + (16 * row_thread_id + col_thread_id) * B_row_inc);  \n");
          }
          else if (!row_major_B && !transpose_B)
          {
            source.append("    __global const "); source.append(numeric_string); source.append(" *bp = B + (b + (16 * row_thread_id + col_thread_id) * B_col_inc * B_internal_rows);  \n");
          }

          // run computations
          source.append("    for(size_t i = 0; i < 16; i++) {  \n");
          if (row_major_B && transpose_B)
          {
            source.append("      "); source.append(numeric_string); source.append(" bv = bp[i * B_col_inc];  \n");
          }
          else if (row_major_B && !transpose_B)
          {
            source.append("      "); source.append(numeric_string); source.append(" bv = bp[i * B_row_inc * B_internal_cols];  \n");
          }
          else if (!row_major_B && transpose_B)
          {
            source.append("      "); source.append(numeric_string); source.append(" bv = bp[i * B_col_inc * B_internal_rows];  \n");
          }
          else if (!row_major_B && !transpose_B)
          {
            source.append("      "); source.append(numeric_string); source.append(" bv = bp[i * B_row_inc];  \n");
          }

          source.append("      for(size_t k = 0; k < 16; k++)   \n");
          source.append("	    cv[k] += ap[k] * bv;  \n");

          source.append("      ap += 16;  \n");
          source.append("    }  \n");

          source.append("    barrier(CLK_LOCAL_MEM_FENCE);  \n");
          source.append("  }  \n");

          // write to C
          if (row_major_C)
          {
              source.append("  int c = C_internal_cols * (C_row_inc * 16 * row_block_id + C_row_start) + 64 * C_col_inc * col_block_id + C_col_start  \n");  //block column index
              source.append("          + C_col_inc * (16 * row_thread_id + col_thread_id);  \n");
          }
          else
          {
              source.append("  int c = C_row_inc * 16 * row_block_id + C_row_start + (64 * C_col_inc * col_block_id + C_col_start) * C_internal_rows  \n");   // block column index
              source.append("          + C_internal_rows * C_col_inc * (16 * row_thread_id + col_thread_id);  \n");
          }

          source.append("  for(size_t i = 0; i < 16; i++) {  \n");

          if (row_major_C)
          {
            source.append("    C[c] = (beta == 0) ? alpha * cv[i] : alpha * cv[i] + beta * C[c];  \n");
            source.append("      c += C_internal_cols * C_row_inc;  \n");
          }
          else
          {
            source.append("    C[c] = (beta == 0) ? alpha * cv[i] : alpha * cv[i] + beta * C[c];  \n");
            source.append("      c += C_row_inc;  \n");
          }

          source.append("  }  \n");
          source.append("}  \n");

        }


        // main kernel class
        /** @brief Main kernel class for the generation of matrix-matrix product kernels C = A * B
          *
          * @param F_A  Row/Column majority tag for A
          * @param F_B  Row/Column majority tag for B
          * @param F_C  Row/Column majority tag for C
          */
        template <class NumericT, typename F_A, typename F_B, typename F_C>
        struct matrix_prod
        {
          static std::string program_name()
          {
            return viennacl::ocl::type_to_string<NumericT>::apply() + "_matrix_prod_" + detail::type_to_string(F_A()) + detail::type_to_string(F_B()) + detail::type_to_string(F_C());
          }

          static void init(viennacl::ocl::context & ctx)
          {
            viennacl::ocl::DOUBLE_PRECISION_CHECKER<NumericT>::apply(ctx);
            std::string numeric_string = viennacl::ocl::type_to_string<NumericT>::apply();
            bool row_major_A = viennacl::is_row_major<F_A>::value;
            bool row_major_B = viennacl::is_row_major<F_B>::value;
            bool row_major_C = viennacl::is_row_major<F_C>::value;


            static std::map<cl_context, bool> init_done;
            if (!init_done[ctx.handle().get()])
            {
              std::string source;
              source.reserve(8192);

              viennacl::ocl::append_double_precision_pragma<NumericT>(ctx, source);

              // only generate for floating points (forces error for integers)
              if (numeric_string == "float" || numeric_string == "double")
              {
                generate_matrix_prod_blas3(source, numeric_string, row_major_A, row_major_B, row_major_C, false, false);
                generate_matrix_prod_blas3(source, numeric_string, row_major_A, row_major_B, row_major_C, false, true);
                generate_matrix_prod_blas3(source, numeric_string, row_major_A, row_major_B, row_major_C, true, false);
                generate_matrix_prod_blas3(source, numeric_string, row_major_A, row_major_B, row_major_C, true, true);

                generate_matrix_prod16_blas3(source, numeric_string, row_major_A, row_major_B, row_major_C, false, false);
                generate_matrix_prod16_blas3(source, numeric_string, row_major_A, row_major_B, row_major_C, false, true);
                generate_matrix_prod16_blas3(source, numeric_string, row_major_A, row_major_B, row_major_C, true, false);
                generate_matrix_prod16_blas3(source, numeric_string, row_major_A, row_major_B, row_major_C, true, true);

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

