#ifndef VIENNACL_LINALG_OPENCL_KERNELS_COMPRESSED_MATRIX_HPP
#define VIENNACL_LINALG_OPENCL_KERNELS_COMPRESSED_MATRIX_HPP

#include "viennacl/tools/tools.hpp"
#include "viennacl/ocl/kernel.hpp"
#include "viennacl/ocl/platform.hpp"
#include "viennacl/ocl/utils.hpp"

#include "viennacl/linalg/opencl/common.hpp"

/** @file viennacl/linalg/opencl/kernels/compressed_matrix.hpp
 *  @brief OpenCL kernel file for compressed_matrix operations */
namespace viennacl
{
  namespace linalg
  {
    namespace opencl
    {
      namespace kernels
      {

        //////////////////////////// Part 1: Kernel generation routines ////////////////////////////////////

        template <typename StringType>
        void generate_compressed_matrix_block_trans_lu_backward(StringType & source, std::string const & numeric_string)
        {
          source.append("__kernel void block_trans_lu_backward( \n");
          source.append("          __global const unsigned int * row_jumper_U,  \n");     //U part (note that U is transposed in memory)
          source.append("          __global const unsigned int * column_indices_U, \n");
          source.append("          __global const "); source.append(numeric_string); source.append(" * elements_U, \n");
          source.append("          __global const "); source.append(numeric_string); source.append(" * diagonal_U, \n");
          source.append("          __global const unsigned int * block_offsets, \n");
          source.append("          __global "); source.append(numeric_string); source.append(" * result, \n");
          source.append("          unsigned int size) \n");
          source.append("{ \n");
          source.append("  unsigned int col_start = block_offsets[2*get_group_id(0)]; \n");
          source.append("  unsigned int col_stop  = block_offsets[2*get_group_id(0)+1]; \n");
          source.append("  unsigned int row_start; \n");
          source.append("  unsigned int row_stop; \n");
          source.append("  "); source.append(numeric_string); source.append(" result_entry = 0; \n");

          source.append("  if (col_start >= col_stop) \n");
          source.append("    return; \n");

            //backward elimination, using U and diagonal_U
          source.append("  for (unsigned int iter = 0; iter < col_stop - col_start; ++iter) \n");
          source.append("  { \n");
          source.append("    unsigned int col = (col_stop - iter) - 1; \n");
          source.append("    result_entry = result[col] / diagonal_U[col]; \n");
          source.append("    row_start = row_jumper_U[col]; \n");
          source.append("    row_stop  = row_jumper_U[col + 1]; \n");
          source.append("    for (unsigned int buffer_index = row_start + get_local_id(0); buffer_index < row_stop; buffer_index += get_local_size(0)) \n");
          source.append("      result[column_indices_U[buffer_index]] -= result_entry * elements_U[buffer_index]; \n");
          source.append("    barrier(CLK_GLOBAL_MEM_FENCE); \n");
          source.append("  } \n");

            //divide result vector by diagonal:
          source.append("  for (unsigned int col = col_start + get_local_id(0); col < col_stop; col += get_local_size(0)) \n");
          source.append("    result[col] /= diagonal_U[col]; \n");
          source.append("} \n");
        }

        template <typename StringType>
        void generate_compressed_matrix_block_trans_unit_lu_forward(StringType & source, std::string const & numeric_string)
        {
          source.append("__kernel void block_trans_unit_lu_forward( \n");
          source.append("          __global const unsigned int * row_jumper_L,  \n");     //L part (note that L is transposed in memory)
          source.append("          __global const unsigned int * column_indices_L, \n");
          source.append("          __global const "); source.append(numeric_string); source.append(" * elements_L, \n");
          source.append("          __global const unsigned int * block_offsets, \n");
          source.append("          __global "); source.append(numeric_string); source.append(" * result, \n");
          source.append("          unsigned int size) \n");
          source.append("{ \n");
          source.append("  unsigned int col_start = block_offsets[2*get_group_id(0)]; \n");
          source.append("  unsigned int col_stop  = block_offsets[2*get_group_id(0)+1]; \n");
          source.append("  unsigned int row_start = row_jumper_L[col_start]; \n");
          source.append("  unsigned int row_stop; \n");
          source.append("  "); source.append(numeric_string); source.append(" result_entry = 0; \n");

          source.append("  if (col_start >= col_stop) \n");
          source.append("    return; \n");

            //forward elimination, using L:
          source.append("  for (unsigned int col = col_start; col < col_stop; ++col) \n");
          source.append("  { \n");
          source.append("    result_entry = result[col]; \n");
          source.append("    row_stop = row_jumper_L[col + 1]; \n");
          source.append("    for (unsigned int buffer_index = row_start + get_local_id(0); buffer_index < row_stop; buffer_index += get_local_size(0)) \n");
          source.append("      result[column_indices_L[buffer_index]] -= result_entry * elements_L[buffer_index]; \n");
          source.append("    row_start = row_stop; \n"); //for next iteration (avoid unnecessary loads from GPU RAM)
          source.append("    barrier(CLK_GLOBAL_MEM_FENCE); \n");
          source.append("  } \n");

          source.append("}; \n");
        }

        namespace detail
        {
          /** @brief Generate kernel for C = A * B with A being a compressed_matrix, B and C dense */
          template <typename StringType>
          void generate_compressed_matrix_dense_matrix_mult(StringType & source, std::string const & numeric_string,
                                                            bool B_transposed, bool B_row_major, bool C_row_major)
          {
            source.append("__kernel void ");
            source.append(viennacl::linalg::opencl::detail::sparse_dense_matmult_kernel_name(B_transposed, B_row_major, C_row_major));
            source.append("( \n");
            source.append("          __global const unsigned int * sp_mat_row_indices, \n");
            source.append("          __global const unsigned int * sp_mat_col_indices, \n");
            source.append("          __global const "); source.append(numeric_string); source.append(" * sp_mat_elements, \n");
            source.append("          __global const "); source.append(numeric_string); source.append(" * d_mat, \n");
            source.append("          unsigned int d_mat_row_start, \n");
            source.append("          unsigned int d_mat_col_start, \n");
            source.append("          unsigned int d_mat_row_inc, \n");
            source.append("          unsigned int d_mat_col_inc, \n");
            source.append("          unsigned int d_mat_row_size, \n");
            source.append("          unsigned int d_mat_col_size, \n");
            source.append("          unsigned int d_mat_internal_rows, \n");
            source.append("          unsigned int d_mat_internal_cols, \n");
            source.append("          __global "); source.append(numeric_string); source.append(" * result, \n");
            source.append("          unsigned int result_row_start, \n");
            source.append("          unsigned int result_col_start, \n");
            source.append("          unsigned int result_row_inc, \n");
            source.append("          unsigned int result_col_inc, \n");
            source.append("          unsigned int result_row_size, \n");
            source.append("          unsigned int result_col_size, \n");
            source.append("          unsigned int result_internal_rows, \n");
            source.append("          unsigned int result_internal_cols) { \n");

              // split work rows (sparse matrix rows) to thread groups
            source.append("  for (unsigned int row = get_group_id(0); row < result_row_size; row += get_num_groups(0)) { \n");

            source.append("    unsigned int row_start = sp_mat_row_indices[row]; \n");
            source.append("    unsigned int row_end = sp_mat_row_indices[row+1]; \n");

                // split result cols between threads in a thread group
            source.append("    for ( unsigned int col = get_local_id(0); col < result_col_size; col += get_local_size(0) ) { \n");

            source.append("      "); source.append(numeric_string); source.append(" r = 0; \n");

            source.append("      for (unsigned int k = row_start; k < row_end; k ++) { \n");

            source.append("        unsigned int j = sp_mat_col_indices[k]; \n");
            source.append("        "); source.append(numeric_string); source.append(" x = sp_mat_elements[k]; \n");

            source.append("        "); source.append(numeric_string);
            if (B_transposed && B_row_major)
              source.append(" y = d_mat[ (d_mat_row_start + col * d_mat_row_inc) * d_mat_internal_cols + d_mat_col_start +   j * d_mat_col_inc ]; \n");
            else if (B_transposed && !B_row_major)
              source.append(" y = d_mat[ (d_mat_row_start + col * d_mat_row_inc)                       + (d_mat_col_start +  j * d_mat_col_inc) * d_mat_internal_rows ]; \n");
            else if (!B_transposed && B_row_major)
              source.append(" y = d_mat[ (d_mat_row_start +   j * d_mat_row_inc) * d_mat_internal_cols + d_mat_col_start + col * d_mat_col_inc ]; \n");
            else
              source.append(" y = d_mat[ (d_mat_row_start +   j * d_mat_row_inc)                       + (d_mat_col_start + col * d_mat_col_inc) * d_mat_internal_rows ]; \n");
            source.append("        r += x * y; \n");
            source.append("      } \n");

            if (C_row_major)
              source.append("      result[ (result_row_start + row * result_row_inc) * result_internal_cols + result_col_start + col * result_col_inc ] = r; \n");
            else
              source.append("      result[ (result_row_start + row * result_row_inc)                        + (result_col_start + col * result_col_inc) * result_internal_rows ] = r; \n");
            source.append("    } \n");
            source.append("  } \n");

            source.append("} \n");

          }
        }
        template <typename StringType>
        void generate_compressed_matrix_dense_matrix_multiplication(StringType & source, std::string const & numeric_string)
        {
          detail::generate_compressed_matrix_dense_matrix_mult(source, numeric_string, false, false, false);
          detail::generate_compressed_matrix_dense_matrix_mult(source, numeric_string, false, false,  true);
          detail::generate_compressed_matrix_dense_matrix_mult(source, numeric_string, false,  true, false);
          detail::generate_compressed_matrix_dense_matrix_mult(source, numeric_string, false,  true,  true);

          detail::generate_compressed_matrix_dense_matrix_mult(source, numeric_string, true, false, false);
          detail::generate_compressed_matrix_dense_matrix_mult(source, numeric_string, true, false,  true);
          detail::generate_compressed_matrix_dense_matrix_mult(source, numeric_string, true,  true, false);
          detail::generate_compressed_matrix_dense_matrix_mult(source, numeric_string, true,  true,  true);
        }

        template <typename StringType>
        void generate_compressed_matrix_jacobi(StringType & source, std::string const & numeric_string)
        {

         source.append(" __kernel void jacobi( \n");
         source.append("  __global const unsigned int * row_indices, \n");
         source.append("  __global const unsigned int * column_indices, \n");
         source.append("  __global const "); source.append(numeric_string); source.append(" * elements, \n");
         source.append("  "); source.append(numeric_string); source.append(" weight, \n");
         source.append("  __global const "); source.append(numeric_string); source.append(" * old_result, \n");
         source.append("  __global "); source.append(numeric_string); source.append(" * new_result, \n");
         source.append("  __global const "); source.append(numeric_string); source.append(" * rhs, \n");
         source.append("  unsigned int size) \n");
         source.append("  { \n");
         source.append("   "); source.append(numeric_string); source.append(" sum, diag=1; \n");
         source.append("   int col; \n");
         source.append("   for (unsigned int i = get_global_id(0); i < size; i += get_global_size(0)) \n");
         source.append("   { \n");
         source.append("     sum = 0; \n");
         source.append("     for (unsigned int j = row_indices[i]; j<row_indices[i+1]; j++) \n");
         source.append("     { \n");
         source.append("       col = column_indices[j]; \n");
         source.append("       if (i == col) \n");
         source.append("   diag = elements[j]; \n");
         source.append("       else \n");
         source.append("   sum += elements[j] * old_result[col]; \n");
         source.append("     } \n");
         source.append("       new_result[i] = weight * (rhs[i]-sum) / diag + (1-weight) * old_result[i]; \n");
         source.append("    } \n");
         source.append("  } \n");

        }

        template <typename StringType>
        void generate_compressed_matrix_lu_backward(StringType & source, std::string const & numeric_string)
        {
          // compute x in Ux = y for incomplete LU factorizations of a sparse matrix in compressed format
          source.append("__kernel void lu_backward( \n");
          source.append("          __global const unsigned int * row_indices, \n");
          source.append("          __global const unsigned int * column_indices, \n");
          source.append("          __global const "); source.append(numeric_string); source.append(" * elements, \n");
          source.append("          __global "); source.append(numeric_string); source.append(" * vector, \n");
          source.append("          unsigned int size) \n");
          source.append("{ \n");
          source.append("  __local unsigned int col_index_buffer[128]; \n");
          source.append("  __local "); source.append(numeric_string); source.append(" element_buffer[128]; \n");
          source.append("  __local "); source.append(numeric_string); source.append(" vector_buffer[128]; \n");

          source.append("  unsigned int nnz = row_indices[size]; \n");
          source.append("  unsigned int current_row = size-1; \n");
          source.append("  unsigned int row_at_window_start = size-1; \n");
          source.append("  "); source.append(numeric_string); source.append(" current_vector_entry = vector[size-1]; \n");
          source.append("  "); source.append(numeric_string); source.append(" diagonal_entry = 0; \n");
          source.append("  unsigned int loop_end = ( (nnz - 1) / get_local_size(0)) * get_local_size(0); \n");
          source.append("  unsigned int next_row = row_indices[size-1]; \n");

          source.append("  unsigned int i = loop_end + get_local_id(0); \n");
          source.append("  while (1) \n");
          source.append("  { \n");
              //load into shared memory (coalesced access):
          source.append("    if (i < nnz) \n");
          source.append("    { \n");
          source.append("      element_buffer[get_local_id(0)] = elements[i]; \n");
          source.append("      unsigned int tmp = column_indices[i]; \n");
          source.append("      col_index_buffer[get_local_id(0)] = tmp; \n");
          source.append("      vector_buffer[get_local_id(0)] = vector[tmp]; \n");
          source.append("    } \n");

          source.append("    barrier(CLK_LOCAL_MEM_FENCE); \n");

              //now a single thread does the remaining work in shared memory:
          source.append("    if (get_local_id(0) == 0) \n");
          source.append("    { \n");
                // traverse through all the loaded data from back to front:
          source.append("      for (unsigned int k2=0; k2<get_local_size(0); ++k2) \n");
          source.append("      { \n");
          source.append("        unsigned int k = (get_local_size(0) - k2) - 1; \n");

          source.append("        if (i+k >= nnz) \n");
          source.append("          continue; \n");

          source.append("        if (col_index_buffer[k] > row_at_window_start) \n"); //use recently computed results
          source.append("          current_vector_entry -= element_buffer[k] * vector_buffer[k]; \n");
          source.append("        else if (col_index_buffer[k] > current_row) \n"); //use buffered data
          source.append("          current_vector_entry -= element_buffer[k] * vector[col_index_buffer[k]]; \n");
          source.append("        else if (col_index_buffer[k] == current_row) \n");
          source.append("          diagonal_entry = element_buffer[k]; \n");

          source.append("        if (i+k == next_row) \n"); //current row is finished. Write back result
          source.append("        { \n");
          source.append("          vector[current_row] = current_vector_entry / diagonal_entry; \n");
          source.append("          if (current_row > 0) //load next row's data \n");
          source.append("          { \n");
          source.append("            --current_row; \n");
          source.append("            next_row = row_indices[current_row]; \n");
          source.append("            current_vector_entry = vector[current_row]; \n");
          source.append("          } \n");
          source.append("        } \n");


          source.append("      } \n"); // for k

          source.append("      row_at_window_start = current_row; \n");
          source.append("    } \n"); // if (get_local_id(0) == 0)

          source.append("    barrier(CLK_GLOBAL_MEM_FENCE); \n");

          source.append("    if (i < get_local_size(0)) \n");
          source.append("      break; \n");

          source.append("    i -= get_local_size(0); \n");
          source.append("  } \n"); //for i
          source.append("} \n");

        }

        template <typename StringType>
        void generate_compressed_matrix_lu_forward(StringType & source, std::string const & numeric_string)
        {

          // compute y in Ly = z for incomplete LU factorizations of a sparse matrix in compressed format
          source.append("__kernel void lu_forward( \n");
          source.append("          __global const unsigned int * row_indices, \n");
          source.append("          __global const unsigned int * column_indices, \n");
          source.append("          __global const "); source.append(numeric_string); source.append(" * elements, \n");
          source.append("          __global "); source.append(numeric_string); source.append(" * vector, \n");
          source.append("          unsigned int size) \n");
          source.append("{ \n");
          source.append("  __local unsigned int col_index_buffer[128]; \n");
          source.append("  __local "); source.append(numeric_string); source.append(" element_buffer[128]; \n");
          source.append("  __local "); source.append(numeric_string); source.append(" vector_buffer[128]; \n");

          source.append("  unsigned int nnz = row_indices[size]; \n");
          source.append("  unsigned int current_row = 0; \n");
          source.append("  unsigned int row_at_window_start = 0; \n");
          source.append("  "); source.append(numeric_string); source.append(" current_vector_entry = vector[0]; \n");
          source.append("  "); source.append(numeric_string); source.append(" diagonal_entry; \n");
          source.append("  unsigned int loop_end = (nnz / get_local_size(0) + 1) * get_local_size(0); \n");
          source.append("  unsigned int next_row = row_indices[1]; \n");

          source.append("  for (unsigned int i = get_local_id(0); i < loop_end; i += get_local_size(0)) \n");
          source.append("  { \n");
              //load into shared memory (coalesced access):
          source.append("    if (i < nnz) \n");
          source.append("    { \n");
          source.append("      element_buffer[get_local_id(0)] = elements[i]; \n");
          source.append("      unsigned int tmp = column_indices[i]; \n");
          source.append("      col_index_buffer[get_local_id(0)] = tmp; \n");
          source.append("      vector_buffer[get_local_id(0)] = vector[tmp]; \n");
          source.append("    } \n");

          source.append("    barrier(CLK_LOCAL_MEM_FENCE); \n");

              //now a single thread does the remaining work in shared memory:
          source.append("    if (get_local_id(0) == 0) \n");
          source.append("    { \n");
                // traverse through all the loaded data:
          source.append("      for (unsigned int k=0; k<get_local_size(0); ++k) \n");
          source.append("      { \n");
          source.append("        if (current_row < size && i+k == next_row) \n"); //current row is finished. Write back result
          source.append("        { \n");
          source.append("          vector[current_row] = current_vector_entry / diagonal_entry; \n");
          source.append("          ++current_row; \n");
          source.append("          if (current_row < size) \n"); //load next row's data
          source.append("          { \n");
          source.append("            next_row = row_indices[current_row+1]; \n");
          source.append("            current_vector_entry = vector[current_row]; \n");
          source.append("          } \n");
          source.append("        } \n");

          source.append("        if (current_row < size && col_index_buffer[k] < current_row) \n"); //substitute
          source.append("        { \n");
          source.append("          if (col_index_buffer[k] < row_at_window_start) \n"); //use recently computed results
          source.append("            current_vector_entry -= element_buffer[k] * vector_buffer[k]; \n");
          source.append("          else if (col_index_buffer[k] < current_row) \n"); //use buffered data
          source.append("            current_vector_entry -= element_buffer[k] * vector[col_index_buffer[k]]; \n");
          source.append("        } \n");
          source.append("        else if (col_index_buffer[k] == current_row) \n");
          source.append("          diagonal_entry = element_buffer[k]; \n");

          source.append("      } \n"); // for k

          source.append("      row_at_window_start = current_row; \n");
          source.append("    } \n"); // if (get_local_id(0) == 0)

          source.append("    barrier(CLK_GLOBAL_MEM_FENCE); \n");
          source.append("  } \n"); //for i
          source.append("} \n");

        }

        template <typename StringType>
        void generate_compressed_matrix_row_info_extractor(StringType & source, std::string const & numeric_string)
        {
          source.append("__kernel void row_info_extractor( \n");
          source.append("          __global const unsigned int * row_indices, \n");
          source.append("          __global const unsigned int * column_indices, \n");
          source.append("          __global const "); source.append(numeric_string); source.append(" * elements, \n");
          source.append("          __global "); source.append(numeric_string); source.append(" * result, \n");
          source.append("          unsigned int size, \n");
          source.append("          unsigned int option \n");
          source.append("          ) \n");
          source.append("{ \n");
          source.append("  for (unsigned int row = get_global_id(0); row < size; row += get_global_size(0)) \n");
          source.append("  { \n");
          source.append("    "); source.append(numeric_string); source.append(" value = 0; \n");
          source.append("    unsigned int row_end = row_indices[row+1]; \n");

          source.append("    switch (option) \n");
          source.append("    { \n");
          source.append("      case 0: \n"); //inf-norm
          source.append("        for (unsigned int i = row_indices[row]; i < row_end; ++i) \n");
          source.append("          value = max(value, fabs(elements[i])); \n");
          source.append("        break; \n");

          source.append("      case 1: \n"); //1-norm
          source.append("        for (unsigned int i = row_indices[row]; i < row_end; ++i) \n");
          source.append("          value += fabs(elements[i]); \n");
          source.append("        break; \n");

          source.append("      case 2: \n"); //2-norm
          source.append("        for (unsigned int i = row_indices[row]; i < row_end; ++i) \n");
          source.append("          value += elements[i] * elements[i]; \n");
          source.append("        value = sqrt(value); \n");
          source.append("        break; \n");

          source.append("      case 3: \n"); //diagonal entry
          source.append("        for (unsigned int i = row_indices[row]; i < row_end; ++i) \n");
          source.append("        { \n");
          source.append("          if (column_indices[i] == row) \n");
          source.append("          { \n");
          source.append("            value = elements[i]; \n");
          source.append("            break; \n");
          source.append("          } \n");
          source.append("        } \n");
          source.append("        break; \n");

          source.append("      default: \n");
          source.append("        break; \n");
          source.append("    } \n");
          source.append("    result[row] = value; \n");
          source.append("  } \n");
          source.append("} \n");

        }

        template <typename StringType>
        void generate_compressed_matrix_trans_lu_backward(StringType & source, std::string const & numeric_string)
        {

          // compute y in Ly = z for incomplete LU factorizations of a sparse matrix in compressed format
          source.append("__kernel void trans_lu_backward( \n");
          source.append("          __global const unsigned int * row_indices, \n");
          source.append("          __global const unsigned int * column_indices, \n");
          source.append("          __global const "); source.append(numeric_string); source.append(" * elements, \n");
          source.append("          __global const "); source.append(numeric_string); source.append(" * diagonal_entries, \n");
          source.append("          __global "); source.append(numeric_string); source.append(" * vector, \n");
          source.append("          unsigned int size) \n");
          source.append("{ \n");
          source.append("  __local unsigned int row_index_lookahead[256]; \n");
          source.append("  __local unsigned int row_index_buffer[256]; \n");

          source.append("  unsigned int row_index; \n");
          source.append("  unsigned int col_index; \n");
          source.append("  "); source.append(numeric_string); source.append(" matrix_entry; \n");
          source.append("  unsigned int nnz = row_indices[size]; \n");
          source.append("  unsigned int row_at_window_start = size; \n");
          source.append("  unsigned int row_at_window_end; \n");
          source.append("  unsigned int loop_end = ( (nnz - 1) / get_local_size(0) + 1) * get_local_size(0); \n");

          source.append("  for (unsigned int i2 = get_local_id(0); i2 < loop_end; i2 += get_local_size(0)) \n");
          source.append("  { \n");
          source.append("    unsigned int i = (nnz - i2) - 1; \n");
          source.append("    col_index    = (i2 < nnz) ? column_indices[i] : 0; \n");
          source.append("    matrix_entry = (i2 < nnz) ? elements[i]       : 0; \n");
          source.append("    row_index_lookahead[get_local_id(0)] = (row_at_window_start >= get_local_id(0)) ? row_indices[row_at_window_start - get_local_id(0)] : 0; \n");

          source.append("    barrier(CLK_LOCAL_MEM_FENCE); \n");

          source.append("    if (i2 < nnz) \n");
          source.append("    { \n");
          source.append("      unsigned int row_index_dec = 0; \n");
          source.append("      while (row_index_lookahead[row_index_dec] > i) \n");
          source.append("        ++row_index_dec; \n");
          source.append("      row_index = row_at_window_start - row_index_dec; \n");
          source.append("      row_index_buffer[get_local_id(0)] = row_index; \n");
          source.append("    } \n");
          source.append("    else \n");
          source.append("    { \n");
          source.append("      row_index = size+1; \n");
          source.append("      row_index_buffer[get_local_id(0)] = 0; \n");
          source.append("    } \n");

          source.append("    barrier(CLK_LOCAL_MEM_FENCE); \n");

          source.append("    row_at_window_start = row_index_buffer[0]; \n");
          source.append("    row_at_window_end   = row_index_buffer[get_local_size(0) - 1]; \n");

              //backward elimination
          source.append("    for (unsigned int row2 = 0; row2 <= (row_at_window_start - row_at_window_end); ++row2) \n");
          source.append("    { \n");
          source.append("      unsigned int row = row_at_window_start - row2; \n");
          source.append("      "); source.append(numeric_string); source.append(" result_entry = vector[row] / diagonal_entries[row]; \n");

          source.append("      if ( (row_index == row) && (col_index < row) ) \n");
          source.append("        vector[col_index] -= result_entry * matrix_entry; \n");

          source.append("      barrier(CLK_GLOBAL_MEM_FENCE); \n");
          source.append("    } \n");

          source.append("    row_at_window_start = row_at_window_end; \n");
          source.append("  } \n");

            // final step: Divide vector by diagonal entries:
          source.append("  for (unsigned int i = get_local_id(0); i < size; i += get_local_size(0)) \n");
          source.append("    vector[i] /= diagonal_entries[i]; \n");
          source.append("} \n");

        }

        template <typename StringType>
        void generate_compressed_matrix_trans_lu_forward(StringType & source, std::string const & numeric_string)
        {

          // compute y in Ly = z for incomplete LU factorizations of a sparse matrix in compressed format
          source.append("__kernel void trans_lu_forward( \n");
          source.append("          __global const unsigned int * row_indices, \n");
          source.append("          __global const unsigned int * column_indices, \n");
          source.append("          __global const "); source.append(numeric_string); source.append(" * elements, \n");
          source.append("          __global const "); source.append(numeric_string); source.append(" * diagonal_entries, \n");
          source.append("          __global "); source.append(numeric_string); source.append(" * vector, \n");
          source.append("          unsigned int size) \n");
          source.append("{ \n");
          source.append("  __local unsigned int row_index_lookahead[256]; \n");
          source.append("  __local unsigned int row_index_buffer[256]; \n");

          source.append("  unsigned int row_index; \n");
          source.append("  unsigned int col_index; \n");
          source.append("  "); source.append(numeric_string); source.append(" matrix_entry; \n");
          source.append("  unsigned int nnz = row_indices[size]; \n");
          source.append("  unsigned int row_at_window_start = 0; \n");
          source.append("  unsigned int row_at_window_end = 0; \n");
          source.append("  unsigned int loop_end = ( (nnz - 1) / get_local_size(0) + 1) * get_local_size(0); \n");

          source.append("  for (unsigned int i = get_local_id(0); i < loop_end; i += get_local_size(0)) \n");
          source.append("  { \n");
          source.append("    col_index    = (i < nnz) ? column_indices[i] : 0; \n");
          source.append("    matrix_entry = (i < nnz) ? elements[i]       : 0; \n");
          source.append("    row_index_lookahead[get_local_id(0)] = (row_at_window_start + get_local_id(0) < size) ? row_indices[row_at_window_start + get_local_id(0)] : size - 1; \n");

          source.append("    barrier(CLK_LOCAL_MEM_FENCE); \n");

          source.append("    if (i < nnz) \n");
          source.append("    { \n");
          source.append("      unsigned int row_index_inc = 0; \n");
          source.append("      while (i >= row_index_lookahead[row_index_inc + 1]) \n");
          source.append("        ++row_index_inc; \n");
          source.append("      row_index = row_at_window_start + row_index_inc; \n");
          source.append("      row_index_buffer[get_local_id(0)] = row_index; \n");
          source.append("    } \n");
          source.append("    else \n");
          source.append("    { \n");
          source.append("      row_index = size+1; \n");
          source.append("      row_index_buffer[get_local_id(0)] = size - 1; \n");
          source.append("    } \n");

          source.append("    barrier(CLK_LOCAL_MEM_FENCE); \n");

          source.append("    row_at_window_start = row_index_buffer[0]; \n");
          source.append("    row_at_window_end   = row_index_buffer[get_local_size(0) - 1]; \n");

              //forward elimination
          source.append("    for (unsigned int row = row_at_window_start; row <= row_at_window_end; ++row) \n");
          source.append("    { \n");
          source.append("      "); source.append(numeric_string); source.append(" result_entry = vector[row] / diagonal_entries[row]; \n");

          source.append("      if ( (row_index == row) && (col_index > row) ) \n");
          source.append("        vector[col_index] -= result_entry * matrix_entry; \n");

          source.append("      barrier(CLK_GLOBAL_MEM_FENCE); \n");
          source.append("    } \n");

          source.append("    row_at_window_start = row_at_window_end; \n");
          source.append("  } \n");

            // final step: Divide vector by diagonal entries:
          source.append("  for (unsigned int i = get_local_id(0); i < size; i += get_local_size(0)) \n");
          source.append("    vector[i] /= diagonal_entries[i]; \n");
          source.append("} \n");

        }

        template <typename StringType>
        void generate_compressed_matrix_trans_unit_lu_backward(StringType & source, std::string const & numeric_string)
        {

          // compute y in Ly = z for incomplete LU factorizations of a sparse matrix in compressed format
          source.append("__kernel void trans_unit_lu_backward( \n");
          source.append("          __global const unsigned int * row_indices, \n");
          source.append("          __global const unsigned int * column_indices, \n");
          source.append("          __global const "); source.append(numeric_string); source.append(" * elements, \n");
          source.append("          __global "); source.append(numeric_string); source.append(" * vector, \n");
          source.append("          unsigned int size) \n");
          source.append("{ \n");
          source.append("  __local unsigned int row_index_lookahead[256]; \n");
          source.append("  __local unsigned int row_index_buffer[256]; \n");

          source.append("  unsigned int row_index; \n");
          source.append("  unsigned int col_index; \n");
          source.append("  "); source.append(numeric_string); source.append(" matrix_entry; \n");
          source.append("  unsigned int nnz = row_indices[size]; \n");
          source.append("  unsigned int row_at_window_start = size; \n");
          source.append("  unsigned int row_at_window_end; \n");
          source.append("  unsigned int loop_end = ( (nnz - 1) / get_local_size(0) + 1) * get_local_size(0); \n");

          source.append("  for (unsigned int i2 = get_local_id(0); i2 < loop_end; i2 += get_local_size(0)) \n");
          source.append("  { \n");
          source.append("    unsigned int i = (nnz - i2) - 1; \n");
          source.append("    col_index    = (i2 < nnz) ? column_indices[i] : 0; \n");
          source.append("    matrix_entry = (i2 < nnz) ? elements[i]       : 0; \n");
          source.append("    row_index_lookahead[get_local_id(0)] = (row_at_window_start >= get_local_id(0)) ? row_indices[row_at_window_start - get_local_id(0)] : 0; \n");

          source.append("    barrier(CLK_LOCAL_MEM_FENCE); \n");

          source.append("    if (i2 < nnz) \n");
          source.append("    { \n");
          source.append("      unsigned int row_index_dec = 0; \n");
          source.append("      while (row_index_lookahead[row_index_dec] > i) \n");
          source.append("        ++row_index_dec; \n");
          source.append("      row_index = row_at_window_start - row_index_dec; \n");
          source.append("      row_index_buffer[get_local_id(0)] = row_index; \n");
          source.append("    } \n");
          source.append("    else \n");
          source.append("    { \n");
          source.append("      row_index = size+1; \n");
          source.append("      row_index_buffer[get_local_id(0)] = 0; \n");
          source.append("    } \n");

          source.append("    barrier(CLK_LOCAL_MEM_FENCE); \n");

          source.append("    row_at_window_start = row_index_buffer[0]; \n");
          source.append("    row_at_window_end   = row_index_buffer[get_local_size(0) - 1]; \n");

              //backward elimination
          source.append("    for (unsigned int row2 = 0; row2 <= (row_at_window_start - row_at_window_end); ++row2) \n");
          source.append("    { \n");
          source.append("      unsigned int row = row_at_window_start - row2; \n");
          source.append("      "); source.append(numeric_string); source.append(" result_entry = vector[row]; \n");

          source.append("      if ( (row_index == row) && (col_index < row) ) \n");
          source.append("        vector[col_index] -= result_entry * matrix_entry; \n");

          source.append("      barrier(CLK_GLOBAL_MEM_FENCE); \n");
          source.append("    } \n");

          source.append("    row_at_window_start = row_at_window_end; \n");
          source.append("  } \n");
          source.append("} \n");

        }


        template <typename StringType>
        void generate_compressed_matrix_trans_unit_lu_forward(StringType & source, std::string const & numeric_string)
        {

          // compute y in Ly = z for incomplete LU factorizations of a sparse matrix in compressed format
          source.append("__kernel void trans_unit_lu_forward( \n");
          source.append("          __global const unsigned int * row_indices, \n");
          source.append("          __global const unsigned int * column_indices, \n");
          source.append("          __global const "); source.append(numeric_string); source.append(" * elements, \n");
          source.append("          __global "); source.append(numeric_string); source.append(" * vector, \n");
          source.append("          unsigned int size) \n");
          source.append("{ \n");
          source.append("  __local unsigned int row_index_lookahead[256]; \n");
          source.append("  __local unsigned int row_index_buffer[256]; \n");

          source.append("  unsigned int row_index; \n");
          source.append("  unsigned int col_index; \n");
          source.append("  "); source.append(numeric_string); source.append(" matrix_entry; \n");
          source.append("  unsigned int nnz = row_indices[size]; \n");
          source.append("  unsigned int row_at_window_start = 0; \n");
          source.append("  unsigned int row_at_window_end = 0; \n");
          source.append("  unsigned int loop_end = ( (nnz - 1) / get_local_size(0) + 1) * get_local_size(0); \n");

          source.append("  for (unsigned int i = get_local_id(0); i < loop_end; i += get_local_size(0)) \n");
          source.append("  { \n");
          source.append("    col_index    = (i < nnz) ? column_indices[i] : 0; \n");
          source.append("    matrix_entry = (i < nnz) ? elements[i]       : 0; \n");
          source.append("    row_index_lookahead[get_local_id(0)] = (row_at_window_start + get_local_id(0) < size) ? row_indices[row_at_window_start + get_local_id(0)] : size - 1; \n");

          source.append("    barrier(CLK_LOCAL_MEM_FENCE); \n");

          source.append("    if (i < nnz) \n");
          source.append("    { \n");
          source.append("      unsigned int row_index_inc = 0; \n");
          source.append("      while (i >= row_index_lookahead[row_index_inc + 1]) \n");
          source.append("        ++row_index_inc; \n");
          source.append("      row_index = row_at_window_start + row_index_inc; \n");
          source.append("      row_index_buffer[get_local_id(0)] = row_index; \n");
          source.append("    } \n");
          source.append("    else \n");
          source.append("    { \n");
          source.append("      row_index = size+1; \n");
          source.append("      row_index_buffer[get_local_id(0)] = size - 1; \n");
          source.append("    } \n");

          source.append("    barrier(CLK_LOCAL_MEM_FENCE); \n");

          source.append("    row_at_window_start = row_index_buffer[0]; \n");
          source.append("    row_at_window_end   = row_index_buffer[get_local_size(0) - 1]; \n");

              //forward elimination
          source.append("    for (unsigned int row = row_at_window_start; row <= row_at_window_end; ++row) \n");
          source.append("    { \n");
          source.append("      "); source.append(numeric_string); source.append(" result_entry = vector[row]; \n");

          source.append("      if ( (row_index == row) && (col_index > row) ) \n");
          source.append("        vector[col_index] -= result_entry * matrix_entry; \n");

          source.append("      barrier(CLK_GLOBAL_MEM_FENCE); \n");
          source.append("    } \n");

          source.append("    row_at_window_start = row_at_window_end; \n");
          source.append("  } \n");
          source.append("} \n");

        }

        template <typename StringType>
        void generate_compressed_matrix_trans_unit_lu_forward_slow(StringType & source, std::string const & numeric_string)
        {

          // compute y in Ly = z for incomplete LU factorizations of a sparse matrix in compressed format
          source.append("__kernel void trans_unit_lu_forward_slow( \n");
          source.append("          __global const unsigned int * row_indices, \n");
          source.append("          __global const unsigned int * column_indices, \n");
          source.append("          __global const "); source.append(numeric_string); source.append(" * elements, \n");
          source.append("          __global "); source.append(numeric_string); source.append(" * vector, \n");
          source.append("          unsigned int size) \n");
          source.append("{ \n");
          source.append("  for (unsigned int row = 0; row < size; ++row) \n");
          source.append("  { \n");
          source.append("    "); source.append(numeric_string); source.append(" result_entry = vector[row]; \n");

          source.append("    unsigned int row_start = row_indices[row]; \n");
          source.append("    unsigned int row_stop  = row_indices[row + 1]; \n");
          source.append("    for (unsigned int entry_index = row_start + get_local_id(0); entry_index < row_stop; entry_index += get_local_size(0)) \n");
          source.append("    { \n");
          source.append("      unsigned int col_index = column_indices[entry_index]; \n");
          source.append("      if (col_index > row) \n");
          source.append("        vector[col_index] -= result_entry * elements[entry_index]; \n");
          source.append("    } \n");

          source.append("    barrier(CLK_GLOBAL_MEM_FENCE); \n");
          source.append("  } \n");
          source.append("} \n");

        }

        template <typename StringType>
        void generate_compressed_matrix_unit_lu_backward(StringType & source, std::string const & numeric_string)
        {

          // compute x in Ux = y for incomplete LU factorizations of a sparse matrix in compressed format
          source.append("__kernel void unit_lu_backward( \n");
          source.append("          __global const unsigned int * row_indices, \n");
          source.append("          __global const unsigned int * column_indices, \n");
          source.append("          __global const "); source.append(numeric_string); source.append(" * elements, \n");
          source.append("          __global "); source.append(numeric_string); source.append(" * vector, \n");
          source.append("          unsigned int size) \n");
          source.append("{ \n");
          source.append("  __local  unsigned int col_index_buffer[128]; \n");
          source.append("  __local  "); source.append(numeric_string); source.append(" element_buffer[128]; \n");
          source.append("  __local  "); source.append(numeric_string); source.append(" vector_buffer[128]; \n");

          source.append("  unsigned int nnz = row_indices[size]; \n");
          source.append("  unsigned int current_row = size-1; \n");
          source.append("  unsigned int row_at_window_start = size-1; \n");
          source.append("  "); source.append(numeric_string); source.append(" current_vector_entry = vector[size-1]; \n");
          source.append("  unsigned int loop_end = ( (nnz - 1) / get_local_size(0)) * get_local_size(0); \n");
          source.append("  unsigned int next_row = row_indices[size-1]; \n");

          source.append("  unsigned int i = loop_end + get_local_id(0); \n");
          source.append("  while (1) \n");
          source.append("  { \n");
              //load into shared memory (coalesced access):
          source.append("    if (i < nnz) \n");
          source.append("    { \n");
          source.append("      element_buffer[get_local_id(0)] = elements[i]; \n");
          source.append("      unsigned int tmp = column_indices[i]; \n");
          source.append("      col_index_buffer[get_local_id(0)] = tmp; \n");
          source.append("      vector_buffer[get_local_id(0)] = vector[tmp]; \n");
          source.append("    } \n");

          source.append("    barrier(CLK_LOCAL_MEM_FENCE); \n");

              //now a single thread does the remaining work in shared memory:
          source.append("    if (get_local_id(0) == 0) \n");
          source.append("    { \n");
              // traverse through all the loaded data from back to front:
          source.append("      for (unsigned int k2=0; k2<get_local_size(0); ++k2) \n");
          source.append("      { \n");
          source.append("        unsigned int k = (get_local_size(0) - k2) - 1; \n");

          source.append("        if (i+k >= nnz) \n");
          source.append("          continue; \n");

          source.append("        if (col_index_buffer[k] > row_at_window_start) \n"); //use recently computed results
          source.append("          current_vector_entry -= element_buffer[k] * vector_buffer[k]; \n");
          source.append("        else if (col_index_buffer[k] > current_row) \n"); //use buffered data
          source.append("          current_vector_entry -= element_buffer[k] * vector[col_index_buffer[k]]; \n");

          source.append("        if (i+k == next_row) \n"); //current row is finished. Write back result
          source.append("        { \n");
          source.append("          vector[current_row] = current_vector_entry; \n");
          source.append("          if (current_row > 0) \n"); //load next row's data
          source.append("          { \n");
          source.append("            --current_row; \n");
          source.append("            next_row = row_indices[current_row]; \n");
          source.append("            current_vector_entry = vector[current_row]; \n");
          source.append("          } \n");
          source.append("        } \n");


          source.append("      } \n"); // for k

          source.append("      row_at_window_start = current_row; \n");
          source.append("    } \n"); // if (get_local_id(0) == 0)

          source.append("    barrier(CLK_GLOBAL_MEM_FENCE); \n");

          source.append("    if (i < get_local_size(0)) \n");
          source.append("      break; \n");

          source.append("    i -= get_local_size(0); \n");
          source.append("  } \n"); //for i
          source.append("} \n");

        }

        template <typename StringType>
        void generate_compressed_matrix_unit_lu_forward(StringType & source, std::string const & numeric_string)
        {

          // compute y in Ly = z for incomplete LU factorizations of a sparse matrix in compressed format
          source.append("__kernel void unit_lu_forward( \n");
          source.append("          __global const unsigned int * row_indices, \n");
          source.append("          __global const unsigned int * column_indices, \n");
          source.append("          __global const "); source.append(numeric_string); source.append(" * elements, \n");
          source.append("          __global "); source.append(numeric_string); source.append(" * vector, \n");
          source.append("          unsigned int size) \n");
          source.append("{ \n");
          source.append("  __local  unsigned int col_index_buffer[128]; \n");
          source.append("  __local  "); source.append(numeric_string); source.append(" element_buffer[128]; \n");
          source.append("  __local  "); source.append(numeric_string); source.append(" vector_buffer[128]; \n");

          source.append("  unsigned int nnz = row_indices[size]; \n");
          source.append("  unsigned int current_row = 0; \n");
          source.append("  unsigned int row_at_window_start = 0; \n");
          source.append("  "); source.append(numeric_string); source.append(" current_vector_entry = vector[0]; \n");
          source.append("  unsigned int loop_end = (nnz / get_local_size(0) + 1) * get_local_size(0); \n");
          source.append("  unsigned int next_row = row_indices[1]; \n");

          source.append("  for (unsigned int i = get_local_id(0); i < loop_end; i += get_local_size(0)) \n");
          source.append("  { \n");
              //load into shared memory (coalesced access):
          source.append("    if (i < nnz) \n");
          source.append("    { \n");
          source.append("      element_buffer[get_local_id(0)] = elements[i]; \n");
          source.append("      unsigned int tmp = column_indices[i]; \n");
          source.append("      col_index_buffer[get_local_id(0)] = tmp; \n");
          source.append("      vector_buffer[get_local_id(0)] = vector[tmp]; \n");
          source.append("    } \n");

          source.append("    barrier(CLK_LOCAL_MEM_FENCE); \n");

              //now a single thread does the remaining work in shared memory:
          source.append("    if (get_local_id(0) == 0) \n");
          source.append("    { \n");
                // traverse through all the loaded data:
          source.append("      for (unsigned int k=0; k<get_local_size(0); ++k) \n");
          source.append("      { \n");
          source.append("        if (i+k == next_row) \n"); //current row is finished. Write back result
          source.append("        { \n");
          source.append("          vector[current_row] = current_vector_entry; \n");
          source.append("          ++current_row; \n");
          source.append("          if (current_row < size) //load next row's data \n");
          source.append("          { \n");
          source.append("            next_row = row_indices[current_row+1]; \n");
          source.append("            current_vector_entry = vector[current_row]; \n");
          source.append("          } \n");
          source.append("        } \n");

          source.append("        if (current_row < size && col_index_buffer[k] < current_row) \n"); //substitute
          source.append("        { \n");
          source.append("          if (col_index_buffer[k] < row_at_window_start) \n"); //use recently computed results
          source.append("            current_vector_entry -= element_buffer[k] * vector_buffer[k]; \n");
          source.append("          else if (col_index_buffer[k] < current_row) \n"); //use buffered data
          source.append("            current_vector_entry -= element_buffer[k] * vector[col_index_buffer[k]]; \n");
          source.append("        } \n");

          source.append("      } \n"); // for k

          source.append("      row_at_window_start = current_row; \n");
          source.append("    } \n"); // if (get_local_id(0) == 0)

          source.append("    barrier(CLK_GLOBAL_MEM_FENCE); \n");
          source.append("  } //for i \n");
          source.append("} \n");

        }

        template <typename StringType>
        void generate_compressed_matrix_vec_mul(StringType & source, std::string const & numeric_string)
        {

          source.append("__kernel void vec_mul( \n");
          source.append("          __global const unsigned int * row_indices, \n");
          source.append("          __global const unsigned int * column_indices, \n");
          source.append("          __global const "); source.append(numeric_string); source.append(" * elements, \n");
          source.append("          __global const "); source.append(numeric_string); source.append(" * x, \n");
          source.append("          uint4 layout_x, \n");
          source.append("          __global "); source.append(numeric_string); source.append(" * result, \n");
          source.append("          uint4 layout_result) \n");
          source.append("{ \n");
          source.append("  for (unsigned int row = get_global_id(0); row < layout_result.z; row += get_global_size(0)) \n");
          source.append("  { \n");
          source.append("    "); source.append(numeric_string); source.append(" dot_prod = 0; \n");
          source.append("    unsigned int row_end = row_indices[row+1]; \n");
          source.append("    for (unsigned int i = row_indices[row]; i < row_end; ++i) \n");
          source.append("      dot_prod += elements[i] * x[column_indices[i] * layout_x.y + layout_x.x]; \n");
          source.append("    result[row * layout_result.y + layout_result.x] = dot_prod; \n");
          source.append("  } \n");
          source.append("} \n");

        }

        template <typename StringType>
        void generate_compressed_matrix_vec_mul4(StringType & source, std::string const & numeric_string)
        {
          source.append("__kernel void vec_mul4( \n");
          source.append("          __global const unsigned int * row_indices, \n");
          source.append("          __global const uint4 * column_indices, \n");
          source.append("          __global const "); source.append(numeric_string); source.append("4 * elements, \n");
          source.append("          __global const "); source.append(numeric_string); source.append(" * x, \n");
          source.append("          uint4 layout_x, \n");
          source.append("          __global "); source.append(numeric_string); source.append(" * result, \n");
          source.append("          uint4 layout_result) \n");
          source.append("{ \n");
          source.append("  "); source.append(numeric_string); source.append(" dot_prod; \n");
          source.append("  unsigned int start, next_stop; \n");
          source.append("  uint4 col_idx; \n");
          source.append("  "); source.append(numeric_string); source.append("4 tmp_vec; \n");
          source.append("  "); source.append(numeric_string); source.append("4 tmp_entries; \n");

          source.append("  for (unsigned int row = get_global_id(0); row < layout_result.z; row += get_global_size(0)) \n");
          source.append("  { \n");
          source.append("    dot_prod = 0; \n");
          source.append("    start = row_indices[row] / 4; \n");
          source.append("    next_stop = row_indices[row+1] / 4; \n");

          source.append("    for (unsigned int i = start; i < next_stop; ++i) \n");
          source.append("    { \n");
          source.append("      col_idx = column_indices[i]; \n");

          source.append("      tmp_entries = elements[i]; \n");
          source.append("      tmp_vec.x = x[col_idx.x * layout_x.y + layout_x.x]; \n");
          source.append("      tmp_vec.y = x[col_idx.y * layout_x.y + layout_x.x]; \n");
          source.append("      tmp_vec.z = x[col_idx.z * layout_x.y + layout_x.x]; \n");
          source.append("      tmp_vec.w = x[col_idx.w * layout_x.y + layout_x.x]; \n");

          source.append("      dot_prod += dot(tmp_entries, tmp_vec); \n");
          source.append("    } \n");
          source.append("    result[row * layout_result.y + layout_result.x] = dot_prod; \n");
          source.append("  } \n");
          source.append("} \n");
        }

        template <typename StringType>
        void generate_compressed_matrix_vec_mul8(StringType & source, std::string const & numeric_string)
        {
          source.append("__kernel void vec_mul8( \n");
          source.append("          __global const unsigned int * row_indices, \n");
          source.append("          __global const uint8 * column_indices, \n");
          source.append("          __global const "); source.append(numeric_string); source.append("8 * elements, \n");
          source.append("          __global const "); source.append(numeric_string); source.append(" * x, \n");
          source.append("          uint4 layout_x, \n");
          source.append("          __global "); source.append(numeric_string); source.append(" * result, \n");
          source.append("          uint4 layout_result) \n");
          source.append("{ \n");
          source.append("  "); source.append(numeric_string); source.append(" dot_prod; \n");
          source.append("  unsigned int start, next_stop; \n");
          source.append("  uint8 col_idx; \n");
          source.append("  "); source.append(numeric_string); source.append("8 tmp_vec; \n");
          source.append("  "); source.append(numeric_string); source.append("8 tmp_entries; \n");

          source.append("  for (unsigned int row = get_global_id(0); row < layout_result.z; row += get_global_size(0)) \n");
          source.append("  { \n");
          source.append("    dot_prod = 0; \n");
          source.append("    start = row_indices[row] / 8; \n");
          source.append("    next_stop = row_indices[row+1] / 8; \n");

          source.append("    for (unsigned int i = start; i < next_stop; ++i) \n");
          source.append("    { \n");
          source.append("      col_idx = column_indices[i]; \n");

          source.append("      tmp_entries = elements[i]; \n");
          source.append("      tmp_vec.s0 = x[col_idx.s0 * layout_x.y + layout_x.x]; \n");
          source.append("      tmp_vec.s1 = x[col_idx.s1 * layout_x.y + layout_x.x]; \n");
          source.append("      tmp_vec.s2 = x[col_idx.s2 * layout_x.y + layout_x.x]; \n");
          source.append("      tmp_vec.s3 = x[col_idx.s3 * layout_x.y + layout_x.x]; \n");
          source.append("      tmp_vec.s4 = x[col_idx.s4 * layout_x.y + layout_x.x]; \n");
          source.append("      tmp_vec.s5 = x[col_idx.s5 * layout_x.y + layout_x.x]; \n");
          source.append("      tmp_vec.s6 = x[col_idx.s6 * layout_x.y + layout_x.x]; \n");
          source.append("      tmp_vec.s7 = x[col_idx.s7 * layout_x.y + layout_x.x]; \n");

          source.append("      dot_prod += dot(tmp_entries.lo, tmp_vec.lo); \n");
          source.append("      dot_prod += dot(tmp_entries.hi, tmp_vec.hi); \n");
          source.append("    } \n");
          source.append("    result[row * layout_result.y + layout_result.x] = dot_prod; \n");
          source.append("  } \n");
          source.append("} \n");
        }

        template <typename StringType>
        void generate_compressed_matrix_vec_mul_cpu(StringType & source, std::string const & numeric_string)
        {
          source.append("__kernel void vec_mul_cpu( \n");
          source.append("          __global const unsigned int * row_indices, \n");
          source.append("          __global const unsigned int * column_indices, \n");
          source.append("          __global const "); source.append(numeric_string); source.append(" * elements, \n");
          source.append("          __global const "); source.append(numeric_string); source.append(" * vector, \n");
          source.append("          __global "); source.append(numeric_string); source.append(" * result, \n");
          source.append("          unsigned int size) \n");
          source.append("{ \n");
          source.append("  unsigned int work_per_item = max((uint) (size / get_global_size(0)), (uint) 1); \n");
          source.append("  unsigned int row_start = get_global_id(0) * work_per_item; \n");
          source.append("  unsigned int row_stop  = min( (uint) ((get_global_id(0) + 1) * work_per_item), (uint) size); \n");
          source.append("  for (unsigned int row = row_start; row < row_stop; ++row) \n");
          source.append("  { \n");
          source.append("    "); source.append(numeric_string); source.append(" dot_prod = ("); source.append(numeric_string); source.append(")0; \n");
          source.append("    unsigned int row_end = row_indices[row+1]; \n");
          source.append("    for (unsigned int i = row_indices[row]; i < row_end; ++i) \n");
          source.append("      dot_prod += elements[i] * vector[column_indices[i]]; \n");
          source.append("    result[row] = dot_prod; \n");
          source.append("  } \n");
          source.append("} \n");

        }


        //////////////////////////// Part 2: Main kernel class ////////////////////////////////////

        // main kernel class
        /** @brief Main kernel class for generating OpenCL kernels for compressed_matrix. */
        template <typename NumericT>
        struct compressed_matrix
        {
          static std::string program_name()
          {
            return viennacl::ocl::type_to_string<NumericT>::apply() + "_compressed_matrix";
          }

          static void init(viennacl::ocl::context & ctx)
          {
            viennacl::ocl::DOUBLE_PRECISION_CHECKER<NumericT>::apply(ctx);
            std::string numeric_string = viennacl::ocl::type_to_string<NumericT>::apply();

            static std::map<cl_context, bool> init_done;
            if (!init_done[ctx.handle().get()])
            {
              std::string source;
              source.reserve(1024);

              viennacl::ocl::append_double_precision_pragma<NumericT>(ctx, source);

              if (numeric_string == "float" || numeric_string == "double")
              {
                generate_compressed_matrix_block_trans_lu_backward(source, numeric_string);
                generate_compressed_matrix_block_trans_unit_lu_forward(source, numeric_string);
                generate_compressed_matrix_jacobi(source, numeric_string);
                generate_compressed_matrix_lu_backward(source, numeric_string);
                generate_compressed_matrix_lu_forward(source, numeric_string);
                generate_compressed_matrix_trans_lu_backward(source, numeric_string);
                generate_compressed_matrix_trans_lu_forward(source, numeric_string);
                generate_compressed_matrix_trans_unit_lu_backward(source, numeric_string);
                generate_compressed_matrix_trans_unit_lu_forward(source, numeric_string);
                generate_compressed_matrix_trans_unit_lu_forward_slow(source, numeric_string);
                generate_compressed_matrix_unit_lu_backward(source, numeric_string);
                generate_compressed_matrix_unit_lu_forward(source, numeric_string);
              }
              generate_compressed_matrix_dense_matrix_multiplication(source, numeric_string);
              generate_compressed_matrix_row_info_extractor(source, numeric_string);
              generate_compressed_matrix_vec_mul(source, numeric_string);
              generate_compressed_matrix_vec_mul4(source, numeric_string);
              generate_compressed_matrix_vec_mul8(source, numeric_string);
              generate_compressed_matrix_vec_mul_cpu(source, numeric_string);

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

