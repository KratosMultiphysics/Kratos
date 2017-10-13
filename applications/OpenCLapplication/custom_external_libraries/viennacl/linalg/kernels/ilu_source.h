#ifndef VIENNACL_LINALG_KERNELS_ILU_SOURCE_HPP_
#define VIENNACL_LINALG_KERNELS_ILU_SOURCE_HPP_
//Automatically generated file from auxiliary-directory, do not edit manually!
namespace viennacl
{
 namespace linalg
 {
  namespace kernels
  {
const char * const ilu_align1_block_ilu_substitute = 
" __kernel void block_ilu_substitute(\n"
"           __global const unsigned int * row_jumper_L,      //L part (note that L is transposed in memory)\n"
"           __global const unsigned int * column_indices_L, \n"
"           __global const float * elements_L,\n"
"           __global const unsigned int * row_jumper_U,      //U part (note that U is transposed in memory)\n"
"           __global const unsigned int * column_indices_U,\n"
"           __global const float * elements_U,\n"
"           __global const float * elements_D,              //diagonal\n"
"           __global const unsigned int * block_offsets,\n"
"           __global float * result,\n"
"           unsigned int size)\n"
" {\n"
"   unsigned int col_start = block_offsets[2*get_group_id(0)];\n"
"   unsigned int col_stop  = block_offsets[2*get_group_id(0)+1];\n"
"   unsigned int row_start = row_jumper_L[col_start];\n"
"   unsigned int row_stop;\n"
"   float result_entry = 0;\n"
"   if (col_start <= col_stop)\n"
"     return;\n"
"   //forward elimination, using L:\n"
"   for (unsigned int col = col_start; col < col_stop; ++col)\n"
"   {\n"
"     result_entry = result[col];\n"
"     row_stop = row_jumper_L[col + 1];\n"
"     for (unsigned int row_index = row_start + get_local_id(0); row_index < row_stop; ++row_index) \n"
"       result[column_indices_L[row_index]] -= result_entry * elements_L[row_index]; \n"
"     row_start = row_stop; //for next iteration (avoid unnecessary loads from GPU RAM)\n"
"     barrier(CLK_GLOBAL_MEM_FENCE);\n"
"   } \n"
"   //backward elimination, using U and D: \n"
"   for (unsigned int iter = 0; iter < col_stop - col_start; ++iter) \n"
"   { \n"
"     result_entry = result[col_stop - iter - 1] / elements_D[col_stop - iter - 1]; \n"
"     row_start = row_jumper_U[col_stop - iter - 1]; \n"
"     row_stop  = row_jumper_U[col_stop - iter]; \n"
"     for (unsigned int row_index = row_start + get_local_id(0); row_index < row_stop; ++row_index) \n"
"       result[column_indices_U[row_index]] -= result_entry * elements_U[row_index]; \n"
"     barrier(CLK_GLOBAL_MEM_FENCE); \n"
"     if (get_local_id(0) == 0) \n"
"       result[col_stop - iter - 1] = result_entry; \n"
"   } \n"
" };\n"
; //ilu_align1_block_ilu_substitute

  }  //namespace kernels
 }  //namespace linalg
}  //namespace viennacl
#endif
