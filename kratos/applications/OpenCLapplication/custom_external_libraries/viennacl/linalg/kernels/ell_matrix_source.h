#ifndef VIENNACL_LINALG_KERNELS_ELL_MATRIX_SOURCE_HPP_
#define VIENNACL_LINALG_KERNELS_ELL_MATRIX_SOURCE_HPP_
//Automatically generated file from auxiliary-directory, do not edit manually!
namespace viennacl
{
 namespace linalg
 {
  namespace kernels
  {
const char * const ell_matrix_align1_vec_mul = 
"__kernel void vec_mul(\n"
"    const __global int* coords,\n"
"    const __global float* elements,\n"
"    const __global const float * vector,\n"
"    __global float * result,\n"
"    const unsigned int row_num,\n"
"    const unsigned int col_num,\n"
"    const unsigned int internal_row_num,\n"
"    const unsigned int items_per_row,\n"
"    const unsigned int aligned_items_per_row\n"
"    )\n"
"{\n"
"    uint glb_id = get_global_id(0);\n"
"    uint glb_sz = get_global_size(0);\n"
"    for(uint row_id = glb_id; row_id < row_num; row_id += glb_sz)\n"
"    {\n"
"        float sum = 0;\n"
"        \n"
"        uint offset = row_id;\n"
"        for(uint item_id = 0; item_id < items_per_row; item_id++, offset += internal_row_num)\n"
"        {\n"
"            float val = elements[offset];\n"
"            if(val != 0.0f)\n"
"            {\n"
"                int col = coords[offset];    \n"
"                sum += (vector[col] * val);\n"
"            }\n"
"            \n"
"        }\n"
"        result[row_id] = sum;\n"
"    }\n"
"}\n"
; //ell_matrix_align1_vec_mul

  }  //namespace kernels
 }  //namespace linalg
}  //namespace viennacl
#endif
