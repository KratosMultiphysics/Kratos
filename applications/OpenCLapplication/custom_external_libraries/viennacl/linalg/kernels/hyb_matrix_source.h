#ifndef VIENNACL_LINALG_KERNELS_HYB_MATRIX_SOURCE_HPP_
#define VIENNACL_LINALG_KERNELS_HYB_MATRIX_SOURCE_HPP_
//Automatically generated file from auxiliary-directory, do not edit manually!
namespace viennacl
{
 namespace linalg
 {
  namespace kernels
  {
const char * const hyb_matrix_align1_vec_mul = 
"__kernel void vec_mul(\n"
"    const __global int* ell_coords,\n"
"    const __global float* ell_elements,\n"
"    const __global uint* csr_rows,\n"
"    const __global uint* csr_cols,\n"
"    const __global float* csr_elements,\n"
"    const __global float * vector,\n"
"    __global float * result,\n"
"    unsigned int row_num,\n"
"    unsigned int internal_row_num,\n"
"    unsigned int items_per_row,\n"
"    unsigned int aligned_items_per_row\n"
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
"            float val = ell_elements[offset];\n"
"            if(val != 0.0f)\n"
"            {\n"
"                int col = ell_coords[offset];    \n"
"                sum += (vector[col] * val);\n"
"            }\n"
"            \n"
"        }\n"
"        uint col_begin = csr_rows[row_id];\n"
"        uint col_end   = csr_rows[row_id + 1];\n"
"        for(uint item_id = col_begin; item_id < col_end; item_id++)\n"
"        {\n"
"            sum += (vector[csr_cols[item_id]] * csr_elements[item_id]);\n"
"        }\n"
"        result[row_id] = sum;\n"
"    }\n"
"}\n"
; //hyb_matrix_align1_vec_mul

  }  //namespace kernels
 }  //namespace linalg
}  //namespace viennacl
#endif
