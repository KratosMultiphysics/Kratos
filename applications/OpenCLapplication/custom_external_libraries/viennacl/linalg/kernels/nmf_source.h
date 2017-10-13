#ifndef VIENNACL_LINALG_KERNELS_NMF_SOURCE_HPP_
#define VIENNACL_LINALG_KERNELS_NMF_SOURCE_HPP_
//Automatically generated file from auxiliary-directory, do not edit manually!
namespace viennacl
{
 namespace linalg
 {
  namespace kernels
  {
const char * const nmf_align1_sub_wise = 
"__kernel void sub_wise(\n"
"          __global const float * matrix1,\n"
"          __global const float * matrix2,\n"
"          __global float * result,\n"
"          unsigned int size)\n"
"{\n"
"  for (unsigned int i = get_global_id(0); i < size; i += get_global_size(0))\n"
"    result[i] = matrix1[i] - matrix2[i];\n"
"}\n"
; //nmf_align1_sub_wise

const char * const nmf_align1_el_wise_mul_div = 
"__kernel void el_wise_mul_div(\n"
"          __global float * matrix1,\n"
"          __global const float * matrix2,\n"
"          __global const float * matrix3,\n"
"          unsigned int size)\n"
"{\n"
"  for (unsigned int i = get_global_id(0); i < size; i += get_global_size(0)) \n"
"  {\n"
"    float val = matrix1[i] * matrix2[i];\n"
"    float divisor = matrix3[i];\n"
"    matrix1[i] = (divisor > 0.00001) ? (val / divisor) : 0;\n"
"  };\n"
"};\n"
; //nmf_align1_el_wise_mul_div

  }  //namespace kernels
 }  //namespace linalg
}  //namespace viennacl
#endif
