#ifndef _VIENNACL_SCALAR_SOURCE_HPP_
#define _VIENNACL_SCALAR_SOURCE_HPP_
//Automatically generated file from aux-directory, do not edit manually!
namespace viennacl
{
 namespace linalg
 {
  namespace kernels
  {
const char * const scalar_align1_inplace_mul = 
" \n"
"__kernel void inplace_mul(\n"
"          __global float * val1,\n"
"          __global const float * val2) \n"
"{ \n"
"  if (get_global_id(0) == 0)\n"
"    *val1 *= *val2;\n"
"}\n"
" \n"
; //scalar_align1_inplace_mul

const char * const scalar_align1_cpu_inplace_mul = 
" \n"
"__kernel void cpu_inplace_mul(\n"
"          __global float * val1,\n"
"          float val2) \n"
"{ \n"
"  if (get_global_id(0) == 0)\n"
"    *val1 *= val2;\n"
"}\n"
" \n"
; //scalar_align1_cpu_inplace_mul

const char * const scalar_align1_sub = 
" \n"
"__kernel void sub(\n"
"          __global const float * val1,\n"
"          __global const float * val2, \n"
"          __global float * result) \n"
"{ \n"
"  if (get_global_id(0) == 0)\n"
"    *result = *val1 - *val2;\n"
"}\n"
" \n"
; //scalar_align1_sub

const char * const scalar_align1_cpu_inplace_sub = 
" \n"
"__kernel void cpu_inplace_sub(\n"
"          __global float * val1,\n"
"          float val2) \n"
"{ \n"
"  if (get_global_id(0) == 0)\n"
"    *val1 -= val2;\n"
"}\n"
" \n"
"\n"
; //scalar_align1_cpu_inplace_sub

const char * const scalar_align1_mul = 
" \n"
"__kernel void mul(\n"
"          __global const float * val1,\n"
"          __global const float * val2, \n"
"          __global float * result) \n"
"{ \n"
"  if (get_global_id(0) == 0)\n"
"    *result = *val1 * *val2;\n"
"}\n"
" \n"
; //scalar_align1_mul

const char * const scalar_align1_inplace_div = 
" \n"
"__kernel void inplace_div(\n"
"          __global float * val1,\n"
"          __global const float * val2) \n"
"{ \n"
"  if (get_global_id(0) == 0)\n"
"    *val1 /= *val2;\n"
"}\n"
" \n"
; //scalar_align1_inplace_div

const char * const scalar_align1_cpu_mul = 
" \n"
"__kernel void cpu_mul(\n"
"          __global const float * val1,\n"
"          float val2, \n"
"          __global float * result) \n"
"{ \n"
"  if (get_global_id(0) == 0)\n"
"    *result = *val1 * val2;\n"
"}\n"
" \n"
; //scalar_align1_cpu_mul

const char * const scalar_align1_divide = 
" \n"
"// note: 'div' seems to produce some name clashes with the OpenCL jit-compiler, thus using 'divide'\n"
"__kernel void divide(\n"
"          __global const float * val1,\n"
"          __global const float * val2, \n"
"          __global float * result) \n"
"{ \n"
"  if (get_global_id(0) == 0)\n"
"    *result = *val1 / *val2;\n"
"}\n"
"\n"
" \n"
; //scalar_align1_divide

const char * const scalar_align1_add = 
" \n"
"__kernel void add(\n"
"          __global const float * val1,\n"
"          __global const float * val2, \n"
"          __global float * result) \n"
"{ \n"
"  if (get_global_id(0) == 0)\n"
"    *result = *val1 + *val2;\n"
"}\n"
" \n"
; //scalar_align1_add

const char * const scalar_align1_cpu_inplace_add = 
" \n"
"__kernel void cpu_inplace_add(\n"
"          __global float * val1,\n"
"          float val2) \n"
"{ \n"
"  if (get_global_id(0) == 0)\n"
"    *val1 += val2;\n"
"}\n"
" \n"
; //scalar_align1_cpu_inplace_add

const char * const scalar_align1_cpu_div = 
" \n"
"__kernel void cpu_div(\n"
"          __global const float * val1,\n"
"          float val2, \n"
"          __global float * result) \n"
"{ \n"
"  if (get_global_id(0) == 0)\n"
"    *result = *val1 / val2;\n"
"}\n"
" \n"
; //scalar_align1_cpu_div

const char * const scalar_align1_cpu_inplace_div = 
" \n"
"__kernel void cpu_inplace_div(\n"
"          __global float * val1,\n"
"          float val2) \n"
"{ \n"
"  if (get_global_id(0) == 0)\n"
"    *val1 /= val2;\n"
"}\n"
"\n"
" \n"
; //scalar_align1_cpu_inplace_div

const char * const scalar_align1_cpu_add = 
" \n"
"__kernel void cpu_add(\n"
"          __global const float * val1,\n"
"          float val2, \n"
"          __global float * result) \n"
"{ \n"
"  if (get_global_id(0) == 0)\n"
"    *result = *val1 + val2;\n"
"}\n"
" \n"
; //scalar_align1_cpu_add

const char * const scalar_align1_inplace_sub = 
" \n"
"__kernel void inplace_sub(\n"
"          __global float * val1,\n"
"          __global const float * val2) \n"
"{ \n"
"  if (get_global_id(0) == 0)\n"
"    *val1 -= *val2;\n"
"}\n"
" \n"
; //scalar_align1_inplace_sub

const char * const scalar_align1_inplace_add = 
" \n"
"__kernel void inplace_add(\n"
"          __global float * val1,\n"
"          __global const float * val2) \n"
"{ \n"
"  if (get_global_id(0) == 0)\n"
"    *val1 += *val2;\n"
"}\n"
" \n"
; //scalar_align1_inplace_add

const char * const scalar_align1_cpu_sub = 
" \n"
"__kernel void cpu_sub(\n"
"          __global const float * val1,\n"
"          float val2, \n"
"          __global float * result) \n"
"{ \n"
"  if (get_global_id(0) == 0)\n"
"    *result = *val1 - val2;\n"
"}\n"
" \n"
; //scalar_align1_cpu_sub

  }  //namespace kernels
 }  //namespace linalg
}  //namespace viennacl
#endif
