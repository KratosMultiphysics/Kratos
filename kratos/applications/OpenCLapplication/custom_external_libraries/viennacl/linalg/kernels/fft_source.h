#ifndef VIENNACL_LINALG_KERNELS_FFT_SOURCE_HPP_
#define VIENNACL_LINALG_KERNELS_FFT_SOURCE_HPP_
//Automatically generated file from auxiliary-directory, do not edit manually!
namespace viennacl
{
 namespace linalg
 {
  namespace kernels
  {
const char * const fft_align1_fft_mult_vec = 
"// elementwise product of two complex vectors\n"
"__kernel void fft_mult_vec(__global const float2* input1,\n"
"                          __global const float2* input2,\n"
"                          __global float2* output,\n"
"                          unsigned int size) {\n"
"    for (unsigned int i = get_global_id(0); i < size; i += get_global_size(0)) {\n"
"        float2 in1 = input1[i];\n"
"        float2 in2 = input2[i];\n"
"        output[i] = (float2)(in1.x * in2.x - in1.y * in2.y, in1.x * in2.y + in1.y * in2.x);\n"
"    }\n"
"}\n"
; //fft_align1_fft_mult_vec

const char * const fft_align1_bluestein_pre = 
"// Preprocessing phase of Bluestein algorithm\n"
"__kernel void bluestein_pre(__global float2* input,\n"
"                            __global float2* A,\n"
"                            __global float2* B,\n"
"                            unsigned int size,\n"
"                            unsigned int ext_size\n"
"                           ) {\n"
"    unsigned int glb_id = get_global_id(0);\n"
"    unsigned int glb_sz = get_global_size(0);\n"
"    unsigned int double_size = size << 1;\n"
"    float sn_a, cs_a;\n"
"    const float NUM_PI = 3.14159265358979323846;\n"
"    for(unsigned int i = glb_id; i < size; i += glb_sz) {\n"
"        unsigned int rm = i * i % (double_size);\n"
"        float angle = (float)rm / size * NUM_PI;\n"
"        sn_a = sincos(-angle, &cs_a);\n"
"        float2 a_i = (float2)(cs_a, sn_a);\n"
"        float2 b_i = (float2)(cs_a, -sn_a);\n"
"        A[i] = (float2)(input[i].x * a_i.x - input[i].y * a_i.y, input[i].x * a_i.y + input[i].y * a_i.x);\n"
"        B[i] = b_i;\n"
"        // very bad instruction, to be fixed\n"
"        if(i) \n"
"          B[ext_size - i] = b_i;\n"
"    }\n"
"}\n"
; //fft_align1_bluestein_pre

const char * const fft_align1_zero2 = 
"// Zero two complex vectors (to avoid kernel launch overhead)\n"
"__kernel void zero2(__global float2* input1,\n"
"                    __global float2* input2,\n"
"                    unsigned int size) {\n"
"    for (unsigned int i = get_global_id(0); i < size; i += get_global_size(0)) {\n"
"        input1[i] = 0;\n"
"        input2[i] = 0;\n"
"    }\n"
"}\n"
; //fft_align1_zero2

const char * const fft_align1_bluestein_post = 
"// Postprocessing phase of Bluestein algorithm\n"
"__kernel void bluestein_post(__global float2* Z,\n"
"                             __global float2* out,\n"
"                             unsigned int size) \n"
"{\n"
"    unsigned int glb_id = get_global_id(0);\n"
"    unsigned int glb_sz = get_global_size(0);\n"
"    unsigned int double_size = size << 1;\n"
"    float sn_a, cs_a;\n"
"    const float NUM_PI = 3.14159265358979323846;\n"
"    for(unsigned int i = glb_id; i < size; i += glb_sz) {\n"
"        unsigned int rm = i * i % (double_size);\n"
"        float angle = (float)rm / size * (-NUM_PI);\n"
"        sn_a = sincos(angle, &cs_a);\n"
"        float2 b_i = (float2)(cs_a, sn_a);\n"
"        out[i] = (float2)(Z[i].x * b_i.x - Z[i].y * b_i.y, Z[i].x * b_i.y + Z[i].y * b_i.x);\n"
"    }\n"
"}\n"
; //fft_align1_bluestein_post

const char * const fft_align1_complex_to_real = 
"__kernel void complex_to_real(__global float2* in,\n"
"                              __global float* out,\n"
"                              unsigned int size) {\n"
"    for (unsigned int i = get_global_id(0); i < size; i += get_global_size(0)) {\n"
"        out[i] = in[i].x;\n"
"    }\n"
"}\n"
; //fft_align1_complex_to_real

const char * const fft_align1_transpose_inplace = 
"// inplace-transpose of a matrix\n"
"__kernel void transpose_inplace(__global float2* input,\n"
"                        unsigned int row_num,\n"
"                        unsigned int col_num) {\n"
"    unsigned int size = row_num * col_num;\n"
"    for(unsigned int i = get_global_id(0); i < size; i+= get_global_size(0)) {\n"
"        unsigned int row = i / col_num;\n"
"        unsigned int col = i - row*col_num;\n"
"        unsigned int new_pos = col * row_num + row;\n"
"        //new_pos = col < row?0:1;\n"
"        //input[i] = new_pos;\n"
"        if(i < new_pos) {\n"
"            float2 val = input[i];\n"
"            input[i] = input[new_pos];\n"
"            input[new_pos] = val;\n"
"        }\n"
"    }\n"
"}\n"
; //fft_align1_transpose_inplace

const char * const fft_align1_transpose = 
"// simplistic matrix transpose function\n"
"__kernel void transpose(__global float2* input,\n"
"                        __global float2* output,\n"
"                        unsigned int row_num,\n"
"                        unsigned int col_num) {\n"
"    unsigned int size = row_num * col_num;\n"
"    for(unsigned int i = get_global_id(0); i < size; i+= get_global_size(0)) {\n"
"        unsigned int row = i / col_num;\n"
"        unsigned int col = i - row*col_num;\n"
"        unsigned int new_pos = col * row_num + row;\n"
"        output[new_pos] = input[i];\n"
"    }\n"
"}\n"
; //fft_align1_transpose

const char * const fft_align1_reverse_inplace = 
"// reverses the entries in a vector\n"
"__kernel void reverse_inplace(__global float* vec, uint size) {\n"
"    for(uint i = get_global_id(0); i < (size >> 1); i+=get_global_size(0)) {\n"
"        float val1 = vec[i];\n"
"        float val2 = vec[size - i - 1];\n"
"        vec[i] = val2;\n"
"        vec[size - i - 1] = val1;\n"
"    }\n"
"}\n"
; //fft_align1_reverse_inplace

const char * const fft_align1_real_to_complex = 
"// embedd a real-valued vector into a complex one\n"
"__kernel void real_to_complex(__global float* in,\n"
"                              __global float2* out,\n"
"                              unsigned int size) {\n"
"    for (unsigned int i = get_global_id(0); i < size; i += get_global_size(0)) {\n"
"        float2 val = 0;\n"
"        val.x = in[i];\n"
"        out[i] = val;\n"
"    }\n"
"}\n"
; //fft_align1_real_to_complex

const char * const fft_align1_fft_div_vec_scalar = 
"// divide a vector by a scalar (to be removed...)\n"
"__kernel void fft_div_vec_scalar(__global float2* input1, unsigned int size, float factor) {\n"
"    for (unsigned int i = get_global_id(0); i < size; i += get_global_size(0)) {\n"
"        input1[i] /= factor;\n"
"    }\n"
"}\n"
; //fft_align1_fft_div_vec_scalar

const char * const fft_align1_vandermonde_prod = 
"// computes the matrix vector product with a Vandermonde matrix\n"
"__kernel void vandermonde_prod(__global float* vander,\n"
"                                __global float* vector,\n"
"                                __global float* result,\n"
"                                uint size) {\n"
"    for(uint i = get_global_id(0); i < size; i+= get_global_size(0)) {\n"
"        float mul = vander[i];\n"
"        float pwr = 1;\n"
"        float val = 0;\n"
"        for(uint j = 0; j < size; j++) {\n"
"            val = val + pwr * vector[j];\n"
"            pwr *= mul;\n"
"        }\n"
"            \n"
"        result[i] = val;\n"
"    }\n"
"}\n"
; //fft_align1_vandermonde_prod

  }  //namespace kernels
 }  //namespace linalg
}  //namespace viennacl
#endif
