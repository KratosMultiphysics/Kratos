#ifndef VIENNACL_LINALG_OPENCL_KERNELS_MATRIX_HPP
#define VIENNACL_LINALG_OPENCL_KERNELS_MATRIX_HPP

#include "viennacl/tools/tools.hpp"
#include "viennacl/ocl/kernel.hpp"
#include "viennacl/ocl/platform.hpp"
#include "viennacl/ocl/utils.hpp"

/** @file viennacl/linalg/opencl/kernels/matrix.hpp
 *  @brief Runtime generation of OpenCL kernels for matrix operations */
namespace viennacl
{
  namespace linalg
  {
    namespace opencl
    {
      namespace kernels
      {

        //////////////////////////// Part 1: Kernel generation routines ////////////////////////////////////

        /** @brief Enumeration for the scalar type in ambm-like operations */
        enum ambm_scalar_type
        {
          VIENNACL_AMBM_NONE = 0, // vector does not exist/contribute
          VIENNACL_AMBM_CPU,
          VIENNACL_AMBM_GPU
        };

        /** @brief Configuration struct for generating OpenCL kernels for linear combinations of matrices */
        struct ambm_config
        {
          ambm_config() : with_stride_and_range(true), is_row_major(true), a(VIENNACL_AMBM_CPU), b(VIENNACL_AMBM_NONE) {}

          bool with_stride_and_range;
          bool is_row_major;
          std::string      assign_op;
          ambm_scalar_type a;
          ambm_scalar_type b;
        };

        // just returns the for-loop
        template <typename StringType>
        void generate_ambm_impl2(StringType & source, ambm_config const & cfg, bool mult_alpha, bool mult_beta)
        {
          if (cfg.is_row_major)
          {
            source.append("  unsigned int row_gid = get_global_id(0) / get_local_size(0);\n");
            source.append("  unsigned int col_gid = get_global_id(0) % get_local_size(0);\n");
            source.append("  for (unsigned int row = row_gid; row < A_size1; row += get_num_groups(0))\n");
            source.append("    for (unsigned int col = col_gid; col < A_size2; col += get_local_size(0))\n");
          }
          else
          {
            source.append("  unsigned int col_gid = get_global_id(0) / get_local_size(0);\n");
            source.append("  unsigned int row_gid = get_global_id(0) % get_local_size(0);\n");
            source.append("  for (unsigned int col = col_gid; col < A_size2; col += get_num_groups(0))\n");
            source.append("    for (unsigned int row = row_gid; row < A_size1; row += get_local_size(0))\n");
          }

          if (cfg.with_stride_and_range)
          {
            if (cfg.is_row_major)
              source.append("      A[(row * A_inc1 + A_start1) * A_internal_size2 + (col * A_inc2 + A_start2)] ");
            else
              source.append("      A[(row * A_inc1 + A_start1) + (col * A_inc2 + A_start2) *  A_internal_size1] ");
            source.append(cfg.assign_op);
            if (cfg.is_row_major)
              source.append(" B[(row * B_inc1 + B_start1) * B_internal_size2 + (col * B_inc2 + B_start2)] ");
            else
              source.append(" B[(row * B_inc1 + B_start1) + (col * B_inc2 + B_start2) * B_internal_size1] ");

            if (mult_alpha)
              source.append("* alpha ");
            else
              source.append("/ alpha ");
            if (cfg.b != VIENNACL_AMBM_NONE)
            {
              if (cfg.is_row_major)
                source.append("+ C[(row * C_inc1 + C_start1) * C_internal_size2 + (col * C_inc2 + C_start2)] ");
              else
                source.append("+ C[(row * C_inc1 + C_start1) + (col * C_inc2 + C_start2) * C_internal_size1] ");
              if (mult_beta)
                source.append("* beta");
              else
                source.append("/ beta");
            }
          }
          else
          {
            if (cfg.is_row_major)
              source.append("    A[row * A_internal_size2 + col] ");
            else
              source.append("    A[row + col * A_internal_size1] ");
            source.append(cfg.assign_op);
            if (cfg.is_row_major)
              source.append(" B[row * B_internal_size2 + col] ");
            else
              source.append(" B[row + col * B_internal_size1] ");

            if (mult_alpha)
              source.append("* alpha ");
            else
              source.append("/ alpha ");
            if (cfg.b != VIENNACL_AMBM_NONE)
            {
              if (cfg.is_row_major)
                source.append("+ C[row * C_internal_size2 + col] ");
              else
                source.append("+ C[row + col * C_internal_size2] ");
              if (mult_beta)
                source.append("* beta");
              else
                source.append("/ beta");
            }
          }
          source.append("; \n");
        }

        template <typename StringType>
        void generate_ambm_impl(StringType & source, std::string const & numeric_string, ambm_config const & cfg)
        {
          source.append("__kernel void am");
          if (cfg.b != VIENNACL_AMBM_NONE)
            source.append("bm");
          if (cfg.assign_op != "=")
            source.append("_m");

          if (cfg.a == VIENNACL_AMBM_CPU)
            source.append("_cpu");
          else if (cfg.a == VIENNACL_AMBM_GPU)
            source.append("_gpu");

          if (cfg.b == VIENNACL_AMBM_CPU)
            source.append("_cpu");
          else if (cfg.b == VIENNACL_AMBM_GPU)
            source.append("_gpu");
          source.append("( \n");
          source.append("  __global "); source.append(numeric_string); source.append(" * A, \n");
          source.append("  unsigned int A_start1, unsigned int A_start2, \n");
          source.append("  unsigned int A_inc1,   unsigned int A_inc2, \n");
          source.append("  unsigned int A_size1,  unsigned int A_size2, \n");
          source.append("  unsigned int A_internal_size1,  unsigned int A_internal_size2, \n");
          if (cfg.a == VIENNACL_AMBM_CPU)
          {
            source.append("  "); source.append(numeric_string); source.append(" fac2, \n");
          }
          else if (cfg.a == VIENNACL_AMBM_GPU)
          {
            source.append("  __global "); source.append(numeric_string); source.append(" * fac2, \n");
          }
          source.append("  unsigned int options2, \n");  // 0: no action, 1: flip sign, 2: take inverse, 3: flip sign and take inverse
          source.append("  __global const "); source.append(numeric_string); source.append(" * B, \n");
          source.append("  unsigned int B_start1, unsigned int B_start2, \n");
          source.append("  unsigned int B_inc1,   unsigned int B_inc2, \n");
          source.append("  unsigned int B_internal_size1,  unsigned int B_internal_size2");

          if (cfg.b != VIENNACL_AMBM_NONE)
          {
            source.append(", \n\n");
            if (cfg.b == VIENNACL_AMBM_CPU)
            {
              source.append("  "); source.append(numeric_string); source.append(" fac3, \n");
            }
            else if (cfg.b == VIENNACL_AMBM_GPU)
            {
              source.append("  __global "); source.append(numeric_string); source.append(" * fac3, \n");
            }
            source.append("  unsigned int options3, \n");  // 0: no action, 1: flip sign, 2: take inverse, 3: flip sign and take inverse
            source.append("  __global const "); source.append(numeric_string); source.append(" * C, \n");
            source.append("  unsigned int C_start1, unsigned int C_start2, \n");
            source.append("  unsigned int C_inc1,   unsigned int C_inc2, \n");
            source.append("  unsigned int C_internal_size1,  unsigned int C_internal_size2 \n");
          }
          source.append(") { \n");

          if (cfg.a == VIENNACL_AMBM_CPU)
          {
            source.append("  "); source.append(numeric_string); source.append(" alpha = fac2; \n");
          }
          else if (cfg.a == VIENNACL_AMBM_GPU)
          {
            source.append("  "); source.append(numeric_string); source.append(" alpha = fac2[0]; \n");
          }
          source.append("  if (options2 & (1 << 0)) \n");
          source.append("    alpha = -alpha; \n");
          source.append(" \n");

          if (cfg.b == VIENNACL_AMBM_CPU)
          {
            source.append("  "); source.append(numeric_string); source.append(" beta = fac3; \n");
          }
          else if (cfg.b == VIENNACL_AMBM_GPU)
          {
            source.append("  "); source.append(numeric_string); source.append(" beta = fac3[0]; \n");
          }
          if (cfg.b != VIENNACL_AMBM_NONE)
          {
            source.append("  if (options3 & (1 << 0)) \n");
            source.append("    beta = -beta; \n");
            source.append(" \n");
          }
          source.append("  if (options2 & (1 << 1)) { \n");
          if (cfg.b != VIENNACL_AMBM_NONE)
          {
            source.append("    if (options3 & (1 << 1)) {\n");
            generate_ambm_impl2(source, cfg, false, false);
            source.append("    } else {\n");
            generate_ambm_impl2(source, cfg, false, true);
            source.append("    } \n");
          }
          else
            generate_ambm_impl2(source, cfg, false, true);
          source.append("  } else { \n");
          if (cfg.b != VIENNACL_AMBM_NONE)
          {
            source.append("    if (options3 & (1 << 1)) {\n");
            generate_ambm_impl2(source, cfg, true, false);
            source.append("    } else {\n");
            generate_ambm_impl2(source, cfg, true, true);
            source.append("    } \n");
          }
          else
            generate_ambm_impl2(source, cfg, true, true);
          source.append("  } \n");
          source.append("} \n");
        }

        template <typename StringType>
        void generate_ambm(StringType & source, std::string const & numeric_string, bool is_row_major)
        {
          ambm_config cfg;
          cfg.assign_op = "=";
          cfg.with_stride_and_range = true;
          cfg.is_row_major = is_row_major;

          // am
          cfg.b = VIENNACL_AMBM_NONE; cfg.a = VIENNACL_AMBM_CPU; generate_ambm_impl(source, numeric_string, cfg);
          cfg.b = VIENNACL_AMBM_NONE; cfg.a = VIENNACL_AMBM_GPU; generate_ambm_impl(source, numeric_string, cfg);

          // ambm
          cfg.a = VIENNACL_AMBM_CPU; cfg.b = VIENNACL_AMBM_CPU; generate_ambm_impl(source, numeric_string, cfg);
          cfg.a = VIENNACL_AMBM_CPU; cfg.b = VIENNACL_AMBM_GPU; generate_ambm_impl(source, numeric_string, cfg);
          cfg.a = VIENNACL_AMBM_GPU; cfg.b = VIENNACL_AMBM_CPU; generate_ambm_impl(source, numeric_string, cfg);
          cfg.a = VIENNACL_AMBM_GPU; cfg.b = VIENNACL_AMBM_GPU; generate_ambm_impl(source, numeric_string, cfg);

          // ambm_m
          cfg.assign_op = "+=";

          cfg.a = VIENNACL_AMBM_CPU; cfg.b = VIENNACL_AMBM_CPU; generate_ambm_impl(source, numeric_string, cfg);
          cfg.a = VIENNACL_AMBM_CPU; cfg.b = VIENNACL_AMBM_GPU; generate_ambm_impl(source, numeric_string, cfg);
          cfg.a = VIENNACL_AMBM_GPU; cfg.b = VIENNACL_AMBM_CPU; generate_ambm_impl(source, numeric_string, cfg);
          cfg.a = VIENNACL_AMBM_GPU; cfg.b = VIENNACL_AMBM_GPU; generate_ambm_impl(source, numeric_string, cfg);
        }

        template <typename StringType>
        void generate_assign_cpu(StringType & source, std::string const & numeric_string, bool is_row_major)
        {
          source.append("__kernel void assign_cpu( \n");
          source.append("  __global "); source.append(numeric_string); source.append(" * A, \n");
          source.append("  unsigned int A_start1, unsigned int A_start2, \n");
          source.append("  unsigned int A_inc1,   unsigned int A_inc2, \n");
          source.append("  unsigned int A_size1,  unsigned int A_size2, \n");
          source.append("  unsigned int A_internal_size1,  unsigned int A_internal_size2, \n");
          source.append("  "); source.append(numeric_string); source.append(" alpha) \n");
          source.append("{ \n");
          if (is_row_major)
          {
            source.append("  unsigned int row_gid = get_global_id(0) / get_local_size(0);\n");
            source.append("  unsigned int col_gid = get_global_id(0) % get_local_size(0);\n");
            source.append("  for (unsigned int row = row_gid; row < A_size1; row += get_num_groups(0))\n");
            source.append("    for (unsigned int col = col_gid; col < A_size2; col += get_local_size(0))\n");
            source.append("      A[(row * A_inc1 + A_start1) * A_internal_size2 + (col * A_inc2 + A_start2)] = alpha; \n");
          }
          else
          {
            source.append("  unsigned int row_gid = get_global_id(0) % get_local_size(0);\n");
            source.append("  unsigned int col_gid = get_global_id(0) / get_local_size(0);\n");
            source.append("  for (unsigned int col = col_gid; col < A_size2; col += get_num_groups(0))\n");
            source.append("    for (unsigned int row = row_gid; row < A_size1; row += get_local_size(0))\n");
            source.append("      A[(row * A_inc1 + A_start1) + (col * A_inc2 + A_start2) *  A_internal_size1] = alpha; \n");
          }
          source.append("} \n");
        }

        template <typename StringType>
        void generate_diagonal_assign_cpu(StringType & source, std::string const & numeric_string, bool is_row_major)
        {
          source.append("__kernel void diagonal_assign_cpu( \n");
          source.append("  __global "); source.append(numeric_string); source.append(" * A, \n");
          source.append("  unsigned int A_start1, unsigned int A_start2, \n");
          source.append("  unsigned int A_inc1,   unsigned int A_inc2, \n");
          source.append("  unsigned int A_size1,  unsigned int A_size2, \n");
          source.append("  unsigned int A_internal_size1,  unsigned int A_internal_size2, \n");
          source.append("  "); source.append(numeric_string); source.append(" alpha) \n");
          source.append("{ \n");
          source.append("  for (unsigned int idx = get_global_id(0); idx < min(A_size1, A_size2); idx += get_global_size(0))\n");
          if (is_row_major)
            source.append("    A[(idx * A_inc1 + A_start1) * A_internal_size2 + (idx * A_inc2 + A_start2)] = alpha; \n");
          else
            source.append("    A[(idx * A_inc1 + A_start1) + (idx * A_inc2 + A_start2) *  A_internal_size1] = alpha; \n");
          source.append("} \n");
        }

        template <typename StringType>
        void generate_element_op(StringType & source, std::string const & numeric_string, bool is_row_major)
        {
          source.append("__kernel void element_op( \n");
          source.append("  __global "); source.append(numeric_string); source.append(" * A, \n");
          source.append("  unsigned int A_start1, unsigned int A_start2, \n");
          source.append("  unsigned int A_inc1,   unsigned int A_inc2, \n");
          source.append("  unsigned int A_size1,  unsigned int A_size2, \n");
          source.append("  unsigned int A_internal_size1,  unsigned int A_internal_size2, \n");
          source.append("  __global "); source.append(numeric_string); source.append(" * B, \n");
          source.append("  unsigned int B_start1, unsigned int B_start2, \n");
          source.append("  unsigned int B_inc1,   unsigned int B_inc2, \n");
          source.append("  unsigned int B_internal_size1,  unsigned int B_internal_size2, \n");
          source.append("  __global "); source.append(numeric_string); source.append(" * C, \n");
          source.append("  unsigned int C_start1, unsigned int C_start2, \n");
          source.append("  unsigned int C_inc1,   unsigned int C_inc2, \n");
          source.append("  unsigned int C_internal_size1,  unsigned int C_internal_size2, \n");
          source.append("  unsigned int op_type) \n"); //0: product, 1: division, 2: pow
          source.append("{ \n");
          if (is_row_major)
          {
            source.append("  unsigned int row_gid = get_global_id(0) / get_local_size(0);\n");
            source.append("  unsigned int col_gid = get_global_id(0) % get_local_size(0);\n");
            source.append("  if (op_type == 2) {");
            if (numeric_string == "float" || numeric_string == "double")
            {
              source.append("    for (unsigned int row = row_gid; row < A_size1; row += get_num_groups(0))\n");
              source.append("      for (unsigned int col = col_gid; col < A_size2; col += get_local_size(0))\n");
              source.append("        A[(row * A_inc1 + A_start1) * A_internal_size2 + (col * A_inc2 + A_start2)] = \n");
              source.append("        pow(B[(row * B_inc1 + B_start1) * B_internal_size2 + (col * B_inc2 + B_start2)], \n");
              source.append("            C[(row * C_inc1 + C_start1) * C_internal_size2 + (col * C_inc2 + C_start2)]); \n");
            }
            source.append("  } else if (op_type == 1) {");
            source.append("    for (unsigned int row = row_gid; row < A_size1; row += get_num_groups(0))\n");
            source.append("      for (unsigned int col = col_gid; col < A_size2; col += get_local_size(0))\n");
            source.append("        A[(row * A_inc1 + A_start1) * A_internal_size2 + (col * A_inc2 + A_start2)] = \n");
            source.append("        B[(row * B_inc1 + B_start1) * B_internal_size2 + (col * B_inc2 + B_start2)] / \n");
            source.append("        C[(row * C_inc1 + C_start1) * C_internal_size2 + (col * C_inc2 + C_start2)]; \n");
            source.append("  } else if (op_type == 0) {");
            source.append("    for (unsigned int row = row_gid; row < A_size1; row += get_num_groups(0))\n");
            source.append("      for (unsigned int col = col_gid; col < A_size2; col += get_local_size(0))\n");
            source.append("        A[(row * A_inc1 + A_start1) * A_internal_size2 + (col * A_inc2 + A_start2)] = \n");
            source.append("        B[(row * B_inc1 + B_start1) * B_internal_size2 + (col * B_inc2 + B_start2)] * \n");
            source.append("        C[(row * C_inc1 + C_start1) * C_internal_size2 + (col * C_inc2 + C_start2)]; \n");
            source.append("  }");
          }
          else
          {
            source.append("  unsigned int row_gid = get_global_id(0) % get_local_size(0);\n");
            source.append("  unsigned int col_gid = get_global_id(0) / get_local_size(0);\n");
            source.append("  if (op_type == 2) {");
            if (numeric_string == "float" || numeric_string == "double")
            {
              source.append("    for (unsigned int col = col_gid; col < A_size2; col += get_num_groups(0))\n");
              source.append("      for (unsigned int row = row_gid; row < A_size1; row += get_local_size(0))\n");
              source.append("        A[(row * A_inc1 + A_start1) + (col * A_inc2 + A_start2) *  A_internal_size1] =  \n");
              source.append("          pow(B[(row * B_inc1 + B_start1) + (col * B_inc2 + B_start2) *  B_internal_size1], \n");
              source.append("              C[(row * C_inc1 + C_start1) + (col * C_inc2 + C_start2) *  C_internal_size1]); \n");
            }
            source.append("  } else if (op_type == 1) {");
            source.append("    for (unsigned int col = col_gid; col < A_size2; col += get_num_groups(0))\n");
            source.append("      for (unsigned int row = row_gid; row < A_size1; row += get_local_size(0))\n");
            source.append("        A[(row * A_inc1 + A_start1) + (col * A_inc2 + A_start2) *  A_internal_size1] =  \n");
            source.append("          B[(row * B_inc1 + B_start1) + (col * B_inc2 + B_start2) *  B_internal_size1] / \n");
            source.append("          C[(row * C_inc1 + C_start1) + (col * C_inc2 + C_start2) *  C_internal_size1]; \n");
            source.append("  } else if (op_type == 0) {");
            source.append("    for (unsigned int col = col_gid; col < A_size2; col += get_num_groups(0))\n");
            source.append("      for (unsigned int row = row_gid; row < A_size1; row += get_local_size(0))\n");
            source.append("        A[(row * A_inc1 + A_start1) + (col * A_inc2 + A_start2) *  A_internal_size1] = \n");
            source.append("          B[(row * B_inc1 + B_start1) + (col * B_inc2 + B_start2) *  B_internal_size1] * \n");
            source.append("          C[(row * C_inc1 + C_start1) + (col * C_inc2 + C_start2) *  C_internal_size1]; \n");
            source.append("  }");
          }
          source.append("} \n");
        }


        template <typename StringType>
        void generate_fft(StringType & source, std::string const & numeric_string, bool is_row_major)
        {
          // naive fourier transform (quadratic complexity, use for reference only)
          source.append("__kernel void fft_direct(__global "); source.append(numeric_string); source.append("2 *input, \n");
          source.append("                         __global "); source.append(numeric_string); source.append("2 *output, \n");
          source.append("                         unsigned int size, \n");
          source.append("                         unsigned int stride, \n");
          source.append("                         unsigned int batch_num, \n");
          source.append("                         "); source.append(numeric_string); source.append(" sign) { \n");
          source.append("    const "); source.append(numeric_string); source.append(" NUM_PI = 3.14159265358979323846; \n");
          source.append(" \n");
          source.append("    for(unsigned int batch_id = 0; batch_id < batch_num; batch_id++) { \n");
          source.append("        for(unsigned int k = get_global_id(0); k < size; k += get_global_size(0)) { \n");
          source.append("            "); source.append(numeric_string); source.append("2 f = 0.0f; \n");
          source.append(" \n");
          source.append("            for(unsigned int n = 0; n < size; n++) { \n");
          source.append("                "); source.append(numeric_string); source.append("2 in = ");
          if (is_row_major)
            source.append("input[batch_id * stride + n]; \n"); //input index here
          else
            source.append("input[n * stride + batch_id]; \n"); //input index here
          source.append(" \n");
          source.append("                "); source.append(numeric_string); source.append(" sn, cs; \n");
          source.append("                "); source.append(numeric_string); source.append(" arg = sign * 2 * NUM_PI * k / size * n; \n");
          source.append("                sn = sincos(arg, &cs); \n");
          source.append(" \n");
          source.append("                "); source.append(numeric_string); source.append("2 ex = ("); source.append(numeric_string); source.append("2)(cs, sn); \n");
          source.append("                f = f + ("); source.append(numeric_string); source.append("2)(in.x * ex.x - in.y * ex.y, in.x * ex.y + in.y * ex.x); \n");
          source.append("            } \n");
          source.append(" \n");
          if (is_row_major)
            source.append("            output[batch_id * stride + k] = f; \n"); // output index here
          else
            source.append("            output[k * stride + batch_id] = f; \n"); // output index here
          source.append("        } \n");
          source.append("    } \n");
          source.append("} \n");

          source.append(" \n"); //////////////////////////////

          source.append("__kernel void fft_radix2(__global "); source.append(numeric_string); source.append("2* input, \n");
          source.append("                         unsigned int s, \n");
          source.append("                         unsigned int bit_size, \n");
          source.append("                         unsigned int size, \n");
          source.append("                         unsigned int stride, \n");
          source.append("                         unsigned int batch_num, \n");
          source.append("                         "); source.append(numeric_string); source.append(" sign) { \n");
          source.append(" \n");
          source.append("    unsigned int ss = 1 << s; \n");
          source.append("    unsigned int half_size = size >> 1; \n");
          source.append(" \n");
          source.append("    "); source.append(numeric_string); source.append(" cs, sn; \n");
          source.append("    const "); source.append(numeric_string); source.append(" NUM_PI = 3.14159265358979323846; \n");
          source.append(" \n");
          source.append("    unsigned int glb_id = get_global_id(0); \n");
          source.append("    unsigned int glb_sz = get_global_size(0); \n");

          source.append("    for(unsigned int batch_id = 0; batch_id < batch_num; batch_id++) { \n");
          source.append("        for(unsigned int tid = glb_id; tid < half_size; tid += glb_sz) { \n");
          source.append("            unsigned int group = (tid & (ss - 1)); \n");
          source.append("            unsigned int pos = ((tid >> s) << (s + 1)) + group; \n");

          if (is_row_major)
          {
            source.append("            unsigned int offset = batch_id * stride + pos; \n");
            source.append("            "); source.append(numeric_string); source.append("2 in1 = input[offset]; \n"); //index
            source.append("            "); source.append(numeric_string); source.append("2 in2 = input[offset + ss]; \n");//index
          }
          else
          {
            source.append("            unsigned int offset = pos * stride + batch_id; \n");
            source.append("            "); source.append(numeric_string); source.append("2 in1 = input[offset]; \n"); //index
            source.append("            "); source.append(numeric_string); source.append("2 in2 = input[offset + ss * stride]; \n");//index
          }

          source.append("            "); source.append(numeric_string); source.append(" arg = group * sign * NUM_PI / ss; \n");

          source.append("            sn = sincos(arg, &cs); \n");

          source.append("            "); source.append(numeric_string); source.append("2 ex = ("); source.append(numeric_string); source.append("2)(cs, sn); \n");

          source.append("            "); source.append(numeric_string); source.append("2 tmp = ("); source.append(numeric_string); source.append("2)(in2.x * ex.x - in2.y * ex.y, in2.x * ex.y + in2.y * ex.x); \n");

          if (is_row_major)
            source.append("            input[offset + ss] = in1 - tmp; \n");//index
          else
            source.append("            input[offset + ss * stride] = in1 - tmp; \n");//index
          source.append("            input[offset] = in1 + tmp; \n");//index
          source.append("        } \n");
          source.append("    } \n");
          source.append("} \n");

          source.append(" \n"); //////////////////////////////

          source.append(" unsigned int get_reorder_num(unsigned int v, unsigned int bit_size) { \n");
          source.append("     v = ((v >> 1) & 0x55555555) | ((v & 0x55555555) << 1); \n");
          source.append("     v = ((v >> 2) & 0x33333333) | ((v & 0x33333333) << 2); \n");
          source.append("     v = ((v >> 4) & 0x0F0F0F0F) | ((v & 0x0F0F0F0F) << 4); \n");
          source.append("     v = ((v >> 8) & 0x00FF00FF) | ((v & 0x00FF00FF) << 8); \n");
          source.append("     v = (v >> 16) | (v << 16); \n");
          source.append("  \n");
          source.append("     v = v >> (32 - bit_size); \n");
          source.append("  \n");
          source.append("     return v; \n");
          source.append(" } \n");

          source.append(" __kernel void fft_radix2_local(__global "); source.append(numeric_string); source.append("2* input, \n");
          source.append("                                 __local "); source.append(numeric_string); source.append("2* lcl_input, \n");
          source.append("                                 unsigned int bit_size, \n");
          source.append("                                 unsigned int size, \n");
          source.append("                                 unsigned int stride, \n");
          source.append("                                 unsigned int batch_num, \n");
          source.append("                                 "); source.append(numeric_string); source.append(" sign) { \n");

          source.append("     unsigned int grp_id = get_group_id(0); \n");
          source.append("     unsigned int grp_num = get_num_groups(0); \n");

          source.append("     unsigned int lcl_sz = get_local_size(0); \n");
          source.append("     unsigned int lcl_id = get_local_id(0); \n");
          source.append("     const "); source.append(numeric_string); source.append(" NUM_PI = 3.14159265358979323846; \n");

          source.append("     for(unsigned int batch_id = grp_id; batch_id < batch_num; batch_id += grp_num) { \n");
                  //unsigned int base_offset = stride * batch_id; \n");
                  //copy chunk of global memory to local \n");
          source.append("         for(unsigned int p = lcl_id; p < size; p += lcl_sz) { \n");
          source.append("             unsigned int v = get_reorder_num(p, bit_size); \n");
          if (is_row_major)
            source.append("             lcl_input[v] = input[batch_id * stride + p]; \n"); //index
          else
            source.append("             lcl_input[v] = input[p * stride + batch_id]; \n"); //index
          source.append("         } \n");

          source.append("         barrier(CLK_LOCAL_MEM_FENCE); \n");

                  //performs Cooley-Tukey FFT on local array
          source.append("         for(unsigned int s = 0; s < bit_size; s++) { \n");
          source.append("             unsigned int ss = 1 << s; \n");

          source.append("             "); source.append(numeric_string); source.append(" cs, sn; \n");

          source.append("             for(unsigned int tid = lcl_id; tid < size; tid += lcl_sz) { \n");
          source.append("                 unsigned int group = (tid & (ss - 1)); \n");
          source.append("                 unsigned int pos = ((tid >> s) << (s + 1)) + group; \n");

          source.append("                 "); source.append(numeric_string); source.append("2 in1 = lcl_input[pos]; \n");
          source.append("                 "); source.append(numeric_string); source.append("2 in2 = lcl_input[pos + ss]; \n");

          source.append("                 "); source.append(numeric_string); source.append(" arg = group * sign * NUM_PI / ss; \n");

          source.append("                 sn = sincos(arg, &cs); \n");
          source.append("                 "); source.append(numeric_string); source.append("2 ex = ("); source.append(numeric_string); source.append("2)(cs, sn); \n");

          source.append("                 "); source.append(numeric_string); source.append("2 tmp = ("); source.append(numeric_string); source.append("2)(in2.x * ex.x - in2.y * ex.y, in2.x * ex.y + in2.y * ex.x); \n");

          source.append("                 lcl_input[pos + ss] = in1 - tmp; \n");
          source.append("                 lcl_input[pos] = in1 + tmp; \n");
          source.append("             } \n");

          source.append("             barrier(CLK_LOCAL_MEM_FENCE); \n");
          source.append("         } \n");

                  //copy local array back to global memory
          source.append("         for(unsigned int p = lcl_id; p < size; p += lcl_sz) { \n");
          if (is_row_major)
            source.append("             input[batch_id * stride + p] = lcl_input[p]; \n");//index
          else
            source.append("             input[p * stride + batch_id] = lcl_input[p]; \n");//index
          source.append("         } \n");
          source.append("     } \n");
          source.append(" } \n");

          source.append(" \n"); //////////////////////////////

          //
          // Performs reordering of input data in bit-reversal order
          // Probably it's better to do in host side,
          //
          source.append("unsigned int get_reorder_num_2(unsigned int v, unsigned int bit_size) { \n");
          source.append("    v = ((v >> 1) & 0x55555555) | ((v & 0x55555555) << 1); \n");
          source.append("    v = ((v >> 2) & 0x33333333) | ((v & 0x33333333) << 2); \n");
          source.append("    v = ((v >> 4) & 0x0F0F0F0F) | ((v & 0x0F0F0F0F) << 4); \n");
          source.append("    v = ((v >> 8) & 0x00FF00FF) | ((v & 0x00FF00FF) << 8); \n");
          source.append("    v = (v >> 16) | (v << 16); \n");

          source.append("    v = v >> (32 - bit_size); \n");

          source.append("    return v; \n");
          source.append("} \n");

          source.append("__kernel void fft_reorder(__global "); source.append(numeric_string); source.append("2* input, \n");
          source.append("                          unsigned int bit_size, \n");
          source.append("                          unsigned int size, \n");
          source.append("                          unsigned int stride, \n");
          source.append("                          int batch_num) { \n");

          source.append("    unsigned int glb_id = get_global_id(0); \n");
          source.append("    unsigned int glb_sz = get_global_size(0); \n");

          source.append("    for(unsigned int batch_id = 0; batch_id < batch_num; batch_id++) { \n");
          source.append("        for(unsigned int i = glb_id; i < size; i += glb_sz) { \n");
          source.append("            unsigned int v = get_reorder_num_2(i, bit_size); \n");

          source.append("            if(i < v) {\n");
          if (is_row_major)
          {
            source.append("                "); source.append(numeric_string); source.append("2 tmp = input[batch_id * stride + i]; \n"); // index
            source.append("                input[batch_id * stride + i] = input[batch_id * stride + v]; \n"); //index
            source.append("                input[batch_id * stride + v] = tmp; \n"); //index
          }
          else
          {
            source.append("                "); source.append(numeric_string); source.append("2 tmp = input[i * stride + batch_id]; \n"); // index
            source.append("                input[i * stride + batch_id] = input[v * stride + batch_id]; \n"); //index
            source.append("                input[v * stride + batch_id] = tmp; \n"); //index
          }
          source.append("            } \n");
          source.append("        } \n");
          source.append("    } \n");
          source.append("} \n");
        }

        template <typename StringType>
        void generate_lu(StringType & source, std::string const & numeric_string, bool is_row_major)
        {
          source.append("__kernel void lu_factorize( \n");
          source.append("          __global "); source.append(numeric_string); source.append(" * matrix, \n");
          source.append("          unsigned int matrix_rows, \n");
          source.append("          unsigned int matrix_cols, \n");
          source.append("          unsigned int matrix_internal_rows, \n");
          source.append("          unsigned int matrix_internal_cols) \n");
          source.append("{ \n");
          source.append("  "); source.append(numeric_string); source.append(" temp; \n");

          if (is_row_major)
          {
            source.append("  unsigned rowi; \n");
            source.append("  unsigned rowk; \n");
            source.append("  for (unsigned int i=1; i<matrix_rows; ++i) \n");
            source.append("  { \n");
            source.append("    rowi = i * matrix_internal_cols; \n");
            source.append("    for (unsigned int k=0; k<i; ++k) \n");
            source.append("    { \n");
            source.append("      rowk = k * matrix_internal_cols; \n");
            source.append("      if (get_global_id(0) == 0) \n");
            source.append("        matrix[rowi + k] /= matrix[rowk + k]; \n");

            source.append("      barrier(CLK_GLOBAL_MEM_FENCE); \n");
            source.append("      temp = matrix[rowi + k]; \n");

            //parallel subtraction:
            source.append("      for (unsigned int j=k+1 + get_global_id(0); j<matrix_rows; j += get_global_size(0)) \n");
            source.append("        matrix[rowi + j] -= temp * matrix[rowk + j]; \n");
          }
          else
          {
            source.append("      for (unsigned int i=1; i<matrix_rows; ++i) \n");
            source.append("      { \n");
            source.append("        for (unsigned int k=0; k<i; ++k) \n");
            source.append("        { \n");

            source.append("          if (get_global_id(0) == 0) \n");
            source.append("            matrix[i + k*matrix_internal_rows] /= matrix[k + k*matrix_internal_rows]; \n");

            source.append("          barrier(CLK_GLOBAL_MEM_FENCE); \n");
            source.append("          temp = matrix[i + k*matrix_internal_rows]; \n");

            //parallel subtraction:
            source.append("          for (unsigned int j=k+1 + get_global_id(0); j<matrix_cols; j += get_global_size(0)) \n");
            source.append("            matrix[i + j*matrix_internal_rows] -= temp * matrix[k + j*matrix_internal_rows]; \n");
          }
          source.append("   }");
          source.append("  }");
          source.append("}");
        }


        template <typename StringType>
        void generate_scaled_rank1_update(StringType & source, std::string const & numeric_string, bool is_row_major, bool alpha_on_cpu)
        {
          source.append("__kernel void scaled_rank1_update_"); alpha_on_cpu ? source.append("cpu") : source.append("gpu"); source.append("( \n");
          source.append("  __global "); source.append(numeric_string); source.append(" * A, \n");
          source.append("  unsigned int A_start1, unsigned int A_start2, \n");
          source.append("  unsigned int A_inc1,   unsigned int A_inc2, \n");
          source.append("  unsigned int A_size1,  unsigned int A_size2, \n");
          source.append("  unsigned int A_internal_size1,  unsigned int A_internal_size2, \n");

          if (alpha_on_cpu) {
            source.append("  "); source.append(numeric_string); source.append(" val, \n");
          } else {
            source.append("  __global const "); source.append(numeric_string); source.append(" *val, \n");
          }
          source.append("  unsigned int options2, \n");

          source.append("  __global const "); source.append(numeric_string); source.append(" * vec1, \n");
          source.append("  unsigned int start1, \n");
          source.append("  unsigned int inc1, \n");
          source.append("  unsigned int size1, \n");

          source.append("  __global const "); source.append(numeric_string); source.append(" * vec2, \n");
          source.append("  unsigned int start2, \n");
          source.append("  unsigned int inc2, \n");
          source.append("  unsigned int size2) \n");
          source.append("{ \n");

          if (alpha_on_cpu) {
            source.append("  "); source.append(numeric_string); source.append(" alpha = val; \n");
          } else {
            source.append("  "); source.append(numeric_string); source.append(" alpha = val[0]; \n");
          }
          source.append("  if (options2 & (1 << 0)) \n");
          source.append("    alpha = -alpha; \n");

          source.append("  unsigned int row_gid = get_global_id(0) / get_local_size(0); \n");
          source.append("  unsigned int col_gid = get_global_id(0) % get_local_size(0); \n");

          source.append("  for (unsigned int row = row_gid; row < A_size1; row += get_num_groups(0)) \n");
          source.append("  { \n");
          source.append("    "); source.append(numeric_string); source.append(" tmp = vec1[row * inc1 + start1];");
          source.append("    tmp = (options2 & (1 << 1)) ? tmp / alpha : tmp * alpha;");
          source.append("    for (unsigned int col = col_gid; col < A_size2; col += get_local_size(0)) \n");
          if (is_row_major)
            source.append("      A[(row * A_inc1 + A_start1) * A_internal_size2 + col * A_inc2 + A_start2] += tmp * vec2[col * inc2 + start2]; \n");
          else
            source.append("      A[(row * A_inc1 + A_start1) + (col * A_inc2 + A_start2) * A_internal_size1] += tmp * vec2[col * inc2 + start2]; \n");
          source.append("  } \n");
          source.append("} \n");
        }

        template <typename StringType>
        void generate_trans_vec_mul(StringType & source, std::string const & numeric_string, bool is_row_major)
        {
          source.append("__kernel void trans_vec_mul( \n");
          source.append("          __global const "); source.append(numeric_string); source.append(" * A, \n");
          source.append("          unsigned int A_row_start, unsigned int A_col_start, \n");
          source.append("          unsigned int A_row_inc, unsigned int A_col_inc, \n");
          source.append("          unsigned int A_row_size, unsigned int A_col_size, \n");
          source.append("          unsigned int A_internal_rows, unsigned int A_internal_cols, \n");
          source.append("          __global const "); source.append(numeric_string); source.append(" * v, \n");
          source.append("          unsigned int v_start, unsigned int v_inc, unsigned int v_size, \n");
          source.append("          __global "); source.append(numeric_string); source.append(" * result, \n");
          source.append("          unsigned int result_start, unsigned int result_inc, unsigned int result_size, \n");
          source.append("          __local "); source.append(numeric_string); source.append(" * work) \n");
          source.append("{ \n");
          if (is_row_major)
          {
            source.append("  for (unsigned int row = get_global_id(0); row < A_col_size; row += get_global_size(0)) \n");
            source.append("  { \n");
            source.append("    "); source.append(numeric_string); source.append(" dot_prod = 0; \n");
            source.append("    for (unsigned int col = 0; col < A_row_size; ++col) \n");
            source.append("      dot_prod += A[(row * A_col_inc + A_col_start) + (col * A_row_inc + A_row_start) * A_internal_cols] * v[v_start + v_inc * col]; \n");
            source.append("    result[row * result_inc + result_start] = dot_prod; \n");
          }
          else
          {
            source.append("  unsigned int row_gid = get_global_id(0) / get_local_size(0); \n");
            source.append("  unsigned int col_gid = get_global_id(0) % get_local_size(0); \n");
            source.append("  unsigned int lid = get_local_id(0); \n");

            source.append("  for (unsigned int row = row_gid; row < A_col_size; row += get_num_groups(0)) \n");
            source.append("  { \n");
            source.append("    "); source.append(numeric_string); source.append(" dot_prod = 0; \n");
            source.append("    for (unsigned int col = col_gid; col < A_row_size; col+=get_local_size(0)) \n");
            source.append("      dot_prod += A[(row * A_col_inc + A_col_start) * A_internal_rows + col * A_row_inc + A_row_start] * v[v_start + v_inc * col]; \n");
            source.append("    work[lid] = dot_prod; \n");

            source.append("    for(unsigned int stride=get_local_size(0)/2 ; stride>0 ; stride>>=1){ \n");
            source.append("      barrier(CLK_LOCAL_MEM_FENCE); \n");
            source.append("      if(lid < stride) \n");
            source.append("        work[lid] += work[lid+stride]; \n");
            source.append("    } \n");

            source.append("    if(lid == 0) \n");
            source.append("      result[row * result_inc + result_start] = work[0]; \n");
          }
          source.append("  } \n");
          source.append("} \n");
        }

        template <typename StringType>
        void generate_triangular_substitute_inplace(StringType & source, std::string const & numeric_string, bool is_row_major)
        {
          source.append("__kernel void triangular_substitute_inplace( \n");
          source.append("          __global "); source.append(numeric_string); source.append(" * A, \n");
          source.append("          unsigned int A_start1, unsigned int A_start2, \n");
          source.append("          unsigned int A_inc1,   unsigned int A_inc2, \n");
          source.append("          unsigned int A_size1,  unsigned int A_size2, \n");
          source.append("          unsigned int A_internal_size1,  unsigned int A_internal_size2, \n");
          source.append("          __global "); source.append(numeric_string); source.append(" * v, \n");
          source.append("          unsigned int v_start, \n");
          source.append("          unsigned int v_inc, \n");
          source.append("          unsigned int v_size, \n");
          source.append("          unsigned int options) \n");
          source.append("{ \n");
          source.append("  "); source.append(numeric_string); source.append(" temp; \n");
          source.append("  unsigned int unit_diagonal_flag  = (options & (1 << 0)); \n");
          source.append("  unsigned int transposed_access_A = (options & (1 << 1)); \n");
          source.append("  unsigned int is_lower_solve      = (options & (1 << 2)); \n");
          source.append("  unsigned int row; \n");
          source.append("  for (unsigned int rows_processed = 0; rows_processed < A_size1; ++rows_processed)  \n");   //Note: A required to be square
          source.append("  { \n");
          source.append("    row = is_lower_solve ? rows_processed : ((A_size1 - rows_processed) - 1); \n");
          source.append("    if (!unit_diagonal_flag) \n");
          source.append("    { \n");
          source.append("      barrier(CLK_GLOBAL_MEM_FENCE); \n");
          source.append("      if (get_global_id(0) == 0) \n");
          if (is_row_major)
            source.append("        v[row * v_inc + v_start] /= A[(row * A_inc1 + A_start1) * A_internal_size2 + (row * A_inc2 + A_start2)]; \n");
          else
            source.append("        v[row * v_inc + v_start] /= A[(row * A_inc1 + A_start1) + (row * A_inc2 + A_start2) * A_internal_size1]; \n");
          source.append("   } \n");

          source.append("    barrier(CLK_GLOBAL_MEM_FENCE); \n");

          source.append("    temp = v[row * v_inc + v_start]; \n");

          source.append("    for (int elim = (is_lower_solve ? (row + get_global_id(0) + 1) : get_global_id(0)); \n");
          source.append("             elim < (is_lower_solve ? A_size1 : row); \n");
          source.append("             elim += get_global_size(0)) \n");
          if (is_row_major)
          {
            source.append("      v[elim * v_inc + v_start] -= temp * A[transposed_access_A ? ((row  * A_inc1 + A_start1) * A_internal_size2 + (elim * A_inc2 + A_start2)) \n");
            source.append("                                                                : ((elim * A_inc1 + A_start1) * A_internal_size2 + (row  * A_inc2 + A_start2))]; \n");
          }
          else
          {
            source.append("      v[elim * v_inc + v_start] -= temp * A[transposed_access_A ? ((row  * A_inc1 + A_start1) + (elim * A_inc2 + A_start2) * A_internal_size1) \n");
            source.append("                                                                : ((elim * A_inc1 + A_start1) + (row  * A_inc2 + A_start2) * A_internal_size1)]; \n");
          }
          source.append("  } \n");
          source.append("} \n");
        }

        template <typename StringType>
        void generate_vec_mul(StringType & source, std::string const & numeric_string, bool is_row_major)
        {
          source.append("__kernel void vec_mul( \n");
          source.append("          __global const "); source.append(numeric_string); source.append(" * A, \n");
          source.append("          unsigned int A_row_start, unsigned int A_col_start, \n");
          source.append("          unsigned int A_row_inc, unsigned int A_col_inc, \n");
          source.append("          unsigned int A_row_size, unsigned int A_col_size, \n");
          source.append("          unsigned int A_internal_rows, unsigned int A_internal_cols, \n");
          source.append("          __global const "); source.append(numeric_string); source.append(" * v, \n");
          source.append("          unsigned int v_start, unsigned int v_inc, unsigned int v_size, \n");
          source.append("          __global "); source.append(numeric_string); source.append(" * result, \n");
          source.append("          unsigned int result_start, unsigned int result_inc, unsigned int result_size, \n");
          source.append("          __local "); source.append(numeric_string); source.append(" * work) \n");
          source.append("{ \n");
          if (is_row_major)
          {
            source.append("  unsigned int row_gid = get_global_id(0) / get_local_size(0); \n");
            source.append("  unsigned int col_gid = get_global_id(0) % get_local_size(0); \n");
            source.append("  unsigned int lid = get_local_id(0); \n");

            source.append("  for (unsigned int row = row_gid; row < A_row_size; row += get_num_groups(0)) \n");
            source.append("  { \n");
            source.append("    "); source.append(numeric_string); source.append(" dot_prod = 0; \n");
            source.append("    for (unsigned int col = col_gid; col < A_col_size; col+=get_local_size(0)) \n");
            source.append("      dot_prod += A[(row * A_row_inc + A_row_start) * A_internal_cols + col * A_col_inc + A_col_start] * v[v_start + v_inc * col]; \n");
            source.append("    work[lid] = dot_prod; \n");

            source.append("    for(unsigned int stride=get_local_size(0)/2 ; stride>0 ; stride>>=1){ \n");
            source.append("      barrier(CLK_LOCAL_MEM_FENCE); \n");
            source.append("      if(lid < stride) \n");
            source.append("        work[lid] += work[lid+stride]; \n");
            source.append("    } \n");

            source.append("    if(lid == 0) \n");
            source.append("      result[row * result_inc + result_start] = work[0]; \n");

          }
          else
          {
            source.append("    for (unsigned int row = get_global_id(0); row < A_row_size; row += get_global_size(0)) \n");
            source.append("    { \n");
            source.append("      "); source.append(numeric_string); source.append(" dot_prod = 0; \n");
            source.append("      for (unsigned int col = 0; col < A_col_size; ++col) \n");
            source.append("        dot_prod += A[(row * A_row_inc + A_row_start) + (col * A_col_inc + A_col_start) * A_internal_rows] * v[v_start + v_inc * col]; \n");
            source.append("      result[row * result_inc + result_start] = dot_prod; \n");
          }
          source.append("  } \n");
          source.append("} \n");
        }

        namespace detail
        {
          inline std::string type_to_string(viennacl::row_major)    { return "row"; }
          inline std::string type_to_string(viennacl::column_major) { return "col"; }
        }

        //////////////////////////// Part 2: Main kernel class ////////////////////////////////////

        // main kernel class
        /** @brief Main kernel class for generating OpenCL kernels for operations on/with dense matrix objects of type viennacl::matrix<>. */
        template <typename NumericT, typename F>
        struct matrix
        {
          static std::string program_name()
          {
            return viennacl::ocl::type_to_string<NumericT>::apply() + "_matrix_" + detail::type_to_string(F());
          }

          static void init(viennacl::ocl::context & ctx)
          {
            viennacl::ocl::DOUBLE_PRECISION_CHECKER<NumericT>::apply(ctx);
            std::string numeric_string = viennacl::ocl::type_to_string<NumericT>::apply();
            bool is_row_major = viennacl::is_row_major<F>::value;

            static std::map<cl_context, bool> init_done;
            if (!init_done[ctx.handle().get()])
            {
              std::string source;
              source.reserve(8192);

              viennacl::ocl::append_double_precision_pragma<NumericT>(ctx, source);

              // fully parametrized kernels:
              generate_ambm(source, numeric_string, is_row_major);

              // kernels with mostly predetermined skeleton:
              generate_assign_cpu(source, numeric_string, is_row_major);
              generate_diagonal_assign_cpu(source, numeric_string, is_row_major);
              generate_element_op(source, numeric_string, is_row_major);
              generate_scaled_rank1_update(source, numeric_string, is_row_major, true);
              generate_scaled_rank1_update(source, numeric_string, is_row_major, false);
              generate_trans_vec_mul(source, numeric_string, is_row_major);
              generate_vec_mul(source, numeric_string, is_row_major);

              if (numeric_string == "float" || numeric_string == "double")
              {
                generate_fft(source, numeric_string, is_row_major);
                generate_lu(source, numeric_string, is_row_major);
                generate_triangular_substitute_inplace(source, numeric_string, is_row_major);
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

