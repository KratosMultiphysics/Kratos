#ifndef VIENNACL_LINALG_OPENCL_KERNELS_ILU_HPP
#define VIENNACL_LINALG_OPENCL_KERNELS_ILU_HPP

#include "viennacl/tools/tools.hpp"
#include "viennacl/ocl/kernel.hpp"
#include "viennacl/ocl/platform.hpp"
#include "viennacl/ocl/utils.hpp"

/** @file viennacl/linalg/opencl/kernels/ilu.hpp
 *  @brief OpenCL kernel file for nonnegative matrix factorization */
namespace viennacl
{
  namespace linalg
  {
    namespace opencl
    {
      namespace kernels
      {
        template <typename StringType>
        void generate_ilu_level_scheduling_substitute(StringType & source, std::string const & numeric_string)
        {
          source.append("__kernel void level_scheduling_substitute( \n");
          source.append("          __global const unsigned int * row_index_array, \n");
          source.append("          __global const unsigned int * row_indices, \n");
          source.append("          __global const unsigned int * column_indices, \n");
          source.append("          __global const "); source.append(numeric_string); source.append(" * elements, \n");
          source.append("          __global "); source.append(numeric_string); source.append(" * vec, \n");
          source.append("          unsigned int size) \n");
          source.append("{ \n");
          source.append("  for (unsigned int row  = get_global_id(0); \n");
          source.append("                    row  < size; \n");
          source.append("                    row += get_global_size(0)) \n");
          source.append("  { \n");
          source.append("    unsigned int eq_row = row_index_array[row]; \n");
          source.append("    "); source.append(numeric_string); source.append(" vec_entry = vec[eq_row]; \n");
          source.append("    unsigned int row_end = row_indices[row+1]; \n");

          source.append("    for (unsigned int j = row_indices[row]; j < row_end; ++j) \n");
          source.append("      vec_entry -= vec[column_indices[j]] * elements[j]; \n");

          source.append("    vec[eq_row] = vec_entry; \n");
          source.append("  } \n");
          source.append("} \n");
        }

        // main kernel class
        /** @brief Main kernel class for generating OpenCL kernels for incomplete LU factorization preconditioners. */
        template <class NumericT>
        struct ilu
        {
          static std::string program_name()
          {
            return viennacl::ocl::type_to_string<NumericT>::apply() + "_ilu";
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

              // only generate for floating points (forces error for integers)
              if (numeric_string == "float" || numeric_string == "double")
              {
                generate_ilu_level_scheduling_substitute(source, numeric_string);
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

