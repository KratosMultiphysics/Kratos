#ifndef VIENNACL_FFT_HPP
#define VIENNACL_FFT_HPP

/* =========================================================================
   Copyright (c) 2010-2014, Institute for Microelectronics,
                            Institute for Analysis and Scientific Computing,
                            TU Wien.
   Portions of this software are copyright by UChicago Argonne, LLC.

                            -----------------
                  ViennaCL - The Vienna Computing Library
                            -----------------

   Project Head:    Karl Rupp                   rupp@iue.tuwien.ac.at

   (A list of authors and contributors can be found in the PDF manual)

   License:         MIT (X11), see file LICENSE in the base directory
============================================================================= */

/** @file viennacl/fft.hpp
    @brief All routines related to the Fast Fourier Transform. Experimental.
*/

#include <viennacl/vector.hpp>
#include <viennacl/matrix.hpp>

#include "viennacl/linalg/opencl/kernels/fft.hpp"

#include <cmath>

#include <stdexcept>

namespace viennacl
{
  namespace detail
  {
    namespace fft
    {
        const vcl_size_t MAX_LOCAL_POINTS_NUM = 512;

        namespace FFT_DATA_ORDER {
            enum DATA_ORDER {
                ROW_MAJOR,
                COL_MAJOR
            };
        }
    }
  }
}

/// @cond
namespace viennacl
{
  namespace detail
  {
    namespace fft
    {

        inline bool is_radix2(vcl_size_t data_size) {
            return !((data_size > 2) && (data_size & (data_size - 1)));

        }

        inline vcl_size_t next_power_2(vcl_size_t n) {
            n = n - 1;

            vcl_size_t power = 1;

            while(power < sizeof(vcl_size_t) * 8) {
                n = n | (n >> power);
                power *= 2;
            }

            return n + 1;
        }

        inline vcl_size_t num_bits(vcl_size_t size)
        {
            vcl_size_t bits_datasize = 0;
            vcl_size_t ds = 1;

            while(ds < size)
            {
                ds = ds << 1;
                bits_datasize++;
            }

            return bits_datasize;
        }


        /**
         * @brief Direct algorithm for computing Fourier transformation.
         *
         * Works on any sizes of data.
         * Serial implementation has o(n^2) complexity
        */
        template<class SCALARTYPE>
        void direct(const viennacl::ocl::handle<cl_mem>& in,
                    const viennacl::ocl::handle<cl_mem>& out,
                    vcl_size_t size,
                    vcl_size_t stride,
                    vcl_size_t batch_num,
                    SCALARTYPE sign = -1.0f,
                    FFT_DATA_ORDER::DATA_ORDER data_order = FFT_DATA_ORDER::ROW_MAJOR
                    )
        {
          viennacl::ocl::context & ctx = const_cast<viennacl::ocl::context &>(in.context());
          viennacl::linalg::opencl::kernels::fft<SCALARTYPE>::init(ctx);

          std::string program_string = viennacl::linalg::opencl::kernels::matrix<SCALARTYPE, row_major>::program_name();
          if (data_order == FFT_DATA_ORDER::COL_MAJOR)
          {
            viennacl::linalg::opencl::kernels::matrix<SCALARTYPE, column_major>::init(ctx);
            program_string = viennacl::linalg::opencl::kernels::matrix<SCALARTYPE, column_major>::program_name();
          }
          else
            viennacl::linalg::opencl::kernels::matrix<SCALARTYPE, row_major>::init(ctx);
          viennacl::ocl::kernel& kernel = ctx.get_kernel(program_string, "fft_direct");
          viennacl::ocl::enqueue(kernel(in, out, static_cast<cl_uint>(size), static_cast<cl_uint>(stride), static_cast<cl_uint>(batch_num), sign));
        }

        /*
        * This function performs reorder of input data. Indexes are sorted in bit-reversal order.
        * Such reordering should be done before in-place FFT.
        */
        template <typename SCALARTYPE>
        void reorder(const viennacl::ocl::handle<cl_mem>& in,
                     vcl_size_t size,
                     vcl_size_t stride,
                     vcl_size_t bits_datasize,
                     vcl_size_t batch_num,
                     FFT_DATA_ORDER::DATA_ORDER data_order = FFT_DATA_ORDER::ROW_MAJOR
                     )
        {
          viennacl::ocl::context & ctx = const_cast<viennacl::ocl::context &>(in.context());
          viennacl::linalg::opencl::kernels::fft<SCALARTYPE>::init(ctx);

          std::string program_string = viennacl::linalg::opencl::kernels::matrix<SCALARTYPE, row_major>::program_name();
          if (data_order == FFT_DATA_ORDER::COL_MAJOR)
          {
            viennacl::linalg::opencl::kernels::matrix<SCALARTYPE, column_major>::init(ctx);
            program_string = viennacl::linalg::opencl::kernels::matrix<SCALARTYPE, column_major>::program_name();
          }
          else
            viennacl::linalg::opencl::kernels::matrix<SCALARTYPE, row_major>::init(ctx);

          viennacl::ocl::kernel& kernel = ctx.get_kernel(program_string, "fft_reorder");
          viennacl::ocl::enqueue(kernel(in,
                                        static_cast<cl_uint>(bits_datasize),
                                        static_cast<cl_uint>(size),
                                        static_cast<cl_uint>(stride),
                                        static_cast<cl_uint>(batch_num)
                                       )
                                );
        }

        /**
         * @brief Radix-2 algorithm for computing Fourier transformation.
         *
         * Works only on power-of-two sizes of data.
         * Serial implementation has o(n * lg n) complexity.
         * This is a Cooley-Tukey algorithm
        */
        template<class SCALARTYPE>
        void radix2(const viennacl::ocl::handle<cl_mem>& in,
                    vcl_size_t size,
                    vcl_size_t stride,
                    vcl_size_t batch_num,
                    SCALARTYPE sign = -1.0f,
                    FFT_DATA_ORDER::DATA_ORDER data_order = FFT_DATA_ORDER::ROW_MAJOR
                    )
        {
          viennacl::ocl::context & ctx = const_cast<viennacl::ocl::context &>(in.context());
          viennacl::linalg::opencl::kernels::fft<SCALARTYPE>::init(ctx);

            assert(batch_num != 0);
            assert(is_radix2(size));

            std::string program_string = viennacl::linalg::opencl::kernels::matrix<SCALARTYPE, row_major>::program_name();
            if (data_order == FFT_DATA_ORDER::COL_MAJOR)
            {
              viennacl::linalg::opencl::kernels::matrix<SCALARTYPE, column_major>::init(ctx);
              program_string = viennacl::linalg::opencl::kernels::matrix<SCALARTYPE, column_major>::program_name();
            }
            else
              viennacl::linalg::opencl::kernels::matrix<SCALARTYPE, row_major>::init(ctx);

            vcl_size_t bits_datasize = num_bits(size);

            if(size <= MAX_LOCAL_POINTS_NUM)
            {
                viennacl::ocl::kernel& kernel = ctx.get_kernel(program_string, "fft_radix2_local");
                viennacl::ocl::enqueue(kernel(in,
                                              viennacl::ocl::local_mem((size * 4) * sizeof(SCALARTYPE)),
                                              static_cast<cl_uint>(bits_datasize),
                                              static_cast<cl_uint>(size),
                                              static_cast<cl_uint>(stride),
                                              static_cast<cl_uint>(batch_num),
                                              sign));
            }
            else
            {
                reorder<SCALARTYPE>(in, size, stride, bits_datasize, batch_num);

                for(vcl_size_t step = 0; step < bits_datasize; step++)
                {
                    viennacl::ocl::kernel& kernel = ctx.get_kernel(program_string, "fft_radix2");
                    viennacl::ocl::enqueue(kernel(in,
                                                  static_cast<cl_uint>(step),
                                                  static_cast<cl_uint>(bits_datasize),
                                                  static_cast<cl_uint>(size),
                                                  static_cast<cl_uint>(stride),
                                                  static_cast<cl_uint>(batch_num),
                                                  sign));
                }

            }
        }

        /**
         * @brief Bluestein's algorithm for computing Fourier transformation.
         *
         * Currently,  Works only for sizes of input data which less than 2^16.
         * Uses a lot of additional memory, but should be fast for any size of data.
         * Serial implementation has something about o(n * lg n) complexity
        */
        template<class SCALARTYPE, unsigned int ALIGNMENT>
        void bluestein(viennacl::vector<SCALARTYPE, ALIGNMENT>& in,
                       viennacl::vector<SCALARTYPE, ALIGNMENT>& out,
                       vcl_size_t /*batch_num*/)
        {
          viennacl::ocl::context & ctx = const_cast<viennacl::ocl::context &>(viennacl::traits::opencl_handle(in).context());
          viennacl::linalg::opencl::kernels::fft<SCALARTYPE>::init(ctx);

          vcl_size_t size = in.size() >> 1;
          vcl_size_t ext_size = next_power_2(2 * size - 1);

          viennacl::vector<SCALARTYPE, ALIGNMENT> A(ext_size << 1);
          viennacl::vector<SCALARTYPE, ALIGNMENT> B(ext_size << 1);

          viennacl::vector<SCALARTYPE, ALIGNMENT> Z(ext_size << 1);

            {
                viennacl::ocl::kernel& kernel = ctx.get_kernel(viennacl::linalg::opencl::kernels::fft<SCALARTYPE>::program_name(), "zero2");
                viennacl::ocl::enqueue(kernel(
                                            A,
                                            B,
                                            static_cast<cl_uint>(ext_size)
                                            ));

            }
            {
                viennacl::ocl::kernel& kernel = ctx.get_kernel(viennacl::linalg::opencl::kernels::fft<SCALARTYPE>::program_name(), "bluestein_pre");
                viennacl::ocl::enqueue(kernel(
                                           in,
                                           A,
                                           B,
                                           static_cast<cl_uint>(size),
                                           static_cast<cl_uint>(ext_size)
                                       ));
            }

            viennacl::linalg::convolve_i(A, B, Z);

            {
                viennacl::ocl::kernel& kernel = ctx.get_kernel(viennacl::linalg::opencl::kernels::fft<SCALARTYPE>::program_name(), "bluestein_post");
                viennacl::ocl::enqueue(kernel(
                                            Z,
                                            out,
                                            static_cast<cl_uint>(size)
                                            ));
            }
        }

        template<class SCALARTYPE, unsigned int ALIGNMENT>
        void multiply(viennacl::vector<SCALARTYPE, ALIGNMENT> const & input1,
                      viennacl::vector<SCALARTYPE, ALIGNMENT> const & input2,
                      viennacl::vector<SCALARTYPE, ALIGNMENT> & output)
        {
          viennacl::ocl::context & ctx = const_cast<viennacl::ocl::context &>(viennacl::traits::opencl_handle(input1).context());
          viennacl::linalg::opencl::kernels::fft<SCALARTYPE>::init(ctx);
          vcl_size_t size = input1.size() >> 1;
          viennacl::ocl::kernel& kernel = ctx.get_kernel(viennacl::linalg::opencl::kernels::fft<SCALARTYPE>::program_name(), "fft_mult_vec");
          viennacl::ocl::enqueue(kernel(input1, input2, output, static_cast<cl_uint>(size)));
        }

        template<class SCALARTYPE, unsigned int ALIGNMENT>
        void normalize(viennacl::vector<SCALARTYPE, ALIGNMENT> & input)
        {
          viennacl::ocl::context & ctx = const_cast<viennacl::ocl::context &>(viennacl::traits::opencl_handle(input).context());
          viennacl::linalg::opencl::kernels::fft<SCALARTYPE>::init(ctx);

          viennacl::ocl::kernel& kernel = ctx.get_kernel(viennacl::linalg::opencl::kernels::fft<SCALARTYPE>::program_name(), "fft_div_vec_scalar");
          vcl_size_t size = input.size() >> 1;
          SCALARTYPE norm_factor = static_cast<SCALARTYPE>(size);
          viennacl::ocl::enqueue(kernel(input, static_cast<cl_uint>(size), norm_factor));
        }

        template<class SCALARTYPE, unsigned int ALIGNMENT>
        void transpose(viennacl::matrix<SCALARTYPE, viennacl::row_major, ALIGNMENT> & input)
        {
          viennacl::ocl::context & ctx = const_cast<viennacl::ocl::context &>(viennacl::traits::opencl_handle(input).context());
          viennacl::linalg::opencl::kernels::fft<SCALARTYPE>::init(ctx);

          viennacl::ocl::kernel& kernel = ctx.get_kernel(viennacl::linalg::opencl::kernels::fft<SCALARTYPE>::program_name(), "transpose_inplace");
          viennacl::ocl::enqueue(kernel(input,
                                        static_cast<cl_uint>(input.internal_size1()),
                                        static_cast<cl_uint>(input.internal_size2()) >> 1));
        }

        template<class SCALARTYPE, unsigned int ALIGNMENT>
        void transpose(viennacl::matrix<SCALARTYPE, viennacl::row_major, ALIGNMENT> const & input,
                       viennacl::matrix<SCALARTYPE, viennacl::row_major, ALIGNMENT> & output)
        {
          viennacl::ocl::context & ctx = const_cast<viennacl::ocl::context &>(viennacl::traits::opencl_handle(input).context());
          viennacl::linalg::opencl::kernels::fft<SCALARTYPE>::init(ctx);

          viennacl::ocl::kernel& kernel = ctx.get_kernel(viennacl::linalg::opencl::kernels::fft<SCALARTYPE>::program_name(), "transpose");
          viennacl::ocl::enqueue(kernel(input,
                                        output,
                                        static_cast<cl_uint>(input.internal_size1()),
                                        static_cast<cl_uint>(input.internal_size2() >> 1))
                                );
        }

        template<class SCALARTYPE>
        void real_to_complex(viennacl::vector_base<SCALARTYPE> const & in,
                             viennacl::vector_base<SCALARTYPE> & out,
                             vcl_size_t size)
        {
          viennacl::ocl::context & ctx = const_cast<viennacl::ocl::context &>(viennacl::traits::opencl_handle(in).context());
          viennacl::linalg::opencl::kernels::fft<SCALARTYPE>::init(ctx);
          viennacl::ocl::kernel & kernel = ctx.get_kernel(viennacl::linalg::opencl::kernels::fft<SCALARTYPE>::program_name(), "real_to_complex");
          viennacl::ocl::enqueue(kernel(in, out, static_cast<cl_uint>(size)));
        }

        template<class SCALARTYPE>
        void complex_to_real(viennacl::vector_base<SCALARTYPE> const & in,
                             viennacl::vector_base<SCALARTYPE>& out,
                             vcl_size_t size)
        {
          viennacl::ocl::context & ctx = const_cast<viennacl::ocl::context &>(viennacl::traits::opencl_handle(in).context());
          viennacl::linalg::opencl::kernels::fft<SCALARTYPE>::init(ctx);
          viennacl::ocl::kernel& kernel = ctx.get_kernel(viennacl::linalg::opencl::kernels::fft<SCALARTYPE>::program_name(), "complex_to_real");
          viennacl::ocl::enqueue(kernel(in, out, static_cast<cl_uint>(size)));
        }

        template<class SCALARTYPE>
        void reverse(viennacl::vector_base<SCALARTYPE>& in)
        {
          viennacl::ocl::context & ctx = const_cast<viennacl::ocl::context &>(viennacl::traits::opencl_handle(in).context());
          viennacl::linalg::opencl::kernels::fft<SCALARTYPE>::init(ctx);
          vcl_size_t size = in.size();
          viennacl::ocl::kernel& kernel = ctx.get_kernel(viennacl::linalg::opencl::kernels::fft<SCALARTYPE>::program_name(), "reverse_inplace");
          viennacl::ocl::enqueue(kernel(in, static_cast<cl_uint>(size)));
        }


    } //namespace fft
  } //namespace detail

  /**
    * @brief Generic inplace version of 1-D Fourier transformation.
    *
    * @param input       Input vector, result will be stored here.
    * @param batch_num   Number of items in batch
    * @param sign        Sign of exponent, default is -1.0
    */
  template<class SCALARTYPE, unsigned int ALIGNMENT>
  void inplace_fft(viennacl::vector<SCALARTYPE, ALIGNMENT>& input,
            vcl_size_t batch_num = 1,
            SCALARTYPE sign = -1.0)
  {
      vcl_size_t size = (input.size() >> 1) / batch_num;

      if(!viennacl::detail::fft::is_radix2(size))
      {
          viennacl::vector<SCALARTYPE, ALIGNMENT> output(input.size());
          viennacl::detail::fft::direct(viennacl::traits::opencl_handle(input),
                                        viennacl::traits::opencl_handle(output),
                                        size,
                                        size,
                                        batch_num,
                                        sign);

          viennacl::copy(output, input);
      } else {
          viennacl::detail::fft::radix2(viennacl::traits::opencl_handle(input), size, size, batch_num, sign);
      }
  }

  /**
    * @brief Generic version of 1-D Fourier transformation.
    *
    * @param input      Input vector.
    * @param output     Output vector.
    * @param batch_num  Number of items in batch.
    * @param sign       Sign of exponent, default is -1.0
    */
  template<class SCALARTYPE, unsigned int ALIGNMENT>
  void fft(viennacl::vector<SCALARTYPE, ALIGNMENT>& input,
            viennacl::vector<SCALARTYPE, ALIGNMENT>& output,
            vcl_size_t batch_num = 1,
            SCALARTYPE sign = -1.0
            )
  {
      vcl_size_t size = (input.size() >> 1) / batch_num;

      if(viennacl::detail::fft::is_radix2(size))
      {
          viennacl::copy(input, output);
          viennacl::detail::fft::radix2(viennacl::traits::opencl_handle(output), size, size, batch_num, sign);
      } else {
          viennacl::detail::fft::direct(viennacl::traits::opencl_handle(input),
                                        viennacl::traits::opencl_handle(output),
                                        size,
                                        size,
                                        batch_num,
                                        sign);
      }
  }

  /**
    * @brief Generic inplace version of 2-D Fourier transformation.
    *
    * @param input       Input matrix, result will be stored here.
    * @param sign        Sign of exponent, default is -1.0
    */
  template<class SCALARTYPE, unsigned int ALIGNMENT>
  void inplace_fft(viennacl::matrix<SCALARTYPE, viennacl::row_major, ALIGNMENT>& input,
            SCALARTYPE sign = -1.0)
  {
      vcl_size_t rows_num = input.size1();
      vcl_size_t cols_num = input.size2() >> 1;

      vcl_size_t cols_int = input.internal_size2() >> 1;

      // batch with rows
      if(viennacl::detail::fft::is_radix2(cols_num))
      {
          viennacl::detail::fft::radix2(viennacl::traits::opencl_handle(input), cols_num, cols_int, rows_num, sign, viennacl::detail::fft::FFT_DATA_ORDER::ROW_MAJOR);
      }
      else
      {
          viennacl::matrix<SCALARTYPE, viennacl::row_major, ALIGNMENT> output(input.size1(), input.size2());

          viennacl::detail::fft::direct(viennacl::traits::opencl_handle(input),
                                        viennacl::traits::opencl_handle(output),
                                        cols_num,
                                        cols_int,
                                        rows_num,
                                        sign,
                                        viennacl::detail::fft::FFT_DATA_ORDER::ROW_MAJOR
                                        );

          input = output;
      }

      // batch with cols
      if (viennacl::detail::fft::is_radix2(rows_num)) {
          viennacl::detail::fft::radix2(viennacl::traits::opencl_handle(input), rows_num, cols_int, cols_num, sign, viennacl::detail::fft::FFT_DATA_ORDER::COL_MAJOR);
      } else {
          viennacl::matrix<SCALARTYPE, viennacl::row_major, ALIGNMENT> output(input.size1(), input.size2());

          viennacl::detail::fft::direct(viennacl::traits::opencl_handle(input),
                                        viennacl::traits::opencl_handle(output),
                                        rows_num,
                                        cols_int,
                                        cols_num,
                                        sign,
                                        viennacl::detail::fft::FFT_DATA_ORDER::COL_MAJOR);

          input = output;
      }

  }

  /**
    * @brief Generic version of 2-D Fourier transformation.
    *
    * @param input      Input vector.
    * @param output     Output vector.
    * @param sign       Sign of exponent, default is -1.0
    */
  template<class SCALARTYPE, unsigned int ALIGNMENT>
  void fft(viennacl::matrix<SCALARTYPE, viennacl::row_major, ALIGNMENT>& input,
            viennacl::matrix<SCALARTYPE, viennacl::row_major, ALIGNMENT>& output,
            SCALARTYPE sign = -1.0)
  {
      vcl_size_t rows_num = input.size1();
      vcl_size_t cols_num = input.size2() >> 1;

      vcl_size_t cols_int = input.internal_size2() >> 1;

      // batch with rows
      if(viennacl::detail::fft::is_radix2(cols_num))
      {
          output = input;
          viennacl::detail::fft::radix2(viennacl::traits::opencl_handle(output), cols_num, cols_int, rows_num, sign, viennacl::detail::fft::FFT_DATA_ORDER::ROW_MAJOR);
      }
      else
      {
          viennacl::detail::fft::direct(viennacl::traits::opencl_handle(input),
                                        viennacl::traits::opencl_handle(output),
                                        cols_num,
                                        cols_int,
                                        rows_num,
                                        sign,
                                        viennacl::detail::fft::FFT_DATA_ORDER::ROW_MAJOR
                                        );
      }

      // batch with cols
      if(viennacl::detail::fft::is_radix2(rows_num))
      {
          viennacl::detail::fft::radix2(viennacl::traits::opencl_handle(output), rows_num, cols_int, cols_num, sign, viennacl::detail::fft::FFT_DATA_ORDER::COL_MAJOR);
      }
      else
      {
          viennacl::matrix<SCALARTYPE, viennacl::row_major, ALIGNMENT> tmp(output.size1(), output.size2());
          tmp = output;

          viennacl::detail::fft::direct(viennacl::traits::opencl_handle(tmp),
                              viennacl::traits::opencl_handle(output),
                              rows_num,
                              cols_int,
                              cols_num,
                              sign,
                              viennacl::detail::fft::FFT_DATA_ORDER::COL_MAJOR);
      }
  }

  /**
    * @brief Generic inplace version of inverse 1-D Fourier transformation.
    *
    * Shorthand function for fft(sign = 1.0)
    *
    * @param input      Input vector.
    * @param batch_num  Number of items in batch.
    * @param sign       Sign of exponent, default is -1.0
    */
  template<class SCALARTYPE, unsigned int ALIGNMENT>
  void inplace_ifft(viennacl::vector<SCALARTYPE, ALIGNMENT>& input,
            vcl_size_t batch_num = 1)
  {
      viennacl::inplace_fft(input, batch_num, SCALARTYPE(1.0));
      viennacl::detail::fft::normalize(input);
  }

  /**
    * @brief Generic version of inverse 1-D Fourier transformation.
    *
    * Shorthand function for fft(sign = 1.0)
    *
    * @param input      Input vector.
    * @param output     Output vector.
    * @param batch_num  Number of items in batch.
    * @param sign       Sign of exponent, default is -1.0
    */
  template<class SCALARTYPE, unsigned int ALIGNMENT>
  void ifft(viennacl::vector<SCALARTYPE, ALIGNMENT>& input,
            viennacl::vector<SCALARTYPE, ALIGNMENT>& output,
            vcl_size_t batch_num = 1
            )
  {
      viennacl::fft(input, output, batch_num, SCALARTYPE(1.0));
      viennacl::detail::fft::normalize(output);
  }

  namespace linalg
  {
    /**
      * @brief 1-D convolution of two vectors.
      *
      * This function does not make any changes to input vectors
      *
      * @param input1     Input vector #1.
      * @param input2     Input vector #2.
      * @param output     Output vector.
      */
    template<class SCALARTYPE, unsigned int ALIGNMENT>
    void convolve(viennacl::vector<SCALARTYPE, ALIGNMENT>& input1,
                  viennacl::vector<SCALARTYPE, ALIGNMENT>& input2,
                  viennacl::vector<SCALARTYPE, ALIGNMENT>& output
                  )
    {
        assert(input1.size() == input2.size());
        assert(input1.size() == output.size());
        //temporal arrays
        viennacl::vector<SCALARTYPE, ALIGNMENT> tmp1(input1.size());
        viennacl::vector<SCALARTYPE, ALIGNMENT> tmp2(input2.size());
        viennacl::vector<SCALARTYPE, ALIGNMENT> tmp3(output.size());

        // align input arrays to equal size
        // FFT of input data
        viennacl::fft(input1, tmp1);
        viennacl::fft(input2, tmp2);

        // multiplication of input data
        viennacl::detail::fft::multiply(tmp1, tmp2, tmp3);
        // inverse FFT of input data
        viennacl::ifft(tmp3, output);
    }

    /**
      * @brief 1-D convolution of two vectors.
      *
      * This function can make changes to input vectors to avoid additional memory allocations.
      *
      * @param input1     Input vector #1.
      * @param input2     Input vector #2.
      * @param output     Output vector.
      */
    template<class SCALARTYPE, unsigned int ALIGNMENT>
    void convolve_i(viennacl::vector<SCALARTYPE, ALIGNMENT>& input1,
                    viennacl::vector<SCALARTYPE, ALIGNMENT>& input2,
                    viennacl::vector<SCALARTYPE, ALIGNMENT>& output
                    )
    {
        assert(input1.size() == input2.size());
        assert(input1.size() == output.size());

        viennacl::inplace_fft(input1);
        viennacl::inplace_fft(input2);

        viennacl::detail::fft::multiply(input1, input2, output);

        viennacl::inplace_ifft(output);
    }
  }
} //namespace linalg

/// @endcond
#endif
