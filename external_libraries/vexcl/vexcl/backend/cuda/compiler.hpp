#ifndef VEXCL_BACKEND_CUDA_COMPILER_HPP
#define VEXCL_BACKEND_CUDA_COMPILER_HPP

/*
The MIT License

Copyright (c) 2012-2018 Denis Demidov <dennis.demidov@gmail.com>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/

/**
 * \file   vexcl/backend/cuda/compiler.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  CUDA source code compilation wrapper.
 */

#include <cstdlib>
#include <cuda.h>

#include <vexcl/backend/common.hpp>
#include <vexcl/detail/backtrace.hpp>

#ifndef VEXCL_NVCC_OPTIONS
#  ifdef NDEBUG
#    define VEXCL_NVCC_OPTIONS "-O3"
#  else
#    define VEXCL_NVCC_OPTIONS "-g"
#  endif
#endif

namespace vex {
namespace backend {
namespace cuda {

/// Create and build a program from source string.
inline vex::backend::program build_sources(
        const command_queue &queue, const std::string &source,
        const std::string &options = ""
        )
{
#ifdef VEXCL_SHOW_KERNELS
    std::cout << source << std::endl;
#else
    if (getenv("VEXCL_SHOW_KERNELS"))
        std::cout << source << std::endl;
#endif

    std::string compile_options = options + " " + get_compile_options(queue);

    queue.context().set_current();

    auto cc = queue.device().compute_capability();
    std::ostringstream ccstr;
    ccstr << std::get<0>(cc) << std::get<1>(cc);

    sha1_hasher sha1;
    sha1.process(source)
        .process(queue.device().name())
        .process(compile_options)
        .process(ccstr.str())
        ;

    std::string hash = static_cast<std::string>(sha1);

    // Write source to a .cu file
    std::string basename = program_binaries_path(hash, true) + "kernel";
    std::string ptxfile  = basename + ".ptx";

    if ( !boost::filesystem::exists(ptxfile) ) {
        std::string cufile = basename + ".cu";

        {
            std::ofstream f(cufile);
            f << source;
        }

        // Compile the source to ptx.
        std::ostringstream cmdline;
        cmdline
            << "nvcc -ptx -arch=sm_" << std::get<0>(cc) << std::get<1>(cc)
            << " " << VEXCL_NVCC_OPTIONS << " " << compile_options
            << " -o " << ptxfile << " " << cufile;
        if (0 != system(cmdline.str().c_str()) ) {
#ifndef VEXCL_SHOW_KERNELS
            std::cerr << source << std::endl;
#endif

            vex::detail::print_backtrace();
            throw std::runtime_error("nvcc invocation failed");
        }
    }

    // Load the compiled ptx.
    CUmodule prg;
    cuda_check( cuModuleLoad(&prg, ptxfile.c_str()) );

    return program(queue.context(), prg);
}

} // namespace cuda
} // namespace backend
} // namespace vex

#endif
