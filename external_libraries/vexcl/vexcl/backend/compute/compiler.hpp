#ifndef VEXCL_BACKEND_COMPUTE_COMPILER_HPP
#define VEXCL_BACKEND_COMPUTE_COMPILER_HPP

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
 * \file   vexcl/backend/compute/compiler.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  OpenCL source code compilation wrapper for Boost.Compute backend.
 */

#include <cstdlib>

#include <boost/thread.hpp>
#include <boost/compute/core.hpp>

#include <vexcl/backend/common.hpp>
#include <vexcl/detail/backtrace.hpp>

namespace vex {
namespace backend {
namespace compute {

/// Create and build a program from source string.
/**
 * If VEXCL_CACHE_KERNELS macro is defined, then program binaries are cached
 * in filesystem and reused in the following runs.
 */
inline boost::compute::program build_sources(
        const boost::compute::command_queue &queue,
        const std::string &source,
        const std::string &options = ""
        )
{
#ifdef VEXCL_SHOW_KERNELS
    std::cout << source << std::endl;
#else
    if (getenv("VEXCL_SHOW_KERNELS"))
        std::cout << source << std::endl;
#endif

    return boost::compute::program::build_with_source(
            source, queue.get_context(),
            options + " " + get_compile_options(queue)
            );
}

} // namespace compute
} // namespace backend
} // namespace vex

#endif
