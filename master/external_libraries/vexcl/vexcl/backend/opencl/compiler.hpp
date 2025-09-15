#ifndef VEXCL_BACKEND_OPENCL_COMPILER_HPP
#define VEXCL_BACKEND_OPENCL_COMPILER_HPP

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
 * \file   vexcl/backend/opencl/compiler.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  OpenCL source code compilation wrapper.
 */

#include <cstdlib>

#include <boost/thread.hpp>

#include <vexcl/backend/common.hpp>
#include <vexcl/detail/backtrace.hpp>

#include <vexcl/backend/opencl/defines.hpp>
#ifdef VEXCL_HAVE_OPENCL_HPP
#  include <CL/opencl.hpp>
#else
#  include <CL/cl2.hpp>
#endif

namespace vex {
namespace backend {
namespace opencl {

/// Saves program binaries for future reuse.
inline void save_program_binaries(
        const std::string &hash, const cl::Program &program
        )
{
    // Prevent writing to the same file by several threads at the same time.
    static boost::mutex mx;
    boost::lock_guard<boost::mutex> lock(mx);

    std::ofstream bfile(program_binaries_path(hash, true) + "kernel", std::ios::binary);
    if (!bfile) return;

    std::vector<size_t> sizes    = program.getInfo<CL_PROGRAM_BINARY_SIZES>();
    std::vector<std::vector<unsigned char>> binaries = program.getInfo<CL_PROGRAM_BINARIES>();

    assert(sizes.size() == 1);

    bfile.write((char*)&sizes[0], sizeof(size_t));
    bfile.write((char*)binaries[0].data(), sizes[0]);
}

/// Tries to read program binaries from file cache.
inline boost::optional<cl::Program> load_program_binaries(
        const std::string &hash, const cl::Context &context,
        const std::vector<cl::Device> &device,
        const std::string &options = ""
        )
{
    std::ifstream bfile(program_binaries_path(hash) + "kernel", std::ios::binary);
    if (!bfile) return boost::optional<cl::Program>();

    size_t n;
    std::vector<std::vector<unsigned char>> buf(1);

    bfile.read((char*)&n, sizeof(size_t));
    buf[0].resize(n);
    bfile.read((char*)buf[0].data(), n);

    cl::Program program(context, device, cl::Program::Binaries(buf));

    try {
        program.build(device, options.c_str());
    } catch(const cl::Error&) {
        std::cerr << "Loading binaries failed:" << std::endl
            << program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(device[0])
            << std::endl;
        return boost::optional<cl::Program>();
    }

    return boost::optional<cl::Program>(program);
}

/// Create and build a program from source string.
/**
 * If VEXCL_CACHE_KERNELS macro is defined, then program binaries are cached
 * in filesystem and reused in the following runs.
 */
inline cl::Program build_sources(
        const cl::CommandQueue &queue, const std::string &source,
        const std::string &options = ""
        )
{
#ifdef VEXCL_SHOW_KERNELS
    std::cout << source << std::endl;
#else
    if (getenv("VEXCL_SHOW_KERNELS"))
        std::cout << source << std::endl;
#endif

    auto context = queue.getInfo<CL_QUEUE_CONTEXT>();
    auto device  = context.getInfo<CL_CONTEXT_DEVICES>();

    std::string compile_options = options + " " + get_compile_options(queue);

#ifdef VEXCL_CACHE_KERNELS
    // Get unique (hopefully) hash string for the kernel.
    std::ostringstream compiler_tag;
    compiler_tag
#if defined(_MSC_VER)
        << "MSC " << _MSC_VER
#elif defined(__clang__)
        << "Clang " << __clang_major__ << "." << __clang_minor__
#elif defined(__GNUC__)
        << "g++ " << __GNUC__ << "." << __GNUC_MINOR__
#else
        << "unknown"
#endif
        ;

    sha1_hasher sha1;
    sha1.process(source)
        .process(cl::Platform(device[0].getInfo<CL_DEVICE_PLATFORM>()).getInfo<CL_PLATFORM_NAME>())
        .process(device[0].getInfo<CL_DEVICE_NAME>())
        .process(compile_options)
        .process(compiler_tag.str())
        ;

    std::string hash = static_cast<std::string>(sha1);

    // Try to get cached program binaries:
    try {
        if (boost::optional<cl::Program> program = load_program_binaries(hash, context, device, compile_options))
            return *program;
    } catch (...) {
        // Shit happens.
        std::cerr << "Failed to load precompiled binaries" << std::endl;
    }
#endif

    // If cache is not available, just compile the sources.
    cl::Program program(context, source);

    try {
        program.build(device, (options + " " + get_compile_options(queue)).c_str());
    } catch(const cl::Error&) {
        std::cerr << source
                  << std::endl
                  << program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(device[0])
                  << std::endl;

        vex::detail::print_backtrace();
        throw;
    }

#ifdef VEXCL_CACHE_KERNELS
    // Save program binaries for future reuse:
    save_program_binaries(hash, program);
#endif

    return program;
}

} // namespace opencl
} // namespace backend
} // namespace vex

#endif
