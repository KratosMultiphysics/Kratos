#ifndef VEXCL_BACKEND_COMPUTE_HPP
#define VEXCL_BACKEND_COMPUTE_HPP

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
 * \file   vexcl/backend/compute.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  OpenCL backend based on Boost.Compute core API.
 */

#ifndef VEXCL_BACKEND_COMPUTE
#  define VEXCL_BACKEND_COMPUTE
#endif

#include <vexcl/backend/compute/error.hpp>
#include <vexcl/backend/compute/context.hpp>
#include <vexcl/backend/compute/filter.hpp>
#include <vexcl/backend/compute/device_vector.hpp>

// Since Boost.Compute is based on OpenCL,
// we can reuse source generator from the OpenCL backend.
#include <vexcl/backend/opencl/source.hpp>

namespace vex {
    namespace backend {
        namespace compute {
            using vex::backend::opencl::standard_kernel_header;
            using vex::backend::opencl::source_generator;
        }
    }
}

#include <vexcl/backend/compute/compiler.hpp>
#include <vexcl/backend/compute/kernel.hpp>
#include <vexcl/backend/compute/event.hpp>

#endif
