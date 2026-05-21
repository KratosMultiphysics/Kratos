#ifndef VEXCL_BACKEND_OPENCL_IMAGE_HPP
#define VEXCL_BACKEND_OPENCL_IMAGE_HPP

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
 * \file   vexcl/backend/opencl/image.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Allow using OpenCL images in vector expressions.
 */

#include <vexcl/operations.hpp>
#include <vexcl/backend/opencl/defines.hpp>
#ifdef VEXCL_HAVE_OPENCL_HPP
#  include <CL/opencl.hpp>
#else
#  include <CL/cl2.hpp>
#endif

namespace vex {

#define VEXCL_ENABLE_CL_IMAGE(dim)                                             \
namespace traits {                                                             \
    template <>                                                                \
    struct is_vector_expr_terminal<cl::Image ## dim ## D> : std::true_type {}; \
}                                                                              \
template <> struct type_name_impl<cl::Image ## dim ## D> {                     \
    static std::string get() { return "image" #dim "d_t"; }                    \
}

VEXCL_ENABLE_CL_IMAGE(1);
VEXCL_ENABLE_CL_IMAGE(2);
VEXCL_ENABLE_CL_IMAGE(3);

#undef VEXCL_ENABLE_CL_IMAGE

} // namespace vex

#endif
