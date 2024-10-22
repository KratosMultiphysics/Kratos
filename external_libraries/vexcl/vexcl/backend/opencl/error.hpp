#ifndef VEXCL_BACKEND_OPENCL_ERROR_HPP
#define VEXCL_BACKEND_OPENCL_ERROR_HPP

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
 * \file   vexcl/backend/opencl/error.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Output OpenCL errors to a std::stream.
 */

#include <iostream>

#include <vexcl/backend/opencl/defines.hpp>
#ifdef VEXCL_HAVE_OPENCL_HPP
#  include <CL/opencl.hpp>
#else
#  include <CL/cl2.hpp>
#endif

namespace vex {
namespace backend {
namespace opencl {

typedef cl::Error error;

} // namespace opencl
} // namespace backend
} // namespace vex

namespace std {

/// Sends description of an OpenCL error to the output stream.
inline std::ostream& operator<<(std::ostream &os, const vex::backend::opencl::error &e) {
    os << e.what() << "(" << e.err() << ": ";

#define VEXCL_CL_ERR2TXT(num, msg) case (num): os << (msg); break

    switch (e.err()) {
        VEXCL_CL_ERR2TXT(  0, "Success");
        VEXCL_CL_ERR2TXT( -1, "Device not found");
        VEXCL_CL_ERR2TXT( -2, "Device not available");
        VEXCL_CL_ERR2TXT( -3, "Compiler not available");
        VEXCL_CL_ERR2TXT( -4, "Mem object allocation failure");
        VEXCL_CL_ERR2TXT( -5, "Out of resources");
        VEXCL_CL_ERR2TXT( -6, "Out of host memory");
        VEXCL_CL_ERR2TXT( -7, "Profiling info not available");
        VEXCL_CL_ERR2TXT( -8, "Mem copy overlap");
        VEXCL_CL_ERR2TXT( -9, "Image format mismatch");
        VEXCL_CL_ERR2TXT(-10, "Image format not supported");
        VEXCL_CL_ERR2TXT(-11, "Build program failure");
        VEXCL_CL_ERR2TXT(-12, "Map failure");
        VEXCL_CL_ERR2TXT(-13, "Misaligned sub buffer offset");
        VEXCL_CL_ERR2TXT(-14, "Exec status error for events in wait list");
        VEXCL_CL_ERR2TXT(-30, "Invalid value");
        VEXCL_CL_ERR2TXT(-31, "Invalid device type");
        VEXCL_CL_ERR2TXT(-32, "Invalid platform");
        VEXCL_CL_ERR2TXT(-33, "Invalid device");
        VEXCL_CL_ERR2TXT(-34, "Invalid context");
        VEXCL_CL_ERR2TXT(-35, "Invalid queue properties");
        VEXCL_CL_ERR2TXT(-36, "Invalid command queue");
        VEXCL_CL_ERR2TXT(-37, "Invalid host ptr");
        VEXCL_CL_ERR2TXT(-38, "Invalid mem object");
        VEXCL_CL_ERR2TXT(-39, "Invalid image format descriptor");
        VEXCL_CL_ERR2TXT(-40, "Invalid image size");
        VEXCL_CL_ERR2TXT(-41, "Invalid sampler");
        VEXCL_CL_ERR2TXT(-42, "Invalid binary");
        VEXCL_CL_ERR2TXT(-43, "Invalid build options");
        VEXCL_CL_ERR2TXT(-44, "Invalid program");
        VEXCL_CL_ERR2TXT(-45, "Invalid program executable");
        VEXCL_CL_ERR2TXT(-46, "Invalid kernel name");
        VEXCL_CL_ERR2TXT(-47, "Invalid kernel definition");
        VEXCL_CL_ERR2TXT(-48, "Invalid kernel");
        VEXCL_CL_ERR2TXT(-49, "Invalid arg index");
        VEXCL_CL_ERR2TXT(-50, "Invalid arg value");
        VEXCL_CL_ERR2TXT(-51, "Invalid arg size");
        VEXCL_CL_ERR2TXT(-52, "Invalid kernel args");
        VEXCL_CL_ERR2TXT(-53, "Invalid work dimension");
        VEXCL_CL_ERR2TXT(-54, "Invalid work group size");
        VEXCL_CL_ERR2TXT(-55, "Invalid work item size");
        VEXCL_CL_ERR2TXT(-56, "Invalid global offset");
        VEXCL_CL_ERR2TXT(-57, "Invalid event wait list");
        VEXCL_CL_ERR2TXT(-58, "Invalid event");
        VEXCL_CL_ERR2TXT(-59, "Invalid operation");
        VEXCL_CL_ERR2TXT(-60, "Invalid gl object");
        VEXCL_CL_ERR2TXT(-61, "Invalid buffer size");
        VEXCL_CL_ERR2TXT(-62, "Invalid mip level");
        VEXCL_CL_ERR2TXT(-63, "Invalid global work size");
        VEXCL_CL_ERR2TXT(-64, "Invalid property");

        default:
            os << "Unknown error";
            break;
    }

#undef VEXCL_CL_ERR2TXT

    return os << ")";
}

}

#endif
