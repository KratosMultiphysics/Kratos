#ifndef VEXCL_BACKEND_CUDA_ERROR_HPP
#define VEXCL_BACKEND_CUDA_ERROR_HPP

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
 * \file   vexcl/backend/cuda/error.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Output CUDA errors to a std::stream.
 */

#include <iostream>
#include <sstream>
#include <stdexcept>

#include <boost/config.hpp>

#ifdef BOOST_NO_NOEXCEPT
#  define noexcept throw()
#endif

#include <cuda.h>

#include <vexcl/detail/backtrace.hpp>

namespace std {

/// Send human-readable representation of CUresult to the output stream.
inline std::ostream& operator<<(std::ostream &os, CUresult rc) {
    os << "CUDA Driver API Error (";
#define VEXCL_CUDA_ERR2TXT(e) case e: os << static_cast<int>(e) << " - " << #e; break
    switch(rc) {
        VEXCL_CUDA_ERR2TXT(CUDA_SUCCESS);
        VEXCL_CUDA_ERR2TXT(CUDA_ERROR_INVALID_VALUE);
        VEXCL_CUDA_ERR2TXT(CUDA_ERROR_OUT_OF_MEMORY);
        VEXCL_CUDA_ERR2TXT(CUDA_ERROR_NOT_INITIALIZED);
        VEXCL_CUDA_ERR2TXT(CUDA_ERROR_DEINITIALIZED);
        VEXCL_CUDA_ERR2TXT(CUDA_ERROR_PROFILER_DISABLED);
        VEXCL_CUDA_ERR2TXT(CUDA_ERROR_PROFILER_NOT_INITIALIZED);
        VEXCL_CUDA_ERR2TXT(CUDA_ERROR_PROFILER_ALREADY_STARTED);
        VEXCL_CUDA_ERR2TXT(CUDA_ERROR_PROFILER_ALREADY_STOPPED);
        VEXCL_CUDA_ERR2TXT(CUDA_ERROR_NO_DEVICE);
        VEXCL_CUDA_ERR2TXT(CUDA_ERROR_INVALID_DEVICE);
        VEXCL_CUDA_ERR2TXT(CUDA_ERROR_INVALID_IMAGE);
        VEXCL_CUDA_ERR2TXT(CUDA_ERROR_INVALID_CONTEXT);
        VEXCL_CUDA_ERR2TXT(CUDA_ERROR_CONTEXT_ALREADY_CURRENT);
        VEXCL_CUDA_ERR2TXT(CUDA_ERROR_MAP_FAILED);
        VEXCL_CUDA_ERR2TXT(CUDA_ERROR_UNMAP_FAILED);
        VEXCL_CUDA_ERR2TXT(CUDA_ERROR_ARRAY_IS_MAPPED);
        VEXCL_CUDA_ERR2TXT(CUDA_ERROR_ALREADY_MAPPED);
        VEXCL_CUDA_ERR2TXT(CUDA_ERROR_NO_BINARY_FOR_GPU);
        VEXCL_CUDA_ERR2TXT(CUDA_ERROR_ALREADY_ACQUIRED);
        VEXCL_CUDA_ERR2TXT(CUDA_ERROR_NOT_MAPPED);
        VEXCL_CUDA_ERR2TXT(CUDA_ERROR_NOT_MAPPED_AS_ARRAY);
        VEXCL_CUDA_ERR2TXT(CUDA_ERROR_NOT_MAPPED_AS_POINTER);
        VEXCL_CUDA_ERR2TXT(CUDA_ERROR_ECC_UNCORRECTABLE);
        VEXCL_CUDA_ERR2TXT(CUDA_ERROR_UNSUPPORTED_LIMIT);
        VEXCL_CUDA_ERR2TXT(CUDA_ERROR_CONTEXT_ALREADY_IN_USE);
        VEXCL_CUDA_ERR2TXT(CUDA_ERROR_PEER_ACCESS_UNSUPPORTED);
        VEXCL_CUDA_ERR2TXT(CUDA_ERROR_INVALID_SOURCE);
        VEXCL_CUDA_ERR2TXT(CUDA_ERROR_FILE_NOT_FOUND);
        VEXCL_CUDA_ERR2TXT(CUDA_ERROR_SHARED_OBJECT_SYMBOL_NOT_FOUND);
        VEXCL_CUDA_ERR2TXT(CUDA_ERROR_SHARED_OBJECT_INIT_FAILED);
        VEXCL_CUDA_ERR2TXT(CUDA_ERROR_OPERATING_SYSTEM);
        VEXCL_CUDA_ERR2TXT(CUDA_ERROR_INVALID_HANDLE);
        VEXCL_CUDA_ERR2TXT(CUDA_ERROR_NOT_FOUND);
        VEXCL_CUDA_ERR2TXT(CUDA_ERROR_NOT_READY);
        VEXCL_CUDA_ERR2TXT(CUDA_ERROR_LAUNCH_FAILED);
        VEXCL_CUDA_ERR2TXT(CUDA_ERROR_LAUNCH_OUT_OF_RESOURCES);
        VEXCL_CUDA_ERR2TXT(CUDA_ERROR_LAUNCH_TIMEOUT);
        VEXCL_CUDA_ERR2TXT(CUDA_ERROR_LAUNCH_INCOMPATIBLE_TEXTURING);
        VEXCL_CUDA_ERR2TXT(CUDA_ERROR_PEER_ACCESS_ALREADY_ENABLED);
        VEXCL_CUDA_ERR2TXT(CUDA_ERROR_PEER_ACCESS_NOT_ENABLED);
        VEXCL_CUDA_ERR2TXT(CUDA_ERROR_PRIMARY_CONTEXT_ACTIVE);
        VEXCL_CUDA_ERR2TXT(CUDA_ERROR_CONTEXT_IS_DESTROYED);
        VEXCL_CUDA_ERR2TXT(CUDA_ERROR_ASSERT);
        VEXCL_CUDA_ERR2TXT(CUDA_ERROR_TOO_MANY_PEERS);
        VEXCL_CUDA_ERR2TXT(CUDA_ERROR_HOST_MEMORY_ALREADY_REGISTERED);
        VEXCL_CUDA_ERR2TXT(CUDA_ERROR_HOST_MEMORY_NOT_REGISTERED);
        VEXCL_CUDA_ERR2TXT(CUDA_ERROR_NOT_PERMITTED);
        VEXCL_CUDA_ERR2TXT(CUDA_ERROR_NOT_SUPPORTED);
        VEXCL_CUDA_ERR2TXT(CUDA_ERROR_UNKNOWN);
        default:
            os << "Unknown error " << static_cast<int>(rc);
    }
#undef VEXCL_CUDA_ERR2TXT
    return os << ")";
}

} // namespace std

namespace vex {
namespace backend {
namespace cuda {

/// CUDA error class to be thrown as exception.
class error : public std::runtime_error {
    public:
        template <class ErrorCode>
        error(ErrorCode code, const char *file, int line)
            : std::runtime_error(get_msg(code, file, line))
        { }
    private:
        template <class ErrorCode>
        static std::string get_msg(ErrorCode code, const char *file, int line) {
            std::ostringstream s;
            s << file << ":" << line << "\n\t" << code;
            return s.str();
        }
};

inline void check(CUresult rc, const char *file, int line) {
    if (rc != CUDA_SUCCESS) {
        vex::detail::print_backtrace();
	throw error(rc, file, line);
    }
}

/// Throws if rc is not CUDA_SUCCESS.
/**
 * Reports offending file and line number on standard error stream.
 */
#define cuda_check(rc) vex::backend::check(rc, __FILE__, __LINE__)

} // namespace cuda
} // namespace backend
} // namespace vex

namespace std {

/// Sends description of a CUDA error to the output stream.
inline std::ostream& operator<<(std::ostream &os, const vex::backend::error &e) {
    return os << e.what();
}

} // namespace std

#endif
