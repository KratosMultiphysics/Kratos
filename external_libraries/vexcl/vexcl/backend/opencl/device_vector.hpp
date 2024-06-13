#ifndef VEXCL_BACKEND_OPENCL_DEVICE_VECTOR_HPP
#define VEXCL_BACKEND_OPENCL_DEVICE_VECTOR_HPP

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
 * \file   vexcl/backend/opencl/device_vector.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  OpenCL device vector.
 */

#include <vexcl/backend/opencl/defines.hpp>
#ifdef VEXCL_HAVE_OPENCL_HPP
#  include <CL/opencl.hpp>
#else
#  include <CL/cl2.hpp>
#endif

namespace vex {
namespace backend {
namespace opencl {

typedef cl_mem_flags mem_flags;

static const mem_flags MEM_READ_ONLY  = CL_MEM_READ_ONLY;
static const mem_flags MEM_WRITE_ONLY = CL_MEM_WRITE_ONLY;
static const mem_flags MEM_READ_WRITE = CL_MEM_READ_WRITE;

template <typename T>
class device_vector {
    public:
        typedef T value_type;
        typedef cl_mem raw_type;

        device_vector() {}

        device_vector(const cl::CommandQueue &q, size_t n,
                const T *host = 0, mem_flags flags = MEM_READ_WRITE)
        {
            if (host && !(flags & CL_MEM_USE_HOST_PTR))
                flags |= CL_MEM_COPY_HOST_PTR;

            if (n)
                buffer = cl::Buffer(q.getInfo<CL_QUEUE_CONTEXT>(), flags,
                        n * sizeof(T), static_cast<void*>(const_cast<T*>(host)));
        }

        device_vector(cl::Buffer buffer) : buffer( std::move(buffer) ) {}

        template <typename U>
        device_vector<U> reinterpret() const {
            return device_vector<U>(buffer);
        }

        void write(const cl::CommandQueue &q, size_t offset, size_t size, const T *host,
                bool blocking = false) const
        {
            if (size)
                q.enqueueWriteBuffer(
                        buffer, blocking ? CL_TRUE : CL_FALSE,
                        sizeof(T) * offset, sizeof(T) * size, host
                        );
        }

        void read(const cl::CommandQueue &q, size_t offset, size_t size, T *host,
                bool blocking = false) const
        {
            if (size)
                q.enqueueReadBuffer(
                        buffer, blocking ? CL_TRUE : CL_FALSE,
                        sizeof(T) * offset, sizeof(T) * size, host
                        );
        }

        size_t size() const {
            if (buffer())
                return buffer.getInfo<CL_MEM_SIZE>() / sizeof(T);
            return 0;
        }

        struct buffer_unmapper {
            const cl::CommandQueue &queue;
            const cl::Buffer       &buffer;

            buffer_unmapper(const cl::CommandQueue &q, const cl::Buffer &b)
                : queue(q), buffer(b)
            {}

            void operator()(T* ptr) const {
                queue.enqueueUnmapMemObject(buffer, ptr);
            }
        };

        typedef std::unique_ptr<T[], buffer_unmapper> mapped_array;

        mapped_array map(const cl::CommandQueue &q) {
            return mapped_array(
                    static_cast<T*>(
                        q.enqueueMapBuffer(
                            buffer, CL_TRUE, CL_MAP_READ | CL_MAP_WRITE,
                            0, size() * sizeof(T)
                            )
                        ),
                    buffer_unmapper(q, buffer)
                    );
        }

        mapped_array map(const cl::CommandQueue &q) const {
            return mapped_array(
                    static_cast<T*>(
                        q.enqueueMapBuffer(
                            buffer, CL_TRUE, CL_MAP_READ,
                            0, size() * sizeof(T)
                            )
                        ),
                    buffer_unmapper(q, buffer)
                    );
        }

        cl_mem raw() const {
            return buffer();
        }

        cl::Buffer raw_buffer() const {
            return buffer;
        }
    private:
        cl::Buffer buffer;
};

} // namespace opencl
} // namespace backend
} // namespace vex

#endif
