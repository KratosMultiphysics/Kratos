#ifndef VEXCL_BACKEND_JIT_DEVICE_VECTOR_HPP
#define VEXCL_BACKEND_JIT_DEVICE_VECTOR_HPP

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
 * \file   vexcl/backend/jit/device_vector.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Wrapper around std::vector for the JIT backend.
 */

#include <vector>
#include <memory>

namespace vex {
namespace backend {
namespace jit {

typedef unsigned mem_flags;

static const mem_flags MEM_READ_ONLY  = 1;
static const mem_flags MEM_WRITE_ONLY = 2;
static const mem_flags MEM_READ_WRITE = 4;

namespace detail {
struct shared_bytes {
    shared_bytes() : size(0) {}

    shared_bytes(size_t n)
        : data(std::shared_ptr<unsigned char>(new unsigned char[n], std::default_delete<unsigned char[]>())),
          size(n)
    {}

    shared_bytes(const shared_bytes &c)
        : data(c.data), size(c.size)
    {}

    template <class T> T* get() {
        return reinterpret_cast<T*>(data.get());
    }

    template <class T> const T* get() const {
        return reinterpret_cast<const T*>(data.get());
    }

    std::shared_ptr<unsigned char> data;
    size_t size;
};
} // namespace detail

template <typename T>
class device_vector {
    public:
        typedef T  value_type;
        typedef T* raw_type;
        typedef detail::shared_bytes buffer_type;

        device_vector() {}

        device_vector(const command_queue &q, size_t n, const T *host = 0, mem_flags = MEM_READ_WRITE)
            : buffer(sizeof(T) * n)
        {
            if (host) std::copy(host, host + n, buffer.get<T>());
        }

        device_vector(buffer_type buffer) : buffer(buffer) {}

        template <typename U>
        device_vector<U> reinterpret() const {
            return device_vector<U>(buffer);
        }

        void write(const command_queue&, size_t offset, size_t size, const T *host, bool /*blocking*/ = false) const
        {
            std::copy(host, host + size, buffer.get<T>() + offset);
        }

        void read(const command_queue&, size_t offset, size_t size, T *host, bool /*blocking*/ = false) const
        {
            std::copy_n(buffer.get<T>() + offset, size, host);
        }

        size_t size() const {
            return buffer.size / sizeof(T);
        }

        typedef T* mapped_array;

        T* map(const command_queue&) {
            return buffer.data ? buffer.get<T>() : nullptr;
        }

        T* map(const command_queue&) const {
            return buffer.data ? buffer.get<T>() : nullptr;
        }

        const T* raw() const {
            return buffer.data ? buffer.get<T>() : nullptr;
        }

        T* raw() {
            return buffer.data ? buffer.get<T>() : nullptr;
        }

        const buffer_type raw_buffer() const {
            return buffer;
        }
    private:
        mutable buffer_type buffer;
};

} // namespace jit
} // namespace backend
} // namespace vex

#endif
