#ifndef VEXCL_BACKEND_OPENCL_SVM_VECTOR_HPP
#define VEXCL_BACKEND_OPENCL_SVM_VECTOR_HPP

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
 * \file   vexcl/backend/opencl/svm_vector.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Shared virtual memory support for OpenCL backend.
 */

#include <memory>

#include <vexcl/backend/opencl/defines.hpp>
#ifdef VEXCL_HAVE_OPENCL_HPP
#  include <CL/opencl.hpp>
#else
#  include <CL/cl2.hpp>
#endif

#include <vexcl/backend/opencl/context.hpp>

namespace vex {

namespace backend { namespace opencl {

typedef cl_map_flags map_flags;

static const map_flags MAP_READ  = CL_MAP_READ;
static const map_flags MAP_WRITE = CL_MAP_WRITE;

} }

template <typename T>
struct svm_ptr {
    T *ptr;

    svm_ptr(T* ptr) : ptr(ptr) {}
};

namespace backend { namespace opencl {

template <typename T>
struct kernel_arg_pusher< svm_ptr<T> > {
    static void set(cl::Kernel &k, unsigned argpos, const svm_ptr<T> &arg) {
        cl_int ret = clSetKernelArgSVMPointer(k(), argpos, arg.ptr);
        if (ret != CL_SUCCESS) throw cl::Error(ret, "clSetKernelArgSVMPointer");
    }
};

} }

/// Shared Virtual Memory wrapper class.
template <typename T>
class svm_vector : public svm_vector_terminal_expression {
    public:
        /// Allocates SVM vector on the given device.
        svm_vector(const cl::CommandQueue &q, size_t n) : n(n), q(q), p(NULL) {
            p = static_cast<T*>(clSVMAlloc(backend::get_context_id(q), CL_MEM_READ_WRITE, n * sizeof(T), 0));
        }

        ~svm_vector() {
            q.finish();
            clSVMFree(backend::get_context_id(q), p);
        }

        /// Returns size of the SVM vector.
        size_t size() const {
            return n;
        }

        svm_ptr<T> get() const {
            return svm_ptr<T>(p);
        }

        /// Returns reference to the command queue associated with the SVM vector.
        const cl::CommandQueue& queue() const {
            return q;
        }

        struct unmapper {
            cl::CommandQueue q;

            unmapper(cl::CommandQueue q) : q(q) {}

            void operator()(T* p) const {
                cl_int ret = clEnqueueSVMUnmap(q(), p, 0, NULL, NULL);
                if (ret != CL_SUCCESS) throw cl::Error(ret, "clEnqueueSVMUnmap");
            }
        };

        typedef std::unique_ptr<T[], unmapper> mapped_pointer;

        /// Returns host pointer ready to be either read or written by the host.
        /**
         * This returns a smart pointer that will be unmapped automatically
         * upon destruction */
        mapped_pointer map(cl_map_flags map_flags = CL_MAP_READ | CL_MAP_WRITE) {
            cl_int ret = clEnqueueSVMMap(q(), CL_TRUE, map_flags, p, n * sizeof(T), 0, NULL, NULL);
            if (ret != CL_SUCCESS) throw cl::Error(ret, "clEnqueueSVMUnmap");
            return mapped_pointer(p, unmapper(q));
        }

        /// Copy assignment operator.
        const svm_vector& operator=(const svm_vector & other) {
            detail::assign_expression<assign::SET>(*this, other);
            return *this;
        }

        VEXCL_ASSIGNMENTS(VEXCL_SVM_ASSIGNMENT)

    private:
        size_t n;
        cl::CommandQueue q;
        T *p;

        svm_vector(const svm_vector&) {}
};

} // namespace vex


#endif
