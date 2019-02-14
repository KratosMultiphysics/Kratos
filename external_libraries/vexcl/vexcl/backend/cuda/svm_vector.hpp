#ifndef VEXCL_BACKEND_CUDA_SVM_VECTOR_HPP
#define VEXCL_BACKEND_CUDA_SVM_VECTOR_HPP

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
 * \file   vexcl/backend/cuda/svm_vector.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Mapping between CUDA UVM and OpenCL SVM concepts.
 */

#include <memory>

#include <cuda.h>

#include <vexcl/backend/cuda/context.hpp>

namespace vex {
namespace backend {
namespace cuda {

typedef int map_flags;

static const map_flags MAP_READ  = 1;
static const map_flags MAP_WRITE = 2;

}
}

template <typename T>
class svm_vector : public svm_vector_terminal_expression {
    public:
        typedef T* mapped_pointer;

        svm_vector(const backend::command_queue &q, size_t n) : n(n), q(q), p(NULL) {
            q.context().set_current();

            CUdeviceptr dptr;
            cuda_check( cuMemAllocManaged(&dptr, n * sizeof(T), CU_MEM_ATTACH_GLOBAL) );
            p = reinterpret_cast<T*>(static_cast<size_t>(dptr));
        }

        ~svm_vector() {
            q.context().set_current();
            cuda_check( cuMemFree(static_cast<CUdeviceptr>(reinterpret_cast<size_t>(p))));
        }

        size_t size() const {
            return n;
        }

        T* get() const {
            return p;
        }

        const backend::command_queue& queue() const {
            return q;
        }

        mapped_pointer map(backend::map_flags) {
            q.finish();
            return p;
        }

        const svm_vector& operator=(const svm_vector & other) {
            detail::assign_expression<assign::SET>(*this, other);
            return *this;
        }

        VEXCL_ASSIGNMENTS(VEXCL_SVM_ASSIGNMENT)
    private:
        size_t n;
        backend::command_queue q;
        T *p;

        svm_vector(const svm_vector&) {}
};

} // namespace vex


#endif
