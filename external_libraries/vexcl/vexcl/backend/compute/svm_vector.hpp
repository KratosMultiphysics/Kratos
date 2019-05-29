#ifndef VEXCL_BACKEND_COMPUTE_SVM_VECTOR_HPP
#define VEXCL_BACKEND_COMPUTE_SVM_VECTOR_HPP

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
 * \file   vexcl/backend/compute/svm_vector.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Shared virtual memory support for Boost.Compute backend.
 */

#include <memory>

#include <boost/compute/command_queue.hpp>
#include <boost/compute/svm.hpp>

namespace vex {

namespace backend { namespace compute {

typedef cl_map_flags map_flags;

static const map_flags MAP_READ  = CL_MAP_READ;
static const map_flags MAP_WRITE = CL_MAP_WRITE;

} }

template <typename T>
class svm_vector : public svm_vector_terminal_expression {
    public:
        svm_vector(const boost::compute::command_queue &q, size_t n)
            : n(n), q(q), p(boost::compute::svm_alloc<T>(q.get_context(), n))
        {}

        ~svm_vector() {
            q.finish();
            boost::compute::svm_free(q.get_context(), p);
        }

        size_t size() const {
            return n;
        }

        boost::compute::svm_ptr<T> get() const {
            return p;
        }

        const boost::compute::command_queue& queue() const {
            return q;
        }

        struct unmapper {
            boost::compute::command_queue &q;

            unmapper(boost::compute::command_queue &q) : q(q) {}

            void operator()(T* p) const {
                q.enqueue_svm_unmap(p);
            }
        };

        typedef std::unique_ptr<T[], unmapper> mapped_pointer;

        mapped_pointer map(cl_map_flags map_flags = CL_MAP_READ | CL_MAP_WRITE) {
            q.enqueue_svm_map(p.get(), n * sizeof(T), map_flags);
            return mapped_pointer(static_cast<T*>(p.get()), unmapper(q));
        }

        const svm_vector& operator=(const svm_vector & other) {
            detail::assign_expression<assign::SET>(*this, other);
            return *this;
        }

        VEXCL_ASSIGNMENTS(VEXCL_SVM_ASSIGNMENT)
    private:
        size_t n;
        boost::compute::command_queue q;
        boost::compute::svm_ptr<T>    p;

        svm_vector(const svm_vector&) {}
};

} // namespace vex

#endif
