#ifndef VEXCL_BACKEND_CUDA_EVENT_HPP
#define VEXCL_BACKEND_CUDA_EVENT_HPP

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
 * \file   vexcl/backend/compute/event.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Bring Boost.Compute events into vex::backend::compute namespace.
 */

#include <vexcl/backend/cuda/context.hpp>

namespace vex {
namespace backend {
namespace cuda {

namespace detail {

template <>
struct deleter_impl<CUevent> {
    static void dispose(CUevent e) {
        cuda_check( cuEventDestroy(e) );
    }
};

} // namespace detail

class event {
    public:
        event(const command_queue &q)
            : q(q), e( create(q), detail::deleter(q.context().raw()) ) { }

        CUevent raw() const { return e.get(); }

        operator CUevent() const { return e.get(); }

        void wait() const {
            cuda_check( cuStreamWaitEvent(q.raw(), e.get(), 0) );
        }

        vex::backend::context context() const {
            return q.context();
        }
    private:
        command_queue q;
        std::shared_ptr<std::remove_pointer<CUevent>::type> e;

        static CUevent create(const command_queue &q) {
            CUevent e;
            q.context().set_current();

            cuda_check( cuEventCreate(&e, CU_EVENT_DEFAULT) );
            cuda_check( cuEventRecord(e, q.raw()) );

            return e;
        }
};

typedef std::vector<event> wait_list;

/// Append event to wait list
inline void wait_list_append(wait_list &dst, const event &e) {
    dst.push_back(e);
}

/// Append wait list to wait list
inline void wait_list_append(wait_list &dst, const wait_list &src) {
    dst.insert(dst.begin(), src.begin(), src.end());
}

/// Get id of the context the event was submitted into
inline context_id get_context_id(const event &e) {
    return e.context().raw();
}

/// Enqueue marker (with wait list) into the queue
inline event enqueue_marker(const command_queue &q, const wait_list &events = wait_list()) {
    context_id ctx = q.context().raw();

    q.context().set_current();

    for(auto e = events.begin(); e != events.end(); ++e)
        if (get_context_id(*e) == ctx)
            cuda_check( cuEventSynchronize(e->raw()) );

    return event(q);
}

/// Enqueue barrier (with wait list) into the queue
inline event enqueue_barrier(command_queue &q,
        const wait_list &events = wait_list())
{
    return enqueue_marker(q, events);
}

/// Wait for events in the list
inline void wait_for_events(const wait_list &events) {
    for(auto e = events.begin(); e != events.end(); ++e)
        e->wait();
}

} // namespace cuda
} // namespace backend
} // namespace vex


#endif
