#ifndef VEXCL_BACKEND_OPENCL_EVENT_HPP
#define VEXCL_BACKEND_OPENCL_EVENT_HPP

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
 * \file   vexcl/backend/opencl/event.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Bring OpenCL events into vex::backend::opencl namespace.
 */

#include <vector>

#include <vexcl/backend/opencl/defines.hpp>
#ifdef VEXCL_HAVE_OPENCL_HPP
#  include <CL/opencl.hpp>
#else
#  include <CL/cl2.hpp>
#endif

#include <vexcl/backend/opencl/filter.hpp>

namespace vex {
namespace backend {
namespace opencl {

typedef cl::Event event;
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
    return e.getInfo<CL_EVENT_CONTEXT>()();
}

#ifdef CL_VERSION_1_2
/// Enqueue marker (with wait list) into the queue
inline event enqueue_marker(command_queue &q,
        const wait_list &events = wait_list())
{
    context_id ctx = get_context_id(q);
    wait_list my_events;

    std::copy_if(events.begin(), events.end(),
            std::back_inserter(my_events),
            [=](const event &e){ return ctx == get_context_id(e); });

    event e;
    q.enqueueMarkerWithWaitList(&my_events, &e);
    return e;
}

/// Enqueue barrier (with wait list) into the queue
inline event enqueue_barrier(command_queue &q,
        const wait_list &events = wait_list())
{

    context_id ctx = get_context_id(q);
    wait_list my_events;

    std::copy_if(events.begin(), events.end(),
            std::back_inserter(my_events),
            [=](const event &e){ return ctx == get_context_id(e); });

    event e;
    q.enqueueBarrierWithWaitList(&my_events, &e);
    return e;
}
#endif

/// Wait for events in the list
inline void wait_for_events(const wait_list &events) {
    std::map<context_id, wait_list> by_context;

    for(auto e = events.begin(); e != events.end(); ++e) {
        context_id ctx = get_context_id(*e);
        by_context[ctx].push_back(*e);
    }

    for(auto c = by_context.begin(); c != by_context.end(); ++c)
        cl::Event::waitForEvents(c->second);
}

} // namespace opencl
} // namespace backend
} // namespace vex

#endif
