#ifndef VEXCL_BACKEND_COMPUTE_EVENT_HPP
#define VEXCL_BACKEND_COMPUTE_EVENT_HPP

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

#include <boost/compute/event.hpp>
#include <boost/compute/wait_list.hpp>
#include <vexcl/backend/compute/context.hpp>

namespace vex {
namespace backend {
namespace compute {

using boost::compute::event;
using boost::compute::wait_list;

/// Append event to wait list
inline void wait_list_append(wait_list &dst, const event &e) {
    dst.insert(e);
}

/// Append wait list to wait list
inline void wait_list_append(wait_list &dst, const wait_list &src) {
    for(size_t i = 0; i < src.size(); ++i)
        dst.insert(src[i]);
}

/// Get id of the context the event was submitted into
inline context_id get_context_id(const event &e) {
    return e.get_info<CL_EVENT_CONTEXT>();
}

#ifdef CL_VERSION_1_2
/// Enqueue marker (with wait list) into the queue
inline event enqueue_marker(command_queue &q,
        const wait_list &events = wait_list())
{
    context_id ctx = get_context_id(q);
    wait_list my_events;

    for(auto e = events.begin(); e != events.end(); ++e)
        if (ctx == get_context_id(*e)) my_events.insert(*e);

    return q.enqueue_marker(my_events);
}

/// Enqueue barrier (with wait list) into the queue
inline event enqueue_barrier(command_queue &q,
        const wait_list &events = wait_list())
{
    context_id ctx = get_context_id(q);
    wait_list my_events;

    for(auto e = events.begin(); e != events.end(); ++e)
        if (ctx == get_context_id(*e)) my_events.insert(*e);

    return q.enqueue_barrier(my_events);
}
#endif

/// Wait for events in the list
inline void wait_for_events(const wait_list &events) {
    std::map<context_id, wait_list> by_context;

    for(auto e = events.begin(); e != events.end(); ++e) {
        context_id ctx = get_context_id(*e);
        by_context[ctx].insert(*e);
    }

    for(auto c = by_context.begin(); c != by_context.end(); ++c)
        c->second.wait();
}

} // namespace compute
} // namespace backend
} // namespace vex

#endif
