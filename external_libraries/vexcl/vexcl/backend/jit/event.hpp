#ifndef VEXCL_BACKEND_JIT_EVENT_HPP
#define VEXCL_BACKEND_JIT_EVENT_HPP

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
 * \file   vexcl/backend/jit/event.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Define dummy events for the JIT backend.
 */

#include <vexcl/backend/jit/context.hpp>

namespace vex {
namespace backend {
namespace jit {

struct event {
    event() {}
    event(const command_queue&) {}
    void wait() const {}
    vex::backend::context context() const {
        return vex::backend::context();
    }
};

struct wait_list {
    template <class... T>
    wait_list(T&&... t) {}
};

/// Append event to wait list
inline void wait_list_append(wait_list&, const event&) { }

/// Append wait list to wait list
inline void wait_list_append(wait_list&, const wait_list&) {}

/// Get id of the context the event was submitted into
inline context_id get_context_id(const event&) {
    return 0;
}

/// Enqueue marker (with wait list) into the queue
inline event enqueue_marker(const command_queue&, const wait_list& = wait_list()) {
    return event();
}

/// Enqueue barrier (with wait list) into the queue
inline event enqueue_barrier(command_queue&, const wait_list& = wait_list()) {
    return event();
}

/// Wait for events in the list
inline void wait_for_events(const wait_list&) {}

} // namespace jit
} // namespace backend
} // namespace vex

#endif
