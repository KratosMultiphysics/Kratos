#ifndef VEXCL_BACKEND_JIT_CONTEXT_HPP
#define VEXCL_BACKEND_JIT_CONTEXT_HPP

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
 * \file   vexcl/backend/jit/context.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  JIT context definitions.
 */

#include <string>
#include <stdexcept>

#include <boost/dll/shared_library.hpp>

namespace vex {
namespace backend {

/// JIT backend with OpenMP support
namespace jit {

struct device {
    device(unsigned id = 0) {}

    std::string name() const {
        return "CPU";
    }

    // Took the constants from Intel OpenCL:
    size_t max_shared_memory_per_block() const { return 32768UL; }
    size_t max_threads_per_block()       const { return 1L; }
};

struct context {};

typedef unsigned command_queue_properties;

struct command_queue {
    void finish() const {}

    vex::backend::context context() const {
        return vex::backend::context();
    }

    vex::backend::device device() const {
        return vex::backend::device();
    }
};

typedef boost::dll::shared_library program;

typedef unsigned device_id;

inline device get_device(const command_queue&) {
    return device();
}

inline device_id get_device_id(const command_queue&) {
    return 0;
}

typedef unsigned context_id;

inline context_id get_context_id(const command_queue&) {
    return 0;
}

inline context get_context(const command_queue&) {
    return context();
}

inline void select_context(const command_queue&) {}

inline command_queue duplicate_queue(const command_queue &q) {
    return command_queue();
}

inline bool is_cpu(const command_queue &q) {
    return true;
}

struct compare_contexts {
    bool operator()(const context&, const context&) const {
        return false;
    }
};

struct compare_queues {
    bool operator()(const command_queue&, const command_queue&) const {
        return false;
    }
};

template<class DevFilter>
std::vector<device> device_list(DevFilter&& filter) {
    std::vector<device> dev;

    device d;
    if (filter(d)) dev.push_back(d);

    return dev;
}

template<class DevFilter>
std::pair< std::vector<context>, std::vector<command_queue> >
queue_list(DevFilter &&filter, unsigned queue_flags = 0)
{
    std::vector<context>       ctx;
    std::vector<command_queue> queue;

    device d;

    if (filter(d)) {
        ctx.push_back(context());
        queue.push_back(command_queue());
    }

    return std::make_pair(ctx, queue);
}

typedef std::exception error;

struct ndrange {
    size_t x, y, z;

    ndrange(size_t x = 1, size_t y = 1, size_t z = 1)
        : x(x), y(y), z(z) {}

    bool operator==(const ndrange &o) const {
        return x == o.x && y == o.y && z == o.z;
    }
};

} // namespace jit
} // namespace backend
} // namespace vex

namespace std {

/// Output device name to stream.
inline std::ostream& operator<<(std::ostream &os, const vex::backend::jit::device &d)
{
    return os << d.name();
}

/// Output device name to stream.
inline std::ostream& operator<<(std::ostream &os, const vex::backend::jit::command_queue &q)
{
    return os << q.device();
}

inline std::ostream& operator<<(std::ostream &os, const vex::backend::jit::error &e) {
    return os << e.what();
}

} // namespace std

#endif
