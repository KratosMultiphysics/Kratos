#ifndef VEXCL_BACKEND_COMPUTE_CONTEXT_HPP
#define VEXCL_BACKEND_COMPUTE_CONTEXT_HPP

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
 * \file   vexcl/backend/compute/context.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Boost.Compute device enumeration and context initialization.
 */

#include <vector>
#include <iostream>

#include <boost/compute/core.hpp>

namespace vex {
namespace backend {

/// The Boost.Compute backend.
namespace compute {

typedef boost::compute::context       context;
typedef boost::compute::device        device;
typedef boost::compute::program       program;
typedef boost::compute::command_queue command_queue;
typedef cl_command_queue_properties   command_queue_properties;
typedef cl_device_id                  device_id;
typedef cl_context                    context_id;

/// Launch grid size.
struct ndrange {
    union {
        size_t dim[3];
        struct {
            size_t x, y, z;
        };
    };

    ndrange(size_t x = 1, size_t y = 1, size_t z = 1)
        : x(x), y(y), z(z) {}
};

/// Binds the specified context to the calling CPU thread.
/**
 * With the OpenCL backend this is an empty stub provided for compatibility
 * with the CUDA backend.
 */
inline void select_context(const command_queue&) { }

/// Returns device associated with the given queue.
inline device get_device(const command_queue &q) {
    return q.get_device();
}

/// Returns id of the device associated with the given queue.
inline device_id get_device_id(const command_queue &q) {
    return q.get_device().get();
}

/// Returns raw context id for the given queue.
inline context_id get_context_id(const command_queue &q) {
    return q.get_context().get();
}

/// Returns context for the given queue.
inline context get_context(const command_queue &q) {
    return q.get_context();
}

/// Compares contexts by raw ids.
struct compare_contexts {
    bool operator()(const context &a, const context &b) const {
        return a.get() < b.get();
    }
};

/// Compares queues by raw ids.
struct compare_queues {
    bool operator()(const command_queue &a, const command_queue &b) const {
        return a.get() < b.get();
    }
};

/// Create command queue on the same context and device as the given one.
inline command_queue duplicate_queue(const command_queue &q) {
    return command_queue(q.get_context(), q.get_device(), q.get_properties());
}

/// Checks if the compute device is CPU.
inline bool is_cpu(const command_queue &q) {
    return q.get_device().get_info<cl_device_type>(CL_DEVICE_TYPE) & CL_DEVICE_TYPE_CPU;
}

/// Select devices by given criteria.
/**
 * \param filter  Device filter functor. Functors may be combined with logical
 *                operators.
 * \returns list of devices satisfying the provided filter.
 *
 * This example selects any GPU which supports double precision arithmetic:
 \code
 auto devices = device_list(
          Filter::Type(CL_DEVICE_TYPE_GPU) && Filter::DoublePrecision
          );
 \endcode
 */
template<class DevFilter>
std::vector<device> device_list(DevFilter&& filter) {
    std::vector<device> dev_list = boost::compute::system::devices();
    std::vector<device> device;

    for(auto d = dev_list.begin(); d != dev_list.end(); d++) {
        if (!d->get_info<cl_bool>(CL_DEVICE_AVAILABLE)) continue;
        if (!filter(*d)) continue;

        device.push_back(*d);
    }

    return device;
}

/// Create command queues on devices by given criteria.
/**
 * \param filter  Device filter functor. Functors may be combined with logical
 *                operators.
 * \param properties Command queue properties.
 *
 * \returns list of queues accociated with selected devices.
 * \see device_list
 */
template<class DevFilter>
std::pair<std::vector<context>, std::vector<command_queue>>
queue_list(DevFilter &&filter, cl_command_queue_properties properties = 0) {
    std::vector<context>       c;
    std::vector<command_queue> q;

    std::vector<device> dev_list = boost::compute::system::devices();

    for(auto d = dev_list.begin(); d != dev_list.end(); d++) {
        if (!d->get_info<cl_bool>(CL_DEVICE_AVAILABLE)) continue;
        if (!filter(*d)) continue;

        try {
            c.push_back(context(*d));
            q.push_back(command_queue(c.back(), *d, properties));
        } catch(...) {
            // Something bad happened. Better skip this device.
        }
    }

    return std::make_pair(c, q);
}

} // namespace compute
} // namespace backend
} // namespace vex

namespace std {

/// Output device name to stream.
inline std::ostream& operator<<(std::ostream &os, const vex::backend::compute::command_queue &q)
{
    boost::compute::device d = q.get_device();
    return os << d.name() << " (" << d.platform().name() << ")";
}

} // namespace std

#endif
