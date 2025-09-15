#ifndef VEXCL_BACKEND_OPENCL_CONTEXT_HPP
#define VEXCL_BACKEND_OPENCL_CONTEXT_HPP

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
 * \file   vexcl/backend/opencl/context.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  OpenCL device enumeration and context initialization.
 */

#include <vector>
#include <iostream>

#include <vexcl/backend/opencl/defines.hpp>
#ifdef VEXCL_HAVE_OPENCL_HPP
#  include <CL/opencl.hpp>
#else
#  include <CL/cl2.hpp>
#endif

namespace vex {
namespace backend {

/// The OpenCL backend.
namespace opencl {

typedef cl::Context                 context;
typedef cl::Device                  device;
typedef cl::Program                 program;
typedef cl::CommandQueue            command_queue;
typedef cl_command_queue_properties command_queue_properties;
typedef cl_device_id                device_id;
typedef cl::NDRange                 ndrange;

/// Binds the specified context to the calling CPU thread.
/**
 * With the OpenCL backend this is an empty stub provided for compatibility
 * with the CUDA backend.
 */
inline void select_context(const command_queue&) {
}

/// Returns device associated with the given queue.
inline device get_device(const command_queue &q) {
    return q.getInfo<CL_QUEUE_DEVICE>();
}

/// Returns id of the device associated with the given queue.
inline device_id get_device_id(const command_queue &q) {
    return q.getInfo<CL_QUEUE_DEVICE>()();
}

typedef cl_context       context_id;
/// Returns raw context id for the given queue.
inline context_id get_context_id(const command_queue &q) {
    return q.getInfo<CL_QUEUE_CONTEXT>()();
}

/// Returns context for the given queue.
inline context get_context(const command_queue &q) {
    return q.getInfo<CL_QUEUE_CONTEXT>();
}

/// Compares contexts by raw ids.
struct compare_contexts {
    bool operator()(const context &a, const context &b) const {
        return a() < b();
    }
};

/// Compares queues by raw ids.
struct compare_queues {
    bool operator()(const command_queue &a, const command_queue &b) const {
        return a() < b();
    }
};

/// Create command queue on the same context and device as the given one.
inline command_queue duplicate_queue(const command_queue &q) {
    return command_queue(
            q.getInfo<CL_QUEUE_CONTEXT>(), q.getInfo<CL_QUEUE_DEVICE>());
}

/// Checks if the compute device is CPU.
inline bool is_cpu(const command_queue &q) {
    cl::Device d = q.getInfo<CL_QUEUE_DEVICE>();
#ifdef _MSC_VER
#  pragma warning(push)
#  pragma warning(disable: 4800)
#endif
    return d.getInfo<CL_DEVICE_TYPE>() & CL_DEVICE_TYPE_CPU;
#ifdef _MSC_VER
#  pragma warning(pop)
#endif
}

/** Returns vector of compute devices satisfying the given criteria without
 * trying to initialize the contexts on the devices.
 */
template<class DevFilter>
std::vector<vex::backend::device> device_list(DevFilter&& filter) {
    std::vector<cl::Device> device;

    std::vector<cl::Platform> platforms;
    cl::Platform::get(&platforms);

    for(auto p = platforms.begin(); p != platforms.end(); p++) {
        std::vector<cl::Device> dev_list;

        p->getDevices(CL_DEVICE_TYPE_ALL, &dev_list);

        for(auto d = dev_list.begin(); d != dev_list.end(); d++) {
            if (!d->getInfo<CL_DEVICE_AVAILABLE>()) continue;
            if (!filter(*d)) continue;

            device.push_back(*d);
        }
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
std::pair<std::vector<cl::Context>, std::vector<command_queue>>
queue_list(DevFilter &&filter, cl_command_queue_properties properties = 0) {
    std::vector<cl::Context>      context;
    std::vector<command_queue> queue;

    std::vector<cl::Platform> platforms;
    cl::Platform::get(&platforms);

    for(auto p = platforms.begin(); p != platforms.end(); p++) {
        std::vector<cl::Device> device;
        std::vector<cl::Device> dev_list;

        p->getDevices(CL_DEVICE_TYPE_ALL, &dev_list);

        for(auto d = dev_list.begin(); d != dev_list.end(); d++) {
            if (!d->getInfo<CL_DEVICE_AVAILABLE>()) continue;
            if (!filter(*d)) continue;

            device.push_back(*d);
        }

        if (device.empty()) continue;

        for(auto d = device.begin(); d != device.end(); d++)
            try {
                context.push_back(cl::Context(std::vector<cl::Device>(1, *d)));
                queue.push_back(command_queue(context.back(), *d, properties));
            } catch(const cl::Error&) {
                // Something bad happened. Better skip this device.
            }
    }

    return std::make_pair(context, queue);
}

} // namespace opencl
} // namespace backend
} // namespace vex

namespace std {

/// Output device name to stream.
inline std::ostream& operator<<(std::ostream &os, const vex::backend::opencl::command_queue &q)
{
    cl::Device   d(q.getInfo<CL_QUEUE_DEVICE>());
    cl::Platform p(d.getInfo<CL_DEVICE_PLATFORM>());

    return os << d.getInfo<CL_DEVICE_NAME>()
              << " (" << p.getInfo<CL_PLATFORM_NAME>() << ")";
}

} // namespace std

#endif
