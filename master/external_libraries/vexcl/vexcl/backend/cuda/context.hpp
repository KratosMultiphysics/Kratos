#ifndef VEXCL_BACKEND_CUDA_CONTEXT_HPP
#define VEXCL_BACKEND_CUDA_CONTEXT_HPP

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
 * \file   vexcl/backend/cuda/context.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  CUDA device enumeration and context initialization.
 */

#include <vector>
#include <tuple>
#include <iostream>
#include <memory>

#include <cuda.h>

namespace vex {
namespace backend {

/// The CUDA backend.
namespace cuda {

inline CUresult do_init() {
    static const CUresult rc = cuInit(0);
    return rc;
}

namespace detail {

template <class Handle>
struct deleter_impl;

template <>
struct deleter_impl<CUcontext> {
    static void dispose(CUcontext context) {
        cuda_check( cuCtxSetCurrent(context) );
        cuda_check( cuCtxSynchronize()       );
        cuda_check( cuCtxDestroy(context)    );
    }
};

template <>
struct deleter_impl<CUmodule> {
    static void dispose(CUmodule module) {
        cuda_check( cuModuleUnload(module) );
    }
};

template <>
struct deleter_impl<CUstream> {
    static void dispose(CUstream stream) {
        cuda_check( cuStreamDestroy(stream) );
    }
};

// Knows how to dispose of various CUDA handles.
struct deleter {
    deleter(CUcontext ctx) : ctx(ctx) {}

    template <class Handle>
    void operator()(Handle handle) const {
        if (ctx) cuda_check( cuCtxSetCurrent(ctx) );
        deleter_impl<Handle>::dispose(handle);
    }

    CUcontext ctx;
};

}

/// Wrapper around CUdevice.
class device {
    public:
        /// Empty constructor.
        device() {}
        
        /// Constructor.
        device(CUdevice d) : d(d) {}

        /// Returns raw CUdevice handle.
        CUdevice raw() const { return d; }

        /// Returns raw CUdevice handle.
        operator CUdevice() const { return d; }

        /// Returns name of the device.
        std::string name() const {
            char name[256];
            cuda_check( cuDeviceGetName(name, 256, d) );
            return name;
        }

        /// Returns device compute capability as a tuple of major and minor version numbers.
        std::tuple<int, int> compute_capability() const {
            int major, minor;

            cuda_check( cuDeviceGetAttribute(&major, CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MAJOR, d) );
            cuda_check( cuDeviceGetAttribute(&minor, CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MINOR, d) );

            return std::make_tuple(major, minor);
        }

        /// Returns of multiprocessors on the device.
        size_t multiprocessor_count() const {
            int n;
            cuda_check( cuDeviceGetAttribute(&n, CU_DEVICE_ATTRIBUTE_MULTIPROCESSOR_COUNT, d) );
            return n;
        }

        /// Returns maximum number of threads per block.
        size_t max_threads_per_block() const {
            int n;
            cuda_check( cuDeviceGetAttribute(&n, CU_DEVICE_ATTRIBUTE_MAX_THREADS_PER_BLOCK, d) );
            return n;
        }

        /// Returns maximum amount of shared memory available to a thread block in bytes.
        size_t max_shared_memory_per_block() const {
            int n;
            cuda_check( cuDeviceGetAttribute(&n, CU_DEVICE_ATTRIBUTE_MAX_SHARED_MEMORY_PER_BLOCK, d) );
            return n;
        }

        size_t warp_size() const {
            int n;
            cuda_check( cuDeviceGetAttribute(&n, CU_DEVICE_ATTRIBUTE_WARP_SIZE, d) );
            return n;
        }
    private:
        CUdevice d;
};

/// Wrapper around CUcontext.
class context {
    public:
        /// Empty constructor.
        context() {}

        /// Creates a CUDA context.
        context(device dev, unsigned flags = 0)
            : c( create(dev, flags), detail::deleter(0) )
        {
            cuda_check( do_init() );
        }

        /// Returns raw CUcontext handle.
        CUcontext raw() const {
            return c.get();
        }

        /// Returns raw CUcontext handle.
        operator CUcontext() const {
            return c.get();
        }

        /// Binds the context to the calling CPU thread.
        void set_current() const {
            cuda_check( cuCtxSetCurrent( c.get() ) );
        }

    private:
        std::shared_ptr<std::remove_pointer<CUcontext>::type> c;

        static CUcontext create(device dev, unsigned flags) {
            CUcontext h = 0;
            cuda_check( cuCtxCreate(&h, flags, dev.raw()) );
            return h;
        }
};

/// Command queue creation flags.
/**
 * typedef'ed to cl_command_queue_properties in the OpenCL backend and to
 * unsigned in the CUDA backend. See OpenCL/CUDA documentation for the possible
 * values.
 */
typedef unsigned command_queue_properties;

/// Command queue.
/** With the CUDA backend, this is a wrapper around CUstream. */
class command_queue {
    public:
        /// Empty constructor.
        command_queue() {}
        
        /// Create command queue for the given context and device.
        command_queue(const vex::backend::context &ctx, vex::backend::device dev, unsigned flags)
            : ctx(ctx), dev(dev), s( create(ctx, flags), detail::deleter(ctx.raw()) ), f(flags)
        { }

        /// Blocks until all previously queued commands in command_queue are issued to the associated device and have completed.
        void finish() const {
            ctx.set_current();
            cuda_check( cuStreamSynchronize( s.get() ) );
        }

        /// Returns the context associated with the command queue.
        vex::backend::context context() const {
            return ctx;
        }

        /// Returns the device associated with the command queue.
        vex::backend::device device() const {
            return dev;
        }

        /// Returns command_queue_properties specified at creation.
        unsigned flags() const {
            return f;
        }

        /// Returns raw CUstream handle for the command queue.
        CUstream raw() const {
            return s.get();
        }

        /// Returns raw CUstream handle for the command queue.
        operator CUstream() const {
            return s.get();
        }

    private:
        vex::backend::context  ctx;
        vex::backend::device   dev;
        std::shared_ptr<std::remove_pointer<CUstream>::type> s;
        unsigned f;

        static CUstream create(const vex::backend::context &ctx, unsigned flags = 0) {
            ctx.set_current();

            CUstream s;
            cuda_check( cuStreamCreate(&s, flags) );

            return s;
        }
};

/// CUmodule wrapper.
struct program {
    public:
        program() {}

        program(vex::backend::context c, CUmodule m)
            : c(c), m(m, detail::deleter(c.raw()))
        {}

        CUmodule raw() const {
            return m.get();
        }

        operator CUmodule() const {
            return m.get();
        }
    private:
        vex::backend::context  c;
        std::shared_ptr<std::remove_pointer<CUmodule>::type> m;
};

/// Binds the specified CUDA context to the calling CPU thread.
inline void select_context(const command_queue &q) {
    q.context().set_current();
}

/// Raw device handle.
typedef CUdevice  device_id;

/// Returns device associated with the given queue.
inline device get_device(const command_queue &q) {
    return q.device();
}

/// Returns id of the device associated with the given queue.
inline device_id get_device_id(const command_queue &q) {
    return q.device().raw();
}

/// Launch grid size.
struct ndrange {
    size_t x, y, z;
    ndrange(size_t x = 1, size_t y = 1, size_t z = 1)
        : x(x), y(y), z(z) {}
};

typedef CUcontext context_id;

/// Returns raw context id for the given queue.
inline context_id get_context_id(const command_queue &q) {
    return q.context().raw();
}

/// Returns context for the given queue.
inline context get_context(const command_queue &q) {
    return q.context();
}

/// Compares contexts by raw ids.
struct compare_contexts {
    bool operator()(const context &a, const context &b) const {
        return a.raw() < b.raw();
    }
};

/// Compares queues by raw ids.
struct compare_queues {
    bool operator()(const command_queue &a, const command_queue &b) const {
        return a.raw() < b.raw();
    }
};

/// Create command queue on the same context and device as the given one.
inline command_queue duplicate_queue(const command_queue &q) {
    return command_queue(q.context(), q.device(), q.flags());
}

/// Checks if the compute device is CPU.
/**
 * Always returns false with the CUDA backend.
 */
inline bool is_cpu(const command_queue&) {
    return false;
}

/// Select devices by given criteria.
/**
 * \param filter  Device filter functor. Functors may be combined with logical
 *                operators.
 * \returns list of devices satisfying the provided filter.
 */
template<class DevFilter>
std::vector<device> device_list(DevFilter&& filter) {
    cuda_check( do_init() );

    std::vector<device> device;

    int ndev;
    cuda_check( cuDeviceGetCount(&ndev) );

    for(int d = 0; d < ndev; ++d) {
        try {
            CUdevice dev;
            cuda_check( cuDeviceGet(&dev, d) );
            if (!filter(dev)) continue;
            device.push_back(dev);
        } catch(const error&) { }
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
std::pair< std::vector<context>, std::vector<command_queue> >
queue_list(DevFilter &&filter, unsigned queue_flags = 0)
{
    cuda_check( do_init() );

    std::vector<context>       ctx;
    std::vector<command_queue> queue;

    int ndev;
    cuda_check( cuDeviceGetCount(&ndev) );

    for(int d = 0; d < ndev; ++d) {
        try {
            CUdevice dev;
            cuda_check( cuDeviceGet(&dev, d) );
            if (!filter(dev)) continue;

            context       c(dev);
            command_queue q(c, dev, queue_flags);

            ctx.push_back(c);
            queue.push_back(q);
        } catch(const error&) { }
    }

    return std::make_pair(ctx, queue);
}

} // namespace cuda
} // namespace backend
} // namespace vex

namespace std {

/// Output device name to stream.
inline std::ostream& operator<<(std::ostream &os, const vex::backend::cuda::device &d)
{
    return os << d.name();
}

/// Output device name to stream.
inline std::ostream& operator<<(std::ostream &os, const vex::backend::cuda::command_queue &q)
{
    return os << q.device();
}

} // namespace std

#endif
