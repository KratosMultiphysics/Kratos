#ifndef VEXCL_BACKEND_COMPUTE_KERNEL_HPP
#define VEXCL_BACKEND_COMPUTE_KERNEL_HPP

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
 * \file   vexcl/backend/compute/kernel.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  An abstraction over Boost.Compute kernel.
 */

#include <functional>

#include <boost/compute/types.hpp>
#include <boost/compute/core.hpp>
#include <boost/compute/memory/local_buffer.hpp>
#include <vexcl/backend/compute/compiler.hpp>

namespace vex {
namespace backend {
namespace compute {

/// An abstraction over OpenCL compute kernel.
class kernel {
    public:
        kernel() : argpos(0), w_size(0), g_size(0) {}

        /// Constructor. Creates a backend::kernel instance from source.
        kernel(const boost::compute::command_queue &queue,
               const std::string &src,
               const std::string &name,
               size_t smem_per_thread = 0,
               const std::string &options = ""
               )
            : argpos(0), K(build_sources(queue, src, options), name)
        {
            config(queue,
                    [smem_per_thread](size_t wgs){ return wgs * smem_per_thread; });
        }

        /// Constructor. Creates a backend::kernel instance from source.
        kernel(const boost::compute::command_queue &queue,
               const std::string &src, const std::string &name,
               std::function<size_t(size_t)> smem,
               const std::string &options = ""
               )
            : argpos(0), K(build_sources(queue, src, options), name)
        {
            config(queue, smem);
        }

        /// Constructor. Extracts a backend::kernel instance from backend::program.
        kernel(const boost::compute::command_queue &queue,
               const boost::compute::program &program,
               const std::string &name,
               size_t smem_per_thread = 0
               )
            : argpos(0), K(program, name)
        {
            config(queue,
                    [smem_per_thread](size_t wgs){ return wgs * smem_per_thread; });
        }

        /// Constructor. Extracts a backend::kernel instance from backend::program.
        kernel(const boost::compute::command_queue &queue,
                const boost::compute::program &program,
               const std::string &name,
               std::function<size_t(size_t)> smem
               )
            : argpos(0), K(program, name)
        {
            config(queue, smem);
        }

        /// Adds an argument to the kernel.
        template <typename T>
        void push_arg(device_vector<T> arg) {
            K.set_arg(argpos++, arg.raw());
        }

        /// Adds an argument to the kernel.
        template <class Arg>
        typename boost::enable_if<boost::compute::is_fundamental<Arg>, void>::type
        push_arg(Arg &&arg) {
            K.set_arg(argpos++, arg);
        }

        /// Adds an argument to the kernel.
        template <class Arg>
        typename boost::disable_if<boost::compute::is_fundamental<Arg>, void>::type
        push_arg(Arg &&arg) {
            K.set_arg(argpos++, sizeof(arg), &arg);
        }

        /// Adds local memory to the kernel.
        void set_smem(size_t smem_per_thread) {
            K.set_arg(
                    argpos++,
                    boost::compute::local_buffer<char>(smem_per_thread * workgroup_size())
                    );
        }

        /// Adds local memory to the kernel.
        template <class F>
        void set_smem(F &&f) {
            K.set_arg(
                    argpos++,
                    boost::compute::local_buffer<char>( f(workgroup_size()) )
                    );
        }

        /// Enqueue the kernel to the specified command queue.
        void operator()(boost::compute::command_queue q) {
            q.enqueue_nd_range_kernel(K, 3, NULL, g_size.dim, w_size.dim);
            argpos = 0;
        }

#ifndef BOOST_NO_VARIADIC_TEMPLATES
        /// Enqueue the kernel to the specified command queue with the given arguments
        template <class Arg1, class... OtherArgs>
        void operator()(const boost::compute::command_queue &q, Arg1 &&arg1, OtherArgs&&... other_args) {
            push_arg(std::forward<Arg1>(arg1));

            (*this)(q, std::forward<OtherArgs>(other_args)...);
        }
#endif

        /// Workgroup size.
        size_t workgroup_size() const {
            return w_size.x * w_size.y * w_size.z;
        }

        /// Standard number of workgroups to launch on a device.
        static inline size_t num_workgroups(const boost::compute::command_queue &q) {
            // This is a simple heuristic-based estimate. More advanced technique may
            // be employed later.
            return 8 * q.get_device().compute_units();
        }

        /// The maximum number of threads per block, beyond which a launch of the kernel would fail.
        size_t max_threads_per_block(const boost::compute::command_queue &q) const {
            return K.get_work_group_info<size_t>(q.get_device(), CL_KERNEL_WORK_GROUP_SIZE);
        }

        /// The size in bytes of shared memory per block available for this kernel.
        size_t max_shared_memory_per_block(const boost::compute::command_queue &q) const {
            boost::compute::device d = q.get_device();

            return d.local_memory_size() - K.get_work_group_info<cl_ulong>(d, CL_KERNEL_LOCAL_MEM_SIZE);
        }

        /// Select best launch configuration for the given shared memory requirements.
        kernel& config(const boost::compute::command_queue &queue, std::function<size_t(size_t)> smem) {
            boost::compute::device dev = queue.get_device();

            size_t ws;

            if ( is_cpu(queue) ) {
                ws = 1;
            } else {
                // Select workgroup size that would fit into the device.
                ws = dev.get_info<std::vector<size_t>>(CL_DEVICE_MAX_WORK_ITEM_SIZES)[0] / 2;

                size_t max_ws   = max_threads_per_block(queue);
                size_t max_smem = max_shared_memory_per_block(queue);

                // Reduce workgroup size until it satisfies resource requirements:
                while( (ws > max_ws) || (smem(ws) > max_smem) )
                    ws /= 2;
            }

            return config(num_workgroups(queue), ws);
        }

        /// Set launch configuration.
        kernel& config(ndrange blocks, ndrange threads, size_t shared_memory = 0) {
            const size_t *b = blocks.dim;
            const size_t *t = threads.dim;

            g_size = ndrange(b[0] * t[0], b[1] * t[1], b[2] * t[2]);
            w_size = threads;

            if (shared_memory)
                K.set_arg(argpos++, boost::compute::local_buffer<char>(shared_memory));

            return *this;
        }

        /// Set launch configuration.
        kernel& config(size_t blocks, size_t threads, size_t shared_memory = 0) {
            return config(ndrange(blocks), ndrange(threads), shared_memory);
        }

        size_t preferred_work_group_size_multiple(const boost::compute::command_queue &q) const {
            return K.get_work_group_info<size_t>(q.get_device(), CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE);
        }

        /// Get reference to the underlying object.
        const boost::compute::kernel& get() const { return K; }
        boost::compute::kernel& get() { return K; }

        /// Reset argument counter.
        void reset() {
            argpos = 0;
        }
    private:
        unsigned argpos;

        boost::compute::kernel K;

        backend::ndrange w_size;
        backend::ndrange g_size;
};

} // namespace compute
} // namespace backend
} // namespace vex

#endif
