#ifndef VEXCL_BACKEND_JIT_KERNEL_HPP
#define VEXCL_BACKEND_JIT_KERNEL_HPP

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
 * \file   vexcl/backend/jit/kernel.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Compute kernel implementation for the JIT backend.
 */

#include <string>
#include <boost/dll/import.hpp>

#ifdef _OPENMP
#  include <omp.h>
#endif

#include <vexcl/util.hpp>
#include <vexcl/backend/jit/compiler.hpp>

namespace vex {
namespace backend {
namespace jit {

namespace detail {

struct kernel_api {
    virtual void execute(
            const ndrange *dim, size_t smem_size, char *prm
            ) const = 0;
};

} // namespace detail

class kernel {
    public:
        kernel() : smem_size(0) {}

        kernel(
                const command_queue &q,
                const std::string &src, const std::string &name,
                size_t smem_per_thread = 0,
                const std::string &options = ""
              )
#if BOOST_VERSION >= 107600
            : K(boost::dll::import_symbol<detail::kernel_api>(build_sources(q, src, options), name)),
#else
            : K(boost::dll::import<detail::kernel_api>(build_sources(q, src, options), name)),
#endif
              grid(num_workgroups(q)), smem_size(smem_per_thread)
        {
            stack.reserve(256);
        }

        kernel(const command_queue &q,
               const std::string &src, const std::string &name,
               std::function<size_t(size_t)> smem,
               const std::string &options = ""
               )
#if BOOST_VERSION >= 107600
            : K(boost::dll::import_symbol<detail::kernel_api>(build_sources(q, src, options), name)),
#else
            : K(boost::dll::import<detail::kernel_api>(build_sources(q, src, options), name)),
#endif
              grid(num_workgroups(q)), smem_size(smem(1))
        {
            stack.reserve(256);
        }

        kernel(const command_queue &q,
               const program &P,
               const std::string &name,
               size_t smem_per_thread = 0
               )
#if BOOST_VERSION >= 107600
            : K(boost::dll::import_symbol<detail::kernel_api>(P, name)),
#else
            : K(boost::dll::import<detail::kernel_api>(P, name)),
#endif
              grid(num_workgroups(q)), smem_size(smem_per_thread)
        {
            stack.reserve(256);
        }

        /// Constructor. Extracts a backend::kernel instance from backend::program.
        kernel(const command_queue &q, const program &P,
               const std::string &name,
               std::function<size_t(size_t)> smem
               )
#if BOOST_VERSION >= 107600
            : K(boost::dll::import_symbol<detail::kernel_api>(P, name)),
#else
            : K(boost::dll::import<detail::kernel_api>(P, name)),
#endif
              grid(num_workgroups(q)), smem_size(smem(1))
        {
            stack.reserve(256);
        }

        template <class Arg>
        void push_arg(const Arg &arg) {
            char *c = (char*)&arg;
            stack.insert(stack.end(), c, c + sizeof(arg));
        }

        template <typename T>
        void push_arg(const device_vector<T> &arg) {
            push_arg(arg.raw());
        }

        void set_smem(size_t smem_per_thread) {
            smem_size = smem_per_thread;
        }

        template <class F>
        void set_smem(F &&f) {
            smem_size = f(1);
        }

        void operator()(const command_queue&) {
            // All parameters have been pushed; time to call the kernel:
            K->execute(&grid, smem_size, stack.data());

            // Reset parameter stack:
            stack.clear();
        }

#ifndef BOOST_NO_VARIADIC_TEMPLATES
        template <class Head, class... Tail>
        void operator()(const command_queue &q, const Head &head, const Tail&... tail) {
            push_arg(head);
            (*this)(q, tail...);
        }
#endif
        size_t workgroup_size() const {
            return 1UL;
        }

        static inline size_t num_workgroups(const command_queue&) {
#ifdef _OPENMP
            return omp_get_num_procs() * 8;
#else
            return 1UL;
#endif
        }

        size_t max_threads_per_block(const command_queue&) const {
            return 1UL;
        }

        size_t max_shared_memory_per_block(const command_queue&) const {
            return 32768UL;
        }

        size_t preferred_work_group_size_multiple(const backend::command_queue &q) const {
            return 1;
        }

        kernel& config(const command_queue &q, std::function<size_t(size_t)> smem) {
            return config(num_workgroups(q), 1);
        }

        kernel& config(ndrange blocks, ndrange threads, size_t shared_memory = 0) {
            precondition(threads == ndrange(), "Maximum workgroup size for the JIT backend is 1");
            grid = blocks;
            if (shared_memory) smem_size = shared_memory;
            return *this;
        }

        kernel& config(size_t blocks, size_t threads) {
            precondition(threads == 1, "Maximum workgroup size for the JIT backend is 1");
            return config(ndrange(blocks), ndrange(threads));
        }

        void reset() {
            stack.clear();
        }
    private:
        boost::shared_ptr<detail::kernel_api> K;
        ndrange grid;
        std::vector<char> stack;
        size_t smem_size;
};

} // namespace jit
} // namespace backend
} // namespace vex

#endif
