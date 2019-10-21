#ifndef VEXCL_GATHER_HPP
#define VEXCL_GATHER_HPP

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
 * \file   vexcl/gather.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Gather scattered points from OpenCL device vector.
 */

#include <vector>
#include <numeric>
#include <cassert>

#include <vexcl/operations.hpp>
#include <vexcl/vector.hpp>
#include <vexcl/vector_view.hpp>

namespace vex {

namespace detail {

template <class T>
class index_partition {
    public:
        index_partition(
                const std::vector<backend::command_queue> &q,
                size_t size, std::vector<size_t> indices
                )
            : queue(q), ptr(q.size() + 1, 0),
              idx(q.size()), val(q.size())
        {
            assert(std::is_sorted(indices.begin(), indices.end()));

            std::vector<size_t> part = partition(size, queue);
            column_owner owner(part);

            for(auto i = indices.begin(); i != indices.end(); ++i) {
                size_t d = owner(*i);
                *i -= part[d];
                ++ptr[d + 1];
            }

            std::partial_sum(ptr.begin(), ptr.end(), ptr.begin());

            for(unsigned d = 0; d < queue.size(); d++) {
                if (size_t n = ptr[d + 1] - ptr[d]) {
                    val[d] = backend::device_vector<T>(queue[d], n, static_cast<const T*>(0));
                    idx[d] = backend::device_vector<size_t>(
                            queue[d], n, &indices[ptr[d]], backend::MEM_READ_ONLY);
                }
            }

            for(unsigned d = 0; d < queue.size(); d++)
                if (ptr[d + 1] - ptr[d]) queue[d].finish();
        }

    protected:
        std::vector<backend::command_queue>           queue;
        std::vector< size_t >                         ptr;
        std::vector< backend::device_vector<size_t> > idx;
        std::vector< backend::device_vector<T> >      val;
};

} // namespace detail

/// Gathers vector elements at specified indices.
template <typename T>
class gather : protected detail::index_partition<T> {
    public:
        /// Constructor.
        /**
         * \param queue   VexCL context.
         * \param indices Indices of elements to be gathered.
         */
        gather(
                const std::vector<backend::command_queue> &q,
                size_t size, std::vector<size_t> indices
              ) : Base(q, size, indices)
        {}

        /// Gather elements of device vector into host vector.
        template <class HostVector>
        void operator()(const vex::vector<T> &src, HostVector &dst) {
            using namespace detail;

            static kernel_cache cache;

            for(unsigned d = 0; d < Base::queue.size(); d++) {
                if (size_t n = Base::ptr[d + 1] - Base::ptr[d]) {
                    vector<T>      v(Base::queue[d], Base::val[d]);
                    vector<T>      s(Base::queue[d], src(d));
                    vector<size_t> i(Base::queue[d], Base::idx[d]);

                    v = permutation(i)(s);

                    Base::val[d].read(Base::queue[d], 0, n, &dst[Base::ptr[d]]);
                }
            }

            for(unsigned d = 0; d < Base::queue.size(); d++)
                if (Base::ptr[d + 1] - Base::ptr[d]) Base::queue[d].finish();
        }
    private:
        typedef detail::index_partition<T> Base;
};

/// Scatters vector elements to specified indices.
template <typename T>
class scatter : protected detail::index_partition<T> {
    public:
        /// Constructor.
        /**
         * \param queue   VexCL context.
         * \param indices Indices of elements to be gathered.
         */
        scatter(
                const std::vector<backend::command_queue> &q,
                size_t size, std::vector<size_t> indices
              ) : Base(q, size, indices)
        {}

        /// Scatter elements of host vector to device vector.
        template <class HostVector>
        void operator()(const HostVector &src, vex::vector<T> &dst) {
            using namespace detail;

            static kernel_cache cache;

            for(unsigned d = 0; d < Base::queue.size(); d++) {
                if (size_t n = Base::ptr[d + 1] - Base::ptr[d]) {
                    Base::val[d].write(Base::queue[d], 0, n, &src[Base::ptr[d]]);

                    vector<T>      v(Base::queue[d], Base::val[d]);
                    vector<T>      s(Base::queue[d], dst(d));
                    vector<size_t> i(Base::queue[d], Base::idx[d]);

                    permutation(i)(s) = v;
                }
            }

            for(unsigned d = 0; d < Base::queue.size(); d++)
                if (Base::ptr[d + 1] - Base::ptr[d]) Base::queue[d].finish();
        }
    private:
        typedef detail::index_partition<T> Base;
};

} // namespace vex

#endif
