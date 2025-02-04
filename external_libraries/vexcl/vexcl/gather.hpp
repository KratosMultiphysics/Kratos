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

class index_partition {
    public:
        index_partition(
                const std::vector<backend::command_queue> &q,
                size_t size, const std::vector<size_t> &indices
                )
            : queue(q), ptr(q.size() + 1, 0)
        {
            if (queue.size() > 1 && !std::is_sorted(indices.begin(), indices.end())) {
                idx.resize(indices.size());
                ord.resize(indices.size());

                std::iota(ord.begin(), ord.end(), 0ul);
                std::sort(ord.begin(), ord.end(), [&indices](size_t i, size_t j){ return indices[i] < indices[j]; });

                std::vector<size_t> I(indices.size());
                for(size_t i = 0; i < indices.size(); ++i)
                    idx[i] = indices[ord[i]];
            } else {
                idx = indices;
            }

            std::vector<size_t> part = partition(size, queue);
            column_owner owner(part);

            for(auto i = idx.begin(); i != idx.end(); ++i) {
                size_t d = owner(*i);
                *i -= part[d];
                ++ptr[d + 1];
            }

            std::partial_sum(ptr.begin(), ptr.end(), ptr.begin());
        }

    protected:
        std::vector<backend::command_queue> queue;
        std::vector< size_t > ptr;
        std::vector< size_t > idx;
        std::vector< size_t > ord;
};

} // namespace detail

/// Gathers vector elements at specified indices.
class gather : protected detail::index_partition {
    public:
        /// Constructor.
        /**
         * \param queue   VexCL context.
         * \param indices Indices of elements to be gathered.
         */
        gather(
                const std::vector<backend::command_queue> &q,
                size_t size, const std::vector<size_t> &indices
              ) : Base(q, size, indices)
        {}

        /// Gather elements of device vector into host vector.
        template <class T, class HostVector>
        void operator()(const vex::vector<T> &src, HostVector &dst) {
            using namespace detail;

            for(unsigned d = 0; d < Base::queue.size(); d++) {
                if (size_t n = Base::ptr[d + 1] - Base::ptr[d]) {
                    auto s = src.map(d);

                    if (Base::ord.empty()) {
                        for(size_t i = Base::ptr[d]; i < Base::ptr[d+1]; ++i)
                            dst[i] = s[Base::idx[i]];
                    } else {
                        for(size_t i = Base::ptr[d]; i < Base::ptr[d+1]; ++i)
                            dst[Base::ord[i]] = s[Base::idx[i]];
                    }
                }
            }
        }
    private:
        typedef detail::index_partition Base;
};

/// Scatters vector elements to specified indices.
class scatter : protected detail::index_partition {
    public:
        /// Constructor.
        /**
         * \param queue   VexCL context.
         * \param indices Indices of elements to be gathered.
         */
        scatter(
                const std::vector<backend::command_queue> &q,
                size_t size, const std::vector<size_t> &indices
              ) : Base(q, size, indices)
        {}

        /// Scatter elements of host vector to device vector.
        template <class HostVector, class T>
        void operator()(const HostVector &src, vex::vector<T> &dst) {
            using namespace detail;

            for(unsigned d = 0; d < Base::queue.size(); d++) {
                if (size_t n = Base::ptr[d + 1] - Base::ptr[d]) {
                    auto v = dst.map(d);
                    if (Base::ord.empty()) {
                        for(size_t i = Base::ptr[d]; i < Base::ptr[d+1]; ++i)
                            v[Base::idx[i]] = src[i];
                    } else {
                        for(size_t i = Base::ptr[d]; i < Base::ptr[d+1]; ++i)
                            v[Base::idx[i]] = src[Base::ord[i]];
                    }
                }
            }
        }
    private:
        typedef detail::index_partition Base;
};

} // namespace vex

#endif
