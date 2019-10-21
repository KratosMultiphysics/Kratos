#ifndef VEXCL_EXTERNAL_BOOST_COMPUTE_HPP
#define VEXCL_EXTERNAL_BOOST_COMPUTE_HPP

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
 * \file   external/boost_compute.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Provides wrappers for some of Boost.Compute (https://github.com/kylelutz/compute) algorithms.
 */

#include <algorithm>
#include <vexcl/vector.hpp>
#include <vexcl/sort.hpp>
#include <boost/compute.hpp>

#if defined(VEXCL_BACKEND_CUDA)
#  error Boost.Compute interoperation is not supported for the CUDA backend!
#endif

#if defined(VEXCL_BACKEND_JIT)
#  error Boost.Compute interoperation is not supported for the JIT backend!
#endif

#if defined(VEXCL_BACKEND_COMPUTE)
#  error The code below is not required for Boost.Compute backend!
#endif

#if !defined(VEXCL_BACKEND_OPENCL)
#  error Unsupported backend!
#endif

namespace vex {

/// Wrapping code for some of Boost.Compute algorithms.
namespace compute {

/// Inclusive scan.
template <typename T>
void inclusive_scan(const vex::vector<T> &src, vex::vector<T> &dst) {
    auto queue = src.queue_list();

    // Scan partitions separately.
    for(unsigned d = 0; d < queue.size(); ++d) {
        if (src.part_size(d)) {
            boost::compute::command_queue q( queue[d]() );

            boost::compute::buffer sbuf( src(d).raw() );
            boost::compute::buffer dbuf( dst(d).raw() );

            boost::compute::inclusive_scan(
                    boost::compute::make_buffer_iterator<T>(sbuf, 0),
                    boost::compute::make_buffer_iterator<T>(sbuf, src.part_size(d)),
                    boost::compute::make_buffer_iterator<T>(dbuf, 0),
                    q
                    );
        }
    }

    // If there are more than one partition,
    // update all of them except for the first.
    if (queue.size() > 1) {
        std::vector<T> tail(queue.size() - 1, T());

        for(unsigned d = 0; d < tail.size(); ++d) {
            if (src.part_size(d))
                tail[d] = dst[src.part_start(d + 1) - 1];
        }

        std::partial_sum(tail.begin(), tail.end(), tail.begin());

        for(unsigned d = 1; d < queue.size(); ++d) {
            if (src.part_size(d)) {
                // Wrap partition into vector for ease of use:
                vex::vector<T> part(queue[d], dst(d));
                part += tail[d - 1];
            }
        }
    }
}

/// Exclusive scan.
template <typename T>
void exclusive_scan(const vex::vector<T> &src, vex::vector<T> &dst) {
    auto queue = src.queue_list();

    std::vector<T> tail;
    /* If there is more than one partition, we need to take a copy the last
     * element in each partition (except the last) as otherwise information
     * about it is lost.
     *
     * This must be captured here rather than later, in case the input and
     * output alias.
     */
    if (queue.size() > 1) {
        tail.resize(queue.size() - 1);
        for (unsigned d = 0; d < tail.size(); ++d) {
            if (src.part_size(d))
                tail[d] = src[src.part_start(d + 1) - 1];
        }
    }

    // Scan partitions separately.
    for(unsigned d = 0; d < queue.size(); ++d) {
        if (src.part_size(d)) {
            boost::compute::command_queue q( queue[d]() );

            boost::compute::buffer sbuf( src(d).raw() );
            boost::compute::buffer dbuf( dst(d).raw() );

            boost::compute::exclusive_scan(
                    boost::compute::make_buffer_iterator<T>(sbuf, 0),
                    boost::compute::make_buffer_iterator<T>(sbuf, src.part_size(d)),
                    boost::compute::make_buffer_iterator<T>(dbuf, 0),
                    q
                    );
        }
    }

    // If there are more than one partition,
    // update all of them except for the first.
    if (queue.size() > 1) {
        T sum{};

        for(unsigned d = 0; d < tail.size(); ++d) {
            if (src.part_size(d)) {
                sum += tail[d];
                sum += dst[src.part_start(d + 1) - 1];
                // Wrap partition into vector for ease of use:
                vex::vector<T> part(queue[d + 1], dst(d + 1));
                part += sum;
            }
        }
    }
}

/// Sort.
/**
 * If there are more than one device in vector's queue list, then all
 * partitions are sorted individually on GPUs and then merged on CPU.
 */
template <typename T>
void sort(vex::vector<T> &x) {
    auto queue = x.queue_list();

    for(unsigned d = 0; d < queue.size(); ++d) {
        if (x.part_size(d)) {
            boost::compute::command_queue q( queue[d]() );
            boost::compute::buffer buf( x(d).raw() );

            boost::compute::sort(
                    boost::compute::make_buffer_iterator<T>(buf, 0),
                    boost::compute::make_buffer_iterator<T>(buf, x.part_size(d)),
                    q
                    );
        }
    }

    // If there are multiple queues, merge the results on the CPU
    if (queue.size() > 1) {
        namespace fusion = boost::fusion;

        auto key_vectors  = fusion::vector_tie(x);
        auto host_vectors = detail::merge(key_vectors, vex::less<T>());
        fusion::for_each( detail::make_zip_view(host_vectors, key_vectors), detail::do_copy() );
    }
}

} // namespace compute
} // namespace vex

#endif
