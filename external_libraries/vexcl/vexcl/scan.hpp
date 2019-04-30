#ifndef VEXCL_SCAN_HPP
#define VEXCL_SCAN_HPP

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
 * \file   vexcl/scan.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Inclusive/Exclusive scan algortihms.

Adopted from Bolt code, see <https://github.com/HSA-Libraries/Bolt>.
The original code came with the following copyright notice:

\verbatim
Copyright 2012 - 2013 Advanced Micro Devices, Inc.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
\endverbatim
*/

#include <string>
#include <functional>

#include <vexcl/backend.hpp>
#include <vexcl/util.hpp>
#include <vexcl/vector.hpp>
#include <vexcl/function.hpp>

namespace vex {

namespace detail {

//---------------------------------------------------------------------------
template <int NT, typename T, typename Oper>
backend::kernel block_inclusive_scan(const backend::command_queue &queue)
{
    static detail::kernel_cache cache;

    auto kernel = cache.find(queue);

    if (kernel == cache.end()) {
        backend::source_generator src(queue);

        Oper::define(src, "oper");

        src.begin_kernel("block_inclusive_scan");
        src.begin_kernel_parameters();
        src.template parameter< size_t              >("n");
        src.template parameter< global_ptr<const T> >("input");
        src.template parameter< T                   >("identity");
        src.template parameter< global_ptr<T>       >("scan_buf1");
        src.template parameter< global_ptr<T>       >("scan_buf2");
        src.template parameter< int                 >("exclusive");
        src.end_kernel_parameters();

        src.new_line() << "size_t l_id  = " << src.local_id(0)   << ";";
        src.new_line() << "size_t g_id  = " << src.global_id(0)  << ";";
        src.new_line() << "size_t block = " << src.group_id(0)   << ";";

        src.new_line() << "size_t offset = 1;";

        {
            std::ostringstream shared;
            shared << "shared[" << 2 * NT << "]";
            src.smem_static_var(type_name<T>(), shared.str());
        }

        // load input into shared memory
        src.new_line()
            << "if(block * " << 2 * NT << " + l_id < n)"
            << " shared[l_id] = input[block * " << 2 * NT << " + l_id];";

        src.new_line()
            << "if(block * " << 2 * NT << " + l_id + " << NT << " < n)"
            << " shared[l_id + " << NT << "] ="
            << " input[block * " << 2 * NT << " + l_id + " << NT << "];";

        // Exclusive case
        src.new_line()
            << "if(exclusive && g_id == 0)"
            << " shared[l_id] = oper(identity, input[0]);";

        src.new_line() << "for (size_t start = " << NT << "; start > 0; start >>= 1, offset *= 2)";
        src.open("{");
        src.new_line().barrier();

        src.new_line() << "if (l_id < start)";
        src.open("{");
        src.new_line() << "size_t temp1 = offset * (2 * l_id + 1) - 1;";
        src.new_line() << "size_t temp2 = offset * (2 * l_id + 2) - 1;";
        src.new_line() << type_name<T>() << " y2 = shared[temp2];";
        src.new_line() << type_name<T>() << " y1 = shared[temp1];";
        src.new_line() << "shared[temp2] = oper(y2, y1);";
        src.close("}");

        src.close("}");
        src.new_line().barrier();

        src.new_line() << "if (l_id == 0)";
        src.open("{");
        src.new_line() << "scan_buf1[ block ] = shared[" << NT * 2 - 1 << "];";
        src.new_line() << "scan_buf2[ block ] = shared[" << NT - 1 << "];";
        src.close("}");
        src.end_kernel();

        kernel = cache.insert(queue, backend::kernel(
                    queue, src.str(), "block_inclusive_scan"));
    }

    return kernel->second;
}

template <int NT, typename T, typename Oper>
backend::kernel intra_block_inclusive_scan(const backend::command_queue &queue)
{
    static detail::kernel_cache cache;

    auto kernel = cache.find(queue);

    if (kernel == cache.end()) {
        backend::source_generator src(queue);

        Oper::define(src, "oper");

        src.begin_kernel("intra_block_inclusive_scan");
        src.begin_kernel_parameters();
        src.template parameter< size_t              >("n");
        src.template parameter< global_ptr<T>       >("post_sum");
        src.template parameter< global_ptr<const T> >("pre_sum");
        src.template parameter< T                   >("identity");
        src.template parameter< uint                >("work_per_thread");
        src.end_kernel_parameters();

        src.new_line() << "size_t l_id   = " << src.local_id(0)   << ";";
        src.new_line() << "size_t g_id   = " << src.global_id(0)  << ";";
        src.new_line() << "size_t map_id = g_id * work_per_thread;";

        {
            std::ostringstream shared;
            shared << "shared[" << NT << "]";
            src.smem_static_var(type_name<T>(), shared.str());
        }

        src.new_line() << "size_t offset;";
        src.new_line() << type_name<T>() << " work_sum;";

        src.new_line() << "if (map_id < n)";
        src.open("{");

        // accumulate zeroth value manually
        src.new_line() << "offset = 0;";
        src.new_line() << "work_sum = pre_sum[map_id];";

        //  Serial accumulation
        src.new_line() << "for( offset = 1; offset < work_per_thread; ++offset )";
        src.open("{");
        src.new_line()
            << "if (map_id + offset < n)"
            << " work_sum = oper( work_sum, pre_sum[map_id + offset] );";
        src.close("}");
        src.close("}");
        src.new_line().barrier();

        src.new_line() << type_name<T>() << " scan_sum = work_sum;";
        src.new_line() << "shared[ l_id ] = work_sum;";

        // scan in shared
        src.new_line() << "for( offset = 1; offset < " << NT << "; offset *= 2 )";
        src.open("{");
        src.new_line().barrier();

        src.new_line()
            << "if (map_id < n && l_id >= offset)"
            << " scan_sum = oper( scan_sum, shared[ l_id - offset ] );";
        src.new_line().barrier();
        src.new_line() << "shared[ l_id ] = scan_sum;";
        src.close("}");
        src.new_line().barrier();

        // write final scan from pre-scan and shared scan
        src.new_line() << "work_sum = pre_sum[map_id];";
        src.new_line() << "if (l_id > 0)";
        src.open("{");
        src.new_line() << "    work_sum = oper(work_sum, shared[l_id - 1]);";
        src.new_line() << "    post_sum[map_id] = work_sum;";
        src.close("}");
        src.new_line() << "else post_sum[map_id] = work_sum;";

        src.new_line() << "for( offset = 1; offset < work_per_thread; ++offset )";
        src.open("{");
        src.new_line().barrier();

        src.new_line() << "if (map_id < n && l_id > 0)";
        src.open("{");
        src.new_line() << type_name<T>() << " y = oper(pre_sum[map_id + offset], work_sum);";
        src.new_line() << "post_sum[ map_id + offset ] = y;";
        src.new_line() << "work_sum = y;";
        src.close("}");
        src.new_line() << "else";
        src.open("{");
        src.new_line() << "post_sum[map_id + offset] = oper(pre_sum[map_id + offset], work_sum);";
        src.new_line() << "work_sum = post_sum[map_id + offset];";
        src.close("}");
        src.close("}");
        src.end_kernel();

        kernel = cache.insert(queue, backend::kernel(
                    queue, src.str(), "intra_block_inclusive_scan"));
    }

    return kernel->second;
}

//---------------------------------------------------------------------------
template <int NT, typename T, typename Oper>
backend::kernel block_addition(
        const backend::command_queue &queue)
{
    static detail::kernel_cache cache;

    auto kernel = cache.find(queue);

    if (kernel == cache.end()) {
        backend::source_generator src(queue);

        Oper::define(src, "oper");

        src.begin_kernel("block_addition");
        src.begin_kernel_parameters();
        src.template parameter< size_t              >("n");
        src.template parameter< global_ptr<const T> >("input");
        src.template parameter< global_ptr<T>       >("output");
        src.template parameter< global_ptr<T>       >("post_sum");
        src.template parameter< global_ptr<T>       >("pre_sum");
        src.template parameter< T                   >("identity");
        src.template parameter< int                 >("exclusive");
        src.end_kernel_parameters();

        src.new_line() << "size_t l_id  = " << src.local_id(0)   << ";";
        src.new_line() << "size_t g_id  = " << src.global_id(0)  << ";";
        src.new_line() << "size_t block = " << src.group_id(0)   << ";";

        src.new_line() << type_name<T>() << " val;";

        {
            std::ostringstream shared;
            shared << "shared[" << NT << "]";
            src.smem_static_var(type_name<T>(), shared.str());
        }

        src.new_line() << "if (g_id < n)";
        src.open("{");
        src.new_line() << "if (exclusive) val = g_id > 0 ? input[g_id - 1] : identity;";
        src.new_line() << "else val = input[g_id];";
        src.close("}");
        src.new_line() << "shared[l_id] = val;";

        src.new_line() << type_name<T>() << " scan_result = val;";
        src.new_line() << type_name<T>() << " post_block_sum, new_result;";
        src.new_line() << type_name<T>() << " y1, y2, sum;";

        src.new_line() << "if(l_id == 0 && g_id < n)";
        src.open("{");
        src.new_line() << "if(block > 0)";
        src.open("{");
        src.new_line() << "if(block % 2 == 0)  post_block_sum = post_sum[ block/2 - 1 ];";
        src.new_line() << "else if(block == 1) post_block_sum = pre_sum[0];";
        src.new_line() << "else";
        src.open("{");
        src.new_line() << "y1 = post_sum[ block/2 - 1 ];";
        src.new_line() << "y2 = pre_sum [ block/2];";
        src.new_line() << "post_block_sum = oper(y1, y2);";
        src.close("}");
        src.new_line() << "new_result = exclusive ? post_block_sum : oper( scan_result, post_block_sum );";
        src.close("}");
        src.new_line() << "else new_result = scan_result;";

        src.new_line() << "shared[ l_id ] = new_result;";
        src.close("}");

        //  Computes a scan within a workgroup
        src.new_line() << "sum = shared[ l_id ];";
        src.new_line() << "for( size_t offset = 1; offset < " << NT << "; offset *= 2 )";
        src.open("{");
        src.new_line().barrier();
        src.new_line() << "if (l_id >= offset) sum = oper( sum, shared[ l_id - offset ] );";
        src.new_line().barrier();
        src.new_line() << "shared[ l_id ] = sum;";
        src.close("}");
        src.new_line().barrier();
        src.new_line() << "if(g_id < n) output[ g_id ] = sum;";

        src.end_kernel();

        kernel = cache.insert(queue, backend::kernel(
                    queue, src.str(), "block_addition"));
    }

    return kernel->second;
}

template <typename T, typename Oper>
void scan(
        backend::command_queue    const &queue,
        backend::device_vector<T> const &input,
        backend::device_vector<T>       &output,
        T init,
        bool exclusive,
        Oper
        )
{
    precondition(
            input.size() == output.size(),
            "Wrong output size in inclusive_scan"
            );

    backend::select_context(queue);

    const int NT_cpu = 1;
    const int NT_gpu = 256;
    const int NT = is_cpu(queue) ? NT_cpu : NT_gpu;
    const int NT2 = 2 * NT;

    int do_exclusive = exclusive ? 1 : 0;

    const size_t count         = input.size();
    const size_t num_blocks    = (count + NT2 - 1) / NT2;
    const size_t scan_buf_size = alignup(num_blocks, NT2);

    backend::device_vector<T> pre_sum1(queue, scan_buf_size);
    backend::device_vector<T> pre_sum2(queue, scan_buf_size);
    backend::device_vector<T> post_sum(queue, scan_buf_size);

    // Kernel0
    auto krn0 = is_cpu(queue) ?
        block_inclusive_scan<NT_cpu, T, Oper>(queue) :
        block_inclusive_scan<NT_gpu, T, Oper>(queue);

    krn0.push_arg(count);
    krn0.push_arg(input);
    krn0.push_arg(init);
    krn0.push_arg(pre_sum1);
    krn0.push_arg(pre_sum2);
    krn0.push_arg(do_exclusive);

    krn0.config(num_blocks, NT);

    krn0(queue);

    // Kernel1
    auto krn1 = is_cpu(queue) ?
        intra_block_inclusive_scan<NT_cpu, T, Oper>(queue) :
        intra_block_inclusive_scan<NT_gpu, T, Oper>(queue);

    uint work_per_thread = std::max<uint>(1U, static_cast<uint>(scan_buf_size / NT));
    krn1.push_arg(num_blocks);
    krn1.push_arg(post_sum);
    krn1.push_arg(pre_sum1);
    krn1.push_arg(init);
    krn1.push_arg(work_per_thread);

    krn1.config(1, NT);

    krn1(queue);

    // Kernel2
    auto krn2 = is_cpu(queue) ?
        block_addition<NT_cpu, T, Oper>(queue) :
        block_addition<NT_gpu, T, Oper>(queue);

    krn2.push_arg(count);
    krn2.push_arg(input);
    krn2.push_arg(output);
    krn2.push_arg(post_sum);
    krn2.push_arg(pre_sum2);
    krn2.push_arg(init);
    krn2.push_arg(do_exclusive);

    krn2.config(num_blocks * 2, NT);

    krn2(queue);
}

} // namespace detail

/// Binary function object class whose call returns the result of adding its two arguments.
template <typename T>
struct plus : std::plus<T> {
    VEX_FUNCTION(T, device, (T, x)(T, y), return x + y;);

    plus() {}
};

/// Inclusive scan.
template <typename T, class Oper>
void inclusive_scan(
        vector<T> const &input,
        vector<T>       &output,
        T init,
        Oper oper
        )
{
    precondition(
            input.nparts() == output.nparts(),
            "Incompatible partitioning"
            );

    auto &queue = input.queue_list();

    for(unsigned d = 0; d < queue.size(); ++d)
        detail::scan(queue[d], input(d), output(d), init, false, oper.device);

    std::vector<T> tail(queue.size() - 1);

    for(unsigned d = 1; d < queue.size(); ++d)
        if (size_t head = output.part_start(d))
            tail[d - 1] = output[head - 1];

    std::partial_sum(tail.begin(), tail.end(), tail.begin(), oper);

    for(unsigned d = 1; d < queue.size(); ++d)
        if (output.part_start(d)) {
            vector<T> part(queue[d], output(d));
            part = oper.device(part, tail[d - 1]);
        }
}

/// Inclusive scan.
template <typename T>
void inclusive_scan(
        vector<T> const &input,
        vector<T>       &output,
        T init = T()
        )
{
    inclusive_scan(input, output, init, plus<T>());
}

/// Exclusive scan.
template <typename T, class Oper>
void exclusive_scan(
        vector<T> const &input,
        vector<T>       &output,
        T init,
        Oper oper
        )
{
    precondition(
            input.nparts() == output.nparts(),
            "Incompatible partitioning"
            );

    auto &queue = input.queue_list();

    std::vector<T> tail(queue.size() - 1);

    for(unsigned d = 1; d < queue.size(); ++d)
        if (size_t head = input.part_start(d))
            tail[d - 1] = input[head - 1];

    for(unsigned d = 0; d < queue.size(); ++d)
        detail::scan(queue[d], input(d), output(d), init, true, oper.device);

    for(unsigned d = 1; d < queue.size(); ++d)
        if (size_t head = output.part_start(d))
            tail[d - 1] = oper(tail[d - 1], output[head - 1]);

    std::partial_sum(tail.begin(), tail.end(), tail.begin(), oper);

    for(unsigned d = 1; d < queue.size(); ++d)
        if (output.part_start(d)) {
            vector<T> part(queue[d], output(d));
            part = oper.device(part, tail[d - 1]);
        }
}

/// Exclusive scan.
template <typename T>
void exclusive_scan(
        vector<T> const &input,
        vector<T>       &output,
        T init = T()
        )
{
    exclusive_scan(input, output, init, plus<T>());
}

} // namespace vex

#endif
