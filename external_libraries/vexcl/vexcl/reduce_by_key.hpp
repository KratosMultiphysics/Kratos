#ifndef VEXCL_REDUCE_BY_KEY_HPP
#define VEXCL_REDUCE_BY_KEY_HPP

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
 * \file   vexcl/reduce_by_key.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Reduce by key algortihm.

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

#include <vexcl/vector.hpp>
#include <vexcl/scan.hpp>
#include <vexcl/detail/fusion.hpp>
#include <vexcl/function.hpp>

namespace vex {
namespace detail {
namespace rbk {

//---------------------------------------------------------------------------
template <typename T, class Comp>
backend::kernel offset_calculation(const backend::command_queue &queue) {
    static detail::kernel_cache cache;

    auto kernel = cache.find(queue);

    if (kernel == cache.end()) {
        backend::source_generator src(queue);

        Comp::define(src, "comp");

        src.begin_kernel("offset_calculation");
        src.begin_kernel_parameters();
        src.template parameter< size_t >("n");

        boost::mpl::for_each<T>(pointer_param<global_ptr, true>(src, "keys"));

        src.template parameter< global_ptr<int> >("offsets");
        src.end_kernel_parameters();

        src.new_line().grid_stride_loop().open("{");
        src.new_line()
            << "if (idx > 0)"
            << " offsets[idx] = !comp(";
        for(int p = 0; p < boost::mpl::size<T>::value; ++p)
            src << (p ? ", " : "") << "keys" << p << "[idx - 1]";
        for(int p = 0; p < boost::mpl::size<T>::value; ++p)
            src << ", keys" << p << "[idx]";
        src << ");";
        src.new_line() << "else offsets[idx] = 0;";
        src.close("}");
        src.end_kernel();

        kernel = cache.insert(queue, backend::kernel(
                    queue, src.str(), "offset_calculation"));
    }

    return kernel->second;
}

//---------------------------------------------------------------------------
template <int NT, typename T, class Oper>
backend::kernel block_scan_by_key(const backend::command_queue &queue) {
    static detail::kernel_cache cache;

    auto kernel = cache.find(queue);

    if (kernel == cache.end()) {
        backend::source_generator src(queue);

        Oper::define(src, "oper");

        src.begin_kernel("block_scan_by_key");
        src.begin_kernel_parameters();
        src.template parameter< size_t                >("n");
        src.template parameter< global_ptr<const int> >("keys");
        src.template parameter< global_ptr<const T>   >("vals");
        src.template parameter< global_ptr<T>         >("output");
        src.template parameter< global_ptr<int>       >("key_buf");
        src.template parameter< global_ptr<T>         >("val_buf");
        src.end_kernel_parameters();

        src.new_line() << "size_t l_id  = " << src.local_id(0)   << ";";
        src.new_line() << "size_t g_id  = " << src.global_id(0)  << ";";
        src.new_line() << "size_t block = " << src.group_id(0)   << ";";

        src.new_line() << "struct Shared";
        src.open("{");
            src.new_line() << "int keys[" << NT << "];";
            src.new_line() << type_name<T>() << " vals[" << NT << "];";
        src.close("};");

        src.smem_static_var("struct Shared", "shared");

        src.new_line() << "int key;";
        src.new_line() << type_name<T>() << " val;";

        src.new_line() << "if (g_id < n)";
        src.open("{");
        src.new_line() << "key = keys[g_id];";
        src.new_line() << "val = vals[g_id];";
        src.new_line() << "shared.keys[l_id] = key;";
        src.new_line() << "shared.vals[l_id] = val;";
        src.close("}");

        // Computes a scan within a workgroup updates vals in lds but not keys
        src.new_line() << type_name<T>() << " sum = val;";
        src.new_line() << "for(size_t offset = 1; offset < " << NT << "; offset *= 2)";
        src.open("{");
        src.new_line().barrier();
        src.new_line() << "if (l_id >= offset && shared.keys[l_id - offset] == key)";
        src.open("{");
        src.new_line() << "sum = oper(sum, shared.vals[l_id - offset]);";
        src.close("}");
        src.new_line().barrier();
        src.new_line() << "shared.vals[l_id] = sum;";
        src.close("}");
        src.new_line().barrier();

        src.new_line() << "if (g_id >= n) return;";

        // Each work item writes out its calculated scan result, relative to the
        // beginning of each work group
        src.new_line() << "int key2 = -1;";
        src.new_line() << "if (g_id < n - 1) key2 = keys[g_id + 1];";
        src.new_line() << "if (key != key2) output[g_id] = sum;";

        src.new_line() << "if (l_id == 0)";
        src.open("{");
        src.new_line() << "key_buf[block] = shared.keys[" << NT - 1 << "];";
        src.new_line() << "val_buf[block] = shared.vals[" << NT - 1 << "];";
        src.close("}");

        src.end_kernel();

        kernel = cache.insert(queue, backend::kernel(
                    queue, src.str(), "block_scan_by_key"));
    }

    return kernel->second;
}

//---------------------------------------------------------------------------
template <int NT, typename T, class Oper>
backend::kernel block_inclusive_scan_by_key(const backend::command_queue &queue)
{
    static detail::kernel_cache cache;

    auto kernel = cache.find(queue);

    if (kernel == cache.end()) {
        backend::source_generator src(queue);

        Oper::define(src, "oper");

        src.begin_kernel("block_inclusive_scan_by_key");
        src.begin_kernel_parameters();
        src.template parameter< size_t                >("n");
        src.template parameter< global_ptr<const int> >("key_sum");
        src.template parameter< global_ptr<const T>   >("pre_sum");
        src.template parameter< global_ptr<T>         >("post_sum");
        src.template parameter< cl_uint               >("work_per_thread");
        src.end_kernel_parameters();

        src.new_line() << "size_t l_id   = " << src.local_id(0)   << ";";
        src.new_line() << "size_t g_id   = " << src.global_id(0)  << ";";
        src.new_line() << "size_t map_id = g_id * work_per_thread;";

        src.new_line() << "struct Shared";
        src.open("{");
            src.new_line() << "int keys[" << NT << "];";
            src.new_line() << type_name<T>() << " vals[" << NT << "];";
        src.close("};");

        src.smem_static_var("struct Shared", "shared");

        src.new_line() << "uint offset;";
        src.new_line() << "int  key;";
        src.new_line() << type_name<T>() << " work_sum;";

        src.new_line() << "if (map_id < n)";
        src.open("{");
        src.new_line() << "int prev_key;";

        // accumulate zeroth value manually
        src.new_line() << "offset   = 0;";
        src.new_line() << "key      = key_sum[map_id];";
        src.new_line() << "work_sum = pre_sum[map_id];";

        src.new_line() << "post_sum[map_id] = work_sum;";

        //  Serial accumulation
        src.new_line() << "for( offset = offset + 1; offset < work_per_thread; ++offset )";
        src.open("{");
        src.new_line() << "prev_key = key;";
        src.new_line() << "key      = key_sum[ map_id + offset ];";

        src.new_line() << "if ( map_id + offset < n )";
        src.open("{");
        src.new_line() << type_name<T>() << " y = pre_sum[ map_id + offset ];";

        src.new_line() << "if ( key == prev_key ) work_sum = oper( work_sum, y );";
        src.new_line() << "else work_sum = y;";

        src.new_line() << "post_sum[ map_id + offset ] = work_sum;";
        src.close("}");
        src.close("}");
        src.close("}");
        src.new_line().barrier();

        // load LDS with register sums
        src.new_line() << "shared.vals[ l_id ] = work_sum;";
        src.new_line() << "shared.keys[ l_id ] = key;";

        // scan in lds
        src.new_line() << type_name<T>() << " scan_sum = work_sum;";

        src.new_line() << "for( offset = 1; offset < " << NT << "; offset *= 2 )";
        src.open("{");
        src.new_line().barrier();

        src.new_line() << "if (map_id < n)";
        src.open("{");
        src.new_line() << "if (l_id >= offset)";
        src.open("{");
        src.new_line() << "int key1 = shared.keys[ l_id ];";
        src.new_line() << "int key2 = shared.keys[ l_id - offset ];";

        src.new_line() << "if ( key1 == key2 ) scan_sum = oper( scan_sum, shared.vals[ l_id - offset ] );";
        src.new_line() << "else scan_sum = shared.vals[ l_id ];";
        src.close("}");

        src.close("}");
        src.new_line().barrier();

        src.new_line() << "shared.vals[ l_id ] = scan_sum;";
        src.close("}");

        src.new_line().barrier();

        // write final scan from pre-scan and lds scan
        src.new_line() << "for( offset = 0; offset < work_per_thread; ++offset )";
        src.open("{");
        src.new_line().barrier(true);

        src.new_line() << "if (map_id < n && l_id > 0)";
        src.open("{");
        src.new_line() << type_name<T>() << " y = post_sum[ map_id + offset ];";
        src.new_line() << "int key1 = key_sum    [ map_id + offset ];";
        src.new_line() << "int key2 = shared.keys[ l_id - 1 ];";

        src.new_line() << "if ( key1 == key2 ) y = oper( y, shared.vals[l_id - 1] );";

        src.new_line() << "post_sum[ map_id + offset ] = y;";
        src.close("}");
        src.close("}");

        src.end_kernel();

        kernel = cache.insert(queue, backend::kernel(
                    queue, src.str(), "block_inclusive_scan_by_key"));
    }

    return kernel->second;
}

//---------------------------------------------------------------------------
template <typename T, class Oper>
backend::kernel block_sum_by_key(const backend::command_queue &queue) {
    static detail::kernel_cache cache;

    auto kernel = cache.find(queue);

    if (kernel == cache.end()) {
        backend::source_generator src(queue);

        Oper::define(src, "oper");

        src.begin_kernel("block_sum_by_key");
        src.begin_kernel_parameters();
        src.template parameter< size_t                >("n");
        src.template parameter< global_ptr<const int> >("key_sum");
        src.template parameter< global_ptr<const T>   >("post_sum");
        src.template parameter< global_ptr<const int> >("keys");
        src.template parameter< global_ptr<T>         >("output");
        src.end_kernel_parameters();

        src.new_line() << "size_t g_id  = " << src.global_id(0)  << ";";
        src.new_line() << "size_t block = " << src.group_id(0)   << ";";

        src.new_line() << "if (g_id >= n) return;";

        // accumulate prefix
        src.new_line() << "int key2 = keys[ g_id ];";
        src.new_line() << "int key1 = (block > 0    ) ? key_sum[ block - 1 ] : key2 - 1;";
        src.new_line() << "int key3 = (g_id  < n - 1) ? keys   [ g_id  + 1 ] : key2 - 1;";

        src.new_line() << "if (block > 0 && key1 == key2 && key2 != key3)";
        src.open("{");
        src.new_line() << type_name<T>() << " scan_result    = output  [ g_id      ];";
        src.new_line() << type_name<T>() << " post_block_sum = post_sum[ block - 1 ];";
        src.new_line() << "output[ g_id ] = oper( scan_result, post_block_sum );";
        src.close("}");

        src.end_kernel();

        kernel = cache.insert(queue, backend::kernel(
                    queue, src.str(), "block_sum_by_key"));
    }

    return kernel->second;
}

//---------------------------------------------------------------------------
template <typename K, typename V>
backend::kernel key_value_mapping(const backend::command_queue &queue) {
    static detail::kernel_cache cache;

    auto kernel = cache.find(queue);

    if (kernel == cache.end()) {
        backend::source_generator src(queue);

        src.begin_kernel("key_value_mapping");
        src.begin_kernel_parameters();
        src.template parameter< size_t >("n");

        boost::mpl::for_each<K>(pointer_param<global_ptr, true>(src, "ikeys"));
        boost::mpl::for_each<K>(pointer_param<global_ptr      >(src, "okeys"));

        src.template parameter< global_ptr<V>       >("ovals");
        src.template parameter< global_ptr<int>     >("offset");
        src.template parameter< global_ptr<const V> >("ivals");
        src.end_kernel_parameters();

        src.new_line().grid_stride_loop().open("{");

        src.new_line() << "int num_sections = offset[n - 1] + 1;";

        src.new_line() << "int off = offset[idx];";
        src.new_line() << "if (idx < (n - 1) && off != offset[idx + 1])";
        src.open("{");
        for(int p = 0; p < boost::mpl::size<K>::value; ++p)
            src.new_line() << "okeys" << p << "[off] = ikeys" << p << "[idx];";
        src.new_line() << "ovals[off] = ivals[idx];";
        src.close("}");

        src.new_line() << "if (idx == (n - 1))";
        src.open("{");
        for(int p = 0; p < boost::mpl::size<K>::value; ++p)
            src.new_line() << "okeys" << p << "[num_sections - 1] = ikeys" << p << "[idx];";
        src.new_line() << "ovals[num_sections - 1] = ivals[idx];";
        src.close("}");

        src.close("}");

        src.end_kernel();

        kernel = cache.insert(queue, backend::kernel(
                    queue, src.str(), "key_value_mapping"));
    }

    return kernel->second;
}

struct do_vex_resize {
    const std::vector<backend::command_queue> &q;
    size_t n;

    do_vex_resize(const std::vector<backend::command_queue> &q, size_t n)
        : q(q), n(n) {}

    template <class V>
    void operator()(V &v) const {
        v.resize(q, n);
    }
};

struct do_push_arg {
    backend::kernel &k;

    do_push_arg(backend::kernel &k) : k(k) {}

    template <class T>
    void operator()(const T &t) const {
        k.push_arg( t(0) );
    }
};

template <typename IKTuple, typename OKTuple, typename V, class Comp, class Oper>
int reduce_by_key_sink(
        IKTuple &&ikeys, vector<V> const &ivals,
        OKTuple &&okeys, vector<V>       &ovals,
        Comp, Oper
        )
{
    namespace fusion = boost::fusion;
    typedef typename extract_value_types<IKTuple>::type K;

    static_assert(
            std::is_same<K, typename extract_value_types<OKTuple>::type>::value,
            "Incompatible input and output key types");

    precondition(
            fusion::at_c<0>(ikeys).nparts() == 1 && ivals.nparts() == 1,
            "reduce_by_key is only supported for single device contexts"
            );

    precondition(fusion::at_c<0>(ikeys).size() == ivals.size(),
            "keys and values should have same size"
            );

    const auto &queue = fusion::at_c<0>(ikeys).queue_list();
    backend::select_context(queue[0]);

    const int NT_cpu = 1;
    const int NT_gpu = 256;
    const int NT = is_cpu(queue[0]) ? NT_cpu : NT_gpu;

    size_t count         = fusion::at_c<0>(ikeys).size();
    size_t num_blocks    = (count + NT - 1) / NT;
    size_t scan_buf_size = alignup(num_blocks, NT);

    backend::device_vector<int> key_sum   (queue[0], scan_buf_size);
    backend::device_vector<V>   pre_sum   (queue[0], scan_buf_size);
    backend::device_vector<V>   post_sum  (queue[0], scan_buf_size);
    backend::device_vector<V>   offset_val(queue[0], count);
    backend::device_vector<int> offset    (queue[0], count);

    /***** Kernel 0 *****/
    auto krn0 = offset_calculation<K, Comp>(queue[0]);

    krn0.push_arg(count);
    boost::fusion::for_each(ikeys, do_push_arg(krn0));
    krn0.push_arg(offset);

    krn0(queue[0]);

    VEX_FUNCTION(int, plus, (int, x)(int, y), return x + y;);
    scan(queue[0], offset, offset, 0, false, plus);

    /***** Kernel 1 *****/
    auto krn1 = is_cpu(queue[0]) ?
        block_scan_by_key<NT_cpu, V, Oper>(queue[0]) :
        block_scan_by_key<NT_gpu, V, Oper>(queue[0]);

    krn1.push_arg(count);
    krn1.push_arg(offset);
    krn1.push_arg(ivals(0));
    krn1.push_arg(offset_val);
    krn1.push_arg(key_sum);
    krn1.push_arg(pre_sum);

    krn1.config(num_blocks, NT);
    krn1(queue[0]);

    /***** Kernel 2 *****/
    uint work_per_thread = std::max<uint>(1U, static_cast<uint>(scan_buf_size / NT));

    auto krn2 = is_cpu(queue[0]) ?
        block_inclusive_scan_by_key<NT_cpu, V, Oper>(queue[0]) :
        block_inclusive_scan_by_key<NT_gpu, V, Oper>(queue[0]);

    krn2.push_arg(num_blocks);
    krn2.push_arg(key_sum);
    krn2.push_arg(pre_sum);
    krn2.push_arg(post_sum);
    krn2.push_arg(work_per_thread);

    krn2.config(1, NT);
    krn2(queue[0]);

    /***** Kernel 3 *****/
    auto krn3 = block_sum_by_key<V, Oper>(queue[0]);

    krn3.push_arg(count);
    krn3.push_arg(key_sum);
    krn3.push_arg(post_sum);
    krn3.push_arg(offset);
    krn3.push_arg(offset_val);

    krn3.config(num_blocks, NT);
    krn3(queue[0]);

    /***** resize okeys and ovals *****/
    int out_elements = 0;
    offset.read(queue[0], count - 1, 1, &out_elements, true);
    ++out_elements;

    boost::fusion::for_each(okeys, do_vex_resize(queue, out_elements));
    ovals.resize(ivals.queue_list(), out_elements);

    /***** Kernel 4 *****/
    auto krn4 = key_value_mapping<K, V>(queue[0]);

    krn4.push_arg(count);
    boost::fusion::for_each(ikeys, do_push_arg(krn4));
    boost::fusion::for_each(okeys, do_push_arg(krn4));
    krn4.push_arg(ovals(0));
    krn4.push_arg(offset);
    krn4.push_arg(offset_val);

    krn4(queue[0]);

    return out_elements;
}

} // namespace rbk
} // namespace detail

/// Reduce by key algorithm.
template <typename IKeys, typename OKeys, typename V, class Comp, class Oper>
int reduce_by_key(
        IKeys &&ikeys, vector<V> const &ivals,
        OKeys &&okeys, vector<V>       &ovals,
        Comp comp, Oper oper
        )
{
    return detail::rbk::reduce_by_key_sink(
            detail::forward_as_sequence(ikeys), ivals,
            detail::forward_as_sequence(okeys), ovals,
            comp, oper);
}

/// Reduce by key algorithm.
template <typename K, typename V>
int reduce_by_key(
        vector<K> const &ikeys,
        vector<V> const &ivals,
        vector<K>       &okeys,
        vector<V>       &ovals
        )
{
    VEX_FUNCTION(bool, equal, (K, x)(K, y), return x == y;);
    VEX_FUNCTION(V, plus, (V, x)(V, y), return x + y;);
    return reduce_by_key(ikeys, ivals, okeys, ovals, equal, plus);
}

}

#endif
