#ifndef VEXCL_SCAN_BY_KEY_HPP
#define VEXCL_SCAN_BY_KEY_HPP

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
 * \file   vexcl/scan_by_key.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Scan by key algortihm.

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
#include <vexcl/detail/fusion.hpp>
#include <vexcl/function.hpp>

namespace vex {
namespace detail {
namespace sbk {

struct gcc46_workaround {
    backend::source_generator &src;
    int wgsz, pos;

    gcc46_workaround(backend::source_generator &src, int wgsz)
        : src(src), wgsz(wgsz), pos(0) {}

    template <class T>
    void operator()(T) {
        src.new_line() << type_name<T>() << " keys" << pos++ << "[" << wgsz << "];";
    }
};

//---------------------------------------------------------------------------
template <int NT, typename K, typename V, class Comp, class Oper, bool exclusive>
backend::kernel block_scan_by_key(const backend::command_queue &queue) {
    static detail::kernel_cache cache;

    auto kernel = cache.find(queue);

    if (kernel == cache.end()) {
        backend::source_generator src(queue);

        Comp::define(src, "comp");
        Oper::define(src, "oper");

        src.begin_kernel( "block_scan_by_key");
        src.begin_kernel_parameters();
        src.template parameter< size_t              >("n");
        src.template parameter< global_ptr<const V> >("ivals");
        src.template parameter< global_ptr<      V> >("ovals1");
        src.template parameter< global_ptr<      V> >("ovals2");

        boost::mpl::for_each<K>(pointer_param<global_ptr, true>(src, "ikeys"));
        boost::mpl::for_each<K>(pointer_param<global_ptr      >(src, "okeys"));

        if (exclusive) src.template parameter<V>("init");

        src.end_kernel_parameters();

        src.new_line() << "size_t g_id   = " << src.global_id(0)  << ";";
        src.new_line() << "size_t l_id   = " << src.local_id(0)   << ";";
        src.new_line() << "size_t block  = " << src.group_id(0)   << ";";
        src.new_line() << "size_t offset = 1;";

        const int    wgsz = NT * 2;
        const size_t nK   = boost::mpl::size<K>::value;

        src.new_line() << "size_t pos = block * " << wgsz << " + l_id;";

        src.new_line() << "struct Shared";
        src.open("{");
        src.new_line() << type_name<V>() << " vals[" << wgsz << "];";

        // gcc 4.6 crashes if the following type iteration is done with lambda.
        // so here it goes:
        boost::mpl::for_each<K>( gcc46_workaround(src, wgsz) );
        src.close("};");
        src.smem_static_var("struct Shared", "shared");

        if (exclusive) {
            src.new_line() << "if (g_id > 0 && pos < n)";
            src.open("{");

            boost::mpl::for_each<K>(
                    type_iterator([&](size_t p, std::string tname) {
                        src.new_line() << tname << " key1" << p << " = ikeys" << p << "[pos];";
                        src.new_line() << tname << " key2" << p << " = ikeys" << p << "[pos - 1];";
                        })
                    );

            src.new_line() << "if (comp(key10";
            for(size_t p = 1; p < nK; ++p) src << ", key1" << p;
            for(size_t p = 0; p < nK; ++p) src << ", key2" << p;
            src << "))";
            src.open("{");
            src.new_line() << "shared.vals[l_id] = ivals[pos];";
            src.close("}");
            src.new_line() << "else";
            src.open("{");
            src.new_line() << "shared.vals[l_id] = oper(init, ivals[pos]);";
            src.close("}");
            for(size_t p = 0; p < nK; ++p)
                src.new_line() << "shared.keys" << p << "[l_id] = ikeys" << p << "[pos];";

            src.close("}");
            src.new_line() << "else";
            src.open("{");

            src.new_line() << "shared.vals[l_id] = oper(init, ivals[0]);";
            for(size_t p = 0; p < nK; ++p)
                src.new_line() << "shared.keys" << p << "[l_id] = ikeys" << p << "[0];";

            src.close("}");

            src.new_line() << "if (pos + " << NT << " < n)";
            src.open("{");

            boost::mpl::for_each<K>(
                    type_iterator([&](size_t p, std::string tname) {
                        src.new_line()
                            << tname << " key1" << p << " = ikeys" << p <<
                            "[pos + " << NT << "];";
                        src.new_line()
                            << tname << " key2" << p << " = ikeys" << p <<
                            "[pos + " << NT << " - 1];";
                        })
                    );

            src.new_line() << "if (comp(key10";
            for(size_t p = 1; p < nK; ++p) src << ", key1" << p;
            for(size_t p = 0; p < nK; ++p) src << ", key2" << p;
            src << "))";
            src.open("{");
            src.new_line()
                << "shared.vals[l_id + " << NT << "] = ivals[pos + " << NT << "];";
            src.close("}");
            src.new_line() << "else";
            src.open("{");
            src.new_line()
                << "shared.vals[l_id + " << NT << "] = oper(init, ivals[pos + " << NT << "]);";
            src.close("}");

            for(size_t p = 0; p < nK; ++p)
                src.new_line()
                    << "shared.keys" << p << "[l_id + " << NT <<
                    "] = ikeys" << p << "[pos + " << NT << "];";

            src.close("}");

        } else { // inclusive
            src.new_line() << "if (pos < n)";
            src.open("{");
            src.new_line() << "shared.vals[l_id] = ivals[pos];";
            for(size_t p = 0; p < nK; ++p)
                src.new_line() << "shared.keys" << p << "[l_id] = ikeys" << p << "[pos];";
            src.close("}");

            src.new_line() << "if (pos + " << NT << " < n)";
            src.open("{");
            src.new_line() << "shared.vals[l_id + " << NT << "] = ivals[pos + " << NT << "];";
            for(size_t p = 0; p < nK; ++p)
                src.new_line() << "shared.keys" << p << "[l_id + " << NT << "] = ikeys" << p << "[pos + " << NT << "];";
            src.close("}");
        }

        src.new_line() << "for(size_t start = " << NT << "; start > 0; start /= 2)";
        src.open("{");
        src.new_line().barrier();
        src.new_line() << "if (l_id < start)";
        src.open("{");

        src.new_line() << "size_t temp1 = offset * (2 * l_id + 1) - 1;";
        src.new_line() << "size_t temp2 = offset * (2 * l_id + 2) - 1;";

        boost::mpl::for_each<K>(
                type_iterator([&](size_t p, std::string tname) {
                    src.new_line() << tname << " key1" << p << " = shared.keys" << p << "[temp1];";
                    src.new_line() << tname << " key2" << p << " = shared.keys" << p << "[temp2];";
                    })
                );

        src.new_line() << "if (comp(key20";
        for(size_t p = 1; p < nK; ++p) src << ", key2" << p;
        for(size_t p = 0; p < nK; ++p) src << ", key1" << p;
        src << "))";
        src.open("{");
        src.new_line() << "shared.vals[temp2] = oper(shared.vals[temp2], shared.vals[temp1]);";
        src.close("}");
        src.close("}");
        src.new_line() << "offset *= 2;";
        src.close("}");

        src.new_line().barrier();
        src.new_line() << "if (l_id == 0)";
        src.open("{");
        for(size_t p = 0; p < nK; ++p)
            src.new_line() << "okeys" << p << "[block] = shared.keys" << p << "[" << wgsz - 1 << "];";
        src.new_line() << "ovals1[block] = shared.vals[" << wgsz - 1 << "];";
        src.new_line() << "ovals2[block] = shared.vals[" << NT - 1 << "];";
        src.close("}");
        src.end_kernel();

        kernel = cache.insert(queue, backend::kernel(
                    queue, src.str(), "block_scan_by_key"));
    }

    return kernel->second;
}

//---------------------------------------------------------------------------
template <int NT, typename K, typename V, class Comp, class Oper>
backend::kernel block_inclusive_scan_by_key(const backend::command_queue &queue)
{
    static detail::kernel_cache cache;

    auto kernel = cache.find(queue);

    if (kernel == cache.end()) {
        backend::source_generator src(queue);

        Comp::define(src, "comp");
        Oper::define(src, "oper");

        const size_t nK = boost::mpl::size<K>::value;

        src.begin_kernel("block_inclusive_scan_by_key");
        src.begin_kernel_parameters();
        src.template parameter< size_t        >("n");
        src.template parameter< global_ptr<V> >("pre_sum");
        src.template parameter< cl_uint       >("work_per_thread");

        boost::mpl::for_each<K>(pointer_param<global_ptr, true>(src, "key_sum"));

        src.end_kernel_parameters();

        src.new_line() << "size_t block  = " << src.group_id(0)  << ";";
        src.new_line() << "size_t g_id   = " << src.global_id(0) << ";";
        src.new_line() << "size_t l_id   = " << src.local_id(0)  << ";";
        src.new_line() << "size_t map_id = g_id * work_per_thread;";

        src.new_line() << "struct Shared";
        src.open("{");
            src.new_line() << type_name<V>() << " vals[" << NT << "];";
            boost::mpl::for_each<K>(
                    type_iterator([&](size_t p, std::string tname) {
                        src.new_line()
                            << tname << " keys" << p << "[" << NT << "];";
                        })
                    );
        src.close("};");
        src.smem_static_var("struct Shared", "shared");

        // do offset of zero manually
        src.new_line() << "uint offset;";
        boost::mpl::for_each<K>(
                type_iterator([&](size_t p, std::string tname) {
                    src.new_line() << tname << " key" << p << ";";
                    })
                );
        src.new_line() << type_name<V>() << " work_sum;";

        src.new_line() << "if (map_id < n)";
        src.open("{");

        boost::mpl::for_each<K>(
                type_iterator([&](size_t p, std::string tname) {
                    src.new_line() << tname << " prev_key" << p << ";";
                    })
                );

        // accumulate zeroth value manually
        src.new_line() << "offset = 0;";
        for(size_t p = 0; p < nK; ++p)
            src.new_line() << "key" << p << " = key_sum" << p << "[map_id];";
        src.new_line() << "work_sum = pre_sum[map_id];";

        // serial accumulation
        src.new_line() << "for(offset = 1; offset < work_per_thread; ++offset)";
        src.open("{");

        for(size_t p = 0; p < nK; ++p) {
            src.new_line() << "prev_key" << p << " = key" << p << ";";
            src.new_line() << "key" << p << " = " << "key_sum" << p << "[map_id + offset];";
        }

        src.new_line() << "if (map_id + offset < n)";
        src.open("{");

        src.new_line() << "if (comp(key0";
        for(size_t p = 1; p < nK; ++p) src << ", key" << p;
        for(size_t p = 0; p < nK; ++p) src << ", prev_key" << p;
        src << ")) work_sum = oper(work_sum, pre_sum[map_id + offset]);";
        src.new_line() << "else work_sum =  pre_sum[map_id + offset];";

        src.new_line() << "pre_sum[map_id + offset] = work_sum;";

        src.close("}");
        src.close("}");
        src.close("}");

        src.new_line().barrier();

        src.new_line() << type_name<V>() << " scan_sum = work_sum;";
        src.new_line() << "shared.vals[l_id] = work_sum;";
        for(size_t p = 0; p < nK; ++p)
            src.new_line() << "shared.keys" << p << "[l_id] = key" << p << ";";

        src.new_line() << "for(offset = 1; offset < " << NT << "; offset *= 2)";
        src.open("{");

        src.new_line().barrier();
        src.new_line() << "if (map_id < n)";
        src.open("{");

        src.new_line() << "if (l_id >= offset)";
        src.open("{");

        boost::mpl::for_each<K>(
                type_iterator([&](size_t p, std::string tname) {
                    src.new_line()
                        << tname << " key1" << p << " = shared.keys" << p
                        << "[l_id];";
                    src.new_line()
                        << tname << " key2" << p << " = shared.keys" << p
                        << "[l_id - offset];";
                    })
                );

        src.new_line() << "if (comp(key10";
        for(size_t p = 1; p < nK; ++p) src << ", key1" << p;
        for(size_t p = 0; p < nK; ++p) src << ", key2" << p;
        src << ")) scan_sum = oper(scan_sum, shared.vals[l_id - offset]);";
        src.new_line() << "else scan_sum = shared.vals[l_id];";

        src.close("}");
        src.close("}");

        src.new_line().barrier();
        src.new_line() << "shared.vals[l_id] = scan_sum;";

        src.close("}");
        src.new_line().barrier();

        // write final scan from pre-scan and shared scan
        src.new_line() << "for(offset = 0; offset < work_per_thread; ++offset)";
        src.open("{");

        src.new_line().barrier(true);

        src.new_line() << "if (map_id < n && l_id > 0)";
        src.open("{");

        src.new_line() << type_name<V>() << " y = pre_sum[map_id + offset];";

        boost::mpl::for_each<K>(
                type_iterator([&](size_t p, std::string tname) {
                    src.new_line()
                        << tname << " key1" << p << " = key_sum" << p
                        << "[map_id + offset];";
                    src.new_line()
                        << tname << " key2" << p << " = shared.keys" << p
                        << "[l_id - 1];";
                    })
                );

        src.new_line() << "if (comp(key10";
        for(size_t p = 1; p < nK; ++p) src << ", key1" << p;
        for(size_t p = 0; p < nK; ++p) src << ", key2" << p;
        src << ")) y = oper(y, shared.vals[l_id - 1]);";
        src.new_line() << "pre_sum[map_id + offset] = y;";
        src.close("}");
        src.close("}");
        src.end_kernel();

        kernel = cache.insert(queue, backend::kernel(
                    queue, src.str(), "block_inclusive_scan_by_key"));
    }

    return kernel->second;
}

//---------------------------------------------------------------------------
template <int NT, typename K, typename V, class Comp, class Oper, bool exclusive>
backend::kernel block_add_by_key(const backend::command_queue &queue) {
    static detail::kernel_cache cache;

    auto kernel = cache.find(queue);
    if (kernel == cache.end()) {
        backend::source_generator src(queue);

        Comp::define(src, "comp");
        Oper::define(src, "oper");

        src.begin_kernel("block_add_by_key");
        src.begin_kernel_parameters();
        src.template parameter<size_t>("n");
        src.template parameter< global_ptr<const V> >("pre_sum");
        src.template parameter< global_ptr<const V> >("pre_sum1");
        src.template parameter< global_ptr<const V> >("ivals");
        src.template parameter< global_ptr<      V> >("ovals");

        boost::mpl::for_each<K>(pointer_param<global_ptr, true>(src, "ikeys"));

        if (exclusive) src.template parameter<V>("init");

        src.end_kernel_parameters();

        src.new_line() << "size_t g_id   = " << src.global_id(0) << ";";
        src.new_line() << "size_t l_id   = " << src.local_id(0)  << ";";
        src.new_line() << "size_t block  = " << src.group_id(0)  << ";";

        src.new_line() << "struct Shared";
        src.open("{");
            src.new_line() << type_name<V>() << " vals[" << NT << "];";
            boost::mpl::for_each<K>(
                    type_iterator([&](size_t p, std::string tname) {
                        src.new_line()
                            << tname << " keys" << p << "[" << NT << "];";
                        })
                    );
        src.close("};");
        src.smem_static_var("struct Shared", "shared");

        const size_t nK = boost::mpl::size<K>::value;

        // if exclusive, load gloId=0 w/ init, and all others shifted-1
        src.new_line() << type_name<V>() << " val;";
        boost::mpl::for_each<K>(
                type_iterator([&](size_t p, std::string tname) {
                    src.new_line() << tname << " key" << p << ";";
                    })
                );

        src.new_line() << "if (g_id < n)";
        src.open("{");

        if (exclusive) {
            src.new_line() << "if (g_id > 0)";
            src.open("{");
            boost::mpl::for_each<K>(
                    type_iterator([&](size_t p, std::string tname) {
                        src.new_line() << tname << " key1" << p << " = key" << p << " = ikeys" << p << "[g_id];";
                        src.new_line() << tname << " key2" << p << " = ikeys" << p << "[g_id-1];";
                        })
                    );

            src.new_line() << "if (comp(key10";
            for(size_t p = 1; p < nK; ++p) src << ", key1" << p;
            for(size_t p = 0; p < nK; ++p) src << ", key2" << p;
            src << ")) val = ivals[g_id - 1];";
            src.new_line() << "else val = init;";

            src.new_line() << "shared.vals[l_id] = val;";
            for(size_t p = 0; p < nK; ++p)
                src.new_line() << "shared.keys" << p << "[l_id] = key" << p << ";";

            src.close("}");
            src.new_line() << "else";
            src.open("{");

            src.new_line() << "val = init;";
            src.new_line() << "shared.vals[l_id] = val;";
            for(size_t p = 0; p < nK; ++p)
                src.new_line() << "shared.keys" << p << "[l_id] = ikeys" << p << "[g_id];";

            src.close("}");
        } else {
            src.new_line() << "shared.vals[l_id] =val = ivals[g_id];";
            for(size_t p = 0; p < nK; ++p)
                src.new_line() << "shared.keys" << p << "[l_id] = key" << p << " = ikeys" << p << "[g_id];";
        }

        src.close("}");

        // Each work item writes out its calculated scan result, relative to
        // the beginning of each work group
        src.new_line() << type_name<V>() << " scan_result = shared.vals[l_id];";
        src.new_line() << type_name<V>() << " post_sum, new_result, sum;";

        boost::mpl::for_each<K>(
                type_iterator([&](size_t p, std::string tname) {
                    src.new_line() << tname
                        << " key1" << p << ", "
                        << " key2" << p << ", "
                        << " key3" << p << ", "
                        << " key4" << p << ";";
                    })
                );

        src.new_line() << "if (l_id == 0 && g_id < n)";
        src.open("{");

        src.new_line() << "if (block > 0)";
        src.open("{");

        for(size_t p = 0; p < nK; ++p) {
            src.new_line() << "key1" << p << " = ikeys" << p << "[g_id];";
            src.new_line() << "key2" << p << " = ikeys" << p << "[block * "<< NT << " - 1];";
        }

        src.new_line() << "if (block % 2 == 0) post_sum = pre_sum[block / 2 - 1];";
        src.new_line() << "else if (block == 1) post_sum = pre_sum1[0];";
        src.new_line() << "else";
        src.open("{");

        for(size_t p = 0; p < nK; ++p) {
            src.new_line() << "key3" << p << " = ikeys" << p << "[block * " << NT << " - 1];";
            src.new_line() << "key4" << p << " = ikeys" << p << "[(block - 1) * " << NT << " - 1];";
        }

        src.new_line() << "if (comp(key30";
        for(size_t p = 1; p < nK; ++p) src << ", key3" << p;
        for(size_t p = 0; p < nK; ++p) src << ", key4" << p;
        src << ")) post_sum = oper(pre_sum[block / 2 - 1], pre_sum1[block / 2]);";
        src.new_line() << "else post_sum = pre_sum1[block / 2];";

        src.close("}");

        if (exclusive) {
            src.new_line() << "if (comp(key10";
            for(size_t p = 1; p < nK; ++p) src << ", key1" << p;
            for(size_t p = 0; p < nK; ++p) src << ", key2" << p;
            src << ")) new_result = post_sum;";
            src.new_line() << "else new_result = init;";
        } else {
            src.new_line() << "if (comp(key10";
            for(size_t p = 1; p < nK; ++p) src << ", key1" << p;
            for(size_t p = 0; p < nK; ++p) src << ", key2" << p;
            src << ")) new_result = oper(scan_result, post_sum);";
            src.new_line() << "else new_result = scan_result;";
        }

        src.close("}");

        src.new_line() << "else new_result = scan_result;";
        src.new_line() << "shared.vals[l_id] = new_result;";

        src.close("}");

        // Computes a scan within a workgroup,
        // updates vals in shared but not keys
        src.new_line() << "sum = shared.vals[l_id];";
        src.new_line() << "for(size_t offset = 1; offset < " << NT << "; offset *= 2)";
        src.open("{");

        src.new_line().barrier();

        src.new_line() << "if (l_id >= offset)";
        src.open("{");

        for(size_t p = 0; p < nK; ++p)
            src.new_line() << "key2" << p << " = shared.keys" << p << "[l_id - offset];";

        src.new_line() << "if (comp(key0";
        for(size_t p = 1; p < nK; ++p) src << ", key" << p;
        for(size_t p = 0; p < nK; ++p) src << ", key2" << p;
        src << ")) sum = oper(sum, shared.vals[l_id - offset]);";

        src.close("}");

        src.new_line().barrier();
        src.new_line() << "shared.vals[l_id] = sum;";

        src.close("}");
        src.new_line().barrier();

        src.new_line() << "if (g_id < n) ovals[g_id] = sum;";

        src.end_kernel();

        kernel = cache.insert(queue, backend::kernel(
                    queue, src.str(), "block_add_by_key"));
    }

    return kernel->second;
}

template <bool exclusive, class KTuple, class V, class Comp, class Oper>
void scan_by_key(
        KTuple &&keys, const vector<V> &ivals, vector<V> &ovals, Comp, Oper, V init
        )
{
    namespace fusion = boost::fusion;
    typedef typename extract_value_types<KTuple>::type K;

    precondition(
            fusion::at_c<0>(keys).nparts() == 1 && ivals.nparts() == 1,
            "scan_by_key is only supported for single device contexts"
            );

    precondition(ivals.size() == ovals.size() && ivals.nparts() == ovals.nparts(),
            "input and output should have same size"
            );


    const auto &queue = fusion::at_c<0>(keys).queue_list()[0];
    backend::select_context(queue);

    const int NT_cpu = 1;
    const int NT_gpu = 256;
    const int NT = is_cpu(queue) ? NT_cpu : NT_gpu;

    size_t count         = fusion::at_c<0>(keys).size();
    size_t num_blocks    = (count + 2 * NT - 1) / (2 * NT);
    size_t scan_buf_size = alignup(num_blocks, NT);

    auto ikeys = fusion::transform(keys, extract_device_vector(0));

    temp_storage<K>           key_sum (queue, scan_buf_size);
    backend::device_vector<V> pre_sum (queue, scan_buf_size);
    backend::device_vector<V> pre_sum1(queue, scan_buf_size);

    /***** Kernel 0 *****/
    auto krn0 = is_cpu(queue) ?
        block_scan_by_key<NT_cpu, K, V, Comp, Oper, exclusive>(queue) :
        block_scan_by_key<NT_gpu, K, V, Comp, Oper, exclusive>(queue);

    krn0.push_arg(count);
    krn0.push_arg(ivals(0));
    krn0.push_arg(pre_sum);
    krn0.push_arg(pre_sum1);

    push_args<boost::mpl::size<K>::value>(krn0, ikeys);
    push_args<boost::mpl::size<K>::value>(krn0, key_sum);

    if (exclusive) krn0.push_arg(init);

    krn0.config(num_blocks, NT);
    krn0(queue);

    /***** Kernel 1 *****/
    auto krn1 = is_cpu(queue) ?
        block_inclusive_scan_by_key<NT_cpu, K, V, Comp, Oper>(queue) :
        block_inclusive_scan_by_key<NT_gpu, K, V, Comp, Oper>(queue);

    uint work_per_thread = std::max<uint>(1U, static_cast<uint>(scan_buf_size / NT));

    krn1.push_arg(scan_buf_size);
    krn1.push_arg(pre_sum);
    krn1.push_arg(work_per_thread);

    push_args<boost::mpl::size<K>::value>(krn1, key_sum);

    krn1.config(1, NT);
    krn1(queue);

    /***** Kernel 2 *****/
    auto krn2 = is_cpu(queue) ?
        block_add_by_key<NT_cpu, K, V, Comp, Oper, exclusive>(queue) :
        block_add_by_key<NT_gpu, K, V, Comp, Oper, exclusive>(queue);

    krn2.push_arg(count);
    krn2.push_arg(pre_sum);
    krn2.push_arg(pre_sum1);
    krn2.push_arg(ivals(0));
    krn2.push_arg(ovals(0));

    push_args<boost::mpl::size<K>::value>(krn2, ikeys);

    if (exclusive) krn2.push_arg(init);

    krn2.config(num_blocks * 2, NT);
    krn2(queue);
}

} // namespace sbk
} // namespace detail

/// Exclusive scan by key.
template<class K, typename V, class Comp, class Oper>
void exclusive_scan_by_key(
        K &&keys, const vector<V> &ivals, vector<V> &ovals,
        Comp comp, Oper oper, V init = V()
        )
{
    return detail::sbk::scan_by_key<true>(
            detail::forward_as_sequence(keys), ivals, ovals, comp, oper, init
            );
}

/// Inclusive scan by key.
template<class K, typename V, class Comp, class Oper>
void inclusive_scan_by_key(
        K &&keys, const vector<V> &ivals, vector<V> &ovals,
        Comp comp, Oper oper, V init = V()
        )
{
    return detail::sbk::scan_by_key<false>(
            detail::forward_as_sequence(keys), ivals, ovals, comp, oper, init
            );
}

/// Exclusive scan by key.
template<typename K, typename V>
void exclusive_scan_by_key(
        const vector<K> &keys, const vector<V> &ivals, vector<V> &ovals,
        V init = V()
        )
{
    VEX_FUNCTION(bool, equal, (K, x)(K, y), return x == y;);
    VEX_FUNCTION(V, plus, (V, x)(V, y), return x + y;);
    return exclusive_scan_by_key(keys, ivals, ovals, equal, plus, init);
}

/// Inclusive scan by key.
template<typename K, typename V>
void inclusive_scan_by_key(
        const vector<K> &keys, const vector<V> &ivals, vector<V> &ovals,
        V init = V()
        )
{
    VEX_FUNCTION(bool, equal, (K, x)(K, y), return x == y;);
    VEX_FUNCTION(V, plus, (V, x)(V, y), return x + y;);
    return inclusive_scan_by_key(keys, ivals, ovals, equal, plus, init);
}

} // namespace vex

#endif
