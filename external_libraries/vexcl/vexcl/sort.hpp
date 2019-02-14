#ifndef VEXCL_SORT_HPP
#define VEXCL_SORT_HPP

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
 * \file   vexcl/sort.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Sorting algorithms.

Adopted from NVIDIA Modern GPU patterns,
see <http://nvlabs.github.io/moderngpu>.
The original code came with the following copyright notice:

\verbatim
Copyright (c) 2013, NVIDIA CORPORATION.  All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the NVIDIA CORPORATION nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL NVIDIA CORPORATION BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
\endverbatim
*/

#include <string>
#include <functional>

#include <vexcl/backend.hpp>
#include <vexcl/util.hpp>
#include <vexcl/vector.hpp>
#include <vexcl/detail/fusion.hpp>
#include <vexcl/function.hpp>

#ifndef VEX_SORT_NT_GPU
#  define VEX_SORT_NT_GPU 256
#endif

namespace vex {
namespace detail {

//---------------------------------------------------------------------------
// Memory transfer functions
//---------------------------------------------------------------------------
template<int NT, int VT, typename T>
std::string global_to_regstr_pred() {
    std::ostringstream s;
    s << "global_to_regstr_pred_" << NT << "_" << VT << "_" << type_name<T>();
    return s.str();
}

template<int NT, int VT, typename T>
void global_to_regstr_pred(backend::source_generator &src) {
    src.begin_function<void>( global_to_regstr_pred<NT,VT,T>() );
    src.begin_function_parameters();
    src.template parameter< int                 >("count");
    src.template parameter< global_ptr<const T> >("data");
    src.template parameter< int                 >("tid");
    src.template parameter< regstr_ptr<T>       >("reg");
    src.end_function_parameters();

    src.new_line() << type_name<int>() << " index;";

    for(int i = 0; i < VT; ++i) {
        src.new_line() << "index = " << NT * i << " + tid;";
        src.new_line() << "if (index < count) reg[" << i << "] = data[index];";
    }
    src.end_function();
}

//---------------------------------------------------------------------------
template<int NT, int VT, typename T>
std::string global_to_regstr() {
    std::ostringstream s;
    s << "global_to_regstr_" << NT << "_" << VT << "_" << type_name<T>();
    return s.str();
}

template<int NT, int VT, typename T>
void global_to_regstr(backend::source_generator &src) {
    global_to_regstr_pred<NT, VT, T>(src);

    src.begin_function<void>( global_to_regstr<NT,VT,T>() );
    src.begin_function_parameters();
    src.template parameter< int                 >("count");
    src.template parameter< global_ptr<const T> >("data");
    src.template parameter< int                 >("tid");
    src.template parameter< regstr_ptr<T>       >("reg");
    src.end_function_parameters();

    src.new_line() << "if (count >= " << NT * VT << ")";
    src.open("{");

    for(int i = 0; i < VT; ++i)
        src.new_line() << "reg[" << i << "] = data[" << NT * i << " + tid];";

    src.close("}") << " else "
        << global_to_regstr_pred<NT, VT, T>() << "(count, data, tid, reg);";

    src.end_function();
}

//---------------------------------------------------------------------------
template<int NT, int VT, typename T>
std::string regstr_to_global() {
    std::ostringstream s;
    s << "regstr_to_global_" << NT << "_" << VT << "_" << type_name<T>();
    return s.str();
}

template<int NT, int VT, typename T>
void regstr_to_global(backend::source_generator &src) {
    src.begin_function<void>( regstr_to_global<NT,VT,T>() );
    src.begin_function_parameters();
    src.template parameter< int                 >("count");
    src.template parameter< regstr_ptr<const T> >("reg");
    src.template parameter< int                 >("tid");
    src.template parameter< global_ptr<T>       >("dest");
    src.end_function_parameters();

    src.new_line() << type_name<int>() << " index;";

    for(int i = 0; i < VT; ++i) {
        src.new_line() << "index = " << NT * i << " + tid;";
        src.new_line() << "if (index < count) dest[index] = reg[" << i << "];";
    }

    src.new_line().barrier();

    src.end_function();
}

//---------------------------------------------------------------------------
template<int NT, int VT, typename T>
std::string shared_to_regstr() {
    std::ostringstream s;
    s << "shared_to_regstr_" << NT << "_" << VT << "_" << type_name<T>();
    return s.str();
}

template<int NT, int VT, typename T>
void shared_to_regstr(backend::source_generator &src) {
    src.begin_function<void>( shared_to_regstr<NT,VT,T>() );
    src.begin_function_parameters();
    src.template parameter< shared_ptr<const T> >("data");
    src.template parameter< int                 >("tid");
    src.template parameter< regstr_ptr<T>       >("reg");
    src.end_function_parameters();

    for(int i = 0; i < VT; ++i)
        src.new_line() << "reg[" << i << "] = data[" << NT * i << " + tid];";

    src.new_line().barrier();

    src.end_function();
}

//---------------------------------------------------------------------------
template<int NT, int VT, typename T>
std::string regstr_to_shared() {
    std::ostringstream s;
    s << "regstr_to_shared_" << NT << "_" << VT << "_" << type_name<T>();
    return s.str();
}

template<int NT, int VT, typename T>
void regstr_to_shared(backend::source_generator &src) {
    src.begin_function<void>( regstr_to_shared<NT,VT,T>() );
    src.begin_function_parameters();
    src.template parameter< regstr_ptr<const T> >("reg");
    src.template parameter< int                 >("tid");
    src.template parameter< shared_ptr<T>       >("dest");
    src.end_function_parameters();

    for(int i = 0; i < VT; ++i)
        src.new_line() << "dest[" << NT * i << " + tid] = reg[" << i << "];";

    src.new_line().barrier();

    src.end_function();
}

//---------------------------------------------------------------------------
template<int NT, int VT, typename T>
std::string global_to_shared() {
    std::ostringstream s;
    s << "global_to_shared_" << NT << "_" << VT << "_" << type_name<T>();
    return s.str();
}

template<int NT, int VT, typename T>
void global_to_shared(backend::source_generator &src) {
    src.begin_function<void>( global_to_shared<NT,VT,T>() );
    src.begin_function_parameters();
    src.template parameter< int                 >("count");
    src.template parameter< global_ptr<const T> >("source");
    src.template parameter< int                 >("tid");
    src.template parameter< shared_ptr<T>       >("dest");
    src.end_function_parameters();

    src.new_line() << type_name<T>() << " reg[" << VT << "];";
    src.new_line() << global_to_regstr<NT, VT, T>() << "(count, source, tid, reg);";
    src.new_line() << regstr_to_shared<NT, VT, T>() << "(reg, tid, dest);";

    src.end_function();
}

//---------------------------------------------------------------------------
template<int NT, int VT, typename T>
std::string shared_to_global() {
    std::ostringstream s;
    s << "shared_to_global_" << NT << "_" << VT << "_" << type_name<T>();
    return s.str();
}

template<int NT, int VT, typename T>
void shared_to_global(backend::source_generator &src) {
    src.begin_function<void>( shared_to_global<NT,VT,T>() );
    src.begin_function_parameters();
    src.template parameter< int                 >("count");
    src.template parameter< shared_ptr<const T> >("source");
    src.template parameter< int                 >("tid");
    src.template parameter< global_ptr<T>       >("dest");
    src.end_function_parameters();

    src.new_line() << type_name<int>() << " index;";

    for(int i = 0; i < VT; ++i) {
        src.new_line() << "index = " << NT * i << " + tid;";
        src.new_line() << "if (index < count) dest[index] = source[index];";
    }

    src.new_line().barrier();
    src.end_function();
}

//---------------------------------------------------------------------------
template<int VT, typename T>
std::string shared_to_thread() {
    std::ostringstream s;
    s << "shared_to_thread_" << VT << "_" << type_name<T>();
    return s.str();
}

template<int VT, typename T>
void shared_to_thread(backend::source_generator &src) {
    src.begin_function<void>( shared_to_thread<VT,T>() );
    src.begin_function_parameters();
    src.template parameter< shared_ptr<const T> >("data");
    src.template parameter< int                 >("tid");
    src.template parameter< regstr_ptr<T>       >("reg");
    src.end_function_parameters();

    for(int i = 0; i < VT; ++i)
        src.new_line() << "reg[" << i << "] = data[" << VT << " * tid + " << i << "];";

    src.new_line().barrier();

    src.end_function();
}

//---------------------------------------------------------------------------
template<int VT, typename T>
std::string thread_to_shared() {
    std::ostringstream s;
    s << "thread_to_shared_" << VT << "_" << type_name<T>();
    return s.str();
}

template<int VT, typename T>
void thread_to_shared(backend::source_generator &src) {
    src.begin_function<void>( thread_to_shared<VT,T>() );
    src.begin_function_parameters();
    src.template parameter< regstr_ptr<const T> >("reg");
    src.template parameter< int                 >("tid");
    src.template parameter< shared_ptr<T>       >("dest");
    src.end_function_parameters();

    for(int i = 0; i < VT; ++i)
        src.new_line() << "dest[" << VT << " * tid + " << i << "] = reg[" << i << "];";

    src.new_line().barrier();

    src.end_function();
}

//---------------------------------------------------------------------------
template<int NT, int VT, typename T>
void transfer_functions(backend::source_generator &src) {
    global_to_regstr<NT, VT, T>(src);
    regstr_to_global<NT, VT, T>(src);

    shared_to_regstr<NT, VT, T>(src);
    regstr_to_shared<NT, VT, T>(src);

    global_to_shared<NT, VT, T>(src);
    shared_to_global<NT, VT, T>(src);

    shared_to_thread<VT, T>(src);
    thread_to_shared<VT, T>(src);
}

//---------------------------------------------------------------------------
// Block merge sort kernel
//---------------------------------------------------------------------------
template <typename T>
std::string swap_function() {
    std::ostringstream s;
    s << "swap";
    print_types<T>(s);
    return s.str();
}

template <typename T>
void swap_function(backend::source_generator &src) {
    src.begin_function<void>(swap_function<T>());
    src.begin_function_parameters();

    boost::mpl::for_each<T>( pointer_param<regstr_ptr>(src, "a") );
    boost::mpl::for_each<T>( pointer_param<regstr_ptr>(src, "b") );

    src.end_function_parameters();

    boost::mpl::for_each<T>( type_iterator([&](size_t pos, std::string tname) {
                src.open("{");
                src.new_line() << tname << " c = *a" << pos << ";";
                src.new_line() << "*a" << pos << " = *b" << pos << ";";
                src.new_line() << "*b" << pos << " = c;";
                src.close("}");
                }) );

    src.end_function();
}

//---------------------------------------------------------------------------
template<int VT, typename K, typename V>
std::string odd_even_transpose_sort() {
    std::ostringstream s;
    s << "odd_even_transpose_sort_" << VT;
    print_types<K>(s);
    print_types<V>(s);
    return s.str();
}

template<int VT, typename K, typename V>
void odd_even_transpose_sort(backend::source_generator &src) {
    swap_function<K>(src);
    if (boost::mpl::size<V>::value && !std::is_same<K, V>::value)
        swap_function<V>(src);

    src.begin_function<void>(odd_even_transpose_sort<VT,K,V>());
    src.begin_function_parameters();
    boost::mpl::for_each<K>( pointer_param<regstr_ptr>(src, "keys") );
    boost::mpl::for_each<V>( pointer_param<regstr_ptr>(src, "vals") );
    src.end_function_parameters();

    for(int I = 0; I < VT; ++I) {
        for(int i = 1 & I; i < VT - 1; i += 2) {
            src.new_line() << "if (comp(";
            for(int p = 0; p < boost::mpl::size<K>::value; ++p)
                src << (p ? ", " : "") << "keys" << p << "[" << i + 1 << "]";
            for(int p = 0; p < boost::mpl::size<K>::value; ++p)
                src << ", keys" << p << "[" << i << "]";
            src << "))";
            src.open("{");
            src.new_line() << swap_function<K>() << "(";
            for(int p = 0; p < boost::mpl::size<K>::value; ++p)
                src << (p ? ", " : "") << "keys" << p << " + " << i;
            for(int p = 0; p < boost::mpl::size<K>::value; ++p)
                src << ", keys" << p << " + " << i + 1;
            src << ");";
            if (boost::mpl::size<V>::value) {
                src.new_line() << swap_function<V>() << "(";
                for(int p = 0; p < boost::mpl::size<V>::value; ++p)
                    src << (p ? ", " : "") << "vals" << p << " + " << i;
                for(int p = 0; p < boost::mpl::size<V>::value; ++p)
                    src << ", vals" << p << " + " << i + 1;
                src << ");";
            }
            src.close("}");
        }
    }

    src.end_function();
}

//---------------------------------------------------------------------------
template<template<class> class Address, typename T>
std::string merge_path() {
    std::ostringstream s;
    s << "merge_path";
    print_types<T>(s);
    return s.str();
}

template<template<class> class Address, typename T>
void merge_path(backend::source_generator &src) {
    src.begin_function<int>(merge_path<Address, T>());
    src.begin_function_parameters();
    src.template parameter< int >("a_count");
    src.template parameter< int >("b_count");
    src.template parameter< int >("diag");

    boost::mpl::for_each<T>( pointer_param<Address, true>(src, "a") );
    boost::mpl::for_each<T>( pointer_param<Address, true>(src, "b") );

    src.end_function_parameters();

    src.new_line() << "int begin = max(0, diag - b_count);";
    src.new_line() << "int end   = min(diag, a_count);";

    src.new_line() << "while (begin < end)";
    src.open("{");
    src.new_line() << "int mid = (begin + end) >> 1;";
    src.new_line() << "if ( !comp(";
    for(int p = 0; p < boost::mpl::size<T>::value; ++p)
        src << (p ? ", " : "") << "b" << p << "[diag - 1 - mid]";
    for(int p = 0; p < boost::mpl::size<T>::value; ++p)
        src << ", a" << p << "[mid]";
    src << ") ) begin = mid + 1;";
    src.new_line() << "else end = mid;";
    src.close("}");

    src.new_line() << "return begin;";

    src.end_function();
}

//---------------------------------------------------------------------------
template<int VT, typename T>
std::string serial_merge() {
    std::ostringstream s;
    s << "serial_merge_" << VT;
    print_types<T>(s);
    return s.str();
}

template<int VT, typename T>
void serial_merge(backend::source_generator &src) {
    src.begin_function<void>(serial_merge<VT, T>());
    src.begin_function_parameters();
    src.template parameter< int                 >("a_begin");
    src.template parameter< int                 >("a_end");
    src.template parameter< int                 >("b_begin");
    src.template parameter< int                 >("b_end");
    src.template parameter< regstr_ptr<int>     >("indices");

    boost::mpl::for_each<T>( pointer_param<shared_ptr, true>(src, "keys_shared") );
    boost::mpl::for_each<T>( pointer_param<regstr_ptr>(src, "results") );

    src.end_function_parameters();

    boost::mpl::for_each<T>( type_iterator([&](size_t pos, std::string tname) {
                src.new_line() << tname << " a_key" << pos << " = keys_shared" << pos << "[a_begin];";
                src.new_line() << tname << " b_key" << pos << " = keys_shared" << pos << "[b_begin];";
                }) );

    src.new_line() << "bool p;";

    for(int i = 0; i < VT; ++i) {
        src.new_line() << "p = (b_begin >= b_end) || ((a_begin < a_end) && !comp(";
        for(int p = 0; p < boost::mpl::size<T>::value; ++p)
            src << (p ? ", " : "") << "b_key" << p;
        for(int p = 0; p < boost::mpl::size<T>::value; ++p)
            src << ", a_key" << p;
        src << "));";

        for(int p = 0; p < boost::mpl::size<T>::value; ++p)
            src.new_line()
                << "results" << p << "[" << i << "] = "
                << "p ? a_key" << p << " : b_key" << p << ";";

        src.new_line() << "indices[" << i << "] = p ? a_begin : b_begin;";

        src.new_line() << "if(p)";
        src.open("{");
        src.new_line() << "++a_begin;";
        for(int p = 0; p < boost::mpl::size<T>::value; ++p)
            src.new_line() << "a_key" << p << " = keys_shared" << p << "[a_begin];";
        src.close("}");
        src.new_line() << "else";
        src.open("{");
        src.new_line() << "++b_begin;";
        for(int p = 0; p < boost::mpl::size<T>::value; ++p)
            src.new_line() << "b_key" << p << " = keys_shared" << p << "[b_begin];";
        src.close("}");
    }

    src.new_line().barrier();

    src.end_function();
}

//---------------------------------------------------------------------------
template<int NT, int VT, typename T>
std::string block_sort_pass() {
    std::ostringstream s;
    s << "block_sort_pass_" << NT << "_" << VT;
    print_types<T>(s);
    return s.str();
}

template<int NT, int VT, typename T>
void block_sort_pass(backend::source_generator &src) {
    merge_path< shared_ptr, T >(src);

    src.begin_function<void>(block_sort_pass<NT, VT, T>());
    src.begin_function_parameters();
    src.template parameter< int                 >("tid");
    src.template parameter< int                 >("count");
    src.template parameter< int                 >("coop");
    src.template parameter< regstr_ptr<int>     >("indices");

    boost::mpl::for_each<T>( pointer_param<shared_ptr, true>(src, "keys_shared") );
    boost::mpl::for_each<T>( pointer_param<regstr_ptr>(src, "keys") );

    src.end_function_parameters();

    src.new_line() << "int list = ~(coop - 1) & tid;";
    src.new_line() << "int diag = min(count, " << VT << " * ((coop - 1) & tid));";
    src.new_line() << "int start = " << VT << " * list;";
    src.new_line() << "int a0 = min(count, start);";
    src.new_line() << "int b0 = min(count, start + " << VT << " * (coop / 2));";
    src.new_line() << "int b1 = min(count, start + " << VT << " * coop);";

    src.new_line() << "int p = "
        << merge_path< shared_ptr, T >()
        << "(b0 - a0, b1 - b0, diag";
    for(int p = 0; p < boost::mpl::size<T>::value; ++p)
        src << ", keys_shared" << p << " + a0";
    for(int p = 0; p < boost::mpl::size<T>::value; ++p)
        src << ", keys_shared" << p << " + b0";
    src << ");";
    src.new_line() << serial_merge<VT, T >() << "(a0 + p, b0, b0 + diag - p, b1, indices";
    for(int p = 0; p < boost::mpl::size<T>::value; ++p)
        src << ", keys_shared" << p;
    for(int p = 0; p < boost::mpl::size<T>::value; ++p)
        src << ", keys" << p;
    src << ");";

    src.end_function();
}

//---------------------------------------------------------------------------
template<int NT, int VT, typename T>
std::string gather() {
    std::ostringstream s;
    s << "gather_" << NT << "_" << VT;
    print_types<T>(s);
    return s.str();
}

template<int NT, int VT, typename T>
void gather(backend::source_generator &src) {
    src.begin_function<void>(gather<NT, VT, T>());
    src.begin_function_parameters();
    src.template parameter< regstr_ptr<const int> >("indices");
    src.template parameter< int                   >("tid");

    boost::mpl::for_each<T>( pointer_param<shared_ptr, true>(src, "data") );
    boost::mpl::for_each<T>( pointer_param<regstr_ptr>(src, "reg") );
    src.end_function_parameters();

    for(int i = 0; i < VT; ++i)
        for(int p = 0; p < boost::mpl::size<T>::value; ++p)
            src.new_line() << "reg" << p << "[" << i << "] = data" << p << "[indices[" << i << "]];";

    src.new_line().barrier();

    src.end_function();
}

//---------------------------------------------------------------------------
template<int NT, int VT, typename K, typename V>
std::string block_sort_loop() {
    std::ostringstream s;
    s << "block_sort_loop_" << NT << "_" << VT;
    print_types<K>(s);
    print_types<V>(s);
    return s.str();
}

template<int VT>
struct call_thread_to_shared {
    backend::source_generator &src;
    const char *tname;
    const char *sname;
    int pos;

    call_thread_to_shared(backend::source_generator &src,
            const char *tname, const char *sname
            ) : src(src), tname(tname), sname(sname), pos(0) {}

    template <typename T>
    void operator()(T) {
        src.new_line() << thread_to_shared<VT, T>()
            << "(" << tname << pos << ", tid, " << sname << pos << ");";
        ++pos;
    }
};

template<int NT, int VT, typename K, typename V>
void block_sort_loop(backend::source_generator &src) {
    block_sort_pass<NT, VT, K >(src);
    if (boost::mpl::size<V>::value) gather<NT, VT, V>(src);

    src.begin_function<void>(block_sort_loop<NT, VT, K, V>());
    src.begin_function_parameters();

    src.template parameter< int >("tid");
    src.template parameter< int >("count");

    boost::mpl::for_each<K>( pointer_param<shared_ptr>(src, "keys_shared") );
    boost::mpl::for_each<V>( pointer_param<regstr_ptr>(src, "thread_vals") );
    boost::mpl::for_each<V>( pointer_param<shared_ptr>(src, "vals_shared") );

    src.end_function_parameters();

    src.new_line() << "int indices[" << VT << "];";

    boost::mpl::for_each<K>( type_iterator([&](size_t pos, std::string tname) {
                src.new_line() << tname << " keys" << pos++ << "[" << VT << "];";
                }) );

    for(int coop = 2; coop <= NT; coop *= 2) {
        src.new_line() << block_sort_pass<NT, VT, K >()
            << "(tid, count, " << coop << ", indices";
        for(int p = 0; p < boost::mpl::size<K>::value; ++p)
            src << ", keys_shared" << p;
        for(int p = 0; p < boost::mpl::size<K>::value; ++p)
            src << ", keys" << p;
        src << ");";

        if (boost::mpl::size<V>::value) {
            // Exchange the values through shared memory.
            boost::mpl::for_each<V>( call_thread_to_shared<VT>(src, "thread_vals", "vals_shared") );

            src.new_line() << gather<NT, VT, V>()
                << "(indices, tid";
            for(int p = 0; p < boost::mpl::size<V>::value; ++p)
                src << ", vals_shared" << p;
            for(int p = 0; p < boost::mpl::size<V>::value; ++p)
                src << ", thread_vals" << p;
            src << ");";
        }

        // Store results in shared memory in sorted order.
        boost::mpl::for_each<K>( call_thread_to_shared<VT>(src, "keys", "keys_shared") );
    }

    src.end_function();
}

//---------------------------------------------------------------------------
template<int NT, int VT, typename K, typename V>
std::string mergesort() {
    std::ostringstream s;
    s << "mergesort_" << NT << "_" << VT;
    print_types<K>(s);
    print_types<V>(s);
    return s.str();
}

template<int NT, int VT, typename K, typename V>
void mergesort(backend::source_generator &src) {
    odd_even_transpose_sort<VT, K, V>(src);
    block_sort_loop<NT, VT, K, V>(src);

    src.begin_function<void>(mergesort<NT, VT, K, V>());
    src.begin_function_parameters();

    src.template parameter< int >("count");
    src.template parameter< int >("tid");

    boost::mpl::for_each<K>( pointer_param<regstr_ptr>(src, "thread_keys") );
    boost::mpl::for_each<K>( pointer_param<shared_ptr>(src, "keys_shared") );
    boost::mpl::for_each<V>( pointer_param<regstr_ptr>(src, "thread_vals") );
    boost::mpl::for_each<V>( pointer_param<shared_ptr>(src, "vals_shared") );

    src.end_function_parameters();

    // Stable sort the keys in the thread.
    src.new_line() << "if(" << VT << " * tid < count) "
        << odd_even_transpose_sort<VT, K, V>()
        << "(";
    for(int p = 0; p < boost::mpl::size<K>::value; ++p)
        src << (p ? ", " : "") << "thread_keys" << p;
    for(int p = 0; p < boost::mpl::size<V>::value; ++p)
        src << ", thread_vals" << p;
    src << ");";

    // Store the locally sorted keys into shared memory.
    boost::mpl::for_each<K>( call_thread_to_shared<VT>(src, "thread_keys", "keys_shared") );

    // Recursively merge lists until the entire CTA is sorted.
    src.new_line() << block_sort_loop<NT, VT, K, V>()
        << "(tid, count";
    for(int p = 0; p < boost::mpl::size<K>::value; ++p)
        src << ", keys_shared" << p;
    for(int p = 0; p < boost::mpl::size<V>::value; ++p)
        src << ", thread_vals" << p;
    for(int p = 0; p < boost::mpl::size<V>::value; ++p)
        src << ", vals_shared" << p;
    src << ");";

    src.end_function();
}

template <int NT, int VT>
struct define_transfer_functions {
    backend::source_generator &src;

    define_transfer_functions(backend::source_generator &src) : src(src) {}

    template <typename T>
    void operator()(T) const {
        transfer_functions<NT, VT, T>(src);
    }
};

template<int NT, int VT>
struct call_global_to_shared {
    backend::source_generator &src;
    const char *gname;
    const char *sname;
    int pos;

    call_global_to_shared(backend::source_generator &src,
            const char *gname, const char *sname
            ) : src(src), gname(gname), sname(sname), pos(0) {}

    template <typename T>
    void operator()(T) {
        src.new_line() << global_to_shared<NT, VT, T>() << "(count2, "
            << gname << pos << " + gid, tid, "
            << sname << pos << ");";
        ++pos;
    }
};

template<int NT, int VT>
struct call_shared_to_global {
    backend::source_generator &src;
    const char *count;
    const char *sname;
    const char *gname;
    const char *gid;
    int pos;

    call_shared_to_global(backend::source_generator &src,
            const char *count, const char *sname, const char *gname, const char *gid
            ) : src(src), count(count), sname(sname), gname(gname), gid(gid), pos(0) {}

    template <typename T>
    void operator()(T) {
        src.new_line() << shared_to_global<NT, VT, T>()
            << "(" << count << ", " << sname << pos << ", tid, " << gname << pos << " + " << gid << ");";
        ++pos;
    }
};

template<int VT>
struct call_shared_to_thread {
    backend::source_generator &src;
    const char *sname;
    const char *tname;
    int pos;

    call_shared_to_thread(backend::source_generator &src,
            const char *sname, const char *tname
            ) : src(src), sname(sname), tname(tname), pos(0) {}

    template <typename T>
    void operator()(T) {
        src.new_line() << shared_to_thread<VT, T>()
            << "(" << sname << pos << ", tid, " << tname << pos << ");";
        ++pos;
    }
};

//---------------------------------------------------------------------------
template <int NT, int VT, typename K, typename V, typename Comp>
backend::kernel& block_sort_kernel(const backend::command_queue &queue) {
    static detail::kernel_cache cache;

    auto kernel = cache.find(queue);

    if (kernel == cache.end()) {
        backend::source_generator src(queue);

        Comp::define(src, "comp");

        boost::mpl::for_each<
            typename boost::mpl::copy<
                typename boost::mpl::copy<
                    V,
                    boost::mpl::back_inserter<K>
                    >::type,
                boost::mpl::inserter<
                    boost::mpl::set<int>,
                    boost::mpl::insert<boost::mpl::_1, boost::mpl::_2>
                    >
                >::type
            >( define_transfer_functions<NT, VT>(src) );

        serial_merge<VT, K >(src);
        mergesort<NT, VT, K, V>(src);

        src.begin_kernel("block_sort");
        src.begin_kernel_parameters();
        src.template parameter< int >("count");

        boost::mpl::for_each<K>( pointer_param<global_ptr, true>(src, "keys_src") );
        boost::mpl::for_each<K>( pointer_param<global_ptr      >(src, "keys_dst") );
        boost::mpl::for_each<V>( pointer_param<global_ptr, true>(src, "vals_src") );
        boost::mpl::for_each<V>( pointer_param<global_ptr      >(src, "vals_dst") );

        src.end_kernel_parameters();

        const int NV = NT * VT;

        src.new_line() << "union Shared";
        src.open("{");

        src.new_line() << "struct";
        src.open("{");
        boost::mpl::for_each<K>( type_iterator([&](size_t pos, std::string tname) {
                    src.new_line() << tname << " keys" << pos << "[" << NT * (VT + 1) << "];";
                    }) );
        src.close("};");

        if (boost::mpl::size<V>::value) {
            src.new_line() << "struct";
            src.open("{");
            boost::mpl::for_each<V>( type_iterator([&](size_t pos, std::string tname) {
                        src.new_line() << tname << " vals" << pos << "[" << NT * VT << "];";
                        }) );
            src.close("};");
        }

        src.close("};");

        src.smem_static_var("union Shared", "shared");

        src.new_line() << "int tid    = " << src.local_id(0) << ";";
        src.new_line() << "int block  = " << src.group_id(0) << ";";
        src.new_line() << "int gid    = " << NV << " * block;";
        src.new_line() << "int count2 = min(" << NV << ", count - gid);";

        // Load the values into thread order.
        boost::mpl::for_each<V>( type_iterator([&](size_t pos, std::string tname) {
                    src.new_line() << tname << " thread_vals" << pos << "[" << VT << "];";
                    }) );

        boost::mpl::for_each<V>( call_global_to_shared<NT, VT>(src, "vals_src", "shared.vals") );
        boost::mpl::for_each<V>( call_shared_to_thread<VT>(src, "shared.vals", "thread_vals") );

        // Load keys into shared memory and transpose into register in thread order.
        boost::mpl::for_each<K>( type_iterator([&](size_t pos, std::string tname) {
                    src.new_line() << tname << " thread_keys" << pos << "[" << VT << "];";
                    }) );

        boost::mpl::for_each<K>( call_global_to_shared<NT, VT>(src, "keys_src", "shared.keys") );
        boost::mpl::for_each<K>( call_shared_to_thread<VT>(src, "shared.keys", "thread_keys") );

        // If we're in the last tile, set the uninitialized keys for the thread with
        // a partial number of keys.
        src.new_line() << "int first = " << VT << " * tid;";
        src.new_line() << "if(first + " << VT << " > count2 && first < count2)";
        src.open("{");

        boost::mpl::for_each<K>( type_iterator([&](size_t pos, std::string tname) {
                    src.new_line() << tname << " max_key" << pos << " = thread_keys" << pos << "[0];";
                    }) );

        for(int i = 1; i < VT; ++i) {
            src.new_line()
                << "if(first + " << i << " < count2 && comp(";
            for(int p = 0; p < boost::mpl::size<K>::value; ++p)
                src << (p ? ", " : "") << "max_key" << p;
            for(int p = 0; p < boost::mpl::size<K>::value; ++p)
                src << ", thread_keys" << p << "[" << i << "]";
            src << ") )";
            src.open("{");
            for(int p = 0; p < boost::mpl::size<K>::value; ++p)
            src.new_line() << "max_key" << p << " = thread_keys" << p << "[" << i << "];";
            src.close("}");
        }

        // Fill in the uninitialized elements with max key.
        for(int i = 0; i < VT; ++i) {
            src.new_line()
                << "if(first + " << i << " >= count2)";
            src.open("{");
            for(int p = 0; p < boost::mpl::size<K>::value; ++p)
                src.new_line() << "thread_keys" << p << "[" << i << "] = max_key" << p << ";";
            src.close("}");
        }

        src.close("}");

        src.new_line() << mergesort<NT, VT, K, V>()
            << "(count2, tid";
        for(int p = 0; p < boost::mpl::size<K>::value; ++p)
            src << ", thread_keys" << p;
        for(int p = 0; p < boost::mpl::size<K>::value; ++p)
            src << ", shared.keys" << p;
        for(int p = 0; p < boost::mpl::size<V>::value; ++p)
            src << ", thread_vals" << p;
        for(int p = 0; p < boost::mpl::size<V>::value; ++p)
            src << ", shared.vals" << p;
        src << ");";

        // Store the sorted keys to global.
        boost::mpl::for_each<K>( call_shared_to_global<NT, VT>(src, "count2", "shared.keys", "keys_dst", "gid") );

        boost::mpl::for_each<V>( call_thread_to_shared<VT>(src, "thread_vals", "shared.vals") );
        boost::mpl::for_each<V>( call_shared_to_global<NT, VT>(src, "count2", "shared.vals", "vals_dst", "gid") );

        src.end_kernel();

        kernel = cache.insert(queue, backend::kernel(
                    queue, src.str(), "block_sort"));
    }

    return kernel->second;
}

//---------------------------------------------------------------------------
// Merge partition kernel
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
inline void find_mergesort_frame(backend::source_generator &src) {
    src.begin_function<cl_int3>("find_mergesort_frame");
    src.begin_function_parameters();
    src.parameter<int>("coop");
    src.parameter<int>("block");
    src.parameter<int>("nv");
    src.end_function_parameters();

    src.new_line() << "int start = ~(coop - 1) & block;";
    src.new_line() << "int size = nv * (coop>> 1);";

    src.new_line() << type_name<cl_int3>() << " frame;";
    src.new_line() << "frame.x = nv * start;";
    src.new_line() << "frame.y = nv * start + size;";
    src.new_line() << "frame.z = size;";
    src.new_line() << "return frame;";

    src.end_function();
}

//---------------------------------------------------------------------------
template <int NT, typename T, typename Comp>
backend::kernel merge_partition_kernel(const backend::command_queue &queue) {
    static detail::kernel_cache cache;

    auto kernel = cache.find(queue);

    if (kernel == cache.end()) {
        backend::source_generator src(queue);

        Comp::define(src, "comp");
        merge_path< global_ptr, T >(src);
        find_mergesort_frame(src);

        src.begin_kernel("merge_partition");
        src.begin_kernel_parameters();
        src.template parameter< int                 >("a_count");
        src.template parameter< int                 >("b_count");
        src.template parameter< int                 >("nv");
        src.template parameter< int                 >("coop");
        src.template parameter< global_ptr<int>     >("mp_global");
        src.template parameter< int                 >("num_searches");

        boost::mpl::for_each<T>( pointer_param<global_ptr, true>(src, "a_global") );
        boost::mpl::for_each<T>( pointer_param<global_ptr, true>(src, "b_global") );

        src.end_kernel_parameters();

        src.new_line() << "int partition = " << src.global_id(0) << ";";
        src.new_line() << "if (partition < num_searches)";
        src.open("{");
        src.new_line() << "int a0 = 0, b0 = 0;";
        src.new_line() << "int gid = nv * partition;";

        src.new_line() << "if(coop)";
        src.open("{");
        src.new_line() << type_name<cl_int3>() << " frame = find_mergesort_frame(coop, partition, nv);";
        src.new_line() << "a0 = frame.x;";
        src.new_line() << "b0 = min(a_count, frame.y);";
        src.new_line() << "b_count = min(a_count, frame.y + frame.z) - b0;";
        src.new_line() << "a_count = min(a_count, frame.x + frame.z) - a0;";

        // Put the cross-diagonal into the coordinate system of the input lists.
        src.new_line() << "gid -= a0;";
        src.close("}");

        src.new_line() << "int mp = " << merge_path< global_ptr, T >()
            << "(a_count, b_count, min(gid, a_count + b_count)";
        for(int p = 0; p < boost::mpl::size<T>::value; ++p)
            src << ", a_global" << p << " + a0";
        for(int p = 0; p < boost::mpl::size<T>::value; ++p)
            src << ", b_global" << p << " + b0";
        src << ");";
        src.new_line() << "mp_global[partition] = mp;";

        src.close("}");

        src.end_kernel();

        kernel = cache.insert(queue, backend::kernel(
                    queue, src.str(), "merge_partition"));
    }

    return kernel->second;
}

//---------------------------------------------------------------------------
template <typename Comp, class KT>
backend::device_vector<int> merge_path_partitions(
        const backend::command_queue &queue,
        const KT &keys,
        int count, int nv, int coop
        )
{
    typedef typename extract_value_types<KT>::type K;

    const int NT_cpu = 1;
    const int NT_gpu = 64;
    const int NT = is_cpu(queue) ? NT_cpu : NT_gpu;

    int num_partitions       = (count + nv - 1) / nv;
    int num_partition_blocks = (num_partitions + NT) / NT;

    backend::device_vector<int> partitions(queue, num_partitions + 1);

    auto merge_partition = is_cpu(queue) ?
        merge_partition_kernel<NT_cpu, K, Comp>(queue) :
        merge_partition_kernel<NT_gpu, K, Comp>(queue);

    int a_count = static_cast<int>(boost::fusion::at_c<0>(keys).size());
    int b_count = 0;

    merge_partition.push_arg(a_count);
    merge_partition.push_arg(b_count);
    merge_partition.push_arg(nv);
    merge_partition.push_arg(coop);
    merge_partition.push_arg(partitions);
    merge_partition.push_arg(num_partitions + 1);

    push_args<boost::mpl::size<K>::value>(merge_partition, keys);
    push_args<boost::mpl::size<K>::value>(merge_partition, keys);

    merge_partition.config(num_partition_blocks, NT);

    merge_partition(queue);

    return partitions;
}

//---------------------------------------------------------------------------
// Merge kernel
//---------------------------------------------------------------------------
inline void find_mergesort_interval(backend::source_generator &src) {
    src.begin_function<cl_int4>("find_mergesort_interval");
    src.begin_function_parameters();
    src.parameter< cl_int3 >("frame");
    src.parameter< int     >("coop");
    src.parameter< int     >("block");
    src.parameter< int     >("nv");
    src.parameter< int     >("count");
    src.parameter< int     >("mp0");
    src.parameter< int     >("mp1");
    src.end_function_parameters();

    // Locate diag from the start of the A sublist.
    src.new_line() << "int diag = nv * block - frame.x;";
    src.new_line() << "int4 interval;";
    src.new_line() << "interval.x = frame.x + mp0;";
    src.new_line() << "interval.y = min(count, frame.x + mp1);";
    src.new_line() << "interval.z = min(count, frame.y + diag - mp0);";
    src.new_line() << "interval.w = min(count, frame.y + diag + nv - mp1);";

    // The end partition of the last block for each merge operation is computed
    // and stored as the begin partition for the subsequent merge. i.e. it is
    // the same partition but in the wrong coordinate system, so its 0 when it
    // should be listSize. Correct that by checking if this is the last block
    // in this merge operation.
    src.new_line() << "if(coop - 1 == ((coop - 1) & block))";
    src.open("{");
    src.new_line() << "interval.y = min(count, frame.x + frame.z);";
    src.new_line() << "interval.w = min(count, frame.y + frame.z);";
    src.close("}");

    src.new_line() << "return interval;";

    src.end_function();
}

//---------------------------------------------------------------------------
inline void compute_merge_range(backend::source_generator &src) {
    find_mergesort_frame(src);
    find_mergesort_interval(src);

    src.begin_function<cl_int4>("compute_merge_range");
    src.begin_function_parameters();
    src.parameter< int >("a_count");
    src.parameter< int >("b_count");
    src.parameter< int >("block");
    src.parameter< int >("coop");
    src.parameter< int >("nv");
    src.parameter< global_ptr<const int> >("mp_global");
    src.end_function_parameters();

    // Load the merge paths computed by the partitioning kernel.
    src.new_line() << "int mp0 = mp_global[block];";
    src.new_line() << "int mp1 = mp_global[block + 1];";
    src.new_line() << "int gid = nv * block;";

    // Compute the ranges of the sources in global memory.
    src.new_line() << "int4 range;";
    src.new_line() << "if(coop)";
    src.open("{");
    src.new_line() << type_name<cl_int3>() << " frame = find_mergesort_frame(coop, block, nv);";
    src.new_line() << "range = find_mergesort_interval(frame, coop, block, nv, a_count, mp0, mp1);";
    src.close("}");
    src.new_line() << "else";
    src.open("{");
    src.new_line() << "range.x = mp0;";
    src.new_line() << "range.y = mp1;";
    src.new_line() << "range.z = gid - range.x;";
    src.new_line() << "range.w = min(a_count + b_count, gid + nv) - range.y;";
    src.close("}");

    src.new_line() << "return range;";

    src.end_function();
}

//---------------------------------------------------------------------------
template<int NT, int vex_VT0, int vex_VT1, typename T>
std::string load2_to_regstr() {
    std::ostringstream s;
    s << "load2_to_regstr_" << NT << "_" << vex_VT0 << "_" << vex_VT1 << "_" << type_name<T>();
    return s.str();
}

template<int NT, int vex_VT0, int vex_VT1, typename T>
void load2_to_regstr(backend::source_generator &src) {
    src.begin_function<void>(load2_to_regstr<NT,vex_VT0,vex_VT1,T>());
    src.begin_function_parameters();
    src.template parameter< global_ptr<const T> >("a_global");
    src.template parameter< int                 >("a_count");
    src.template parameter< global_ptr<const T> >("b_global");
    src.template parameter< int                 >("b_count");
    src.template parameter< int                 >("tid");
    src.template parameter< regstr_ptr<T>       >("reg");
    src.end_function_parameters();

    src.new_line() << "b_global -= a_count;";
    src.new_line() << "int total = a_count + b_count;";
    src.new_line() << "int index;";
    src.new_line() << "if (total >= " << NT * vex_VT0 << ")";
    src.open("{");

    for(int i = 0; i < vex_VT0; ++i) {
        src.new_line() << "index = " << NT * i << " + tid;";
        src.new_line() << "if (index < a_count) reg[" << i << "] = a_global[index];";
        src.new_line() << "else reg[" << i << "] = b_global[index];";
    }

    src.close("}");
    src.new_line() << "else";
    src.open("{");

    for(int i = 0; i < vex_VT0; ++i) {
        src.new_line() << "index = " << NT * i << " + tid;";
        src.new_line() << "if (index < a_count) reg[" << i << "] = a_global[index];";
        src.new_line() << "else if (index < total) reg[" << i << "] = b_global[index];";
    }

    src.close("}");

    for(int i = vex_VT0; i < vex_VT1; ++i) {
        src.new_line() << "index = " << NT * i << " + tid;";
        src.new_line() << "if (index < a_count) reg[" << i << "] = a_global[index];";
        src.new_line() << "else if (index < total) reg[" << i << "] = b_global[index];";
    }

    src.end_function();
}

//---------------------------------------------------------------------------
template<int NT, int vex_VT0, int vex_VT1, typename T>
std::string load2_to_shared() {
    std::ostringstream s;
    s << "load2_to_shared_" << NT << "_" << vex_VT0 << "_" << vex_VT1 << "_" << type_name<T>();
    return s.str();
}

template<int NT, int vex_VT0, int vex_VT1, typename T>
void load2_to_shared(backend::source_generator &src) {
    load2_to_regstr<NT, vex_VT0, vex_VT1, T>(src);

    src.begin_function<void>(load2_to_shared<NT,vex_VT0,vex_VT1,T>());
    src.begin_function_parameters();
    src.template parameter< global_ptr<const T> >("a_global");
    src.template parameter< int                 >("a_count");
    src.template parameter< global_ptr<const T> >("b_global");
    src.template parameter< int                 >("b_count");
    src.template parameter< int                 >("tid");
    src.template parameter< shared_ptr<T>       >("shared");
    src.end_function_parameters();

    src.new_line() << type_name<T>() << " reg[" << vex_VT1 << "];";
    src.new_line() << load2_to_regstr<NT, vex_VT0, vex_VT1, T>()
        << "(a_global, a_count, b_global, b_count, tid, reg);";
    src.new_line() << regstr_to_shared<NT, vex_VT1, T>()
        << "(reg, tid, shared);";

    src.end_function();
}

//---------------------------------------------------------------------------
template <int NT, int VT, typename T>
std::string merge_keys_indices() {
    std::ostringstream s;
    s << "merge_keys_indices_" << NT << "_" << VT;
    print_types<T>(s);
    return s.str();
}


template <int NT, int VT>
struct define_load2_to_shared {
    backend::source_generator &src;

    define_load2_to_shared(backend::source_generator &src) : src(src) {}

    template <typename T>
    void operator()(T) const {
        load2_to_shared<NT, VT, VT, T>(src);
    }
};

template <int NT, int VT>
struct call_load2_to_shared {
    backend::source_generator &src;
    int pos;

    call_load2_to_shared(backend::source_generator &src) : src(src), pos(0) {}

    template <typename T>
    void operator()(T) {
        src.new_line() << load2_to_shared<NT, VT, VT, T>()
            << "(a_global" << pos << " + a0, a_count, b_global" << pos
            << " + b0, b_count, tid, keys_shared" << pos << ");";

        ++pos;
    }
};

template <int NT, int VT, typename T>
void merge_keys_indices(backend::source_generator &src) {
    serial_merge<VT, T >(src);
    merge_path< shared_ptr, T >(src);

    boost::mpl::for_each<
        typename boost::mpl::copy<
            T,
            boost::mpl::inserter<
                boost::mpl::set0<>,
                boost::mpl::insert<boost::mpl::_1, boost::mpl::_2>
                >
            >::type
        >( define_load2_to_shared<NT, VT>(src) );

    src.begin_function<void>(merge_keys_indices<NT, VT, T>());
    src.begin_function_parameters();
    src.template parameter< int                 >("a_count");
    src.template parameter< int                 >("b_count");
    src.template parameter< cl_int4             >("range");
    src.template parameter< int                 >("tid");
    src.template parameter< regstr_ptr<int>     >("indices");

    boost::mpl::for_each<T>( pointer_param<global_ptr, true>(src, "a_global") );
    boost::mpl::for_each<T>( pointer_param<global_ptr, true>(src, "b_global") );
    boost::mpl::for_each<T>( pointer_param<shared_ptr      >(src, "keys_shared") );
    boost::mpl::for_each<T>( pointer_param<regstr_ptr      >(src, "results") );

    src.end_function_parameters();

    src.new_line() << "int a0 = range.x;";
    src.new_line() << "int a1 = range.y;";
    src.new_line() << "int b0 = range.z;";
    src.new_line() << "int b1 = range.w;";

    // Use the input intervals from the ranges between the merge path
    // intersections.
    src.new_line() << "a_count = a1 - a0;";
    src.new_line() << "b_count = b1 - b0;";

    // Load the data into shared memory.
    boost::mpl::for_each<T>( call_load2_to_shared<NT, VT>(src) );

    // Run a merge path to find the start of the serial merge for each
    // thread.
    src.new_line() << "int diag = " << VT << " * tid;";
    src.new_line() << "int mp = " << merge_path< shared_ptr, T >()
        << "(a_count, b_count, diag";
    for(int p = 0; p < boost::mpl::size<T>::value; ++p)
        src << ", keys_shared" << p;
    for(int p = 0; p < boost::mpl::size<T>::value; ++p)
        src << ", keys_shared" << p << " + a_count";
    src << ");";

    // Compute the ranges of the sources in shared memory.
    src.new_line() << "int a0tid = mp;";
    src.new_line() << "int a1tid = a_count;";
    src.new_line() << "int b0tid = a_count + diag - mp;";
    src.new_line() << "int b1tid = a_count + b_count;";

    // Serial merge into register.
    src.new_line() << serial_merge<VT, T >()
        << "(a0tid, a1tid, b0tid, b1tid, indices";
    for(int p = 0; p < boost::mpl::size<T>::value; ++p)
        src << ", keys_shared" << p;
    for(int p = 0; p < boost::mpl::size<T>::value; ++p)
        src << ", results" << p;
    src << ");";

    src.end_function();
}

//---------------------------------------------------------------------------
template <int NT, int VT, typename T>
std::string transfer_merge_values_regstr() {
    std::ostringstream s;
    s << "transfer_merge_values_regstr_" << NT << "_" << VT;
    print_types<T>(s);
    return s.str();
}

template <int NT, int VT, typename T>
void transfer_merge_values_regstr(backend::source_generator &src) {
    src.begin_function<void>(transfer_merge_values_regstr<NT, VT, T>());
    src.begin_function_parameters();
    src.template parameter< int                   >("count");
    src.template parameter< int                   >("b_start");
    src.template parameter< regstr_ptr<const int> >("indices");
    src.template parameter< int                   >("tid");

    boost::mpl::for_each<T>( pointer_param<global_ptr, true>(src, "a_global") );
    boost::mpl::for_each<T>( pointer_param<global_ptr, true>(src, "b_global") );
    boost::mpl::for_each<T>( pointer_param<regstr_ptr      >(src, "reg") );

    src.end_function_parameters();

    for(int p = 0; p < boost::mpl::size<T>::value; ++p)
        src.new_line() << "b_global" << p << " -= b_start;";
    src.new_line() << "if(count >= " << NT * VT << ")";
    src.open("{");

    for(int i = 0; i < VT; ++i) {
        src.new_line() << "if (indices[" << i << "] < b_start)";
        src.open("{");
        for(int p = 0; p < boost::mpl::size<T>::value; ++p)
            src.new_line() << "reg" << p << "[" << i << "] = a_global" << p << "[indices[" << i << "]];";
        src.close("}");
        src.new_line() << "else";
        src.open("{");
        for(int p = 0; p < boost::mpl::size<T>::value; ++p)
            src.new_line() << "reg" << p << "[" << i << "] = b_global" << p << "[indices[" << i << "]];";
        src.close("}");
    }

    src.close("}");
    src.new_line() << "else";
    src.open("{");

    src.new_line() << "int index;";

    for(int i = 0; i < VT; ++i) {
        src.new_line() << "index = " << NT * i << " + tid;";
        src.new_line() << "if(index < count)";
        src.open("{");
        src.new_line() << "if (indices[" << i << "] < b_start)";
        src.open("{");
        for(int p = 0; p < boost::mpl::size<T>::value; ++p)
            src.new_line() << "reg" << p << "[" << i << "] = a_global" << p << "[indices[" << i << "]];";
        src.close("}");
        src.new_line() << "else";
        src.open("{");
        for(int p = 0; p < boost::mpl::size<T>::value; ++p)
            src.new_line() << "reg" << p << "[" << i << "] = b_global" << p << "[indices[" << i << "]];";
        src.close("}");
        src.close("}");
    }
    src.close("}");

    src.new_line().barrier();

    src.end_function();
}

//---------------------------------------------------------------------------
template <int NT, int VT, typename T>
std::string transfer_merge_values_shared() {
    std::ostringstream s;
    s << "transfer_merge_values_shared_" << NT << "_" << VT;
    print_types<T>(s);
    return s.str();
}

template <int NT, int VT>
struct call_regstr_to_global {
    backend::source_generator &src;
    int pos;

    call_regstr_to_global(backend::source_generator &src) : src(src), pos(0) {}

    template <typename T>
    void operator()(T) {
        src.new_line() << regstr_to_global<NT, VT, T>()
            << "(count, reg" << pos << ", tid, dest_global" << pos << ");";
        ++pos;
    }
};

template <int NT, int VT, typename T>
void transfer_merge_values_shared(backend::source_generator &src) {
    transfer_merge_values_regstr<NT, VT, T >(src);

    src.begin_function<void>(transfer_merge_values_shared<NT, VT, T>());
    src.begin_function_parameters();
    src.template parameter< int                   >("count");
    src.template parameter< int                   >("b_start");
    src.template parameter< shared_ptr<const int> >("indices_shared");
    src.template parameter< int                   >("tid");

    boost::mpl::for_each<T>( pointer_param<global_ptr, true>(src, "a_global") );
    boost::mpl::for_each<T>( pointer_param<global_ptr, true>(src, "b_global") );
    boost::mpl::for_each<T>( pointer_param<global_ptr      >(src, "dest_global") );

    src.end_function_parameters();

    src.new_line() << "int indices[" << VT << "];";
    src.new_line() << shared_to_regstr<NT, VT, int>()
        << "(indices_shared, tid, indices);";

    boost::mpl::for_each<T>( type_iterator([&](size_t pos, std::string tname) {
                src.new_line() << tname << " reg" << pos++ << "[" << VT << "];";
                }) );

    src.new_line() << transfer_merge_values_regstr<NT, VT, T >()
        << "(count, b_start, indices, tid";
    for(int p = 0; p < boost::mpl::size<T>::value; ++p)
        src << ", a_global" << p;
    for(int p = 0; p < boost::mpl::size<T>::value; ++p)
        src << ", b_global" << p;
    for(int p = 0; p < boost::mpl::size<T>::value; ++p)
        src << ", reg" << p;
    src << ");";

    boost::mpl::for_each<T>( call_regstr_to_global<NT, VT>(src) );

    src.end_function();
}

//---------------------------------------------------------------------------
template<int NT, int VT, typename K, typename V>
std::string device_merge() {
    std::ostringstream s;
    s << "device_merge_" << NT << "_" << VT;
    print_types<K>(s);
    print_types<V>(s);
    return s.str();
}

template<int NT, int VT, typename K, typename V>
void device_merge(backend::source_generator &src) {
    merge_keys_indices<NT, VT, K >(src);
    if (boost::mpl::size<V>::value)
        transfer_merge_values_shared<NT, VT, V >(src);

    src.begin_function<void>(device_merge<NT, VT, K, V>());
    src.begin_function_parameters();
    src.template parameter< int >("a_count");
    src.template parameter< int >("b_count");

    boost::mpl::for_each<K>( pointer_param<global_ptr, true>(src, "a_keys_global"));
    boost::mpl::for_each<K>( pointer_param<global_ptr, true>(src, "b_keys_global"));
    boost::mpl::for_each<K>( pointer_param<global_ptr      >(src, "keys_global"));
    boost::mpl::for_each<K>( pointer_param<shared_ptr      >(src, "keys_shared"));

    boost::mpl::for_each<V>( pointer_param<global_ptr, true>(src, "a_vals_global"));
    boost::mpl::for_each<V>( pointer_param<global_ptr, true>(src, "b_vals_global"));
    boost::mpl::for_each<V>( pointer_param<global_ptr      >(src, "vals_global"));

    src.template parameter< int             >("tid");
    src.template parameter< int             >("block");
    src.template parameter< cl_int4         >("range");
    src.template parameter< shared_ptr<int> >("indices_shared");
    src.end_function_parameters();

    boost::mpl::for_each<K>( type_iterator([&](size_t pos, std::string tname) {
                src.new_line() << tname << " results" << pos++ << "[" << VT << "];";
                }) );

    src.new_line() << "int indices[" << VT << "];";

    src.new_line() << merge_keys_indices<NT, VT, K >()
        << "(a_count, b_count, range, tid, indices";

    for(int p = 0; p < boost::mpl::size<K>::value; ++p)
        src << ", a_keys_global" << p;
    for(int p = 0; p < boost::mpl::size<K>::value; ++p)
        src << ", b_keys_global" << p;
    for(int p = 0; p < boost::mpl::size<K>::value; ++p)
        src << ", keys_shared" << p;
    for(int p = 0; p < boost::mpl::size<K>::value; ++p)
        src << ", results" << p;

    src << ");";

    // Store merge results back to shared memory.
    boost::mpl::for_each<K>( call_thread_to_shared<VT>(src, "results", "keys_shared") );

    // Store merged keys to global memory.
    src.new_line() << "a_count = range.y - range.x;";
    src.new_line() << "b_count = range.w - range.z;";

    {
        std::ostringstream s;
        s << NT * VT << " * block";
        boost::mpl::for_each<K>( call_shared_to_global<NT, VT>(src, "a_count + b_count", "keys_shared", "keys_global", s.str().c_str() ) );
    }

    // Copy the values.
    if (boost::mpl::size<V>::value) {
        src.new_line() << thread_to_shared<VT, int>()
            << "(indices, tid, indices_shared);";

        src.new_line() << transfer_merge_values_shared<NT, VT, V >()
            << "(a_count + b_count, a_count, indices_shared, tid";

        for(int p = 0; p < boost::mpl::size<V>::value; ++p)
            src << ", a_vals_global" << p << " + range.x";
        for(int p = 0; p < boost::mpl::size<V>::value; ++p)
            src << ", b_vals_global" << p << " + range.z";
        for(int p = 0; p < boost::mpl::size<V>::value; ++p)
            src << ", vals_global" << p << " + " << NT * VT << " * block";
        src << ");";
    }

    src.end_function();
}

//---------------------------------------------------------------------------
template <int NT, int VT, typename K, typename V, typename Comp>
backend::kernel merge_kernel(const backend::command_queue &queue) {
    static detail::kernel_cache cache;

    auto kernel = cache.find(queue);

    if (kernel == cache.end()) {
        backend::source_generator src(queue);

        Comp::define(src, "comp");
        compute_merge_range(src);

        boost::mpl::for_each<
            typename boost::mpl::copy<
                typename boost::mpl::copy<
                    V,
                    boost::mpl::back_inserter<K>
                    >::type,
                boost::mpl::inserter<
                    boost::mpl::set<int>,
                    boost::mpl::insert<boost::mpl::_1, boost::mpl::_2>
                    >
                >::type
            >( define_transfer_functions<NT, VT>(src) );

        device_merge<NT, VT, K, V >(src);

        src.begin_kernel("merge");
        src.begin_kernel_parameters();
        src.template parameter< int >("a_count");
        src.template parameter< int >("b_count");

        boost::mpl::for_each<K>( pointer_param<global_ptr, true>(src, "a_keys_global"));
        boost::mpl::for_each<K>( pointer_param<global_ptr, true>(src, "b_keys_global"));
        boost::mpl::for_each<K>( pointer_param<global_ptr      >(src, "keys_global"));

        boost::mpl::for_each<V>( pointer_param<global_ptr, true>(src, "a_vals_global"));
        boost::mpl::for_each<V>( pointer_param<global_ptr, true>(src, "b_vals_global"));
        boost::mpl::for_each<V>( pointer_param<global_ptr      >(src, "vals_global"));

        src.template parameter< global_ptr<const int> >("mp_global");
        src.template parameter< int                   >("coop");
        src.end_kernel_parameters();

        const int NV = NT * VT;

        src.new_line() << "union Shared";
        src.open("{");

        src.new_line() << "struct";
        src.open("{");
        boost::mpl::for_each<K>( type_iterator([&](size_t pos, std::string tname) {
                    src.new_line() << tname << " keys" << pos++ << "[" << NT * (VT + 1) << "];";
                    }) );
        src.close("};");

        src.new_line() << "int indices[" << NV << "];";
        src.close("};");

        src.smem_static_var("union Shared", "shared");

        src.new_line() << "int tid    = " << src.local_id(0) << ";";
        src.new_line() << "int block  = " << src.group_id(0) << ";";

        src.new_line() << "int4 range = compute_merge_range("
            "a_count, b_count, block, coop, " << NV << ", mp_global);";

        src.new_line() << device_merge<NT, VT, K, V >()
            << "(a_count, b_count";

        for(int p = 0; p < boost::mpl::size<K>::value; ++p)
            src << ", a_keys_global" << p;
        for(int p = 0; p < boost::mpl::size<K>::value; ++p)
            src << ", b_keys_global" << p;
        for(int p = 0; p < boost::mpl::size<K>::value; ++p)
            src << ", keys_global" << p;
        for(int p = 0; p < boost::mpl::size<K>::value; ++p)
            src << ", shared.keys" << p;

        for(int p = 0; p < boost::mpl::size<V>::value; ++p)
            src << ", a_vals_global" << p;
        for(int p = 0; p < boost::mpl::size<V>::value; ++p)
            src << ", b_vals_global" << p;
        for(int p = 0; p < boost::mpl::size<V>::value; ++p)
            src << ", vals_global" << p;

        src << ", tid, block, range, shared.indices);";

        src.end_kernel();

        kernel = cache.insert(queue, backend::kernel(
                    queue, src.str(), "merge"));
    }

    return kernel->second;
}

//---------------------------------------------------------------------------
inline int clz(int x) {
    for(int i = 31; i >= 0; --i)
        if((1 << i) & x) return 31 - i;
    return 32;
}

//---------------------------------------------------------------------------
inline int find_log2(int x, bool round_up = false) {
    int a = 31 - clz(x);

    if(round_up) {
        bool is_pow_2 = (0 == (x & (x - 1)));
        a += !is_pow_2;
    }

    return a;
}

/// Sorts single partition of a vector.
template <class KT, class Comp>
void sort(const backend::command_queue &queue, KT &keys, Comp) {
    typedef typename extract_value_types<KT>::type K;
    using boost::fusion::at_c;

    typedef
        typename boost::mpl::accumulate<
            K,
            boost::mpl::int_<0>,
            boost::mpl::plus<boost::mpl::_1, boost::mpl::sizeof_<boost::mpl::_2> >
            >::type
        sizeof_keys;

    backend::select_context(queue);

    const int NT_cpu = 1;
    const int NT_gpu = VEX_SORT_NT_GPU;
    const int NT = is_cpu(queue) ? NT_cpu : NT_gpu;
    const int VT = (sizeof_keys::value > 4) ? 7 : 11;
    const int NV = NT * VT;

    const int count = static_cast<int>(at_c<0>(keys).size());
    const int num_blocks = (count + NV - 1) / NV;
    const int num_passes = detail::find_log2(num_blocks, true);

    temp_storage<K> tmp(queue, count);

    auto block_sort = is_cpu(queue) ?
        detail::block_sort_kernel<NT_cpu, VT, K, boost::mpl::vector<>, Comp>(queue) :
        detail::block_sort_kernel<NT_gpu, VT, K, boost::mpl::vector<>, Comp>(queue);

    block_sort.push_arg(count);

    push_args<boost::mpl::size<K>::value>(block_sort, keys);
    if (1 & num_passes)
        push_args<boost::mpl::size<K>::value>(block_sort, tmp);
    else
        push_args<boost::mpl::size<K>::value>(block_sort, keys);

    block_sort.config(num_blocks, NT);

    block_sort(queue);

    if (1 & num_passes) tmp.swap(keys);

    auto merge = is_cpu(queue) ?
        detail::merge_kernel<NT_cpu, VT, K, boost::mpl::vector<>, Comp>(queue) :
        detail::merge_kernel<NT_gpu, VT, K, boost::mpl::vector<>, Comp>(queue);

    for(int pass = 0; pass < num_passes; ++pass) {
        int coop = 2 << pass;

        auto partitions = detail::merge_path_partitions<Comp>(queue, keys, count, NV, coop);

        merge.push_arg(count);
        merge.push_arg(0);

        push_args<boost::mpl::size<K>::value>(merge, keys);
        push_args<boost::mpl::size<K>::value>(merge, keys);
        push_args<boost::mpl::size<K>::value>(merge, tmp);
        merge.push_arg(partitions);
        merge.push_arg(coop);

        merge.config(num_blocks, NT);
        merge(queue);

        tmp.swap(keys);
    }
}

/// Sorts single partition of a vector.
template <class KTup, class VTup, class Comp>
void sort_by_key(const backend::command_queue &queue, KTup &&keys, VTup &&vals, Comp) {
    typedef typename extract_value_types<KTup>::type K;
    typedef typename extract_value_types<VTup>::type V;

    using boost::fusion::at_c;

    precondition(at_c<0>(keys).size() == at_c<0>(vals).size(),
            "keys and values should have same size"
            );

    backend::select_context(queue);

    const int NT_cpu = 1;
    const int NT_gpu = VEX_SORT_NT_GPU;
    const int NT = is_cpu(queue) ? NT_cpu : NT_gpu;
    const int VT = (sizeof(K) > 4) ? 7 : 11;
    const int NV = NT * VT;

    const int count = static_cast<int>(at_c<0>(keys).size());
    const int num_blocks = (count + NV - 1) / NV;
    const int num_passes = detail::find_log2(num_blocks, true);

    temp_storage<K> keys_tmp(queue, count);
    temp_storage<V> vals_tmp(queue, count);

    auto block_sort = is_cpu(queue) ?
        detail::block_sort_kernel<NT_cpu, VT, K, V, Comp>(queue) :
        detail::block_sort_kernel<NT_gpu, VT, K, V, Comp>(queue);

    block_sort.push_arg(count);

    push_args<boost::mpl::size<K>::value>(block_sort, keys);
    if (1 & num_passes)
        push_args<boost::mpl::size<K>::value>(block_sort, keys_tmp);
    else
        push_args<boost::mpl::size<K>::value>(block_sort, keys);

    push_args<boost::mpl::size<V>::value>(block_sort, vals);
    if (1 & num_passes)
        push_args<boost::mpl::size<V>::value>(block_sort, vals_tmp);
    else
        push_args<boost::mpl::size<V>::value>(block_sort, vals);

    block_sort.config(num_blocks, NT);

    block_sort(queue);

    if (1 & num_passes) {
        keys_tmp.swap(keys);
        vals_tmp.swap(vals);
    }

    auto merge = is_cpu(queue) ?
        detail::merge_kernel<NT_cpu, VT, K, V, Comp>(queue) :
        detail::merge_kernel<NT_gpu, VT, K, V, Comp>(queue);

    for(int pass = 0; pass < num_passes; ++pass) {
        int coop = 2 << pass;

        auto partitions = detail::merge_path_partitions<Comp>(queue, keys, count, NV, coop);

        merge.push_arg(count);
        merge.push_arg(0);

        push_args<boost::mpl::size<K>::value>(merge, keys);
        push_args<boost::mpl::size<K>::value>(merge, keys);
        push_args<boost::mpl::size<K>::value>(merge, keys_tmp);

        push_args<boost::mpl::size<V>::value>(merge, vals);
        push_args<boost::mpl::size<V>::value>(merge, vals);
        push_args<boost::mpl::size<V>::value>(merge, vals_tmp);

        merge.push_arg(partitions);
        merge.push_arg(coop);

        merge.config(num_blocks, NT);
        merge(queue);

        keys_tmp.swap(keys);
        vals_tmp.swap(vals);
    }
}

template <class S1, class S2>
boost::fusion::zip_view< boost::fusion::vector<S1&, S2&> >
make_zip_view(S1 &s1, S2 &s2) {
    typedef boost::fusion::vector<S1&, S2&> Z;
    return boost::fusion::zip_view<Z>( Z(s1, s2));
}

struct do_resize {
    size_t n;
    do_resize(size_t n) : n(n) {}
    template <class V>
    void operator()(V &v) const {
        v.resize(n);
    }
};

struct do_copy {
    template <class T>
    void operator()(T t) const {
        using boost::fusion::at_c;
        vex::copy(at_c<0>(t), at_c<1>(t));
    }
};

struct copy_element {
    size_t dst, src;
    copy_element(size_t dst, size_t src) : dst(dst), src(src) {}

    template <class T>
    void operator()(T t) const {
        using boost::fusion::at_c;
        at_c<0>(t)[dst] = at_c<1>(t)[src];
    }
};

struct do_index {
    size_t pos;
    do_index(size_t pos) : pos(pos) {}

    template <class T> struct result;

    template <class This, class T>
    struct result< This(T) > {
        typedef typename std::decay<T>::type::value_type type;
    };

    template <class T>
    typename result<do_index(T)>::type operator()(const T &t) const {
        return t[pos];
    }
};

/// Merges partially sorted vector partitions into host vector
template <typename KTuple, class Comp>
typename boost::fusion::result_of::as_vector<
    typename boost::mpl::transform<
        typename extract_value_types<KTuple>::type,
        std::vector<boost::mpl::_1>
    >::type
>::type
merge(const KTuple &keys, Comp comp) {
    namespace fusion = boost::fusion;

    typedef typename extract_value_types<KTuple>::type K;

    const auto  &queue = fusion::at_c<0>(keys).queue_list();
    const size_t count = fusion::at_c<0>(keys).size();

    typedef typename fusion::result_of::as_vector<
        typename boost::mpl::transform< K, std::vector<boost::mpl::_1> >::type
    >::type host_vectors;

    host_vectors dst;

    fusion::for_each(dst, do_resize(count) );
    fusion::for_each( make_zip_view(keys, dst), do_copy() );

    if (queue.size() > 1) {
        host_vectors src;
        fusion::for_each(src, do_resize(count) );

        fusion::swap(src, dst);

        std::vector<size_t> begin(queue.size());
        std::vector<size_t> end  (queue.size());

        for(unsigned d = 0; d < queue.size(); ++d) {
            begin[d] = fusion::at_c<0>(keys).part_start(d);
            end  [d] = fusion::at_c<0>(keys).part_start(d + 1);
        }

        for(size_t pos = 0; pos < count; ++pos) {
            int winner = -1;
            for(unsigned d = 0; d < queue.size(); ++d) {
                if (begin[d] == end[d]) continue;

                if (winner < 0) {
                    winner = d;
                    continue;
                }

                auto curr = fusion::transform(src, do_index(begin[d]));
                auto best = fusion::transform(src, do_index(begin[winner]));

                if (fusion::invoke(comp, fusion::join(curr, best)))
                    winner = d;
            }

            fusion::for_each(make_zip_view(dst, src), copy_element(pos, begin[winner]++));
        }
    }

    return dst;
}

/// Merges partially sorted vector partitions into host vector
template <typename KTuple, typename VTuple, class Comp>
typename boost::fusion::result_of::as_vector<
    typename boost::mpl::transform<
        typename extract_value_types<
            typename boost::fusion::result_of::as_vector<
                typename boost::fusion::result_of::join<KTuple, VTuple>::type
            >::type
        >::type,
        std::vector<boost::mpl::_1>
    >::type
>::type
merge(const KTuple &keys, const VTuple &vals, Comp comp) {
    namespace fusion = boost::fusion;

    typedef typename extract_value_types<KTuple>::type K;
    typedef typename extract_value_types<VTuple>::type V;

    const auto  &queue = fusion::at_c<0>(keys).queue_list();
    const size_t count = fusion::at_c<0>(keys).size();

    typedef typename fusion::result_of::as_vector<
        typename boost::mpl::transform< K, std::vector<boost::mpl::_1> >::type
    >::type host_keys;

    typedef typename fusion::result_of::as_vector<
        typename boost::mpl::transform< V, std::vector<boost::mpl::_1> >::type
    >::type host_vals;

    host_keys dst_keys;
    host_vals dst_vals;

    fusion::for_each(dst_keys, do_resize(count) );
    fusion::for_each(dst_vals, do_resize(count) );

    fusion::for_each( make_zip_view(keys, dst_keys), do_copy() );
    fusion::for_each( make_zip_view(vals, dst_vals), do_copy() );

    if (queue.size() > 1) {
        host_keys src_keys;
        host_vals src_vals;

        fusion::for_each(src_keys, do_resize(count) );
        fusion::for_each(src_vals, do_resize(count) );

        fusion::swap(src_keys, dst_keys);
        fusion::swap(src_vals, dst_vals);

        std::vector<size_t> begin(queue.size()), end(queue.size());

        for(unsigned d = 0; d < queue.size(); ++d) {
            begin[d] = fusion::at_c<0>(keys).part_start(d);
            end  [d] = fusion::at_c<0>(keys).part_start(d + 1);
        }

        for(size_t pos = 0; pos < count; ++pos) {
            int winner = -1;
            for(unsigned d = 0; d < queue.size(); ++d) {
                if (begin[d] == end[d]) continue;

                if (winner < 0) {
                    winner = d;
                    continue;
                }

                auto curr = fusion::transform(src_keys, do_index(begin[d]));
                auto best = fusion::transform(src_keys, do_index(begin[winner]));

                if (fusion::invoke(comp, fusion::join(curr, best)))
                    winner = d;
            }

            fusion::for_each(make_zip_view(dst_keys, src_keys), copy_element(pos, begin[winner]));
            fusion::for_each(make_zip_view(dst_vals, src_vals), copy_element(pos, begin[winner]));

            ++begin[winner];
        }
    }

    return fusion::as_vector(fusion::join(dst_keys, dst_vals));
}

template <class K, class Comp>
void sort_sink(K &&keys, Comp comp) {
    namespace fusion = boost::fusion;

    const auto &queue = boost::fusion::at_c<0>(keys).queue_list();

    for(unsigned d = 0; d < queue.size(); ++d)
        if (fusion::at_c<0>(keys).part_size(d)) {
            auto part = fusion::transform(keys, extract_device_vector(d));
            sort(queue[d], part, comp.device);
        }

    if (queue.size() <= 1) return;

    // Vector partitions have been sorted on compute devices.
    // Now we need to merge them on a CPU. This is a linear time operation,
    // so total performance should be good enough.
    auto host_vectors = merge(keys, comp);
    fusion::for_each( make_zip_view(host_vectors, keys), do_copy() );
}

template <class K, class V, class Comp>
void sort_by_key_sink(K &&keys, V &&vals, Comp comp) {
    namespace fusion = boost::fusion;

    precondition(
            fusion::at_c<0>(keys).nparts() == fusion::at_c<0>(vals).nparts(),
            "Keys and values span different devices"
            );

    const auto &queue = fusion::at_c<0>(keys).queue_list();

    for(unsigned d = 0; d < queue.size(); ++d)
        if (fusion::at_c<0>(keys).part_size(d)) {
            auto kpart = fusion::transform(keys, extract_device_vector(d));
            auto vpart = fusion::transform(vals, extract_device_vector(d));
            sort_by_key(queue[d], kpart, vpart, comp.device);
        }

    if (queue.size() <= 1) return;

    // Vector partitions have been sorted on compute devices.
    // Now we need to merge them on a CPU. This is a linear time operation,
    // so total performance should be good enough.
    auto host_vectors = merge(keys, vals, comp);
    auto dev_vectors  = fusion::join(keys, vals);
    fusion::for_each( make_zip_view(host_vectors, dev_vectors), do_copy() );
}

} // namespace detail

/// Function object class for less-than inequality comparison.
/**
 * The need for host-side and device-side parts comes from the fact that
 * vectors are partially sorted on device and then final merge step is done on
 * host.
 */
template <typename T>
struct less : std::less<T> {
    VEX_FUNCTION(bool, device, (T, x)(T, y), return x < y;);

    less() {}
};

/// Function object class for less-than-or-equal inequality comparison.
template <typename T>
struct less_equal : std::less_equal<T> {
    VEX_FUNCTION(bool, device, (T, x)(T, y), return x <= y;);

    less_equal() {}
};

/// Function object class for greater-than inequality comparison.
template <typename T>
struct greater : std::greater<T> {
    VEX_FUNCTION(bool, device, (T, x)(T, y), return x > y;);

    greater() {}
};

/// Function object class for greater-than-or-equal inequality comparison.
template <typename T>
struct greater_equal : std::greater_equal<T> {
    VEX_FUNCTION(bool, device, (T, x)(T, y), return x >= y;);

    greater_equal() {}
};

/// Sorts the vector into ascending order.
template <class K, class Comp>
void sort(K &&keys, Comp comp) {
    detail::sort_sink(detail::forward_as_sequence(keys), comp);
}

/// Sorts the elements in keys and values into ascending key order.
template <class K>
void sort(vector<K> &keys) {
    sort(keys, less<K>());
}

/// Sorts the elements in keys and values into ascending key order.
template <class K, class V, class Comp>
void sort_by_key(K &&keys, V &&vals, Comp comp) {
    detail::sort_by_key_sink(
            detail::forward_as_sequence(keys),
            detail::forward_as_sequence(vals),
            comp);
}

/// Sorts the elements in keys and values into ascending key order.
template <class K, class V>
void sort_by_key(vector<K> &keys, vector<V> &vals) {
    sort_by_key(keys, vals, less<K>());
}

} // namespace vex

#endif
