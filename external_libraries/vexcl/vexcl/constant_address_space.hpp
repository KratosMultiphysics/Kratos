#ifndef VEXCL_CONSTANT_ADDRESS_SPACE_HPP
#define VEXCL_CONSTANT_ADDRESS_SPACE_HPP
/*
The MIT License

Copyright (c) 2012-2018 Denis Demidov <dennis.demidov@gmail.com>
Copyright (c) 2015 Boris Smidt <smidtboris1@gmail.com>

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
 * \file   vexcl/constant_address_space.hpp
 * \author Boris Smidt <smidtboris1@gmail.com>
 * \brief  Functions enabling use of constant cache in OpenCL backend.
 */

#if defined(VEXCL_BACKEND_CUDA)
#  error Constant address space is not supported by CUDA backend!
#endif

#include <vexcl/operations.hpp>

namespace vex {

// Will be used as an identificator for constant vector expressions.
struct constant_vector_terminal {};

// Introduce the new kind of expression to Boost.Proto.
typedef vector_expression<
    typename boost::proto::terminal< constant_vector_terminal >::type
    > constant_vector_expression;

// The actual terminal.
// This struct will hold a reference to a vector that needs to be cached in
// constant cache.
template <class T>
struct constant_vector : constant_vector_expression {
    typedef T value_type;
    const vector<T> &v;
    constant_vector(const vector<T> &v) : v(v) {}
};

// Code generation helpers.
namespace traits {

// Tell proto it is ok to use these terminals in vexcl expressions.
template <>
struct is_vector_expr_terminal< constant_vector_terminal > : std::true_type {};

// Tell proto that it is not required to unpack the terminal with
// proto::value() function.
template <>
struct proto_terminal_is_value< constant_vector_terminal > : std::true_type {};

// Declare kernel parameter for the constant vector.
// This is the only specification that is different from vector<T>.
template <typename T>
struct kernel_param_declaration< constant_vector<T> > {
    static void get(backend::source_generator &src,
            const constant_vector<T>&,
            const backend::command_queue&, const std::string &prm_name,
            detail::kernel_generator_state_ptr)
    {
        src.parameter< constant_ptr<T> >( prm_name);
    }
};

// The following specifications just delegate the work to vector<T>:
template <typename T>
struct partial_vector_expr< constant_vector<T> > {
    static void get(backend::source_generator &src,
            const constant_vector<T> &rv,
            const backend::command_queue &q, const std::string &prm_name,
            detail::kernel_generator_state_ptr state)
    {
        get_partial_vector_expr(src, rv.v, q, prm_name, state);
    }
};

template <typename T>
struct kernel_arg_setter< constant_vector<T> > {
    static void set(const constant_vector<T> &rv,
            backend::kernel &kernel, unsigned device, size_t index_offset,
            detail::kernel_generator_state_ptr state)
    {
        set_kernel_args(rv.v, kernel, device, index_offset, state);
    }
};

template <class T>
struct expression_properties< constant_vector<T> > {
    static void get(const constant_vector<T> &rv,
            std::vector<backend::command_queue> &queue_list,
            std::vector<size_t> &partition,
            size_t &size
            )
    {
        extract_expression_properties(rv.v, queue_list, partition, size);
    }
};

} // namespace traits

/// Uses constant cache for access to the wrapped vector.
/**
\rst
.. note::
    Only available for OpenCL-based backends.
\endrst
 */
template <class T>
inline constant_vector<T> constant(const vector<T> &v) {
    return constant_vector<T>(v);
}

} // namespace vex

#endif // VEXCL_CONSTANT_ADDRESS_SPACE_HPP
