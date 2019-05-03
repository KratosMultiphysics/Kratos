#ifndef VEXCL_SVM_VECTOR_HPP
#define VEXCL_SVM_VECTOR_HPP

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
 * \file   vexcl/svm_vector.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Shared virtual memory support.
 */

#include <vexcl/operations.hpp>

namespace vex {

struct svm_vector_terminal {};

typedef vector_expression<
    typename boost::proto::terminal< svm_vector_terminal >::type
    > svm_vector_terminal_expression;

namespace traits {

// Hold SVM vector terminals by reference:
template <class T>
struct hold_terminal_by_reference< T,
        typename std::enable_if<
            boost::proto::matches<
                typename boost::proto::result_of::as_expr< T >::type,
                boost::proto::terminal< svm_vector_terminal >
            >::value
        >::type
    >
    : std::true_type
{ };

} // namespace traits
} // namespace vex

#define VEXCL_SVM_ASSIGNMENT(op, op_type)                                      \
  /** Expression assignment operator. */                                       \
  template <class Expr>                                                        \
  auto operator op(const Expr & expr) ->                                       \
      typename std::enable_if<                                                 \
          boost::proto::matches<                                               \
              typename boost::proto::result_of::as_expr<Expr>::type,           \
              vector_expr_grammar>::value,                                     \
          const svm_vector &>::type                                            \
  {                                                                            \
    detail::assign_expression<op_type>(*this, expr);                           \
    return *this;                                                              \
  }

#if defined(VEXCL_BACKEND_OPENCL)
#  include <vexcl/backend/opencl/svm_vector.hpp>
#elif defined(VEXCL_BACKEND_COMPUTE)
#  include <vexcl/backend/compute/svm_vector.hpp>
#elif defined(VEXCL_BACKEND_CUDA)
#  include <vexcl/backend/cuda/svm_vector.hpp>
#else
#  error Unsupported backend!
#endif

#undef VEXCL_SVM_ASSIGNMENT

namespace vex {
namespace traits {

template <> struct is_vector_expr_terminal< svm_vector_terminal > : std::true_type {};
template <> struct proto_terminal_is_value< svm_vector_terminal > : std::true_type {};

template <typename T>
struct kernel_param_declaration< svm_vector<T> >
{
    static void get(backend::source_generator &src,
            const svm_vector<T>&,
            const backend::command_queue&, const std::string &prm_name,
            detail::kernel_generator_state_ptr)
    {
        src.parameter< global_ptr<T> >(prm_name);
    }
};

template <typename T>
struct partial_vector_expr< svm_vector<T> > {
    static void get(backend::source_generator &src,
            const svm_vector<T>&,
            const backend::command_queue&, const std::string &prm_name,
            detail::kernel_generator_state_ptr)
    {
        src << prm_name << "[idx]";
    }
};

template <typename T>
struct kernel_arg_setter< svm_vector<T> >
{
    static void set(const svm_vector<T> &term,
            backend::kernel &kernel, unsigned/*part*/, size_t/*index_offset*/,
            detail::kernel_generator_state_ptr)
    {
        kernel.push_arg(term.get());
    }
};

template <typename T>
struct expression_properties< svm_vector<T> >
{
    static void get(const svm_vector<T> &term,
            std::vector<backend::command_queue> &queue_list,
            std::vector<size_t> &partition,
            size_t &size
            )
    {
        queue_list.clear();
        queue_list.push_back(term.queue());

        size = term.size();

        partition.clear();
        partition.push_back(0);
        partition.push_back(size);
    }
};

} // namespace traits

// --- raw_pointer(svm_vector<T>) ---
template <typename T>
struct svm_vector_pointer {
    typedef T value_type;

    const svm_vector<T> &v;

    svm_vector_pointer(const svm_vector<T> &v) : v(v) {}
};

namespace traits {

template <typename T>
struct is_vector_expr_terminal< svm_vector_pointer<T> > : std::true_type {};

template <typename T>
struct kernel_param_declaration< svm_vector_pointer<T> >
{
    static void get(backend::source_generator &src,
            const svm_vector_pointer<T>&,
            const backend::command_queue&, const std::string &prm_name,
            detail::kernel_generator_state_ptr)
    {
        src.parameter< global_ptr<T> >(prm_name);
    }
};

template <typename T>
struct kernel_arg_setter< svm_vector_pointer<T> >
{
    static void set(const svm_vector_pointer<T> &term,
            backend::kernel &kernel, unsigned/*part*/, size_t/*index_offset*/,
            detail::kernel_generator_state_ptr)
    {
        kernel.push_arg(term.v.get());
    }
};

template <typename T>
struct expression_properties< svm_vector_pointer<T> >
{
    static void get(const svm_vector_pointer<T> &term,
            std::vector<backend::command_queue> &queue_list,
            std::vector<size_t> &partition,
            size_t &size
            )
    {
        detail::get_expression_properties prop;
        detail::extract_terminals()(boost::proto::as_child(term.v), prop);

        queue_list = prop.queue;
        partition  = prop.part;
        size       = prop.size;
    }
};

} // namespace traits

/// Cast vex::svm_vector to a raw pointer.
template <typename T>
inline auto
raw_pointer(const svm_vector<T> &v) ->
    typename boost::proto::result_of::as_expr<svm_vector_pointer<T>, vector_domain>::type
{
    return boost::proto::as_expr<vector_domain>(svm_vector_pointer<T>(v));
}

} // namespace vex



#endif
