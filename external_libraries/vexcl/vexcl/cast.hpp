#ifndef VEXCL_CAST_HPP
#define VEXCL_CAST_HPP

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
 * \file   vexcl/cast.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Cast an expression to a type.
 */

#include <vexcl/operations.hpp>

namespace vex {

struct cast_terminal {};

typedef vector_expression<
    typename boost::proto::terminal< cast_terminal >::type
    > cast_terminal_expression;

template <typename T, class Expr>
struct casted_expession : public cast_terminal_expression
{
    typedef T value_type;
    const Expr expr;

    casted_expession(const Expr &expr) : expr(expr) {}
};

/// Cast an expression to a given type.
template <typename T, class Expr>
auto cast(const Expr &expr) ->
    typename std::enable_if<
        boost::proto::matches<
                typename boost::proto::result_of::as_expr< Expr >::type,
                vector_expr_grammar
        >::value,
        casted_expession<T, typename boost::proto::result_of::as_child<const Expr, vector_domain>::type
        >
    >::type
{
    return casted_expession<
                T,
                typename boost::proto::result_of::as_child<const Expr, vector_domain>::type
            >(boost::proto::as_child<vector_domain>(expr));
}

namespace traits {

template <>
struct is_vector_expr_terminal< cast_terminal > : std::true_type {};

template <>
struct proto_terminal_is_value< cast_terminal > : std::true_type {};

template <typename T, class Expr>
struct terminal_preamble< casted_expession<T, Expr> > {
    static void get(backend::source_generator &src,
            const casted_expession<T, Expr> &term,
            const backend::command_queue &queue, const std::string &prm_name,
            detail::kernel_generator_state_ptr state)
    {
        detail::output_terminal_preamble termpream(src, queue, prm_name, state);
        boost::proto::eval(boost::proto::as_child(term.expr), termpream);
    }
};

template <typename T, class Expr>
struct kernel_param_declaration< casted_expession<T, Expr> > {
    static void get(backend::source_generator &src,
            const casted_expession<T, Expr> &term,
            const backend::command_queue &queue, const std::string &prm_name,
            detail::kernel_generator_state_ptr state)
    {
        detail::declare_expression_parameter declare(src, queue, prm_name, state);
        detail::extract_terminals()(boost::proto::as_child(term.expr),  declare);
    }
};

template <typename T, class Expr>
struct local_terminal_init< casted_expession<T, Expr> > {
    static void get(backend::source_generator &src,
            const casted_expession<T, Expr> &term,
            const backend::command_queue &queue, const std::string &prm_name,
            detail::kernel_generator_state_ptr state)
    {
        detail::output_local_preamble init_ctx(src, queue, prm_name, state);
        boost::proto::eval(boost::proto::as_child(term.expr), init_ctx);
    }
};

template <typename T, class Expr>
struct partial_vector_expr< casted_expession<T, Expr> > {
    static void get(backend::source_generator &src,
            const casted_expession<T, Expr> &term,
            const backend::command_queue &queue, const std::string &prm_name,
            detail::kernel_generator_state_ptr state)
    {
        detail::vector_expr_context expr_ctx(src, queue, prm_name, state);
        boost::proto::eval(boost::proto::as_child(term.expr), expr_ctx);
    }
};

template <typename T, class Expr>
struct kernel_arg_setter< casted_expession<T, Expr> > {
    static void set(const casted_expession<T, Expr> &term,
            backend::kernel &kernel, unsigned device, size_t index_offset,
            detail::kernel_generator_state_ptr state)
    {
        detail::set_expression_argument setarg(kernel, device, index_offset, state);
        detail::extract_terminals()( boost::proto::as_child(term.expr), setarg);
    }
};

template <typename T, class Expr>
struct expression_properties< casted_expession<T, Expr> > {
    static void get(const casted_expession<T, Expr> &term,
            std::vector<backend::command_queue> &queue_list,
            std::vector<size_t> &partition,
            size_t &size
            )
    {
        detail::get_expression_properties prop;
        detail::extract_terminals()(boost::proto::as_child(term.expr), prop);

        queue_list = prop.queue;
        partition  = prop.part;
        size       = prop.size;
    }
};

} // namespace traits

#if !defined(VEXCL_BACKEND_CUDA)

#define VEXCL_CONVERT_FUNCTIONS(to)                                            \
  struct convert_##to##_func : builtin_function {                              \
    static const char *name() { return "convert_" #to; }                       \
  };                                                                           \
  template <typename Arg>                                                      \
  auto convert_##to(const Arg &arg)                                            \
      ->casted_expession<cl_##to,                                              \
             typename boost::proto::result_of::make_expr<                      \
                 boost::proto::tag::function, convert_##to##_func,             \
                 const Arg &>::type> const                                     \
  {                                                                            \
    return cast<cl_##to>(boost::proto::make_expr<boost::proto::tag::function>( \
        convert_##to##_func(), boost::ref(arg)));                              \
  }                                                                            \
  struct as_##to##_func : builtin_function {                                   \
    static const char *name() { return "as_" #to; }                            \
  };                                                                           \
  template <typename Arg>                                                      \
  auto as_##to(const Arg &arg) ->                                              \
      casted_expession<cl_##to, typename boost::proto::result_of::make_expr<   \
                                boost::proto::tag::function, as_##to##_func,   \
                                const Arg &>::type> const                      \
  {                                                                            \
    return cast<cl_##to>(boost::proto::make_expr<boost::proto::tag::function>( \
        as_##to##_func(), boost::ref(arg)));                                   \
  }

#define VEXCL_TYPES(name)                                                      \
  VEXCL_CONVERT_FUNCTIONS(name)                                                \
  VEXCL_CONVERT_FUNCTIONS(name##2)                                             \
  VEXCL_CONVERT_FUNCTIONS(name##4)                                             \
  VEXCL_CONVERT_FUNCTIONS(name##8)                                             \
  VEXCL_CONVERT_FUNCTIONS(name##16)

VEXCL_TYPES(float)
VEXCL_TYPES(double)
VEXCL_TYPES(char)
VEXCL_TYPES(uchar)
VEXCL_TYPES(short)
VEXCL_TYPES(ushort)
VEXCL_TYPES(int)
VEXCL_TYPES(uint)
VEXCL_TYPES(long)
VEXCL_TYPES(ulong)

#undef VEXCL_TYPES
#undef VEXCL_CONVERT_FUNCTION

#endif

} // namespace vex

#endif
