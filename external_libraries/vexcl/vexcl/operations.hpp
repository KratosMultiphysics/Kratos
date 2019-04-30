#ifndef VEXCL_OPERATIONS_HPP
#define VEXCL_OPERATIONS_HPP

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
 * \file   vexcl/operations.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Templates used for expression tree traversal and kernel generation.
 */

#include <array>
#include <tuple>
#include <deque>
#include <set>
#include <memory>

#include <boost/proto/proto.hpp>
#include <boost/mpl/max.hpp>
#include <boost/any.hpp>

#include <vexcl/backend.hpp>
#include <vexcl/types.hpp>
#include <vexcl/util.hpp>
#include <vexcl/cache.hpp>

// Workaround for gcc bug http://gcc.gnu.org/bugzilla/show_bug.cgi?id=35722
#if defined(BOOST_NO_VARIADIC_TEMPLATES) || (defined(__GNUC__) && !defined(__clang__) && __GNUC__ == 4 && __GNUC_MINOR__ == 6)
#  include <boost/preprocessor/repetition.hpp>
#  ifndef VEXCL_MAX_ARITY
#    define VEXCL_MAX_ARITY BOOST_PROTO_MAX_ARITY
#  endif
#endif

/// Vector expression template library for OpenCL.
namespace vex {

//---------------------------------------------------------------------------
// Assignment operators.
//---------------------------------------------------------------------------
namespace assign {

#define VEXCL_ASSIGN_OP(name, op)                                              \
  struct name {                                                                \
    static std::string string() { return #op; }                                \
  };

    VEXCL_ASSIGN_OP(SET, =)
    VEXCL_ASSIGN_OP(ADD, +=)
    VEXCL_ASSIGN_OP(SUB, -=)
    VEXCL_ASSIGN_OP(MUL, *=)
    VEXCL_ASSIGN_OP(DIV, /=)
    VEXCL_ASSIGN_OP(MOD, %=)
    VEXCL_ASSIGN_OP(AND, &=)
    VEXCL_ASSIGN_OP(OR,  |=)
    VEXCL_ASSIGN_OP(XOR, ^=)
    VEXCL_ASSIGN_OP(LSH, <<=)
    VEXCL_ASSIGN_OP(RSH, >>=)

#undef VEXCL_ASSIGN_OP

#define VEXCL_ASSIGNMENTS(ASSIGNMENT_MACRO)                                    \
    ASSIGNMENT_MACRO(=,   assign::SET)                                         \
    ASSIGNMENT_MACRO(+=,  assign::ADD)                                         \
    ASSIGNMENT_MACRO(-=,  assign::SUB)                                         \
    ASSIGNMENT_MACRO(*=,  assign::MUL)                                         \
    ASSIGNMENT_MACRO(/=,  assign::DIV)                                         \
    ASSIGNMENT_MACRO(%=,  assign::MOD)                                         \
    ASSIGNMENT_MACRO(&=,  assign::AND)                                         \
    ASSIGNMENT_MACRO(|=,  assign::OR)                                          \
    ASSIGNMENT_MACRO(^=,  assign::XOR)                                         \
    ASSIGNMENT_MACRO(<<=, assign::LSH)                                         \
    ASSIGNMENT_MACRO(>>=, assign::RSH)
}

namespace detail {

// Used as a state parameter in kernel generation functions.
typedef std::map<std::string, boost::any> kernel_generator_state;
typedef std::shared_ptr<kernel_generator_state> kernel_generator_state_ptr;

inline kernel_generator_state_ptr empty_state() {
    return std::make_shared<kernel_generator_state>();
}

} // namespace detail

namespace traits {

// Terminals allowed in vector expressions
template <class Term, class Enable = void>
struct is_vector_expr_terminal : std::false_type { };

template <class T>
struct is_vector_expr_terminal< T,
    typename std::enable_if< is_cl_native< T >::value >::type
    > : std::true_type
{ };

// Hold everything by value inside proto expressions unless explicitly
// specified otherwise.
template <class T, class Enable = void>
struct hold_terminal_by_reference : std::false_type {};

// Value type of a terminal
template <class T, class Enable = void>
struct value_type { typedef T type; };

// If a terminal has typedef'ed value_type, then use it:
template <class T>
struct value_type<T,
    typename std::enable_if<
        std::is_same<typename T::value_type, typename T::value_type>::value
    >::type>
{
    typedef typename T::value_type type;
};


//---------------------------------------------------------------------------
// Kernel source generation
//---------------------------------------------------------------------------

// Some terminals need preamble (e.g. struct declaration or helper function).
// But most of them do not:
template <class T>
struct terminal_preamble {
    static void get(backend::source_generator&, const T&,
            const backend::command_queue&, const std::string &/*prm_name*/,
            detail::kernel_generator_state_ptr)
    { }
};

template <class T>
void get_terminal_preamble(backend::source_generator &src,
        const T &term, const backend::command_queue &queue, const std::string &prm_name,
        detail::kernel_generator_state_ptr state)
{
    terminal_preamble<
        typename std::decay<T>::type
    >::get(src, term, queue, prm_name, state);
}

// How to declare OpenCL kernel parameters for a terminal:
template <class Term, class Enable = void>
struct kernel_param_declaration {
    static void get(backend::source_generator &src, const Term&,
            const backend::command_queue&, const std::string &prm_name,
            detail::kernel_generator_state_ptr)
    {
        src.parameter<Term>(prm_name);
    }
};

template <class T>
void get_kernel_param_declaration(backend::source_generator &src,
        const T &term, const backend::command_queue &queue, const std::string &prm_name,
        detail::kernel_generator_state_ptr state)
{
    kernel_param_declaration<
        typename std::decay<T>::type
    >::get(src, term, queue, prm_name, state);
}

// Local terminal initialization (e.g. temporary declaration)
template <class Term, class Enable = void>
struct local_terminal_init {
    static void get(backend::source_generator&, const Term&,
            const backend::command_queue&, const std::string &/*prm_name*/,
            detail::kernel_generator_state_ptr)
    { }
};

template <class T>
void get_local_terminal_init(backend::source_generator &src,
        const T &term, const backend::command_queue &queue, const std::string &prm_name,
        detail::kernel_generator_state_ptr state)
{
    local_terminal_init<
        typename std::decay<T>::type
    >::get(src, term, queue, prm_name, state);
}

// Partial expression for a terminal:
template <class Term, class Enable = void>
struct partial_vector_expr {
    static void get(backend::source_generator &src, const Term&,
            const backend::command_queue&, const std::string &prm_name,
            detail::kernel_generator_state_ptr)
    {
        src << prm_name;
    }
};

template <class T>
void get_partial_vector_expr(backend::source_generator &src,
        const T &term, const backend::command_queue &queue, const std::string &prm_name,
        detail::kernel_generator_state_ptr state)
{
    partial_vector_expr<
        typename std::decay<T>::type
    >::get(src, term, queue, prm_name, state);
}

// How to set OpenCL kernel arguments for a terminal:
template <class Term, class Enable = void>
struct kernel_arg_setter {
    static void set(const Term &term,
            backend::kernel &kernel, unsigned/*device*/, size_t/*index_offset*/,
            detail::kernel_generator_state_ptr)
    {
        kernel.push_arg(term);
    }
};

template <class T>
void set_kernel_args(
        const T &term, backend::kernel &kernel, unsigned device, size_t index_offset,
        detail::kernel_generator_state_ptr state)
{
    kernel_arg_setter<
        typename std::decay<T>::type
    >::set(term, kernel, device, index_offset, state);
}

// How to deduce queue list, partitioning and size from a terminal:
template <class T, class Enable = void>
struct expression_properties {
    static void get(const T &/*term*/,
            std::vector<backend::command_queue> &/*queue_list*/,
            std::vector<size_t> &/*partition*/,
            size_t &/*size*/
            )
    { }
};

template <class T>
void extract_expression_properties(
        const T &term,
        std::vector<backend::command_queue> &queue_list,
        std::vector<size_t> &partition,
        size_t &size
        )
{
    expression_properties<
        typename std::decay<T>::type
    >::get(term, queue_list, partition, size);
}

//---------------------------------------------------------------------------
// Scalars and helper types/functions used in multivector expressions
//---------------------------------------------------------------------------
template <class T, class Enable = void>
struct is_multiscalar : std::false_type {};

// Arithmetic scalars
template <class T>
struct is_multiscalar< T,
    typename std::enable_if< is_cl_native<T>::value >::type
    > : std::true_type
{};

// Number of components in a multivector expression terminal.
template <class T>
struct number_of_components : boost::mpl::size_t<1> {};

// Type of I-th component of a multivector expression terminal.
template <size_t I, class T, class Enable = void>
struct component { typedef T type; };

#ifndef BOOST_NO_VARIADIC_TEMPLATES
template <typename... T>
struct And : std::true_type {};

template <typename Head, typename... Tail>
struct And<Head, Tail...>
    : std::conditional<Head::value, And<Tail...>, std::false_type>::type
{};
#endif

// std::tuple<...>
#ifndef BOOST_NO_VARIADIC_TEMPLATES
template <class... Args>
struct is_multiscalar<std::tuple<Args...>,
    typename std::enable_if<And< is_cl_native<Args>... >::type::value >::type
    > : std::true_type
{};

template <class... Args>
struct number_of_components< std::tuple<Args...> >
    : boost::mpl::size_t<sizeof...(Args)>
{};

template <size_t I, class... Args>
struct component< I, std::tuple<Args...> >
    : std::tuple_element< I, std::tuple<Args...> >
{};

#else

#define VEXCL_TUPLE_IS_MS(z, n, unused)                                        \
  template <BOOST_PP_ENUM_PARAMS(n, class Arg)>                                \
  struct is_multiscalar<                                                       \
      std::tuple<BOOST_PP_ENUM_PARAMS(n, Arg)> > : std::true_type {            \
  };

#define VEXCL_TUPLE_COMP(z, n, unused)                                         \
  template <size_t I, BOOST_PP_ENUM_PARAMS(n, class Arg)>                      \
  struct component<                                                            \
      I, std::tuple<BOOST_PP_ENUM_PARAMS(n, Arg)> > : std::tuple_element<      \
      I, std::tuple<BOOST_PP_ENUM_PARAMS(n, Arg)> > {                          \
  };

BOOST_PP_REPEAT_FROM_TO(1, VEXCL_MAX_ARITY, VEXCL_TUPLE_IS_MS, ~)
BOOST_PP_REPEAT_FROM_TO(1, VEXCL_MAX_ARITY, VEXCL_TUPLE_COMP, ~)

#undef VEXCL_TUPLE_IS_MS
#undef VEXCL_TUPLE_COMP

#endif

// std::array<T,N>
template <class T, size_t N>
struct is_multiscalar< std::array<T, N>,
    typename std::enable_if< is_cl_native<T>::value >::type
    > : std::true_type
{};

template <class T, size_t N>
struct number_of_components< std::array<T, N> > : boost::mpl::size_t<N> {};

template <size_t I, class T, size_t N>
struct component< I, std::array<T, N> > { typedef T type; };

// C-style arrays
template <class T, size_t N>
struct is_multiscalar< T[N],
    typename std::enable_if< is_cl_native<T>::value >::type
    > : std::true_type
{};

template <class T, size_t N>
struct number_of_components< T[N] > : boost::mpl::size_t<N> {};

template <size_t I, class T, size_t N>
struct component< I, T[N] > { typedef T type; };

// Terminals allowed in multivector expressions.
template <class Term, class Enable = void>
struct is_multivector_expr_terminal : std::false_type { };

template <class T>
struct is_multivector_expr_terminal< T,
    typename std::enable_if< is_multiscalar<T>::value >::type
    >
    : std::true_type
{ };

// Extract component directly from terminal rather than from value(terminal):
template <class T, class Enable = void>
struct proto_terminal_is_value : std::false_type { };

template <class T>
struct terminal_is_value :
    boost::proto::matches<
            typename boost::proto::result_of::as_expr<T>::type,
            boost::proto::and_<
                boost::proto::terminal< boost::proto::_ >,
                boost::proto::if_< proto_terminal_is_value< boost::proto::_value >() >
            >
    >
{};

/* Type trait to determine if an expression is scalable.
 *
 * The expression should have a type `value_type` and a field `scale` of that
 * type, this enables operator* and operator/.
 */
template <class T> struct is_scalable : std::false_type {};

} // namespace traits

//---------------------------------------------------------------------------
// Extracting components from multivector expression terminals
//---------------------------------------------------------------------------
template <size_t I, typename T, size_t N>
inline T& get(T t[N]) {
    static_assert(I < N, "Component number out of bounds");
    return t[I];
}

template <size_t I, typename T, size_t N>
inline const T& get(const T t[N]) {
    static_assert(I < N, "Component number out of bounds");
    return t[I];
}

template <size_t I, typename T>
inline T& get(T &t) {
    return t;
}

// Scalable expressions may be multiplied by a scalar:
template <class T>
typename std::enable_if<vex::traits::is_scalable<T>::value, T>::type
operator*(const T &expr, const typename T::value_type &factor) {
    T scaled_expr(expr);
    scaled_expr.scale *= factor;
    return scaled_expr;
}

// Scalable expressions may be multiplied by a scalar:
template <class T>
typename std::enable_if<vex::traits::is_scalable<T>::value, T>::type
operator*(const typename T::value_type &factor, const T &expr) {
    return expr * factor;
}

// Scalable expressions may be divided by a scalar:
template <class T> typename std::enable_if<vex::traits::is_scalable<T>::value, T>::type
operator/(const T &expr, const typename T::value_type &factor) {
    T scaled_expr(expr);
    scaled_expr.scale /= factor;
    return scaled_expr;
}

//---------------------------------------------------------------------------
// Standard grammar (no terminals)
//---------------------------------------------------------------------------
struct builtin_function {};
struct user_function {};

#define VEXCL_BUILTIN_OPERATIONS(grammar)                                      \
    boost::proto::or_<                                                         \
        boost::proto::unary_plus< grammar >,                                   \
        boost::proto::negate< grammar >,                                       \
        boost::proto::logical_not< grammar >,                                  \
        boost::proto::pre_inc< grammar >,                                      \
        boost::proto::pre_dec< grammar >,                                      \
        boost::proto::post_inc< grammar >,                                     \
        boost::proto::post_dec< grammar >                                      \
    >,                                                                         \
    boost::proto::or_<                                                         \
        boost::proto::or_<                                                     \
            boost::proto::plus          < grammar, grammar >,                  \
            boost::proto::minus         < grammar, grammar >,                  \
            boost::proto::multiplies    < grammar, grammar >,                  \
            boost::proto::divides       < grammar, grammar >,                  \
            boost::proto::modulus       < grammar, grammar >,                  \
            boost::proto::shift_left    < grammar, grammar >,                  \
            boost::proto::shift_right   < grammar, grammar >                   \
        >,                                                                     \
        boost::proto::or_<                                                     \
            boost::proto::less          < grammar, grammar >,                  \
            boost::proto::greater       < grammar, grammar >,                  \
            boost::proto::less_equal    < grammar, grammar >,                  \
            boost::proto::greater_equal < grammar, grammar >,                  \
            boost::proto::equal_to      < grammar, grammar >,                  \
            boost::proto::not_equal_to  < grammar, grammar >                   \
        >,                                                                     \
        boost::proto::or_<                                                     \
            boost::proto::logical_and   < grammar, grammar >,                  \
            boost::proto::logical_or    < grammar, grammar >                   \
        >,                                                                     \
        boost::proto::or_<                                                     \
            boost::proto::bitwise_and   < grammar, grammar >,                  \
            boost::proto::bitwise_or    < grammar, grammar >,                  \
            boost::proto::bitwise_xor   < grammar, grammar >                   \
        >,                                                                     \
        boost::proto::or_<                                                     \
            boost::proto::address_of < grammar >,                              \
            boost::proto::dereference< grammar >,                              \
            boost::proto::subscript  < grammar, grammar >                      \
        >,                                                                     \
        boost::proto::or_<                                                     \
            boost::proto::if_else_< grammar, grammar, grammar >                \
        >                                                                      \
    >,                                                                         \
    boost::proto::function<                                                    \
        boost::proto::terminal<                                                \
            boost::proto::convertible_to<builtin_function>                     \
        >,                                                                     \
        boost::proto::vararg<grammar>                                          \
    >

#define VEXCL_USER_FUNCTIONS(grammar)                                          \
  boost::proto::function<                                                      \
      boost::proto::terminal<boost::proto::convertible_to<user_function> >,    \
      boost::proto::vararg<grammar>                                            \
  >

#define VEXCL_VECTOR_EXPR_EXTRACTOR(name, VG, AG, FG)                          \
struct name                                                                    \
    : boost::proto::or_<                                                       \
          boost::proto::when<VG, boost::proto::_>,                             \
          boost::proto::when<                                                  \
             boost::proto::plus< FG, AG >,                                     \
             name(boost::proto::_left)                                         \
          >,                                                                   \
          boost::proto::when<                                                  \
             boost::proto::plus< AG, FG >,                                     \
             name(boost::proto::_right)                                        \
          >,                                                                   \
          boost::proto::when<                                                  \
             boost::proto::minus< FG, AG >,                                    \
             name(boost::proto::_left)                                         \
          >,                                                                   \
          boost::proto::when<                                                  \
             boost::proto::minus< AG, FG >,                                    \
             boost::proto::_make_negate( name(boost::proto::_right) )          \
          >,                                                                   \
          boost::proto::when<                                                  \
             boost::proto::binary_expr<boost::proto::_, name, name >           \
          >                                                                    \
      >                                                                        \
{}

#define VEXCL_ADDITIVE_EXPR_EXTRACTOR(name, VG, AG, FG)                        \
struct name                                                                    \
    : boost::proto::or_<                                                       \
          boost::proto::when<AG, boost::proto::_>                              \
        , boost::proto::when<                                                  \
             boost::proto::plus<FG, VG >,                                      \
             name(boost::proto::_left)                                         \
          >                                                                    \
        , boost::proto::when<                                                  \
             boost::proto::plus< VG, FG >,                                     \
             name(boost::proto::_right)                                        \
          >                                                                    \
        , boost::proto::when<                                                  \
             boost::proto::minus<FG, VG >,                                     \
             name(boost::proto::_left)                                         \
          >                                                                    \
        , boost::proto::when<                                                  \
             boost::proto::minus<VG, FG >,                                     \
             boost::proto::_make_negate( name(boost::proto::_right) )          \
          >                                                                    \
        , boost::proto::when<                                                  \
             boost::proto::binary_expr<boost::proto::_, name, name >           \
          >                                                                    \
      >                                                                        \
{}

//---------------------------------------------------------------------------
// User-defined functions
//---------------------------------------------------------------------------
template <class C, class T>
struct UserFunction {};

// Workaround for gcc bug http://gcc.gnu.org/bugzilla/show_bug.cgi?id=35722
#if !defined(BOOST_NO_VARIADIC_TEMPLATES) && ((!defined(__GNUC__) || (__GNUC__ > 4 || __GNUC__ == 4 && __GNUC_MINOR__ > 6)) || defined(__clang__))

template<class Impl, class RetType, class... ArgType>
struct UserFunction<Impl, RetType(ArgType...)> : user_function
{
    UserFunction() {}

    typedef RetType value_type;

    template <class... Arg>
    typename boost::proto::result_of::make_expr<
        boost::proto::tag::function,
        Impl,
        const Arg&...
    >::type const
    operator()(const Arg&... arg) const {
        static_assert(sizeof...(Arg) == sizeof...(ArgType),
                "Wrong number of arguments for a user-defined function!"
                );
        return boost::proto::make_expr<boost::proto::tag::function>(
                Impl(), boost::ref(arg)...
                );
    }

    static std::string preamble() { return ""; }

    static std::string name() {
        return "user_function";
    }

    static void define(backend::source_generator &src) {
        define(src, Impl::name());
    }

    static void define(backend::source_generator &src, const std::string &name)
    {
        src << Impl::preamble();
        src.begin_function<RetType>(name);
        src.begin_function_parameters();
        show_arg<ArgType...>(src, 1);
        src.end_function_parameters();
        src.new_line() << Impl::body();
        src.end_function();
    }

    template <class Head>
    static void show_arg(backend::source_generator &src, unsigned pos) {
        src.parameter<Head>("prm") << pos;
    }

    template <class Head, class... Tail>
    static typename std::enable_if<(sizeof...(Tail) > 0), void>::type
    show_arg(backend::source_generator &src, unsigned pos) {
        show_arg<Tail...>(src.parameter<Head>("prm") << pos, pos + 1);
    }
};

#else

#define VEXCL_PRINT_PRM_DEF(z, n, data) src.parameter<ArgType##n>("prm") << n + 1;
#define VEXCL_PRINT_BOOST_REF(z, n, data) boost::ref(arg##n)

#define VEXCL_USER_FUNCTION(z, n, data)                                        \
  template <class Impl, class RetType, BOOST_PP_ENUM_PARAMS(n, class ArgType)> \
  struct UserFunction<                                                         \
      Impl, RetType(BOOST_PP_ENUM_PARAMS(n, ArgType))> : user_function {       \
    UserFunction() {}                                                          \
    static std::string name() {                                                \
        return "user_function";                                                \
    }                                                                          \
    static void define(backend::source_generator &src) {                       \
        define(src, Impl::name());                                             \
    }                                                                          \
    typedef RetType value_type;                                                \
    template <BOOST_PP_ENUM_PARAMS(n, class Arg)>                              \
    typename boost::proto::result_of::make_expr<                               \
        boost::proto::tag::function, Impl,                                     \
        BOOST_PP_ENUM_BINARY_PARAMS(                                           \
                n, const Arg,& BOOST_PP_INTERCEPT)>::type const                \
    operator()(                                                                \
        BOOST_PP_ENUM_BINARY_PARAMS(n, const Arg, &arg)) const {               \
      return boost::proto::make_expr<boost::proto::tag::function>(             \
          Impl(), BOOST_PP_ENUM(n, VEXCL_PRINT_BOOST_REF, ~));                 \
    }                                                                          \
    static std::string preamble() { return ""; }                               \
    static void define(backend::source_generator &src,                         \
                       const std::string &name) {                              \
      src << Impl::preamble();                                                 \
      src.begin_function<RetType>(name);                                       \
      src.begin_function_parameters();                                         \
      BOOST_PP_REPEAT(n, VEXCL_PRINT_PRM_DEF, n)                               \
      src.end_function_parameters();                                           \
      src.new_line() << Impl::body();                                          \
      src.end_function();                                                      \
    }                                                                          \
  };

BOOST_PP_REPEAT_FROM_TO(1, VEXCL_MAX_ARITY, VEXCL_USER_FUNCTION, ~)

#undef VEXCL_PRINT_PRM_DEF
#undef VEXCL_PRINT_BOOST_REF
#undef VEXCL_USER_FUNCTION

#endif

//---------------------------------------------------------------------------
// Vector grammar
//---------------------------------------------------------------------------
// Grammar for vector expressions that may be processed with single kernel:
struct vector_expr_grammar
    : boost::proto::or_<
          boost::proto::and_<
              boost::proto::terminal< boost::proto::_ >,
              boost::proto::if_< traits::is_vector_expr_terminal< boost::proto::_value >() >
          >,
          VEXCL_BUILTIN_OPERATIONS(vector_expr_grammar),
          VEXCL_USER_FUNCTIONS(vector_expr_grammar)
      >
{};

// Grammar for additive expressions
// (each additive term requires separate kernel):
struct additive_vector_transform {};

struct additive_vector_transform_grammar
    : boost::proto::or_<
        boost::proto::terminal< additive_vector_transform >,
        boost::proto::plus<
            additive_vector_transform_grammar,
            additive_vector_transform_grammar
        >,
        boost::proto::minus<
            additive_vector_transform_grammar,
            additive_vector_transform_grammar
        >,
        boost::proto::negate<
            additive_vector_transform_grammar
        >
      >
{};

// Joined grammar for vector expressions and additive expressions.
struct vector_full_grammar
    : boost::proto::or_<
        vector_expr_grammar,
        boost::proto::terminal< additive_vector_transform >,
        boost::proto::plus< vector_full_grammar, vector_full_grammar >,
        boost::proto::minus< vector_full_grammar, vector_full_grammar >,
        boost::proto::negate< vector_full_grammar >
      >
{};


// Boost.Proto domain for vector expressions.
template <class Expr>
struct vector_expression;

struct vector_domain
    : boost::proto::domain<
        boost::proto::generator<vector_expression>,
        vector_full_grammar
        >
{
    // Store everything by value inside expressions...
    template <typename T, class Enable = void>
    struct as_child : proto_base_domain::as_expr<T>
    {};

    // ... except for terminals that explicitly request storage by reference:
    template <typename T>
    struct as_child< T,
        typename std::enable_if<
                traits::hold_terminal_by_reference<T>::value
            >::type
        > : proto_base_domain::as_child< T >
    {};
};

template <class Expr>
struct vector_expression
    : boost::proto::extends< Expr, vector_expression<Expr>, vector_domain>
{
    vector_expression(const Expr &expr = Expr())
        : boost::proto::extends< Expr, vector_expression<Expr>, vector_domain>(expr) {}
};

template <class M, class V>
struct additive_operator
    : vector_expression< boost::proto::terminal< additive_vector_transform >::type >
{
    typedef typename V::value_type value_type;

    const M &A;
    const V &x;

    typename cl_scalar_of<value_type>::type scale;

    additive_operator(const M &A, const V &x) : A(A), x(x), scale(1) {}

    template<bool negate, bool append>
    void apply(V &y) const {
        A.apply(x, y, negate ? -scale : scale, append);
    }
};


namespace traits {

template <class M, class V>
struct is_scalable< additive_operator<M, V> > : std::true_type {};

} // namespace traits

//---------------------------------------------------------------------------
// Multivector grammar
//---------------------------------------------------------------------------
struct multivector_expr_grammar
    : boost::proto::or_<
          boost::proto::and_<
              boost::proto::terminal< boost::proto::_ >,
              boost::proto::if_< traits::is_multivector_expr_terminal< boost::proto::_value >() >
          >,
          VEXCL_BUILTIN_OPERATIONS(multivector_expr_grammar),
          VEXCL_USER_FUNCTIONS(multivector_expr_grammar)
      >
{};

struct additive_multivector_transform {};

struct additive_multivector_transform_grammar
    : boost::proto::or_<
        boost::proto::terminal< additive_multivector_transform >,
        boost::proto::plus<
            additive_multivector_transform_grammar,
            additive_multivector_transform_grammar
        >,
        boost::proto::minus<
            additive_multivector_transform_grammar,
            additive_multivector_transform_grammar
        >,
        boost::proto::negate<
            additive_multivector_transform_grammar
        >
      >
{};

struct multivector_full_grammar
    : boost::proto::or_<
        multivector_expr_grammar,
        boost::proto::terminal< additive_multivector_transform >,
        boost::proto::plus< multivector_full_grammar, multivector_full_grammar >,
        boost::proto::minus< multivector_full_grammar, multivector_full_grammar >,
        boost::proto::negate< multivector_full_grammar >
      >
{};

template <class Expr>
struct multivector_expression;

struct multivector_domain
    : boost::proto::domain<
        boost::proto::generator<multivector_expression>,
        multivector_full_grammar
      >
{
    // Store everything by value inside expressions...
    template <typename T, class Enable = void>
    struct as_child : proto_base_domain::as_expr<T>
    {};

    // ... except for terminals that explicitly request storage by reference:
    template <typename T>
    struct as_child< T,
        typename std::enable_if<
                traits::hold_terminal_by_reference<T>::value
            >::type
        > : proto_base_domain::as_child< T >
    {};
};

template <class Expr>
struct multivector_expression
    : boost::proto::extends< Expr, multivector_expression<Expr>, multivector_domain>
{
    multivector_expression(const Expr &expr = Expr())
        : boost::proto::extends< Expr, multivector_expression<Expr>, multivector_domain>(expr) {}
};

template <class M, class V>
struct multiadditive_operator
    : multivector_expression<
        boost::proto::terminal< additive_multivector_transform >::type
        >
{
    typedef typename V::sub_value_type value_type;

    const M &A;
    const V &x;

    typename cl_scalar_of<value_type>::type scale;

    multiadditive_operator(const M &A, const V &x) : A(A), x(x), scale(1) {}

    template <bool negate, bool append>
    void apply(V &y) const {
        for(size_t i = 0; i < traits::number_of_components<V>::value; i++)
            A.apply(x(i), y(i), negate ? -scale : scale, append);
    }
};

namespace traits {

template <class M, class V>
struct is_scalable< multiadditive_operator<M, V> > : std::true_type {};

} // namespace traits

namespace traits {

// Number of components in a multivector expression.
struct multiex_dimension :
    boost::proto::or_ <
        boost::proto::when <
            boost::proto::terminal< boost::proto::_ >,
            traits::number_of_components<boost::proto::_>()
        > ,
        boost::proto::when <
            boost::proto::nary_expr<boost::proto::_, boost::proto::vararg<boost::proto::_> >,
            boost::proto::fold<boost::proto::_,
                boost::mpl::size_t<0>(),
                boost::mpl::max<multiex_dimension, boost::proto::_state>()>()
        >
    >
{};

template <class Expr, class Enable = void>
struct get_dimension {};

template <class Expr>
struct get_dimension<Expr, typename std::enable_if<
        boost::proto::matches<
            typename boost::proto::result_of::as_expr<Expr>::type,
            multivector_expr_grammar
        >::value &&
        !is_tuple<typename std::decay<Expr>::type>::value
    >::type>
{
    const static size_t value = std::result_of<traits::multiex_dimension(Expr)>::type::value;
};

template <class Expr>
struct get_dimension<Expr, typename std::enable_if<
        is_tuple<typename std::decay<Expr>::type>::value
    >::type>
{
    const static size_t value = std::tuple_size<typename std::decay<Expr>::type>::value;
};

} // namespace traits

//---------------------------------------------------------------------------
// Expression Transforms and evaluation contexts
//---------------------------------------------------------------------------
namespace detail {

// Helper functor for use with boost::fusion::for_each
template <class Context>
struct do_eval {
    Context &ctx;

    do_eval(Context &ctx) : ctx(ctx) {}

    template <class Expr>
    void operator()(const Expr &expr) const {
        boost::proto::eval(expr, ctx);
    }
};

// Generic terminal processing functor
struct process_terminal
    : boost::proto::transform < process_terminal >
{
    template<typename Expr, typename Unused1, typename Unused2>
    struct impl : boost::proto::transform_impl<Expr, Unused1, Unused2>
    {
        typedef typename impl::expr_param result_type;

        result_type operator ()(
              typename impl::expr_param term
            , typename impl::state_param process
            , typename impl::data_param) const
        {
            process(term);
            return term;
        }
    };
};

// Extract (and process) terminals from a vector expression.
struct extract_terminals
    : boost::proto::or_ <
        boost::proto::when <
            boost::proto::terminal<boost::proto::_>,
            process_terminal
        > ,
        boost::proto::function<
            boost::proto::_,
            boost::proto::vararg< extract_terminals >
        > ,
        boost::proto::when <
            boost::proto::nary_expr<
                boost::proto::_,
                boost::proto::vararg< extract_terminals >
            >
        >
    >
{};

// Base class for stateful expression evaluation contexts .
struct expression_context {
    backend::source_generator &src;
    const backend::command_queue &queue;
    mutable int prm_idx;
    int fun_idx;
    std::string prefix;
    kernel_generator_state_ptr state;

    expression_context(
            backend::source_generator &src, const backend::command_queue &queue,
            const std::string &prefix, kernel_generator_state_ptr state
            )
        : src(src), queue(queue), prm_idx(0), fun_idx(0),
          prefix(prefix), state(state)
    {}
};

// Outputs kernel preamble.
struct output_terminal_preamble : public expression_context {

    output_terminal_preamble(
            backend::source_generator &src, const backend::command_queue &queue,
            const std::string &prefix, kernel_generator_state_ptr state
            )
        : expression_context(src, queue, prefix, state)
    {}

    // Any expression except user function or terminal is only interesting
    // for its children:
    template <typename Expr, typename Tag = typename Expr::proto_tag>
    struct eval {
        typedef void result_type;

        void operator()(const Expr &expr, output_terminal_preamble &ctx) const
        {
            boost::fusion::for_each( expr,
                    do_eval<output_terminal_preamble>(ctx));
        }
    };

    // Function is either builtin (not interesting) or user-defined:
    template <typename Expr>
    struct eval<Expr, boost::proto::tag::function> {
        typedef void result_type;

        // Builtin function is only interesting for its children:
        template <class FunCall>
        typename std::enable_if<
            std::is_base_of<
                builtin_function,
                typename boost::proto::result_of::value<
                    typename boost::proto::result_of::child_c<FunCall,0>::type
                >::type
            >::value,
        void
        >::type
        operator()(const FunCall &expr, output_terminal_preamble &ctx) const
        {
            boost::fusion::for_each(
                    boost::fusion::pop_front(expr),
                    do_eval<output_terminal_preamble>(ctx)
                    );
        }

        // User-defined function needs to be defined.
        // Then look at its children:
        template <class FunCall>
        typename std::enable_if<
            std::is_base_of<
                user_function,
                typename boost::proto::result_of::value<
                    typename boost::proto::result_of::child_c<FunCall,0>::type
                >::type
            >::value,
        void
        >::type
        operator()(const FunCall &expr, output_terminal_preamble &ctx) const {
            // Output function definition (once) and continue with parameters.
            auto s = ctx.state->find("user_functions");
            if (s == ctx.state->end()) {
                s = ctx.state->insert(std::make_pair(
                            std::string("user_functions"),
                            boost::any( std::set<std::string>() )
                            )).first;
            }
            auto &seen = boost::any_cast< std::set<std::string>& >(s->second);

            typedef typename boost::proto::result_of::value<
                typename boost::proto::result_of::child_c<FunCall,0>::type
            >::type fun;

            std::string fname = fun::name();

            if (seen.find(fname) == seen.end()) {
                seen.insert(fname);
                fun::define(ctx.src);
            }

            boost::fusion::for_each(
                    boost::fusion::pop_front(expr),
                    do_eval<output_terminal_preamble>(ctx)
                    );
        }
    };

    // Some terminals have preambles too:
    template <typename T>
    struct eval<T, boost::proto::tag::terminal> {
        typedef void result_type;

        template <class Term>
        typename std::enable_if<traits::terminal_is_value<Term>::value, void>::type
        operator()(const Term &term, output_terminal_preamble &ctx) const
        {
            std::ostringstream prm_name;
            prm_name << ctx.prefix << "_" << ++ctx.prm_idx;

            traits::get_terminal_preamble(ctx.src, term, ctx.queue,
                    prm_name.str(), ctx.state);
        }

        template <class Term>
        typename std::enable_if<!traits::terminal_is_value<Term>::value, void>::type
        operator()(const Term &term, output_terminal_preamble &ctx) const
        {
            std::ostringstream prm_name;
            prm_name << ctx.prefix << "_" << ++ctx.prm_idx;

            traits::get_terminal_preamble(ctx.src, boost::proto::value(term),
                    ctx.queue, prm_name.str(), ctx.state);
        }
    };
};

// Performs local initialization (such as declaring and initializing temporary values).
struct output_local_preamble : public expression_context {

    output_local_preamble(
            backend::source_generator &src, const backend::command_queue &queue,
            const std::string &prefix,
            kernel_generator_state_ptr state
            )
        : expression_context(src, queue, prefix, state)
    {}

    // Any expression except user function or terminal is only interesting
    // for its children:
    template <typename Expr, typename Tag = typename Expr::proto_tag>
    struct eval {
        typedef void result_type;

        void operator()(const Expr &expr, output_local_preamble &ctx) const
        {
            boost::fusion::for_each( expr,
                    do_eval<output_local_preamble>(ctx));
        }
    };

    // Functions are only interesting for their parameters:
    template <typename Expr>
    struct eval<Expr, boost::proto::tag::function> {
        typedef void result_type;

        // Builtin function is only interesting for its children:
        template <class FunCall>
        void operator()(const FunCall &expr, output_local_preamble &ctx) const
        {
            boost::fusion::for_each(
                    boost::fusion::pop_front(expr),
                    do_eval<output_local_preamble>(ctx)
                    );
        }
    };

    // Some terminals need to be initialized:
    template <typename T>
    struct eval<T, boost::proto::tag::terminal> {
        typedef void result_type;

        template <class Term>
        typename std::enable_if<traits::terminal_is_value<Term>::value, void>::type
        operator()(const Term &term, output_local_preamble &ctx) const
        {
            std::ostringstream prm_name;
            prm_name << ctx.prefix << "_" << ++ctx.prm_idx;

            traits::get_local_terminal_init(ctx.src, term, ctx.queue,
                    prm_name.str(), ctx.state);
        }

        template <class Term>
        typename std::enable_if<!traits::terminal_is_value<Term>::value, void>::type
        operator()(const Term &term, output_local_preamble &ctx) const
        {
            std::ostringstream prm_name;
            prm_name << ctx.prefix << "_" << ++ctx.prm_idx;

            traits::get_local_terminal_init(ctx.src, boost::proto::value(term),
                    ctx.queue, prm_name.str(), ctx.state);
        }
    };
};

// Builds textual representation for a vector expression.
struct vector_expr_context : public expression_context {

    vector_expr_context(
            backend::source_generator &src, const backend::command_queue &queue,
            const std::string &prefix,
            kernel_generator_state_ptr state
            )
        : expression_context(src, queue, prefix, state)
    {}

    template <typename Expr, typename Tag = typename Expr::proto_tag>
    struct eval {};

#define VEXCL_BINARY_OPERATION(the_tag, the_op)                                \
  template <typename Expr> struct eval<Expr, boost::proto::tag::the_tag> {     \
    typedef void result_type;                                                  \
    void operator()(const Expr &expr, vector_expr_context &ctx) const {        \
      ctx.src << "( ";                                                         \
      boost::proto::eval(boost::proto::left(expr), ctx);                       \
      ctx.src << " " #the_op " ";                                              \
      boost::proto::eval(boost::proto::right(expr), ctx);                      \
      ctx.src << " )";                                                         \
    }                                                                          \
  }

    VEXCL_BINARY_OPERATION(plus,          +);
    VEXCL_BINARY_OPERATION(minus,         -);
    VEXCL_BINARY_OPERATION(multiplies,    *);
    VEXCL_BINARY_OPERATION(divides,       /);
    VEXCL_BINARY_OPERATION(modulus,       %);
    VEXCL_BINARY_OPERATION(shift_left,   <<);
    VEXCL_BINARY_OPERATION(shift_right,  >>);
    VEXCL_BINARY_OPERATION(less,          <);
    VEXCL_BINARY_OPERATION(greater,       >);
    VEXCL_BINARY_OPERATION(less_equal,   <=);
    VEXCL_BINARY_OPERATION(greater_equal,>=);
    VEXCL_BINARY_OPERATION(equal_to,     ==);
    VEXCL_BINARY_OPERATION(not_equal_to, !=);
    VEXCL_BINARY_OPERATION(logical_and,  &&);
    VEXCL_BINARY_OPERATION(logical_or,   ||);
    VEXCL_BINARY_OPERATION(bitwise_and,   &);
    VEXCL_BINARY_OPERATION(bitwise_or,    |);
    VEXCL_BINARY_OPERATION(bitwise_xor,   ^);

#undef VEXCL_BINARY_OPERATION

#define VEXCL_UNARY_PRE_OPERATION(the_tag, the_op)                             \
  template <typename Expr> struct eval<Expr, boost::proto::tag::the_tag> {     \
    typedef void result_type;                                                  \
    void operator()(const Expr &expr, vector_expr_context &ctx) const {        \
      ctx.src << "( " #the_op "( ";                                            \
      boost::proto::eval(boost::proto::child(expr), ctx);                      \
      ctx.src << " ) )";                                                       \
    }                                                                          \
  }

    VEXCL_UNARY_PRE_OPERATION(unary_plus,   +);
    VEXCL_UNARY_PRE_OPERATION(negate,       -);
    VEXCL_UNARY_PRE_OPERATION(logical_not,  !);
    VEXCL_UNARY_PRE_OPERATION(pre_inc,     ++);
    VEXCL_UNARY_PRE_OPERATION(pre_dec,     --);
    VEXCL_UNARY_PRE_OPERATION(address_of,   &);
    VEXCL_UNARY_PRE_OPERATION(dereference,  *);

#undef VEXCL_UNARY_PRE_OPERATION

#define VEXCL_UNARY_POST_OPERATION(the_tag, the_op)                            \
  template <typename Expr> struct eval<Expr, boost::proto::tag::the_tag> {     \
    typedef void result_type;                                                  \
    void operator()(const Expr &expr, vector_expr_context &ctx) const {        \
      ctx.src << "( ( ";                                                       \
      boost::proto::eval(boost::proto::child(expr), ctx);                      \
      ctx.src << " )" #the_op " )";                                            \
    }                                                                          \
  }

    VEXCL_UNARY_POST_OPERATION(post_inc, ++);
    VEXCL_UNARY_POST_OPERATION(post_dec, --);

#undef VEXCL_UNARY_POST_OPERATION

    template <typename Expr>
    struct eval<Expr, boost::proto::tag::if_else_> {
        typedef void result_type;
        void operator()(const Expr &expr, vector_expr_context &ctx) const {
            ctx.src << "( ";
            boost::proto::eval(boost::proto::child_c<0>(expr), ctx);
            ctx.src << " ? ";
            boost::proto::eval(boost::proto::child_c<1>(expr), ctx);
            ctx.src << " : ";
            boost::proto::eval(boost::proto::child_c<2>(expr), ctx);
            ctx.src << " )";
        }
    };

    template <typename Expr>
    struct eval<Expr, boost::proto::tag::subscript> {
        typedef void result_type;
        void operator()(const Expr &expr, vector_expr_context &ctx) const {
            ctx.src << "( ( ";
            boost::proto::eval(boost::proto::child_c<0>(expr), ctx);
            ctx.src << " )[ ";
            boost::proto::eval(boost::proto::child_c<1>(expr), ctx);
            ctx.src << " ] )";
        }
    };

    template <typename Expr>
    struct eval<Expr, boost::proto::tag::function> {
        typedef void result_type;

        struct do_eval {
            mutable int pos;
            vector_expr_context &ctx;

            do_eval(vector_expr_context &ctx) : pos(0), ctx(ctx) {}

            template <typename Arg>
            void operator()(const Arg &arg) const {
                if (pos++) ctx.src << ", ";
                boost::proto::eval(arg, ctx);
            }
        };

        template <class FunCall>
        void operator()(const FunCall &expr, vector_expr_context &ctx) const {
            ctx.src << boost::proto::value(boost::proto::child_c<0>(expr)).name() << "( ";
            boost::fusion::for_each(
                    boost::fusion::pop_front(expr), do_eval(ctx)
                    );
            ctx.src << " )";
        }
    };

    template <typename Expr>
    struct eval<Expr, boost::proto::tag::terminal> {
        typedef void result_type;

        template <typename Term>
        typename std::enable_if<traits::terminal_is_value<Term>::value, void>::type
        operator()(const Term &term, vector_expr_context &ctx) const {
            std::ostringstream prm_name;
            prm_name << ctx.prefix << "_" << ++ctx.prm_idx;

            traits::get_partial_vector_expr(ctx.src, term, ctx.queue,
                    prm_name.str(), ctx.state);
        }

        template <typename Term>
        typename std::enable_if<!traits::terminal_is_value<Term>::value, void>::type
        operator()(const Term &term, vector_expr_context &ctx) const {
            std::ostringstream prm_name;
            prm_name << ctx.prefix << "_" << ++ctx.prm_idx;

            traits::get_partial_vector_expr(ctx.src, boost::proto::value(term),
                    ctx.queue, prm_name.str(), ctx.state);
        }
    };
};

struct declare_expression_parameter : expression_context {

    declare_expression_parameter(backend::source_generator &src,
            const backend::command_queue &queue, const std::string &prefix,
            kernel_generator_state_ptr state
            )
        : expression_context(src, queue, prefix, state)
    {}

    template <typename Term>
    typename std::enable_if<traits::terminal_is_value<Term>::value, void>::type
    operator()(const Term &term) const {
        std::ostringstream prm_name;
        prm_name << prefix << "_" << ++prm_idx;

        traits::get_kernel_param_declaration(src, term, queue,
                prm_name.str(), state);
    }

    template <typename Term>
    typename std::enable_if<!traits::terminal_is_value<Term>::value, void>::type
    operator()(const Term &term) const {
        std::ostringstream prm_name;
        prm_name << prefix << "_" << ++prm_idx;

        traits::get_kernel_param_declaration(src, boost::proto::value(term),
                queue, prm_name.str(), state);
    }
};

struct set_expression_argument {
    backend::kernel &krn;
    unsigned part;
    size_t part_start;
    kernel_generator_state_ptr state;

    set_expression_argument(backend::kernel &krn, unsigned part, size_t part_start,
            kernel_generator_state_ptr state
            )
        : krn(krn), part(part), part_start(part_start), state(state)
    {}

    template <typename Term>
    typename std::enable_if<traits::terminal_is_value<Term>::value, void>::type
    operator()(const Term &term) const {
        traits::set_kernel_args(term, krn, part, part_start, state);
    }

    template <typename Term>
    typename std::enable_if<!traits::terminal_is_value<Term>::value, void>::type
    operator()(const Term &term) const {
        traits::set_kernel_args(boost::proto::value(term), krn, part, part_start, state);
    }
};

struct get_expression_properties {
    mutable std::vector<backend::command_queue> queue;
    mutable std::vector<size_t> part;
    mutable size_t size;

    get_expression_properties() : size(0) {}

    size_t part_start(unsigned d) const {
        return part.empty() ? 0 : part[d];
    }

    size_t part_size(unsigned d) const {
        return part.empty() ? 0 : part[d + 1] - part[d];
    }

    template <typename Term>
    typename std::enable_if<traits::terminal_is_value<Term>::value, void>::type
    operator()(const Term &term) const {
        get(term);
    }

    template <typename Term>
    typename std::enable_if<!traits::terminal_is_value<Term>::value, void>::type
    operator()(const Term &term) const {
        get(boost::proto::value(term));
    }

    template <typename Term>
    void get(const Term &term) const {
        if (queue.empty())
            traits::extract_expression_properties(term, queue, part, size);
#if (VEXCL_CHECK_SIZES > 0)
        else {
            std::vector<backend::command_queue> q;
            std::vector<size_t> p;
            size_t s = 0;

            traits::extract_expression_properties(term, q, p, s);

            precondition(
                    q.empty() || queue.empty() || q.size() == queue.size(),
                    "Incompatible queue lists");

            precondition(
                    s == 0 || size == 0 || s == size,
                    "Incompatible expression sizes");
        }
#endif
    }
};

//---------------------------------------------------------------------------
VEXCL_VECTOR_EXPR_EXTRACTOR(extract_vector_expressions,
        vector_expr_grammar,
        additive_vector_transform_grammar,
        vector_full_grammar
        );

VEXCL_ADDITIVE_EXPR_EXTRACTOR(extract_additive_vector_transforms,
        vector_expr_grammar,
        additive_vector_transform_grammar,
        vector_full_grammar
        );

struct simplify_additive_transform
    : boost::proto::or_<
          boost::proto::terminal< boost::proto::_ >,
          boost::proto::when<
             boost::proto::negate< boost::proto::terminal< boost::proto::_ > >,
             boost::proto::_
          >,
          boost::proto::when<
             boost::proto::negate< boost::proto::negate< boost::proto::_ > >,
             simplify_additive_transform(boost::proto::_child(boost::proto::_child))
          >,
          boost::proto::plus< simplify_additive_transform, simplify_additive_transform >,
          boost::proto::when<
            boost::proto::minus< boost::proto::_, boost::proto::_ >,
            boost::proto::_make_plus(
                    simplify_additive_transform(boost::proto::_left),
                    simplify_additive_transform(boost::proto::_make_negate(boost::proto::_right))
                    )
          >,
          boost::proto::when<
             boost::proto::negate< boost::proto::plus<boost::proto::_, boost::proto::_> >,
             boost::proto::_make_plus(
                     simplify_additive_transform(boost::proto::_make_negate(boost::proto::_left(boost::proto::_child))),
                     simplify_additive_transform(boost::proto::_make_negate(boost::proto::_right(boost::proto::_child)))
                     )
          >,
          boost::proto::when<
             boost::proto::negate< boost::proto::minus<boost::proto::_, boost::proto::_> >,
             boost::proto::_make_plus(
                     simplify_additive_transform(boost::proto::_make_negate(boost::proto::_left(boost::proto::_child))),
                     simplify_additive_transform(boost::proto::_right(boost::proto::_child))
                     )
          >
      >
{};

template <bool append, class Vector>
struct additive_applicator {
    Vector &dest;

    additive_applicator(Vector &dest) : dest(dest) {}

    template <typename Expr>
    typename std::enable_if<
        boost::proto::matches<
            typename boost::proto::result_of::as_expr<Expr>::type,
            boost::proto::terminal<boost::proto::_>
        >::value,
        void
    >::type
    operator()(const Expr &expr) const {
        expr.template apply</*negate=*/false, append>(dest);
    }

    template <typename Expr>
    typename std::enable_if<
        boost::proto::matches<
            typename boost::proto::result_of::as_expr<Expr>::type,
            boost::proto::negate<boost::proto::_>
        >::value,
        void
    >::type
    operator()(const Expr &expr) const {
        boost::proto::child(expr).template apply</*negate=*/true, append>(dest);
    }
};

template <bool append, class Vector, class Expr>
typename std::enable_if<
    boost::proto::matches<
        typename boost::proto::result_of::as_expr<Expr>::type,
        boost::proto::terminal<boost::proto::_>
    >::value ||
    boost::proto::matches<
        typename boost::proto::result_of::as_expr<Expr>::type,
        boost::proto::negate<boost::proto::_>
    >::value,
    void
>::type apply_additive_transform(Vector &dest, const Expr &expr) {
    (additive_applicator<append, Vector>(dest))(expr);
}

template <bool append, class Vector, class Expr>
typename std::enable_if<
    !boost::proto::matches<
        typename boost::proto::result_of::as_expr<Expr>::type,
        boost::proto::terminal<boost::proto::_>
    >::value &&
    !boost::proto::matches<
        typename boost::proto::result_of::as_expr<Expr>::type,
        boost::proto::negate<boost::proto::_>
    >::value,
    void
>::type apply_additive_transform(Vector &dest, const Expr &expr) {
    auto flat_expr = boost::proto::flatten(expr);

    (additive_applicator<append, Vector>(dest))(boost::fusion::front(flat_expr));

    boost::fusion::for_each(boost::fusion::pop_front(flat_expr),
            additive_applicator</*append=*/true, Vector>(dest)
            );
}

//---------------------------------------------------------------------------
// Multiexpression component extractor
//---------------------------------------------------------------------------
template <size_t C>
struct extract_component : boost::proto::callable {
    template <class T>
    struct result;

    template <class This, class T>
    struct result< This(T) > {
        typedef const typename traits::component< C,
                typename std::decay<T>::type
            >::type type;
    };

    template <class T>
    typename result<extract_component(const T&)>::type
    operator()(const T &t) const {
        using namespace std;
        return get<C>(t);
    }
};

template <size_t C>
struct extract_subexpression
    : boost::proto::or_ <
        boost::proto::when <
            boost::proto::and_<
                boost::proto::terminal< boost::proto::_ >,
                boost::proto::if_< traits::proto_terminal_is_value< boost::proto::_value >() >
            >,
            extract_component<C>( boost::proto::_ )
        > ,
        boost::proto::when <
            boost::proto::terminal<boost::proto::_>,
            boost::proto::_make_terminal (
                    extract_component<C>(
                        boost::proto::_value(boost::proto::_)
                        )
                    )
        > ,
        boost::proto::function<
            boost::proto::_,
            boost::proto::vararg< extract_subexpression<C> >
        > ,
        boost::proto::when <
            boost::proto::nary_expr<
                boost::proto::_,
                boost::proto::vararg< extract_subexpression<C> >
            >
        >
    >
{};

template <size_t C>
struct subexpression {
    template<class Expr>
    static auto get(const Expr &expr) ->
        typename std::enable_if<
            boost::proto::matches<
                typename boost::proto::result_of::as_expr<Expr>::type,
                multivector_expr_grammar
            >::value,
            decltype(extract_subexpression<C>()(boost::proto::as_child(expr)))
        >::type
    {
        return extract_subexpression<C>()(boost::proto::as_child(expr));
    }

    // If expression does not match multivector_expr_grammar, assume its a
    // tuple of vector expressions.
    template<class Expr>
    static auto get(const Expr &expr) ->
        typename std::enable_if<
            !boost::proto::matches<
                typename boost::proto::result_of::as_expr<Expr>::type,
                multivector_expr_grammar
            >::value,
            decltype(std::get<C>(expr))
        >::type
    {
        return std::get<C>(expr);
    }
};

VEXCL_VECTOR_EXPR_EXTRACTOR(extract_multivector_expressions,
        multivector_expr_grammar,
        additive_multivector_transform_grammar,
        multivector_full_grammar
        );

VEXCL_ADDITIVE_EXPR_EXTRACTOR(extract_additive_multivector_transforms,
        multivector_expr_grammar,
        additive_multivector_transform_grammar,
        multivector_full_grammar
        );

//---------------------------------------------------------------------------
// Expression result type deduction
//---------------------------------------------------------------------------

// Proxy for value_type<>
struct get_value_type : boost::proto::callable {
    template <class T> struct result;

    template <class This, class T>
    struct result< This(T) > {
        typedef
            typename traits::value_type< typename std::decay<T>::type >::type
            type;
    };
};

// Proxy for std::common_type<>
struct common_type : boost::proto::callable {
    template <class T, class Enable = void> struct result;

    template <class This, class T1, class T2>
    struct result< This(T1, T2) >
    {
        typedef typename std::decay<T1>::type D1;
        typedef typename std::decay<T2>::type D2;

        static_assert(
                cl_vector_length<D1>::value == 1 ||
                cl_vector_length<D2>::value == 1 ||
                cl_vector_length<D1>::value == cl_vector_length<D2>::value,
                "Operations with vectors of different lengths are not supported"
                );

        typedef
            typename cl_vector_of<
                typename std::common_type<
                    typename cl_scalar_of<D1>::type,
                    typename cl_scalar_of<D2>::type
                >::type,
                boost::mpl::max<
                    boost::mpl::size_t<cl_vector_length<D1>::value>,
                    boost::mpl::size_t<cl_vector_length<D2>::value>
                >::type::value
            >::type type;
    };
};


struct deduce_value_type
    : boost::proto::or_<
        // Terminals are passed to value_type<>
        boost::proto::when <
            boost::proto::and_<
                boost::proto::terminal< boost::proto::_ >,
                boost::proto::if_< traits::proto_terminal_is_value< boost::proto::_value >() >
            >,
            get_value_type( boost::proto::_ )
        > ,
        boost::proto::when <
            boost::proto::terminal< boost::proto::_ >,
            get_value_type( boost::proto::_value )
        >,
        // Result of logical operations is bool for scalars and long for vector
        // types. Lets keep it simple and return long.
        boost::proto::when <
            boost::proto::or_<
                boost::proto::or_<
                    boost::proto::less          < boost::proto::_, boost::proto::_ >,
                    boost::proto::greater       < boost::proto::_, boost::proto::_ >,
                    boost::proto::less_equal    < boost::proto::_, boost::proto::_ >,
                    boost::proto::greater_equal < boost::proto::_, boost::proto::_ >,
                    boost::proto::equal_to      < boost::proto::_, boost::proto::_ >,
                    boost::proto::not_equal_to  < boost::proto::_, boost::proto::_ >
                >,
                boost::proto::or_<
                    boost::proto::logical_and   < boost::proto::_, boost::proto::_ >,
                    boost::proto::logical_or    < boost::proto::_, boost::proto::_ >,
                    boost::proto::logical_not   < boost::proto::_ >
                >
            >,
            cl_long()
        >,
        boost::proto::when <
            boost::proto::if_else_< boost::proto::_, boost::proto::_, boost::proto::_ >,
            common_type( deduce_value_type(boost::proto::_child1), deduce_value_type(boost::proto::_child2) )
        >,
        // We assume that type of builtin function is the common type of its
        // arguments (TODO: this could be wrong for some functions).
        boost::proto::when <
            boost::proto::function<
                boost::proto::terminal<
                    boost::proto::convertible_to< builtin_function >
                >,
                boost::proto::vararg< boost::proto::_ >
                >,
            boost::proto::fold<
                boost::proto::functional::pop_front( boost::proto::_ ),
                char(),
                common_type(deduce_value_type, boost::proto::_state)
            >()
        >,
        // User-defined functions know their return type
        boost::proto::when <
            boost::proto::function<
                boost::proto::terminal<
                    boost::proto::convertible_to< user_function >
                >,
                boost::proto::vararg< boost::proto::_ >
                >,
            get_value_type( boost::proto::functional::value(boost::proto::_child0) )
        >,
        // Fold the operands of nary epxressions with std::common_type<>
        boost::proto::when <
            boost::proto::nary_expr<boost::proto::_, boost::proto::vararg<boost::proto::_> >,
            boost::proto::fold<
                boost::proto::_,
                char(),
                common_type(deduce_value_type, boost::proto::_state)
            >()
        >
      >
{};

// Hide the ugly type deduction details in an easy to use metafunction:
template <class Expr>
struct return_type {
    typedef
        typename std::decay<
                typename boost::result_of<
                    deduce_value_type(
                        typename boost::proto::result_of::as_expr<
                            typename std::decay<Expr>::type
                        >::type
                        )
                >::type
            >::type
        type;
};


//---------------------------------------------------------------------------
// Assign expression to lhs
//---------------------------------------------------------------------------
template <class OP, class LHS, class RHS>
void assign_expression(LHS &lhs, const RHS &rhs,
        const std::vector<backend::command_queue> &queue,
        const std::vector<size_t> &part
        )
{
#if (VEXCL_CHECK_SIZES > 0)
    {
        get_expression_properties prop;
        extract_terminals()(boost::proto::as_child(lhs), prop);
        extract_terminals()(boost::proto::as_child(rhs), prop);

        precondition(
                prop.queue.empty() || prop.queue.size() == queue.size(),
                "Incompatible queue lists"
                );

        precondition(
                prop.size == 0 || prop.size == part.back(),
                "Incompatible expression sizes"
                );
    }
#endif
    static kernel_cache cache;

    for(unsigned d = 0; d < queue.size(); d++) {
        auto kernel = cache.find(queue[d]);

        backend::select_context(queue[d]);

        if (kernel == cache.end()) {
            backend::source_generator source(queue[d]);

            output_terminal_preamble termpream(source, queue[d], "prm", empty_state());

            boost::proto::eval(boost::proto::as_child(lhs), termpream);
            boost::proto::eval(boost::proto::as_child(rhs), termpream);

            source.begin_kernel("vexcl_vector_kernel");
            source.begin_kernel_parameters();
            source.parameter<size_t>("n");

            declare_expression_parameter declare(source, queue[d], "prm", empty_state());

            extract_terminals()(boost::proto::as_child(lhs), declare);
            extract_terminals()(boost::proto::as_child(rhs), declare);

            source.end_kernel_parameters();
            source.grid_stride_loop().open("{");

            output_local_preamble loc_init(source, queue[d], "prm", empty_state());
            boost::proto::eval(boost::proto::as_child(lhs), loc_init);
            boost::proto::eval(boost::proto::as_child(rhs), loc_init);

            vector_expr_context expr_ctx(source, queue[d], "prm", empty_state());

            source.new_line();
            boost::proto::eval(boost::proto::as_child(lhs), expr_ctx);
            source << " " << OP::string() << " ";
            boost::proto::eval(boost::proto::as_child(rhs), expr_ctx);

            source << ";";
            source.close("}").end_kernel();

            kernel = cache.insert(queue[d], backend::kernel(
                        queue[d], source.str(), "vexcl_vector_kernel"));
        }

        if (size_t psize = part[d + 1] - part[d]) {
            kernel->second.push_arg(psize);

            set_expression_argument setarg(kernel->second, d, part[d], empty_state());

            extract_terminals()( boost::proto::as_child(lhs), setarg);
            extract_terminals()( boost::proto::as_child(rhs), setarg);

            kernel->second(queue[d]);
        }
    }
}

template <class OP, class LHS, class RHS>
void assign_expression(LHS &lhs, const RHS &rhs) {
    get_expression_properties prop;
    extract_terminals()(boost::proto::as_child(lhs), prop);

    precondition(!prop.queue.empty() && !prop.part.empty(),
            "Can not determine expression size and queue list"
            );

    assign_expression<OP>(lhs, rhs, prop.queue, prop.part);
}

// Static for loop
template <size_t Begin, size_t End>
class static_for {
    public:
        template <class Func>
        static void loop(Func &&f) {
            iterate<Begin>(f);
        }

    private:
        template <size_t I, class Func>
        static typename std::enable_if<(I < End)>::type
        iterate(Func &&f) {
            f.template apply<I>();
            iterate<I + 1>(f);
        }

        template <size_t I, class Func>
        static typename std::enable_if<(I >= End)>::type
        iterate(Func&&)
        { }
};

template <class OP, class LHS, class RHS>
struct subexpression_assigner {
    const LHS &lhs;
    const RHS &rhs;
    const std::vector<backend::command_queue> &queue;
    const std::vector<size_t> &part;

    subexpression_assigner(LHS &lhs, const RHS &rhs,
            const std::vector<backend::command_queue> &queue,
            const std::vector<size_t> &part
            )
        : lhs(lhs), rhs(rhs), queue(queue), part(part) {}

    template <size_t I>
    void apply() const {
        detail::assign_expression<OP>(
                subexpression<I>::get(lhs),
                subexpression<I>::get(rhs),
                queue, part);
    }
};

template <class LHS, class RHS>
struct preamble_constructor {
    const LHS &lhs;
    const RHS &rhs;

    kernel_generator_state_ptr state;
    mutable detail::output_terminal_preamble lhs_ctx;
    mutable detail::output_terminal_preamble rhs_ctx;

    preamble_constructor(const LHS &lhs, const RHS &rhs,
            backend::source_generator &source, const backend::command_queue &queue
            )
        : lhs(lhs), rhs(rhs), state(empty_state()),
          lhs_ctx(source, queue, "lhs", state),
          rhs_ctx(source, queue, "rhs", state)
    { }

    template <size_t I>
    void apply() const {
        boost::proto::eval(subexpression<I>::get(lhs), lhs_ctx);
        boost::proto::eval(subexpression<I>::get(rhs), rhs_ctx);
    }
};

template <class LHS, class RHS>
struct parameter_declarator {
    const LHS &lhs;
    const RHS &rhs;

    kernel_generator_state_ptr state;
    mutable detail::declare_expression_parameter lhs_ctx;
    mutable detail::declare_expression_parameter rhs_ctx;

    parameter_declarator(const LHS &lhs, const RHS &rhs,
            backend::source_generator &source, const backend::command_queue &queue)
        : lhs(lhs), rhs(rhs), state(empty_state()),
          lhs_ctx(source, queue, "lhs", state),
          rhs_ctx(source, queue, "rhs", state)
    { }

    template <size_t I>
    void apply() const {
        extract_terminals()(subexpression<I>::get(lhs), lhs_ctx);
        extract_terminals()(subexpression<I>::get(rhs), rhs_ctx);
    }
};

template <class LHS, class RHS>
struct expression_init {
    const LHS &lhs;
    const RHS &rhs;

    backend::source_generator &source;

    kernel_generator_state_ptr state;
    mutable detail::output_local_preamble rhs_pre;
    mutable detail::vector_expr_context   rhs_ctx;

    expression_init(const LHS &lhs, const RHS &rhs,
            backend::source_generator &source, const backend::command_queue &queue)
        : lhs(lhs), rhs(rhs), source(source), state(empty_state()),
          rhs_pre(source, queue, "rhs", state),
          rhs_ctx(source, queue, "rhs", state)
    { }

    template <size_t I>
    void apply() const {
        boost::proto::eval(subexpression<I>::get(rhs), rhs_pre);

        typedef
            typename return_type<decltype(subexpression<I>::get(lhs))>::type
            RT;

        source.new_line() << type_name<RT>() << " buf_" << I + 1 << " = ";

        boost::proto::eval(subexpression<I>::get(rhs), rhs_ctx);
        source << ";";
    }
};

template <class OP, class LHS>
struct expression_finalize {
    const LHS &lhs;

    backend::source_generator &source;

    kernel_generator_state_ptr state;
    mutable detail::output_local_preamble lhs_pre;
    mutable detail::vector_expr_context   lhs_ctx;

    expression_finalize(const LHS &lhs,
            backend::source_generator &source, const backend::command_queue &queue)
        : lhs(lhs), source(source), state(empty_state()),
          lhs_pre(source, queue, "lhs", state),
          lhs_ctx(source, queue, "lhs", state)
    { }

    template <size_t I>
    void apply() const {
        boost::proto::eval(subexpression<I>::get(lhs), lhs_pre);
        source.new_line();
        boost::proto::eval(subexpression<I>::get(lhs), lhs_ctx);
        source << " " << OP::string() << " buf_" << I + 1 << ";";
    }
};

template <class LHS, class RHS>
struct kernel_arg_setter {
    const LHS &lhs;
    const RHS &rhs;

    mutable detail::set_expression_argument ctx;

    kernel_arg_setter(const LHS &lhs, const RHS &rhs,
            backend::kernel &krn, unsigned part, size_t offset)
        : lhs(lhs), rhs(rhs), ctx(krn, part, offset, empty_state())
    { }

    template <size_t I>
    void apply() const {
        detail::extract_terminals()(subexpression<I>::get(lhs), ctx);
        detail::extract_terminals()(subexpression<I>::get(rhs), ctx);
    }
};

template <class OP, class LHS, class RHS>
void assign_multiexpression( LHS &lhs, const RHS &rhs,
        const std::vector<backend::command_queue> &queue,
        const std::vector<size_t> &part
        )
{

#if (VEXCL_CHECK_SIZES > 0)
    {
        get_expression_properties prop;
        extract_terminals()(subexpression<0>::get(lhs), prop);
        extract_terminals()(subexpression<0>::get(rhs), prop);

        precondition(
                prop.queue.empty() || prop.queue.size() == queue.size(),
                "Incompatible queue lists"
                );

        precondition(
                prop.size == 0 || prop.size == part.back(),
                "Incompatible expression sizes"
                );
    }
#endif

    typedef traits::get_dimension<LHS> N;

    static kernel_cache cache;

    // 1. If any device in context is CPU, then do not fuse the kernel,
    //    but assign components individually (this works better with CPU
    //    caches).
    // 2. If dimension of the multiexpression is 1, then assign_expression()
    //    would work better as well (no need to spend registers on temp
    //    variables).
    if (
#ifdef VEXCL_SPLIT_MULTIEXPRESSIONS
            1 ||
#endif
            (N::value == 1) ||
            std::any_of(queue.begin(), queue.end(),
                [](const backend::command_queue &q) { return backend::is_cpu(q); })
       )
    {
        static_for<0, N::value>::loop(
                subexpression_assigner<OP, LHS, RHS>(lhs, rhs, queue, part)
                );
        return;
    }

    for(unsigned d = 0; d < queue.size(); d++) {
        auto kernel = cache.find(queue[d]);

        backend::select_context(queue[d]);

        if (kernel == cache.end()) {
            backend::source_generator source(queue[d]);

            static_for<0, N::value>::loop(
                    preamble_constructor<LHS, RHS>(lhs, rhs, source, queue[d])
                    );

            source.begin_kernel("vexcl_multivector_kernel");
            source.begin_kernel_parameters();
            source.parameter<size_t>("n");

            static_for<0, N::value>::loop(
                    parameter_declarator<LHS, RHS>(lhs, rhs, source, queue[d])
                    );

            source.end_kernel_parameters();
            source.grid_stride_loop().open("{");

            static_for<0, N::value>::loop(expression_init<LHS, RHS>(lhs, rhs, source, queue[d]));
            static_for<0, N::value>::loop(expression_finalize<OP, LHS>(lhs, source, queue[d]));

            source.close("}").end_kernel();

            kernel = cache.insert(queue[d], backend::kernel(
                        queue[d], source.str(), "vexcl_multivector_kernel") );
        }

        if (size_t psize = part[d + 1] - part[d]) {
            kernel->second.push_arg(psize);

            static_for<0, N::value>::loop(
                    kernel_arg_setter<LHS, RHS>(lhs, rhs, kernel->second, d, part[d])
                    );

            kernel->second(queue[d]);
        }
    }
}

template <class OP, class LHS, class RHS>
void assign_multiexpression( LHS &lhs, const RHS &rhs) {
    get_expression_properties prop;
    extract_terminals()(subexpression<0>::get(lhs), prop);

    precondition(!prop.queue.empty() && !prop.part.empty(),
            "Can not determine expression size and queue list"
            );

    assign_multiexpression<OP>(lhs, rhs, prop.queue, prop.part);
}

} // namespace detail

/// Assignable tuple of expressions
template <class LHS>
struct expression_tuple {
    static const size_t NDIM = std::tuple_size<LHS>::value;

    const LHS lhs;

    expression_tuple(const LHS &lhs) : lhs(lhs) {}

#define VEXCL_ASSIGNMENT(cop, op)                                              \
  /** \brief Multiexpression assignment.                                       \
   \details All operations are delegated to components of the multivector.     \
   */                                                                          \
  template <class RHS>                                                         \
  auto operator cop(const RHS & rhs) const ->                                  \
      typename std::enable_if <                                                \
          boost::proto::matches<                                               \
              typename boost::proto::result_of::as_expr<RHS>::type,            \
              multivector_expr_grammar>::value || is_tuple<RHS>::value,        \
      const expression_tuple & > ::type                                        \
  {                                                                            \
    detail::assign_multiexpression<op>(lhs, rhs);                              \
    return *this;                                                              \
  }

    VEXCL_ASSIGNMENT(=,   assign::SET)
    VEXCL_ASSIGNMENT(+=,  assign::ADD)
    VEXCL_ASSIGNMENT(-=,  assign::SUB)
    VEXCL_ASSIGNMENT(*=,  assign::MUL)
    VEXCL_ASSIGNMENT(/=,  assign::DIV)
    VEXCL_ASSIGNMENT(%=,  assign::MOD)
    VEXCL_ASSIGNMENT(&=,  assign::AND)
    VEXCL_ASSIGNMENT(|=,  assign::OR)
    VEXCL_ASSIGNMENT(^=,  assign::XOR)
    VEXCL_ASSIGNMENT(<<=, assign::LSH)
    VEXCL_ASSIGNMENT(>>=, assign::RSH)

#undef VEXCL_ASSIGNMENT
};

#ifndef BOOST_NO_VARIADIC_TEMPLATES
/// Ties several vector expressions into a writeable tuple.
/**
 The following example results in a single kernel:

 \code
 vex::vector<double> x(ctx, 1024);
 vex::vector<double> y(ctx, 1024);

 vex::tie(x,y) = std::make_tuple( x + y, y - x );
 \endcode

 This is functionally equivalent to the following code, 
 but does not use temporary vectors and is more efficient:

 \code
 tmp_x = x + y;
 tmp_y = y - x;
 x = tmp_x;
 y = tmp_y;
 \endcode
 */
template<class... Expr>
auto tie(const Expr&... expr) ->
    expression_tuple< std::tuple<const Expr&...> >
{
    return expression_tuple< std::tuple<const Expr&...> >( std::tie(expr...) );
}
#else

#define VEXCL_TIE_VECTORS(z, n, data)                                          \
  template <BOOST_PP_ENUM_PARAMS(n, class Expr)>                               \
  expression_tuple<std::tuple<                                                 \
      BOOST_PP_ENUM_BINARY_PARAMS(n, const Expr, &BOOST_PP_INTERCEPT)> >       \
  tie(BOOST_PP_ENUM_BINARY_PARAMS(n, const Expr, &expr)) {                     \
    return expression_tuple<std::tuple<                                        \
        BOOST_PP_ENUM_BINARY_PARAMS(n, const Expr, &BOOST_PP_INTERCEPT)> >(    \
        std::tie(BOOST_PP_ENUM_PARAMS(n, expr)));                              \
  }

BOOST_PP_REPEAT_FROM_TO(1, VEXCL_MAX_ARITY, VEXCL_TIE_VECTORS, ~)

#undef VEXCL_TIE_VECTORS

#endif

/// Meta-filter for VexCL vector expressions
template <class T>
struct is_vector_expression
    : std::is_same< typename boost::proto::domain_of<T>::type, vector_domain >
{};

/// Meta-filter for VexCL multivector expressions
template <class T>
struct is_multivector_expression
    : boost::mpl::or_<
        boost::proto::matches<
          typename boost::proto::result_of::as_expr<T>::type,
          multivector_full_grammar
          >,
        is_tuple<T>
        >
{};

/// Helper function for getting command queue and size info from a vector expression
template <class Expr>
typename std::enable_if<
    !is_multivector_expression<Expr>::value,
    std::tuple<
        std::vector<backend::command_queue>,
        size_t
        >
    >::type
expression_properties(const Expr &expr) {
    detail::get_expression_properties prop;
    detail::extract_terminals()(boost::proto::as_child(expr), prop);
    return std::make_tuple(prop.queue, prop.size);
}

/// Helper function for getting command queue and size info from a vector multiexpression
template <class Expr>
typename std::enable_if<
    is_multivector_expression<Expr>::value,
    std::tuple<
        std::vector<backend::command_queue>,
        size_t
        >
    >::type
expression_properties(const Expr &expr) {
    detail::get_expression_properties prop;
    detail::extract_terminals()(detail::subexpression<0>::get(expr), prop);
    return std::make_tuple(prop.queue, prop.size);
}

} // namespace vex;


#endif
