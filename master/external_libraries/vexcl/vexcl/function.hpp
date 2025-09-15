#ifndef VEXCL_FUNCTION_HPP
#define VEXCL_FUNCTION_HPP

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
 * \file   vexcl/function.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  User-defined device functions.
 */

#include <boost/preprocessor/enum.hpp>
#include <boost/preprocessor/enum_params.hpp>
#include <boost/preprocessor/seq.hpp>
#include <boost/preprocessor/cat.hpp>
#include <boost/preprocessor/tuple.hpp>
#include <boost/preprocessor/stringize.hpp>

#include <vexcl/operations.hpp>

namespace vex {

/// Converts an unquoted text into a string literal.
#define VEX_STRINGIZE_SOURCE(...) #__VA_ARGS__

//---------------------------------------------------------------------------
// VEX_FUNCTION (v1) macros
//---------------------------------------------------------------------------
/// Macro to declare a user function type.
/**
 \code
 VEX_FUNCTION_V1_TYPE(pow3_t, double(double), "", "return pow(prm1, 3.0);");
 pow3_t pow3;
 output = pow3(input);
 \endcode

 \deprecated

 \note This version of the macro uses function call signature in order to
 define the function paramaters. Parameters are named automatically (prm1,
 prm2, ...), which reduces readability of the code. Use of VEX_FUNCTION is
 recommended instead.

 \note Should be used in case same function is used in several places (to
 save on OpenCL kernel recompilations). Otherwise VEX_FUNCTION should
 be used locally.
 */
#define VEX_FUNCTION_V1_TYPE(fname, signature, preamble_str, body_str)         \
  struct vex_function_##fname                                                  \
    : vex::UserFunction<vex_function_##fname, signature>                       \
  {                                                                            \
    vex_function_##fname() {}                                                  \
    static std::string name() { return #fname; }                               \
    static std::string preamble() { return preamble_str; }                     \
    static std::string body() { return body_str; }                             \
  }

/// Macro to declare a user function.
/**
 \code
 VEX_FUNCTION_V1(pow3, double(double), "return pow(prm1, 3.0);");
 output = pow3(input);
 \endcode

 \deprecated

 \note This version of the macro uses function call signature in order to
 define the function paramaters. Parameters are named automatically (prm1,
 prm2, ...), which reduces readability of the code. Use of VEX_FUNCTION is
 recommended instead.
 */
#define VEX_FUNCTION_V1(name, signature, body)                                 \
  VEX_FUNCTION_V1_TYPE(name, signature, "", body) const name


/// Macro to declare a user function with preamble.
/**
 * The preamble may be used to define helper functions or macros.
 \code
 VEX_FUNCTION_V1_WITH_PREAMBLE(one, double(double),
         "double sin2(double x) { return pow(sin(x), 2.0); }\n"
         "double cos2(double x) { return pow(cos(x), 2.0); }\n",
         "return sin2(prm1) + cos2(prm1);"
         );
 y = one(x);
 \endcode

 \deprecated

 \note This version of the macro uses function call signature in order to
 define the function paramaters. Parameters are named automatically (prm1,
 prm2, ...), which reduces readability of the code. Use of VEX_FUNCTION is
 recommended instead.
 */
#define VEX_FUNCTION_V1_WITH_PREAMBLE(name, signature, preamble, body)         \
  VEX_FUNCTION_V1_TYPE(name, signature, preamble, body) const name



//---------------------------------------------------------------------------
// VEX_FUNCTION (v2) macros
//---------------------------------------------------------------------------
#define VEXCL_FUNCTION_ARG_TYPE(s, data, arg) BOOST_PP_TUPLE_ELEM(2, 0, arg)

#define VEXCL_FUNCTION_ARG_TYPES(args) \
    BOOST_PP_SEQ_ENUM(BOOST_PP_SEQ_TRANSFORM(VEXCL_FUNCTION_ARG_TYPE, ~, args))

#define VEXCL_FUNCTION_NTH_ARG_TYPE(n, args)                                   \
    BOOST_PP_TUPLE_ELEM(2, 0, BOOST_PP_SEQ_ELEM(n, args))

#define VEXCL_FUNCTION_NTH_ARG_NAME(n, args)                                   \
    BOOST_PP_STRINGIZE(BOOST_PP_TUPLE_ELEM(2, 1, BOOST_PP_SEQ_ELEM(n, args)))

#define VEXCL_FUNCTION_DEF_ARG(z, n, args)                                     \
    src.parameter<VEXCL_FUNCTION_NTH_ARG_TYPE(n, args)>(                       \
            VEXCL_FUNCTION_NTH_ARG_NAME(n, args));

#define VEXCL_FUNCTION_DEFINE_DEP(z, data, dep)                                \
    {                                                                          \
        typedef decltype(dep) dep_type;                                        \
        dep_type::define(src);                                                 \
    }

#define VEX_FUNCTION_SINK(rtype, func_name, nargs, args, deps, body)           \
struct vex_function_##func_name                                                \
    : vex::UserFunction<                                                       \
        vex_function_##func_name,                                              \
        rtype( VEXCL_FUNCTION_ARG_TYPES(args) )                                \
      >                                                                        \
{                                                                              \
    vex_function_##func_name() {}                                              \
    static std::string name() { return #func_name; }                           \
    static void define(vex::backend::source_generator &src) {                  \
        define(src, name());                                                   \
    }                                                                          \
    static void define(vex::backend::source_generator &src,                    \
            const std::string &fname)                                          \
    {                                                                          \
        BOOST_PP_SEQ_FOR_EACH(VEXCL_FUNCTION_DEFINE_DEP, ~, deps)              \
        src.begin_function< rtype >(fname);                                    \
        src.begin_function_parameters();                                       \
        BOOST_PP_REPEAT(nargs, VEXCL_FUNCTION_DEF_ARG, args)                   \
        src.end_function_parameters();                                         \
        src.new_line() << body;                                                \
        src.end_function();                                                    \
    }                                                                          \
} const func_name

#define VEXCL_FUNCTION_MAKE_SEQ_0(...) ((__VA_ARGS__)) VEXCL_FUNCTION_MAKE_SEQ_1
#define VEXCL_FUNCTION_MAKE_SEQ_1(...) ((__VA_ARGS__)) VEXCL_FUNCTION_MAKE_SEQ_0
#define VEXCL_FUNCTION_MAKE_SEQ_0_END
#define VEXCL_FUNCTION_MAKE_SEQ_1_END

#define VEXCL_FUNCTION_MAKE_SEQ(args)                                          \
    BOOST_PP_CAT(VEXCL_FUNCTION_MAKE_SEQ_0 args,_END)

#define VEXCL_DUAL_FUNCTOR_ARG(s, data, arg)                                   \
    BOOST_PP_TUPLE_ELEM(2, 0, arg) BOOST_PP_TUPLE_ELEM(2, 1, arg)

#define VEXCL_DUAL_FUNCTOR_ARGS(args) \
    BOOST_PP_SEQ_ENUM(BOOST_PP_SEQ_TRANSFORM(VEXCL_DUAL_FUNCTOR_ARG, ~, args))

#define VEX_DUAL_FUNCTOR_SINK(rtype, nargs, args, ...)                         \
rtype operator()(VEXCL_DUAL_FUNCTOR_ARGS(args)) const {                        \
    __VA_ARGS__                                                                \
}

/// Creates a user-defined function with dependencies.
/** The body of the function is passed as a string literal or a static string
 * expression.
 */
#define VEX_FUNCTION_SD(return_type, name, arguments, dependencies, body)      \
    VEX_FUNCTION_SINK(return_type, name,                                       \
            BOOST_PP_SEQ_SIZE(VEXCL_FUNCTION_MAKE_SEQ(arguments)),             \
            VEXCL_FUNCTION_MAKE_SEQ(arguments), dependencies, body)

/// Creates a user-defined function with dependencies.
/** The body of the function is passed as a string literal or a static string
 * expression.
 */
#define VEX_FUNCTION_DS VEX_FUNCTION_SD

/// Creates a user-defined function.
/** The body of the function is passed as a string literal or a static string
 * expression.
 */
#define VEX_FUNCTION_S(return_type, name, arguments, body)                     \
    VEX_FUNCTION_SD(return_type, name, arguments, , body)

/// Creates a user-defined function with dependencies.
/** The body of the function is specified as unquoted C source at the end of
 * the macro. The source will be stringized with VEX_STRINGIZE_SOURCE macro.
 */
#define VEX_FUNCTION_D(return_type, name, arguments, dependencies, ...)        \
    VEX_FUNCTION_SD(return_type, name, arguments, dependencies, VEX_STRINGIZE_SOURCE(__VA_ARGS__) )


/// Creates a user-defined function.
/**
 The body of the function is specified as unquoted C source at the end of the
 macro. The source will be stringized with VEX_STRINGIZE_SOURCE macro.
 */
#define VEX_FUNCTION(return_type, name, arguments, ...)                        \
    VEX_FUNCTION_S(return_type, name, arguments, VEX_STRINGIZE_SOURCE(__VA_ARGS__))

/// Defines both device and host versions of a function call operator.
/**
 The intended use is the creation of comparison and reduction functors for
 use with scan/sort/reduce algorithms.

 Example:
 \code
 template <typename T>
 struct less {
     VEX_DUAL_FUNCTOR(bool, (T, a)(T, b),
         return a < b;
         )
 };
 \endcode
 */
#define VEX_DUAL_FUNCTOR(type, args, ...) \
    VEX_FUNCTION(type, device, args, __VA_ARGS__);                             \
    VEX_DUAL_FUNCTOR_SINK(type,                                                \
            BOOST_PP_SEQ_SIZE(VEXCL_FUNCTION_MAKE_SEQ(args)),                  \
            VEXCL_FUNCTION_MAKE_SEQ(args), __VA_ARGS__)

//---------------------------------------------------------------------------
// Builtin functions
//---------------------------------------------------------------------------
#define VEXCL_BUILTIN_PRINT_BOOST_REF(z, n, data) boost::ref(arg##n)

/// Define builtin function.
#define VEX_BUILTIN_FUNCTION(nargs, func)                                      \
    struct func##_func : vex::builtin_function {                               \
        static const char *name() { return #func; }                            \
    };                                                                         \
    template <BOOST_PP_ENUM_PARAMS(nargs, class Arg)>                          \
    typename boost::proto::result_of::make_expr<                               \
        boost::proto::tag::function, func##_func,                              \
        BOOST_PP_ENUM_BINARY_PARAMS(nargs, const Arg,                          \
                                    &BOOST_PP_INTERCEPT)>::type const          \
    func(BOOST_PP_ENUM_BINARY_PARAMS(nargs, const Arg, &arg)) {                \
        return boost::proto::make_expr<boost::proto::tag::function>(           \
            func##_func(), BOOST_PP_ENUM(                                      \
                nargs, VEXCL_BUILTIN_PRINT_BOOST_REF, ~));                     \
    }

#define VEX_BUILTIN_FUNCTION_ALIAS(nargs, alias, func)                         \
    struct func##_alias : vex::builtin_function {                              \
        static const char *name() { return #func; }                            \
    };                                                                         \
    template <BOOST_PP_ENUM_PARAMS(nargs, class Arg)>                          \
    typename boost::proto::result_of::make_expr<                               \
        boost::proto::tag::function, func##_alias,                             \
        BOOST_PP_ENUM_BINARY_PARAMS(nargs, const Arg,                          \
                                    &BOOST_PP_INTERCEPT)>::type const          \
    alias(BOOST_PP_ENUM_BINARY_PARAMS(nargs, const Arg, &arg)) {               \
        return boost::proto::make_expr<boost::proto::tag::function>(           \
            func##_alias(), BOOST_PP_ENUM(                                     \
                nargs, VEXCL_BUILTIN_PRINT_BOOST_REF, ~));                     \
    }

/// \defgroup builtins Builtin device functions
/** @{ */
VEX_BUILTIN_FUNCTION( 2, abs_diff )
VEX_BUILTIN_FUNCTION( 1, acos )
VEX_BUILTIN_FUNCTION( 1, acosh )
VEX_BUILTIN_FUNCTION( 1, acospi )
VEX_BUILTIN_FUNCTION( 2, add_sat )
VEX_BUILTIN_FUNCTION( 1, all )
VEX_BUILTIN_FUNCTION( 1, any )
VEX_BUILTIN_FUNCTION( 1, asin )
VEX_BUILTIN_FUNCTION( 1, asinh )
VEX_BUILTIN_FUNCTION( 1, asinpi )
VEX_BUILTIN_FUNCTION( 1, atan )
VEX_BUILTIN_FUNCTION( 2, atan2 )
VEX_BUILTIN_FUNCTION( 2, atan2pi )
VEX_BUILTIN_FUNCTION( 1, atanh )
VEX_BUILTIN_FUNCTION( 1, atanpi )
VEX_BUILTIN_FUNCTION( 3, bitselect )
VEX_BUILTIN_FUNCTION( 1, cbrt )
VEX_BUILTIN_FUNCTION( 1, ceil )
VEX_BUILTIN_FUNCTION( 3, clamp )
VEX_BUILTIN_FUNCTION( 1, clz )
VEX_BUILTIN_FUNCTION( 2, copysign )
VEX_BUILTIN_FUNCTION( 1, cos )
VEX_BUILTIN_FUNCTION( 1, cosh )
VEX_BUILTIN_FUNCTION( 1, cospi )
VEX_BUILTIN_FUNCTION( 2, cross )
VEX_BUILTIN_FUNCTION( 1, degrees )
VEX_BUILTIN_FUNCTION( 2, distance )
VEX_BUILTIN_FUNCTION( 2, dot )
VEX_BUILTIN_FUNCTION( 1, erf )
VEX_BUILTIN_FUNCTION( 1, erfc )
VEX_BUILTIN_FUNCTION( 1, exp )
VEX_BUILTIN_FUNCTION( 1, exp10 )
VEX_BUILTIN_FUNCTION( 1, exp2 )
VEX_BUILTIN_FUNCTION( 1, expm1 )
VEX_BUILTIN_FUNCTION( 1, fabs )
VEX_BUILTIN_FUNCTION( 2, fast_distance )
VEX_BUILTIN_FUNCTION( 1, fast_length )
VEX_BUILTIN_FUNCTION( 1, fast_normalize )
VEX_BUILTIN_FUNCTION( 2, fdim )
VEX_BUILTIN_FUNCTION( 1, floor )
VEX_BUILTIN_FUNCTION( 3, fma )
VEX_BUILTIN_FUNCTION( 2, fmax )
VEX_BUILTIN_FUNCTION( 2, fmin )
VEX_BUILTIN_FUNCTION( 2, fmod )
VEX_BUILTIN_FUNCTION( 2, fract )
VEX_BUILTIN_FUNCTION( 2, frexp )
VEX_BUILTIN_FUNCTION( 2, hadd )
VEX_BUILTIN_FUNCTION( 2, hypot )
VEX_BUILTIN_FUNCTION( 1, ilogb )
VEX_BUILTIN_FUNCTION( 2, isequal )
VEX_BUILTIN_FUNCTION( 1, isfinite )
VEX_BUILTIN_FUNCTION( 2, isgreater )
VEX_BUILTIN_FUNCTION( 2, isgreaterequal )
VEX_BUILTIN_FUNCTION( 1, isinf )
VEX_BUILTIN_FUNCTION( 2, isless )
VEX_BUILTIN_FUNCTION( 2, islessequal )
VEX_BUILTIN_FUNCTION( 2, islessgreater )
VEX_BUILTIN_FUNCTION( 1, isnan )
VEX_BUILTIN_FUNCTION( 1, isnormal )
VEX_BUILTIN_FUNCTION( 2, isnotequal )
VEX_BUILTIN_FUNCTION( 2, isordered )
VEX_BUILTIN_FUNCTION( 2, isunordered )
VEX_BUILTIN_FUNCTION( 2, ldexp )
VEX_BUILTIN_FUNCTION( 1, length )
VEX_BUILTIN_FUNCTION( 1, lgamma )
VEX_BUILTIN_FUNCTION( 2, lgamma_r )
VEX_BUILTIN_FUNCTION( 1, log )
VEX_BUILTIN_FUNCTION( 1, log10 )
VEX_BUILTIN_FUNCTION( 1, log1p )
VEX_BUILTIN_FUNCTION( 1, log2 )
VEX_BUILTIN_FUNCTION( 1, logb )
VEX_BUILTIN_FUNCTION( 3, mad )
VEX_BUILTIN_FUNCTION( 3, mad24 )
VEX_BUILTIN_FUNCTION( 3, mad_hi )
VEX_BUILTIN_FUNCTION( 3, mad_sat )
VEX_BUILTIN_FUNCTION( 2, max )
VEX_BUILTIN_FUNCTION( 2, maxmag )
VEX_BUILTIN_FUNCTION( 2, min )
VEX_BUILTIN_FUNCTION( 2, minmag )
VEX_BUILTIN_FUNCTION( 3, mix )
VEX_BUILTIN_FUNCTION( 2, modf )
VEX_BUILTIN_FUNCTION( 2, mul_hi )
VEX_BUILTIN_FUNCTION( 1, nan )
VEX_BUILTIN_FUNCTION( 2, nextafter )
VEX_BUILTIN_FUNCTION( 1, normalize )
#if defined(VEXCL_BACKEND_CUDA)
VEX_BUILTIN_FUNCTION( 1, __popc )
VEX_BUILTIN_FUNCTION( 1, __popcll )
#else
VEX_BUILTIN_FUNCTION( 1, popcount )
#endif
VEX_BUILTIN_FUNCTION( 2, pow )
VEX_BUILTIN_FUNCTION( 2, pown )
VEX_BUILTIN_FUNCTION( 2, powr )
VEX_BUILTIN_FUNCTION( 1, radians )
VEX_BUILTIN_FUNCTION( 2, remainder )
VEX_BUILTIN_FUNCTION( 3, remquo )
VEX_BUILTIN_FUNCTION( 2, rhadd )
VEX_BUILTIN_FUNCTION( 1, rint )
VEX_BUILTIN_FUNCTION( 2, rootn )
VEX_BUILTIN_FUNCTION( 2, rotate )
VEX_BUILTIN_FUNCTION( 1, round )
VEX_BUILTIN_FUNCTION( 1, rsqrt )
VEX_BUILTIN_FUNCTION( 3, select )
VEX_BUILTIN_FUNCTION( 2, shuffle )
VEX_BUILTIN_FUNCTION( 3, shuffle2 )
VEX_BUILTIN_FUNCTION( 1, sign )
VEX_BUILTIN_FUNCTION( 1, signbit )
VEX_BUILTIN_FUNCTION( 1, sin )
VEX_BUILTIN_FUNCTION( 2, sincos )
VEX_BUILTIN_FUNCTION( 1, sinh )
VEX_BUILTIN_FUNCTION( 1, sinpi )
VEX_BUILTIN_FUNCTION( 3, smoothstep )
VEX_BUILTIN_FUNCTION( 1, sqrt )
VEX_BUILTIN_FUNCTION( 2, step )
VEX_BUILTIN_FUNCTION( 2, sub_sat )
VEX_BUILTIN_FUNCTION( 1, tan )
VEX_BUILTIN_FUNCTION( 1, tanh )
VEX_BUILTIN_FUNCTION( 1, tanpi )
VEX_BUILTIN_FUNCTION( 1, tgamma )
VEX_BUILTIN_FUNCTION( 1, trunc )
VEX_BUILTIN_FUNCTION( 2, upsample )

// Atomic functions
#if defined(VEXCL_BACKEND_CUDA)

VEX_BUILTIN_FUNCTION( 2, atomicAdd  )
VEX_BUILTIN_FUNCTION( 2, atomicSub  )
VEX_BUILTIN_FUNCTION( 2, atomicExch )
VEX_BUILTIN_FUNCTION( 2, atomicMin  )
VEX_BUILTIN_FUNCTION( 2, atomicMax  )
VEX_BUILTIN_FUNCTION( 2, atomicInc  )
VEX_BUILTIN_FUNCTION( 2, atomicDec  )
VEX_BUILTIN_FUNCTION( 3, atomicCAS  )
VEX_BUILTIN_FUNCTION( 2, atomicAnd  )
VEX_BUILTIN_FUNCTION( 2, atomicOr   )
VEX_BUILTIN_FUNCTION( 2, atomicXor  )

// Also provide aliases for OpenCL-style functions
VEX_BUILTIN_FUNCTION_ALIAS(2, atomic_add,     atomicAdd  )
VEX_BUILTIN_FUNCTION_ALIAS(2, atomic_sub,     atomicSub  )
VEX_BUILTIN_FUNCTION_ALIAS(2, atomic_xchg,    atomicExch )
VEX_BUILTIN_FUNCTION_ALIAS(2, atomic_min,     atomicMin  )
VEX_BUILTIN_FUNCTION_ALIAS(2, atomic_max,     atomicMax  )
VEX_BUILTIN_FUNCTION_ALIAS(3, atomic_cmpxchg, atomicCAS  )
VEX_BUILTIN_FUNCTION_ALIAS(2, atomic_and,     atomicAnd  )
VEX_BUILTIN_FUNCTION_ALIAS(2, atomic_or,      atomicOr   )
VEX_BUILTIN_FUNCTION_ALIAS(2, atomic_xor,     atomicXor  )

#else

VEX_BUILTIN_FUNCTION(2, atomic_add     )
VEX_BUILTIN_FUNCTION(2, atomic_sub     )
VEX_BUILTIN_FUNCTION(2, atomic_xchg    )
VEX_BUILTIN_FUNCTION(2, atomic_min     )
VEX_BUILTIN_FUNCTION(2, atomic_max     )
VEX_BUILTIN_FUNCTION(1, atomic_inc     )
VEX_BUILTIN_FUNCTION(1, atomic_dec     )
VEX_BUILTIN_FUNCTION(3, atomic_cmpxchg )
VEX_BUILTIN_FUNCTION(2, atomic_and     )
VEX_BUILTIN_FUNCTION(2, atomic_or      )
VEX_BUILTIN_FUNCTION(2, atomic_xor     )

#endif


// Special case: abs() overloaded with floating point arguments should call
// fabs in the OpenCL code
struct abs_func : builtin_function {
    static const char* name() {
        return "abs";
    }
};

namespace detail {
    template <class Expr> struct return_type;
}

template <typename Arg>
auto abs(const Arg &arg) ->
    typename std::enable_if<
        std::is_integral<
            typename cl_scalar_of<
                typename detail::return_type<Arg>::type
            >::type
        >::value,
        typename boost::proto::result_of::make_expr<
            boost::proto::tag::function,
            abs_func,
            const Arg&
        >::type const
    >::type
{
    return boost::proto::make_expr<boost::proto::tag::function>(
            abs_func(),
            boost::ref(arg)
            );
}

template <typename Arg>
auto abs(const Arg &arg) ->
    typename std::enable_if<
        !std::is_integral<
            typename cl_scalar_of<
                typename detail::return_type<Arg>::type
            >::type
        >::value,
        typename boost::proto::result_of::make_expr<
            boost::proto::tag::function,
            fabs_func,
            const Arg&
        >::type const
    >::type
{
    return boost::proto::make_expr<boost::proto::tag::function>(
            fabs_func(),
            boost::ref(arg)
            );
}

/** @} */

} // namespace vex

#endif
