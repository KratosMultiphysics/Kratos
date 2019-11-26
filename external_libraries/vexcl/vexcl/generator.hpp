#ifndef VEXCL_GENERATOR_HPP
#define VEXCL_GENERATOR_HPP

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
 * \file   generator.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  OpenCL kernel generator.
 */

#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <stdexcept>
#include <memory>

#include <boost/proto/proto.hpp>
#include <boost/function_types/parameter_types.hpp>
#include <boost/function_types/result_type.hpp>
#include <boost/function_types/function_arity.hpp>

#include <boost/fusion/include/for_each.hpp>
#include <boost/fusion/include/vector_tie.hpp>

#include <vexcl/util.hpp>
#include <vexcl/operations.hpp>
#include <vexcl/function.hpp>
#include <vexcl/vector.hpp>

#include <boost/preprocessor/repetition.hpp>
#ifndef VEXCL_MAX_ARITY
#  define VEXCL_MAX_ARITY BOOST_PROTO_MAX_ARITY
#endif

/// Vector expression template library for OpenCL.
namespace vex {

template <typename T> class symbolic;

/// Sends name of the symbolic variable to output stream.
template <typename T>
std::ostream& operator<<(std::ostream &os, const symbolic<T> &sym);

/// Kernel generation interface.
namespace generator {

//---------------------------------------------------------------------------
// The recorder class. Holds static output stream for kernel recording and
// static variable index (used in variable names).
//---------------------------------------------------------------------------
template <bool dummy = true>
class recorder {
    static_assert(dummy, "dummy parameter should be true");

    public:
        static void set(std::ostream &s) {
            os = &s;

            // Reset preamble and state.
            preamble.reset(new backend::source_generator);
            state = vex::detail::empty_state();
        }

        static std::ostream& get() {
            return os ? *os : std::cout;
        }

        static backend::source_generator& get_preamble() {
            return *preamble;
        }

        static vex::detail::kernel_generator_state_ptr get_state() {
            return state;
        }

        static size_t var_id() {
            return ++index;
        }
    private:
        static size_t index;
        static std::ostream *os;
        static std::unique_ptr<backend::source_generator> preamble;
        static vex::detail::kernel_generator_state_ptr state;
};

template <bool dummy>
size_t recorder<dummy>::index = 0;

template <bool dummy>
std::ostream *recorder<dummy>::os = 0;

template <bool dummy>
std::unique_ptr<backend::source_generator> recorder<dummy>::preamble;

template <bool dummy>
vex::detail::kernel_generator_state_ptr recorder<dummy>::state;

inline size_t var_id() {
    return recorder<>::var_id();
}

inline std::ostream& get_recorder() {
    return recorder<>::get();
}

inline backend::source_generator& get_preamble() {
    return recorder<>::get_preamble();
}

inline vex::detail::kernel_generator_state_ptr get_state() {
    return recorder<>::get_state();
}

/// Set output stream for the kernel recorder.
inline void set_recorder(std::ostream &os) {
    recorder<>::set(os);
}

//---------------------------------------------------------------------------
// Setting up boost::proto.
//---------------------------------------------------------------------------
struct variable {};

// --- The grammar ----------------------------------------------------------
struct symbolic_grammar
    : boost::proto::or_<
          boost::proto::or_<
              boost::proto::terminal< variable >,
              boost::proto::and_<
                  boost::proto::terminal< boost::proto::_ >,
                  boost::proto::if_< is_cl_native< boost::proto::_value >() >
              >
          >,
          VEXCL_BUILTIN_OPERATIONS(symbolic_grammar),
          VEXCL_USER_FUNCTIONS(symbolic_grammar)
      >
{};

template <class Expr>
struct symbolic_expr;

struct symbolic_domain
    : boost::proto::domain< boost::proto::generator< symbolic_expr >, symbolic_grammar >
{
    // Store everything by value inside expressions...
    template <typename T, class Enable = void>
    struct as_child : proto_base_domain::as_expr<T> {};

    // ... except for symbolic variables:
    template <typename T>
    struct as_child< T,
        typename std::enable_if<
            boost::proto::matches<
                typename boost::proto::result_of::as_expr< T >::type,
                boost::proto::terminal< variable >
                >::value
            >::type
        > : proto_base_domain::as_child< T > {};
};

template <class Expr>
struct symbolic_expr
    : boost::proto::extends< Expr, symbolic_expr< Expr >, symbolic_domain >
{
    typedef boost::proto::extends< Expr, symbolic_expr< Expr >, symbolic_domain > base_type;

    symbolic_expr(const Expr &expr = Expr()) : base_type(expr) {}
};

//---------------------------------------------------------------------------
struct index_expr
    : public generator::symbolic_expr< boost::proto::terminal< generator::variable >::type >
{};

inline auto index()
    -> boost::proto::result_of::as_expr<index_expr, symbolic_domain>::type const
{
    return boost::proto::as_expr<symbolic_domain>(index_expr());
}

//---------------------------------------------------------------------------
namespace detail {

struct symbolic_context {
    template <typename Expr, typename Tag = typename Expr::proto_tag>
    struct eval {};

#define VEXCL_BINARY_OPERATION(bin_tag, bin_op)                                \
  template <typename Expr> struct eval<Expr, boost::proto::tag::bin_tag> {     \
    typedef void result_type;                                                  \
    void operator()(const Expr &expr, symbolic_context &ctx) const {           \
      get_recorder() << "( ";                                                  \
      boost::proto::eval(boost::proto::left(expr), ctx);                       \
      get_recorder() << " " #bin_op " ";                                       \
      boost::proto::eval(boost::proto::right(expr), ctx);                      \
      get_recorder() << " )";                                                  \
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
    void operator()(const Expr &expr, symbolic_context &ctx) const {           \
      get_recorder() << "( " #the_op "( ";                                     \
      boost::proto::eval(boost::proto::child(expr), ctx);                      \
      get_recorder() << " ) )";                                                \
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
    void operator()(const Expr &expr, symbolic_context &ctx) const {           \
      get_recorder() << "( ( ";                                                \
      boost::proto::eval(boost::proto::child(expr), ctx);                      \
      get_recorder() << " )" #the_op " )";                                     \
    }                                                                          \
  }

    VEXCL_UNARY_POST_OPERATION(post_inc, ++);
    VEXCL_UNARY_POST_OPERATION(post_dec, --);

#undef VEXCL_UNARY_POST_OPERATION

    template <typename Expr>
    struct eval<Expr, boost::proto::tag::if_else_> {
        typedef void result_type;
        void operator()(const Expr &expr, symbolic_context &ctx) const {
            get_recorder() << "( ";
            boost::proto::eval(boost::proto::child_c<0>(expr), ctx);
            get_recorder() << " ? ";
            boost::proto::eval(boost::proto::child_c<1>(expr), ctx);
            get_recorder() << " : ";
            boost::proto::eval(boost::proto::child_c<2>(expr), ctx);
            get_recorder() << " )";
        }
    };

    template <class Expr>
    struct eval<Expr, boost::proto::tag::function> {
        typedef void result_type;

        struct display {
            mutable int pos;
            symbolic_context &ctx;

            display(symbolic_context &ctx) : pos(0), ctx(ctx) {}

            template <class Arg>
            void operator()(const Arg &arg) const {
                if (pos++) get_recorder() << ", ";
                boost::proto::eval(arg, ctx);
            }
        };

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
        operator()(const FunCall &expr, symbolic_context &ctx) const {
            get_recorder() << boost::proto::value(boost::proto::child_c<0>(expr)).name() << "( ";

            boost::fusion::for_each(
                    boost::fusion::pop_front(expr),
                    display(ctx)
                    );

            get_recorder() << " )";
        }

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
        operator()(const FunCall &expr, symbolic_context &ctx) const {
            typedef typename boost::proto::result_of::value<
                typename boost::proto::result_of::child_c<FunCall,0>::type
            >::type fun;

            // Output function definition (once).
            auto s = get_state()->find("user_functions");
            if (s == get_state()->end()) {
                s = get_state()->insert(std::make_pair(
                            std::string("user_functions"),
                            boost::any( std::set<std::string>() )
                            )).first;
            }
            auto &seen = boost::any_cast< std::set<std::string>& >(s->second);

            std::string fname = fun::name();

            if (seen.find(fname) == seen.end()) {
                seen.insert(fname);
                fun::define(get_preamble());
            }

            get_recorder() << fun::name() << "( ";

            boost::fusion::for_each(
                    boost::fusion::pop_front(expr),
                    display(ctx)
                    );

            get_recorder() << " )";
        }
    };

    template <typename Expr>
    struct eval<Expr, boost::proto::tag::terminal> {
        typedef void result_type;

        template <typename Term>
        void operator()(const Term &term, symbolic_context &) const {
            get_recorder() << std::scientific << std::setprecision(12)
                << boost::proto::value(term);
        }

        template <typename T>
        void operator()(const symbolic<T> &v, symbolic_context &) const {
            get_recorder() << v;
        }

        void operator()(const index_expr&, symbolic_context&) const {
            get_recorder() << "idx";
        }
    };
};

} // namespace detail

} // namespace generator

//---------------------------------------------------------------------------
// The symbolic class.
//---------------------------------------------------------------------------
/// Symbolic variable
template <typename T>
class symbolic
    : public generator::symbolic_expr< boost::proto::terminal< generator::variable >::type >
{
    public:
        typedef T value_type;

        /// Scope/Type of the symbolic variable.
        enum scope_type {
            LocalVar        = 0, ///< Local variable.
            VectorParameter = 1, ///< Vector kernel parameter.
            ScalarParameter = 2  ///< Scalar kernel parameter.
        };

        /// Constness of vector parameter.
        enum constness_type {
            NonConst = 0,   ///< Parameter should be written back at kernel exit.
            Const = 1       ///< Parameter is readonly.
        };

        /// Default constructor. Results in a local variable declaration.
        symbolic() : num(generator::var_id()), scope(LocalVar), constness(NonConst)
        {
            generator::get_recorder() << "\t\t" << type_name<T>() << " " << *this << " = " << T() << ";\n";
        }

        /// Constructor.
        explicit symbolic(scope_type scope, constness_type constness = NonConst)
            : num(generator::var_id()), scope(scope), constness(constness)
        {
            if (scope == LocalVar) {
                generator::get_recorder() << "\t\t" << type_name<T>() << " " << *this << ";\n";
            }
        }

        /// Copy constructor.
        symbolic(const symbolic &expr)
            : num(generator::var_id()), scope(LocalVar), constness(NonConst)
        {
            generator::get_recorder() << "\t\t" << type_name<T>() << " " << *this << " = ";
            record(expr);
            generator::get_recorder() << ";\n";
        }

        /// Expression constructor. Results in a local variable declaration initialized by the expression.
        template <class Expr>
        symbolic(const Expr &expr)
            : num(generator::var_id()), scope(LocalVar), constness(NonConst)
        {
            generator::get_recorder() << "\t\t" << type_name<T>() << " " << *this << " = ";
            record(expr);
            generator::get_recorder() << ";\n";
        }

        /// Assignment operator. Results in the assignment expression written to the recorder.
        const symbolic& operator=(const symbolic &c) const {
            generator::get_recorder() << "\t\t" << *this << " = ";
            record(c);
            generator::get_recorder() << ";\n";
            return *this;
        }

#define VEXCL_ASSIGNMENT(cop, op)                                              \
  /** Assignment operator.
   Results in the assignment expression written to the recorder. */            \
  template <class Expr> const symbolic &operator cop(const Expr & expr) {      \
    generator::get_recorder() << "\t\t" << *this << " " #cop " ";              \
    record(expr);                                                              \
    generator::get_recorder() << ";\n";                                        \
    return *this;                                                              \
  }

        VEXCL_ASSIGNMENTS(VEXCL_ASSIGNMENT)

#undef VEXCL_ASSIGNMENT

        size_t id() const {
            return num;
        }

        // Initialize local variable at kernel enter.
        std::string init() const {
            std::ostringstream s;

            if (scope == VectorParameter) {
                s << "\t\t" << type_name<T>() << " " << *this
                    << " = p_" << *this << "[idx];\n";
            } else if (scope == ScalarParameter) {
                s << "\t\t" << type_name<T>() << " " << *this
                    << " = p_" << *this << ";\n";
            }

            return s.str();
        }

        // Write local variable to parameter at kernel exit.
        std::string write() const {
            std::ostringstream s;

            if (scope == VectorParameter && constness == NonConst)
                s << "\t\tp_" << *this << "[idx] = " << *this << ";\n";

            return s.str();
        }

        // Returns parameter type and name as strings.
        std::tuple<std::string, std::string> prmdecl() const {
            std::ostringstream name;
            name << "p_" << *this;

            std::string prm_type;

            if (scope == VectorParameter) {
                if (constness == Const)
                    prm_type = type_name< global_ptr<const T> >();
                else
                    prm_type = type_name< global_ptr<T> >();
            } else {
                prm_type = type_name<T>();
            }

            return std::make_tuple(prm_type, name.str());
        }
    private:
        size_t         num;
        scope_type     scope;
        constness_type constness;

        template <class Expr>
        static void record(const Expr &expr) {
            generator::detail::symbolic_context ctx;
            boost::proto::eval(boost::proto::as_child(expr), ctx);
        }
};

template <typename T>
std::ostream& operator<<(std::ostream &os, const symbolic<T> &sym) {
    return os << "var" << sym.id();
}

namespace generator {

/// Autogenerated kernel.
class kernel {
    public:
        kernel(
                const std::vector<backend::command_queue> &queue,
                const std::string &name
              ) : queue(queue), name(name), psize(queue.size(), 0)
        {
            prm_read.reset(new std::ostringstream);
            prm_save.reset(new std::ostringstream);
        }

        template <class SymVar>
        void add_param(const SymVar &var) {
            prm_decl.push_back(var.prmdecl());
            *prm_read << var.init();
            *prm_save << var.write();
        }

        void build(const std::string &body) {
            for(auto q = queue.begin(); q != queue.end(); q++) {
                backend::source_generator source(*q);

                source << get_preamble().str();

                source.begin_kernel(name);
                source.begin_kernel_parameters();

                for(auto p = prm_decl.begin(); p != prm_decl.end(); ++p)
                    source.parameter(std::get<0>(*p), std::get<1>(*p));

                source.parameter<size_t>("n");

                source.end_kernel_parameters();
                source.grid_stride_loop().open("{");

                source.new_line() << prm_read->str() << body << prm_save->str();

                source.close("}");
                source.end_kernel();

                backend::select_context(*q);
                cache.insert(std::make_pair(
                            backend::get_context_id(*q),
                            backend::kernel(*q, source.str(), name.c_str())
                            ));
            }
        }

        template <class T>
        void push_arg(const T &v) {
            for(unsigned d = 0; d < queue.size(); d++) {
                cache.find(backend::get_context_id(queue[d]))->second.push_arg(v);
            }
        }

        template <class T>
        void push_arg(const vector<T> &v) {
            for(unsigned d = 0; d < queue.size(); d++) {
                cache.find(backend::get_context_id(queue[d]))->second.push_arg(v(d));
                psize[d] = std::max(psize[d], v.part_size(d));
            }
        }

        template <class T>
        void push_arg(const std::vector<T> &args) {
            for(unsigned d = 0; d < queue.size(); d++) {
                cache.find(backend::get_context_id(queue[d]))->second.push_arg(args[d]);
            }
        }

        void operator()() {
            for(unsigned d = 0; d < queue.size(); d++) {
                auto &K = cache.find(backend::get_context_id(queue[d]))->second;

                if (psize[d]) {
                    K.push_arg(psize[d]);
                    K(queue[d]);

                    psize[d] = 0;
                } else {
                    K.reset();
                }
            }
        }

#ifndef BOOST_NO_VARIADIC_TEMPLATES
        /// Launches the kernel with the provided parameters.
        template <class Head, class... Tail>
        void operator()(const Head &head, const Tail&... tail) {
            push_arg(head);
            (*this)(tail...);
        }
#else

#define VEXCL_FUNCALL_OPERATOR(z, n, data)                                     \
  template <BOOST_PP_ENUM_PARAMS(n, class Param)>                              \
  void operator()(BOOST_PP_ENUM_BINARY_PARAMS(n, const Param, &param)) {       \
    boost::fusion::for_each(                                                   \
            boost::fusion::vector_tie(BOOST_PP_ENUM_PARAMS(n, param)),         \
            push_args(*this)                                                   \
            );                                                                 \
    (*this)();                                                                 \
  }

BOOST_PP_REPEAT_FROM_TO(1, VEXCL_MAX_ARITY, VEXCL_FUNCALL_OPERATOR, ~)

#undef VEXCL_FUNCALL_OPERATOR

#endif

        static void add_params(kernel &) {}

        template <class Head, class... Tail>
        static void add_params(kernel &K, const Head &head, const Tail&... tail) {
            K.add_param(head);
            add_params(K, tail...);
        }
    private:
        std::vector<backend::command_queue> queue;
        std::string name;
        std::vector<size_t> psize;
        std::vector< std::tuple<std::string, std::string> > prm_decl;
        std::unique_ptr<std::ostringstream> prm_read, prm_save;
        std::map<vex::backend::context_id, vex::backend::kernel> cache;

        struct push_args {
            kernel &K;
            push_args(kernel &K) : K(K) {}

            template <class T>
            void operator()(const T &p) const {
                K.push_arg(p);
            }
        };
};

/// Function body generator.
class Function {
    public:
        template <class Ret, class ArgTuple>
        Function(const std::string &body, const Ret &ret, const ArgTuple &arg)
        {
            boost::fusion::for_each(arg, read_params(source));

            source << body;

            source << "\t\treturn " << ret << ";\n";
        }

        std::string get() const {
            return source.str();
        }
    private:
        std::ostringstream source;

        struct read_params {
            std::ostream &os;
            mutable int prm_idx;

            read_params(std::ostream &os) : os(os), prm_idx(0) {}

            template <class T>
            void operator()(const T &v) const {
                os << "\t\t" << type_name<typename T::value_type>() << " "
                   << v << " = prm" << ++prm_idx << ";\n";
            }
        };
};

#ifndef BOOST_NO_VARIADIC_TEMPLATES
/// Builds kernel from the recorded expression sequence and the symbolic parameter list.
/** The symbolic variables passed to the function should have participated in
 * the recorded algorithm and will be converted to the generated kernel
 * arguments.
 */
template <class... Args>
kernel build_kernel(
        const std::vector<backend::command_queue> &queue,
        const std::string &name, const std::string& body, const Args&... args
        )
{
    kernel K(queue, name);
    kernel::add_params(K, args...);
    K.build(body);
    return K;
}

/// Builds function body from the recorded expression.
/** The symbolic variables passed to the function should have participated in
 * the recorded algorithm and will be converted to the output value and the
 * input arguments of the generated function.
 */
template <class Ret, class... Args>
std::string make_function(std::string body, const Ret &ret, const Args&... args) {
    return Function(body, ret, boost::fusion::vector_tie(args...)).get();
}
#else

#define VEXCL_BUILD_KERNEL(z, n, data)                                         \
  template <BOOST_PP_ENUM_PARAMS(n, class Arg)>                                \
  kernel build_kernel(const std::vector<backend::command_queue> & queue,    \
                         const std::string & name, const std::string & body,   \
                         BOOST_PP_ENUM_BINARY_PARAMS(n, const Arg, &arg)) {    \
    kernel K(queue, name);                                                     \
    boost::fusion::for_each(                                                   \
            boost::fusion::vector_tie(BOOST_PP_ENUM_PARAMS(n, arg)),           \
            detail::kernel_add_param(K));                                      \
    K.build(body);                                                             \
    return K;                                                                  \
  }

#define VEXCL_MAKE_FUNCTION(z, n, data)                                        \
  template <class Ret, BOOST_PP_ENUM_PARAMS(n, class Arg)>                     \
  std::string make_function(std::string body, const Ret &ret,                  \
                            BOOST_PP_ENUM_BINARY_PARAMS(n, const Arg, &arg)) { \
    return Function(body, ret, boost::fusion::vector_tie(                      \
                                   BOOST_PP_ENUM_PARAMS(n, arg))).get();       \
  }

BOOST_PP_REPEAT_FROM_TO(1, VEXCL_MAX_ARITY, VEXCL_BUILD_KERNEL, ~)
BOOST_PP_REPEAT_FROM_TO(1, VEXCL_MAX_ARITY, VEXCL_MAKE_FUNCTION, ~)

#undef VEXCL_BUILD_KERNEL
#undef VEXCL_MAKE_FUNCTION

#endif

// UserFunction implementation from a generic functor
template <class Signature, class Functor>
struct FunctorAdapter : UserFunction<FunctorAdapter<Signature, Functor>, Signature>
{
    static std::string name_string;
    static std::string body_string;

    FunctorAdapter(Functor &&f, std::string fname) {
        using boost::function_types::function_arity;

        name_string = fname;
        body_string = get_body(std::forward<Functor>(f),
                boost::mpl::size_t< function_arity<Signature>::value >() );
    }

    // Empty constructor. Used in UserFunction::operator(). Hopefuly the body
    // string is already constructed by the time the constructor is called.
    FunctorAdapter() {}

    static std::string name() { return name_string; }
    static std::string body() { return body_string; }

#define VEXCL_PRINT_PRM(z, n, data)                                            \
  typedef symbolic<                                                            \
      typename boost::mpl::at<params, boost::mpl::int_<n> >::type> Prm##n;     \
  Prm##n prm##n(Prm##n::ScalarParameter);                                      \
  source << "\t\t" << type_name<typename Prm##n::value_type>() << " "          \
         << prm##n << " = prm" << n + 1 << ";\n";

#define VEXCL_BODY_GETTER(z, n, data)                                          \
  static std::string get_body(Functor && f, boost::mpl::size_t<n>) {           \
    typedef typename boost::function_types::result_type<Signature>::type       \
        result;                                                                \
    typedef typename boost::function_types::parameter_types<Signature>::type   \
        params;                                                                \
    std::ostringstream source;                                                 \
    set_recorder(source);                                                      \
    BOOST_PP_REPEAT(n, VEXCL_PRINT_PRM, ~) symbolic<result> ret =              \
        f(BOOST_PP_ENUM_PARAMS(n, prm));                                       \
    source << "\t\treturn " << ret << ";\n";                                   \
    return source.str();                                                       \
  }

    BOOST_PP_REPEAT_FROM_TO(1, VEXCL_MAX_ARITY, VEXCL_BODY_GETTER, ~)

#undef VEXCL_BODY_GETTER
#undef VEXCL_PRINT_PRM
};

template <class Signature, class Functor>
std::string FunctorAdapter<Signature, Functor>::name_string;

template <class Signature, class Functor>
std::string FunctorAdapter<Signature, Functor>::body_string;

inline size_t get_gen_fun_id() {
    static size_t id = 0;
    return id++;
}

/// Generates a user-defined function from a generic functor.
/**
 * Takes the function signature as template parameter and a generic functor as
 * a single argument.
 * Returns user-defined function ready to be used in vector expressions.
 */
template <class Signature, class Functor>
auto make_function(Functor &&f) ->
    FunctorAdapter<Signature, Functor>
{
    size_t id = get_gen_fun_id();
    std::ostringstream name;
    name << "generated_function_" << id;
    return FunctorAdapter<Signature, Functor>(std::forward<Functor>(f), name.str());
}

} // namespace generator;

} // namespace vex;

#endif
