#ifndef VEXCL_REDUCTOR_HPP
#define VEXCL_REDUCTOR_HPP

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
 * \file   vexcl/reductor.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Vector expression reduction.
 */

#include <vector>
#include <array>
#include <string>
#include <sstream>
#include <numeric>
#include <limits>

#include <vexcl/vector.hpp>
#include <vexcl/operations.hpp>

namespace vex {

/// Summation.
struct SUM {
    // In order to define a reduction kind for vex::Reductor, one should define
    // a struct like the following:
    template <class T>
    struct impl {
        typedef T result_type;

        // Initial value for the operation.
        static T initial() {
            return T();
        }

        // Device-side reduction function.
        struct device_in : UserFunction<device_in, T(T, T)> {
            static std::string name() { return "SUM_" + type_name<T>(); }
            static std::string body() { return "return prm1 + prm2;"; }
        };

        typedef device_in device_out;

        // Host-side reduction function.
        T operator()(const T &a, const T &b) const {
            return a + b;
        }
    };
};

/// Compensated summation.
/**
 * Reduces the numerical error in the result with <a
 * href="http://en.wikipedia.org/wiki/Kahan_summation_algorithm">
 * Kahan summation algorithm</a>.
 */
struct SUM_Kahan : SUM {};

/// Maximum element.
struct MAX {
    template <class T>
    struct impl {
        typedef T result_type;

        // Initial value for the operation.
        static T initial() {
            return std::numeric_limits<T>::lowest();
        }

        struct device_in : UserFunction<device_in, T(T, T)> {
            static std::string name() { return "MAX_" + type_name<T>(); }
            static std::string body() { return "return prm1 > prm2 ? prm1 : prm2;"; }
        };

        typedef device_in device_out;

        T operator()(const T &a, const T &b) const {
            return a > b ? a : b;
        }
    };
};

/// Minimum element.
struct MIN {
    template <class T>
    struct impl {
        typedef T result_type;

        // Initial value for the operation.
        static T initial() {
            return std::numeric_limits<T>::max();
        }

        struct device_in : UserFunction<device_in, T(T, T)> {
            static std::string name() { return "MIN_" + type_name<T>(); }
            static std::string body() { return "return prm1 < prm2 ? prm1 : prm2;"; }
        };

        typedef device_in device_out;

        T operator()(const T &a, const T &b) const {
            return a < b ? a : b;
        }
    };
};

#ifndef BOOST_NO_VARIADIC_TEMPLATES
/// Combines several reduce operations.
template <class... R>
struct CombineReductors {
    template <class T>
    struct impl {
        static_assert(is_cl_scalar<T>::value, "Only scalar reductors may be combined");
        typedef typename cl_vector_of<T, cl_fit_vec_size<sizeof...(R)>::value>::type result_type;

        static result_type initial() {
            result_type r;
            assign_initial<R...>(r.s);
            return r;
        }

        struct device_in : UserFunction<device_in, result_type(result_type, T)> {
            static std::string name() { return "CombinedReductor_in_" + type_name<T>(); }

            static void define(backend::source_generator &src, const std::string &fname = name()) {
                define_scalar_fun<R...>(src);

                src.begin_function<result_type>(fname);
                src.begin_function_parameters();
                src.template parameter<result_type>("a");
                src.template parameter<T>("b");
                src.end_function_parameters();
                src.new_line() << type_name<result_type>() << " r = { ";

                int pos = 0;
                partial_reduce<R...>(src, pos);

                src << "};";
                src.new_line() << "return r;";
                src.end_function();
            }

            template <class Head>
            static void define_scalar_fun(backend::source_generator &src) {
                Head::template impl<T>::device_in::define(src);
            }

            template <class Head, class... Tail>
            static
            typename std::enable_if<(sizeof...(Tail) > 0), void>::type
            define_scalar_fun(backend::source_generator &src) {
                Head::template impl<T>::device_in::define(src);
                define_scalar_fun<Tail...>(src);
            }

            template <class Head>
            static void partial_reduce(backend::source_generator &src, int &pos) {
                src << Head::template impl<T>::device_in::name() << "(a." << vcmp(pos) << ", b) ";
                pos++;
            }

            template <class Head, class... Tail>
            static
            typename std::enable_if<(sizeof...(Tail) > 0), void>::type
            partial_reduce(backend::source_generator &src, int &pos) {
                src << Head::template impl<T>::device_in::name() << "(a." << vcmp(pos) << ", b), ";
                pos++;
                partial_reduce<Tail...>(src, pos);
            }
        };

        struct device_out : UserFunction<device_out, result_type(result_type, T)> {
            static std::string name() { return "CombinedReductor_out_" + type_name<T>(); }

            static void define(backend::source_generator &src, const std::string &fname = name()) {
                src.begin_function<result_type>(fname);
                src.begin_function_parameters();
                src.template parameter<result_type>("a");
                src.template parameter<result_type>("b");
                src.end_function_parameters();
                src.new_line() << type_name<result_type>() << " r = { ";

                int pos = 0;
                partial_reduce<R...>(src, pos);

                src << "};";
                src.new_line() << "return r;";
                src.end_function();
            }

            template <class Head>
            static void partial_reduce(backend::source_generator &src, int &pos) {
                src << Head::template impl<T>::device_out::name()
                    << "(a." << vcmp(pos) << ", b." << vcmp(pos) << ") ";
                pos++;
            }

            template <class Head, class... Tail>
            static
            typename std::enable_if<(sizeof...(Tail) > 0), void>::type
            partial_reduce(backend::source_generator &src, int &pos) {
                src << Head::template impl<T>::device_out::name()
                    << "(a." << vcmp(pos) << ", b." << vcmp(pos) << "), ";
                pos++;
                partial_reduce<Tail...>(src, pos);
            }
        };

        result_type operator()(result_type a, result_type b) const {
            result_type r;
            partial_reduce<R...>(r.s, a.s, b.s);
            return r;
        }

        private:
            template <class Head>
            static void assign_initial(T *p) {
                *p++ = Head::template impl<T>::initial();
            }

            template <class Head, class... Tail>
            static
            typename std::enable_if<(sizeof...(Tail) > 0), void>::type
            assign_initial(T *p) {
                *p++ = Head::template impl<T>::initial();
                assign_initial<Tail...>(p);
            }

            template <class Head>
            void partial_reduce(T *ptr, const T *a, const T *b) const {
                *ptr++ = typename Head::template impl<T>()(*a++, *b++);
            }

            template <class Head, class... Tail>
            typename std::enable_if<(sizeof...(Tail) > 0), void>::type
            partial_reduce(T *ptr, const T *a, const T *b) const {
                *ptr++ = typename Head::template impl<T>()(*a++, *b++);
                partial_reduce<Tail...>(ptr, a, b);
            }

            static const char* vcmp(int i) {
                static const char *cmp[] = {
                    "x",   "y",   "z",   "w",
                    // CUDA does not provide vector types of width more than 4,
                    // so we can assume the rest is OpenCL-specific:
                    "s4", "s5", "s6", "s7",
                    "s8", "s9", "sa", "sb",
                    "sc", "sd", "se", "sf"
                };

                return cmp[i];
            }
    };
};

/// Combined MIN and MAX operation.
typedef CombineReductors<MIN, MAX> MIN_MAX;
#endif

/// Parallel reduction of arbitrary expression.
/**
 * Reduction uses small temporary buffer on each device present in the queue
 * parameter. One Reductor class for each reduction kind is enough per thread
 * of execution.
 */
template <typename ScalarType, class RDC = SUM>
class Reductor {
    public:
        typedef typename RDC::template impl<ScalarType>::result_type result_type;

        /// Constructor.
        Reductor(const std::vector<backend::command_queue> &queue
#ifndef VEXCL_NO_STATIC_CONTEXT_CONSTRUCTORS
                = current_context().queue()
#endif
                ) : queue(queue) {}

        /// Compute reduction of a vector expression.
        template <class Expr>
        auto operator()(const Expr &expr) const ->
            typename std::enable_if<
                boost::proto::matches<Expr, vector_expr_grammar>::value,
                result_type
            >::type
        {
            using namespace detail;

            static kernel_cache cache;

            auto &data_cache = get_data_cache();

            get_expression_properties prop;
            extract_terminals()(expr, prop);

            result_type initial = RDC::template impl<ScalarType>::initial();

            // If expression is of zero size, then there is nothing to do. Hurray!
            if (prop.size == 0) return initial;

            // Sometimes the expression only knows its size:
            if (prop.size && prop.part.empty())
                prop.part = vex::partition(prop.size, queue);

            for(unsigned d = 0; d < queue.size(); ++d) {
                auto kernel = cache.find(queue[d]);

                backend::select_context(queue[d]);

                if (kernel == cache.end()) {
                    backend::source_generator source(queue[d]);

                    output_terminal_preamble termpream(source, queue[d], "prm", empty_state());
                    boost::proto::eval(boost::proto::as_child(expr),  termpream);

                    typedef typename RDC::template impl<ScalarType>::device_in  fun_in;
                    typedef typename RDC::template impl<ScalarType>::device_out fun_out;
                    boost::proto::eval(boost::proto::as_child( fun_in()( result_type(), ScalarType()) ), termpream);
                    boost::proto::eval(boost::proto::as_child( fun_out()( result_type(), result_type()) ), termpream);

                    source.begin_kernel("vexcl_reductor_kernel");
                    source.begin_kernel_parameters();
                    source.template parameter<size_t>("n");

                    extract_terminals()( expr, declare_expression_parameter(source, queue[d], "prm", empty_state()) );

                    source.template parameter< global_ptr<result_type> >("g_odata");

                    if (!backend::is_cpu(queue[d]))
                        source.template smem_parameter<result_type>();

                    source.end_kernel_parameters();

                    local_sum<Expr, RDC>::get(queue[d], expr, source);

                    if ( backend::is_cpu(queue[d]) ) {
                        source.new_line() << "g_odata[" << source.group_id(0) << "] = mySum;";
                        source.end_kernel();

                        kernel = cache.insert(queue[d], backend::kernel(
                                    queue[d], source.str(), "vexcl_reductor_kernel"));
                    } else {
                        source.smem_declaration<result_type>();
                        source.new_line() << type_name< shared_ptr<result_type> >() << " sdata = smem;";

                        source.new_line() << "size_t tid = " << source.local_id(0) << ";";
                        source.new_line() << "size_t block_size = " << source.local_size(0) << ";";

                        source.new_line() << "sdata[tid] = mySum;";
                        source.new_line().barrier();
                        for(unsigned bs = 512; bs > 0; bs /= 2) {
                            source.new_line() << "if (block_size >= " << bs * 2 << ")";
                            source.open("{").new_line() << "if (tid < " << bs << ") "
                                "{ sdata[tid] = mySum = " << fun_out::name() << "(mySum, sdata[tid + " << bs << "]); }";
                            source.new_line().barrier().close("}");
                        }
                        source.new_line() << "if (tid == 0) g_odata[" << source.group_id(0) << "] = sdata[0];";
                        source.end_kernel();

                        kernel = cache.insert(queue[d], backend::kernel(
                                    queue[d], source.str(), "vexcl_reductor_kernel",
                                    sizeof(ScalarType)));
                    }
                }

                if (size_t psize = prop.part_size(d)) {
                    auto data = data_cache.find(queue[d]);
                    if (data == data_cache.end())
                        data = data_cache.insert(queue[d], reductor_data(queue[d]));

                    kernel->second.push_arg(psize);

                    extract_terminals()(
                            expr,
                            set_expression_argument(kernel->second, d, prop.part_start(d), empty_state())
                            );

                    kernel->second.push_arg(data->second.dbuf);

                    if (!backend::is_cpu(queue[d]))
                        kernel->second.set_smem(
                                [](size_t wgs){
                                return wgs * sizeof(result_type);
                                });

                    kernel->second(queue[d]);
                }
            }

            for(unsigned d = 0; d < queue.size(); d++) {
                if (prop.part_size(d)) {
                    auto data = data_cache.find(queue[d]);

                    data->second.dbuf.read(queue[d], 0, data->second.hbuf.size(), data->second.hbuf.data());
                }
            }

            result_type result = initial;
            typename RDC::template impl<ScalarType> rdc;
            for(unsigned d = 0; d < queue.size(); d++) {
                if (prop.part_size(d)) {
                    auto data = data_cache.find(queue[d]);

                    queue[d].finish();

                    for(auto h = data->second.hbuf.cbegin(),
                             e = data->second.hbuf.cend();
                             h != e; ++h
                       )
                    {
                        result = rdc(result, *h);
                    }
                }
            }

            return result;
        }

        /// Compute reduction of a multivector expression.
        template <class Expr>
#ifdef DOXYGEN
        std::array<result_type, N>
#else
        typename std::enable_if<
            boost::proto::matches<Expr, multivector_expr_grammar>::value &&
            !boost::proto::matches<Expr, vector_expr_grammar>::value,
            std::array<result_type, std::result_of<traits::multiex_dimension(Expr)>::type::value>
        >::type
#endif
        operator()(const Expr &expr) const {
            const size_t dim = std::result_of<traits::multiex_dimension(Expr)>::type::value;
            std::array<result_type, dim> result;

            assign_subexpressions<0, dim, Expr>(result, expr);

            return result;
        }
    private:
        mutable std::vector<backend::command_queue> queue;

        struct reductor_data {
            std::vector<result_type>            hbuf;
            backend::device_vector<result_type> dbuf;

            reductor_data(const backend::command_queue &q)
                : hbuf(backend::kernel::num_workgroups(q)),
                  dbuf(q, backend::kernel::num_workgroups(q))
            { }
        };

        typedef
            detail::object_cache<detail::index_by_queue, reductor_data>
            reductor_data_cache;

        static reductor_data_cache& get_data_cache() {
            static reductor_data_cache cache;
            return cache;
        }

        template <size_t I, size_t N, class Expr>
        typename std::enable_if<I == N, void>::type
        assign_subexpressions(std::array<result_type, N> &, const Expr &) const
        { }

        template <size_t I, size_t N, class Expr>
        typename std::enable_if<(I < N), void>::type
        assign_subexpressions(std::array<result_type, N> &result, const Expr &expr) const
        {
            result[I] = (*this)(detail::extract_subexpression<I>()(expr));

            assign_subexpressions<I + 1, N, Expr>(result, expr);
        }

        template <typename result_type>
        static typename std::enable_if<cl_vector_length<result_type>::value == 1, void>::type
        initial_value(backend::source_generator &src, const result_type &initial) {
            src.new_line() << type_name<result_type>() << " mySum = " << initial << ";";
        }

        template <typename result_type>
        static typename std::enable_if<(cl_vector_length<result_type>::value > 1), void>::type
        initial_value(backend::source_generator &src, const result_type &initial) {
            src.new_line() << type_name<result_type>() << " mySum = {" << initial.s[0];
            for(unsigned i = 1; i < cl_vector_length<result_type>::value; ++i)
                src << ", " << initial.s[i];
            src << "};";
        }

        template <class Expr, class OP>
        struct local_sum {
            static void get(const backend::command_queue &q, const Expr &expr,
                    backend::source_generator &source)
            {
                using namespace detail;

                typedef typename OP::template impl<ScalarType>::device_in fun;
                result_type initial = OP::template impl<ScalarType>::initial();

                initial_value(source, initial);

                source.grid_stride_loop().open("{");

                output_local_preamble loc_init(source, q, "prm", empty_state());
                boost::proto::eval(expr, loc_init);
                source.new_line() << "mySum = " << fun::name() << "(mySum, ";
                vector_expr_context expr_ctx(source, q, "prm", empty_state());
                boost::proto::eval(expr, expr_ctx);
                source << ");";

                source.close("}");
            }
        };

        // http://en.wikipedia.org/wiki/Kahan_summation_algorithm
        template <class Expr>
        struct local_sum<Expr, SUM_Kahan> {
            static void get(const backend::command_queue &q, const Expr &expr,
                    backend::source_generator &source)
            {
                using namespace detail;

                source.new_line()
                    << type_name<result_type>() << " mySum = ("
                    << type_name<result_type>() << ")0, c = ("
                    << type_name<result_type>() << ")0;";
                source.grid_stride_loop().open("{");

                output_local_preamble loc_init(source, q, "prm", empty_state());
                boost::proto::eval(expr, loc_init);

                source.new_line() << type_name<result_type>() << " y = (";
                vector_expr_context expr_ctx(source, q, "prm", empty_state());
                boost::proto::eval(expr, expr_ctx);
                source << ") - c;";

                source.new_line() << type_name<result_type>() << " t = mySum + y;";
                source.new_line() << "c = (t - mySum) - y;";
                source.new_line() << "mySum = t;";

                source.close("}");
            }
        };
};

/// Returns an instance of vex::Reductor<T,R>
/**
 * \deprecated
 */
template <typename T, class R>
vex::Reductor<T, R> get_reductor(const std::vector<backend::command_queue> &queue)
{
    return vex::Reductor<T, R>(queue);
}

} // namespace vex

#endif
