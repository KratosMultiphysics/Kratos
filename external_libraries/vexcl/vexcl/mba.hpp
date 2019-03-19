#ifndef VEXCL_MBA_HPP
#define VEXCL_MBA_HPP

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
 * \file   vexcl/mba.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Scattered data interpolation with multilevel B-Splines.
 */

#include <vector>
#include <array>
#include <sstream>
#include <memory>
#include <algorithm>
#include <numeric>
#include <type_traits>
#include <cassert>

#include <boost/tuple/tuple.hpp>
#include <boost/fusion/adapted/boost_tuple.hpp>

#include <vexcl/operations.hpp>

// Include boost.preprocessor header if variadic templates are not available.
// Also include it if we use gcc v4.6.
// This is required due to bug http://gcc.gnu.org/bugzilla/show_bug.cgi?id=35722
#if defined(BOOST_NO_VARIADIC_TEMPLATES) || (defined(__GNUC__) && !defined(__clang__) && __GNUC__ == 4 && __GNUC_MINOR__ == 6)
#  include <boost/preprocessor/repetition.hpp>
#  ifndef VEXCL_MAX_ARITY
#    define VEXCL_MAX_ARITY BOOST_PROTO_MAX_ARITY
#  endif
#endif
namespace vex {

struct mba_terminal {};

typedef vector_expression<
    typename boost::proto::terminal< mba_terminal >::type
    > mba_terminal_expression;

template <class MBA, class ExprTuple>
struct mba_interp : public mba_terminal_expression {
    typedef typename MBA::value_type value_type;

    const MBA      &cloud;
    const ExprTuple coord;

    mba_interp(const MBA &cloud, const ExprTuple coord)
        : cloud(cloud), coord(coord) {}
};

namespace detail {
    // Compile time value of N^M.
    template <size_t N, size_t M>
    struct power : std::integral_constant<size_t, N * power<N, M-1>::value> {};

    template <size_t N>
    struct power<N, 0> : std::integral_constant<size_t, 1> {};

    // Nested loop counter of compile-time size (M loops of size N).
    template <size_t N, size_t M>
    class scounter {
        public:
            scounter() : idx(0) {
                std::fill(i.begin(), i.end(), static_cast<size_t>(0));
            }

            size_t operator[](size_t d) const {
                return i[d];
            }

            scounter& operator++() {
                for(size_t d = M; d--; ) {
                    if (++i[d] < N) break;
                    i[d] = 0;
                }

                ++idx;

                return *this;
            }

            operator size_t() const {
                return idx;
            }

            bool valid() const {
                return idx < power<N, M>::value;
            }
        private:
            size_t idx;
            std::array<size_t, M> i;
    };

    // Nested loop counter of run-time size (M loops of given sizes).
    template <size_t M>
    class dcounter {
        public:
            dcounter(const std::array<size_t, M> &N)
                : idx(0),
                  size(std::accumulate(N.begin(), N.end(),
                            static_cast<size_t>(1), std::multiplies<size_t>())),
                  N(N)
            {
                std::fill(i.begin(), i.end(), static_cast<size_t>(0));
            }

            size_t operator[](size_t d) const {
                return i[d];
            }

            dcounter& operator++() {
                for(size_t d = M; d--; ) {
                    if (++i[d] < N[d]) break;
                    i[d] = 0;
                }

                ++idx;

                return *this;
            }

            operator size_t() const {
                return idx;
            }

            bool valid() const {
                return idx < size;
            }
        private:
            size_t idx, size;
            std::array<size_t, M> N, i;
    };
} // namespace detail

/// Scattered data interpolation with multilevel B-Splines.
template <size_t NDIM, typename real = double>
class mba {
    public:
        typedef real value_type;
        typedef std::array<real,   NDIM> point;
        typedef std::array<size_t, NDIM> index;

        static const size_t ndim = NDIM;

        std::vector< backend::command_queue >       queue;
        std::vector< backend::device_vector<real> > phi;
        point xmin, hinv;
        index n, stride;

        /** Creates the approximation functor.
         * `cmin` and `cmax` specify the domain boundaries, `coo` and `val`
         * contain coordinates and values of the data points. `grid` is the
         * initial control grid size. The approximation hierarchy will have at
         * most `levels` and will stop when the desired approximation precision
         * `tol` will be reached.
         */
        mba(
                const std::vector<backend::command_queue> &queue,
                const point &cmin, const point &cmax,
                const std::vector<point> &coo, std::vector<real> val,
                std::array<size_t, NDIM> grid, size_t levels = 8, real tol = 1e-8
           ) : queue(queue)
        {
            init(cmin, cmax, coo.begin(), coo.end(), val.begin(), grid, levels, tol);
        }

        /** Creates the approximation functor.
         * `cmin` and `cmax` specify the domain boundaries. Coordinates and
         * values of the data points are passed as iterator ranges. `grid` is
         * the initial control grid size. The approximation hierarchy will have
         * at most `levels` and will stop when the desired approximation
         * precision `tol` will be reached.
         */
        template <class CooIter, class ValIter>
        mba(
                const std::vector<backend::command_queue> &queue,
                const point &cmin, const point &cmax,
                CooIter coo_begin, CooIter coo_end, ValIter val_begin,
                std::array<size_t, NDIM> grid, size_t levels = 8, real tol = 1e-8
           ) : queue(queue)
        {
            init(cmin, cmax, coo_begin, coo_end, val_begin, grid, levels, tol);
        }

#if !defined(BOOST_NO_VARIADIC_TEMPLATES) && ((!defined(__GNUC__) || (__GNUC__ > 4 || __GNUC__ == 4 && __GNUC_MINOR__ > 6)) || defined(__clang__))
        /// Provide interpolated values at given coordinates.
        template <class... Expr>
        auto operator()(const Expr&... expr) const ->
            mba_interp< mba, boost::tuple<const Expr&...> >
        {
            static_assert(sizeof...(Expr) == NDIM, "Wrong number of parameters");
            return mba_interp< mba, boost::tuple<const Expr&...> >(*this, boost::tie(expr...));
        }
#else

#define VEXCL_FUNCALL_OPERATOR(z, n, data)                                     \
  template <BOOST_PP_ENUM_PARAMS(n, class Expr)>                               \
  mba_interp<mba, boost::tuple<BOOST_PP_ENUM_BINARY_PARAMS(                    \
                      n, const Expr, &BOOST_PP_INTERCEPT)> >                   \
  operator()(BOOST_PP_ENUM_BINARY_PARAMS(n, const Expr, &expr)) {              \
    return mba_interp<mba, boost::tuple<BOOST_PP_ENUM_BINARY_PARAMS(           \
                               n, const Expr, &BOOST_PP_INTERCEPT)> >(         \
        *this, boost::tie(BOOST_PP_ENUM_PARAMS(n, expr)));                     \
  }

BOOST_PP_REPEAT_FROM_TO(1, 10, VEXCL_FUNCALL_OPERATOR, ~)

#undef VEXCL_FUNCALL_OPERATOR
#endif
    private:
        template <class CooIter, class ValIter>
        void init(
                const point &cmin, const point &cmax,
                CooIter coo_begin, CooIter coo_end, ValIter val_begin,
                std::array<size_t, NDIM> grid, size_t levels = 8, real tol = 1e-8
                )
        {
            for(size_t k = 0; k < NDIM; ++k)
                assert(grid[k] > 1);

            double res0 = std::accumulate(
                    val_begin, val_begin + (coo_end - coo_begin),
                    static_cast<real>(0),
                    [](real sum, real v) { return sum + v * v; }
                    );

            std::unique_ptr<lattice> psi(
                    new lattice(cmin, cmax, grid, coo_begin, coo_end, val_begin)
                    );
            double res = psi->update_data(coo_begin, coo_end, val_begin);
#ifdef VEXCL_MBA_VERBOSE
            std::cout << "level  0: res = " << std::scientific << res << std::endl;
#endif

            for (size_t k = 1; (res > res0 * tol) && (k < levels); ++k) {
                for(size_t d = 0; d < NDIM; ++d) grid[d] = 2 * grid[d] - 1;

                std::unique_ptr<lattice> f(
                        new lattice(cmin, cmax, grid, coo_begin, coo_end, val_begin)
                        );
                res = f->update_data(coo_begin, coo_end, val_begin);
#ifdef VEXCL_MBA_VERBOSE
                std::cout << "level " << k << std::scientific << ": res = " << res << std::endl;
#endif

                f->append_refined(*psi);
                psi = std::move(f);
            }

            xmin   = psi->xmin;
            hinv   = psi->hinv;
            n      = psi->n;
            stride = psi->stride;

            phi.reserve(queue.size());

            for(auto q = queue.begin(); q != queue.end(); ++q)
                phi.push_back( backend::device_vector<real>(
                            *q, psi->phi.size(), psi->phi.data(), backend::MEM_READ_ONLY
                            ) );
        }

        // Control lattice.
        struct lattice {
            point xmin, hinv;
            index n, stride;
            std::vector<real> phi;

            template <class CooIter, class ValIter>
            lattice(
                    const point &cmin, const point &cmax, std::array<size_t, NDIM> grid,
                    CooIter coo_begin, CooIter coo_end, ValIter val_begin
                   ) : xmin(cmin), n(grid)
            {
                for(size_t d = 0; d < NDIM; ++d) {
                    hinv[d] = (grid[d] - 1) / (cmax[d] - cmin[d]);
                    xmin[d] -= 1 / hinv[d];
                    n[d]    += 2;
                }

                stride[NDIM - 1] = 1;
                for(size_t d = NDIM - 1; d--; )
                    stride[d] = stride[d + 1] * n[d + 1];

                std::vector<real> delta(n[0] * stride[0], 0.0);
                std::vector<real> omega(n[0] * stride[0], 0.0);

                auto p = coo_begin;
                auto v = val_begin;
                for(; p != coo_end; ++p, ++v) {
                    if (!contained(cmin, cmax, *p)) continue;

                    index i;
                    point s;

                    for(size_t d = 0; d < NDIM; ++d) {
                        real u = ((*p)[d] - xmin[d]) * hinv[d];
                        i[d] = static_cast<size_t>(std::floor(u) - 1);
                        s[d] = u - std::floor(u);
                    }

                    std::array<real, detail::power<4, NDIM>::value> w;
                    real sw2 = 0;

                    for(detail::scounter<4, NDIM> d; d.valid(); ++d) {
                        real buf = 1;
                        for(size_t k = 0; k < NDIM; ++k)
                            buf *= B(d[k], s[k]);

                        w[d] = buf;
                        sw2 += buf * buf;
                    }

                    for(detail::scounter<4, NDIM> d; d.valid(); ++d) {
                        real phi = (*v) * w[d] / sw2;

                        size_t idx = 0;
                        for(size_t k = 0; k < NDIM; ++k) {
                            assert(i[k] + d[k] < n[k]);

                            idx += (i[k] + d[k]) * stride[k];
                        }

                        real w2 = w[d] * w[d];

                        assert(idx < delta.size());

                        delta[idx] += w2 * phi;
                        omega[idx] += w2;
                    }
                }

                phi.resize(omega.size());

                for(auto w = omega.begin(), d = delta.begin(), f = phi.begin();
                        w != omega.end();
                        ++w, ++d, ++f
                   )
                {
                    if (std::fabs(*w) < 1e-32)
                        *f = 0;
                    else
                        *f = (*d) / (*w);
                }
            }

            // Get interpolated value at given position.
            real operator()(const point &p) const {
                index i;
                point s;

                for(size_t d = 0; d < NDIM; ++d) {
                    real u = (p[d] - xmin[d]) * hinv[d];
                    i[d] = static_cast<size_t>(std::floor(u) - 1);
                    s[d] = u - std::floor(u);
                }

                real f = 0;

                for(detail::scounter<4, NDIM> d; d.valid(); ++d) {
                    real w = 1;
                    for(size_t k = 0; k < NDIM; ++k)
                        w *= B(d[k], s[k]);

                    f += w * get(i, d);
                }

                return f;
            }

            // Subtract interpolated values from data points.
            template <class CooIter, class ValIter>
            real update_data(
                    CooIter coo_begin, CooIter coo_end, ValIter val_begin
                    ) const
            {
                auto c = coo_begin;
                auto v = val_begin;

                real res = 0;

                for(; c != coo_end; ++c, ++v) {
                    *v -= (*this)(*c);

                    res += (*v) * (*v);
                }

                return res;
            }

            // Refine r and append it to the current control lattice.
            void append_refined(const lattice &r) {
                static const std::array<real, 5> s = {{
                    0.125, 0.500, 0.750, 0.500, 0.125
                }};

                for(detail::dcounter<NDIM> i(r.n); i.valid(); ++i) {
                    real f = r.phi[i];
                    for(detail::scounter<5, NDIM> d; d.valid(); ++d) {
                        index j;
                        bool skip = false;
                        size_t idx = 0;
                        for(size_t k = 0; k < NDIM; ++k) {
                            j[k] = 2 * i[k] + d[k] - 3;
                            if (j[k] >= n[k]) { skip = true; break; }

                            idx += j[k] * stride[k];
                        }

                        if (skip) continue;

                        real c = 1;
                        for(size_t k = 0; k < NDIM; ++k) c *= s[d[k]];

                        phi[idx] += f * c;
                    }
                }
            }

            private:
                // Value of k-th B-Spline at t.
                static inline real B(size_t k, real t) {
                    assert(0 <= t && t < 1);
                    assert(k < 4);

                    switch (k) {
                        case 0:
                            return (t * (t * (-t + 3) - 3) + 1) / 6;
                        case 1:
                            return (t * t * (3 * t - 6) + 4) / 6;
                        case 2:
                            return (t * (t * (-3 * t + 3) + 3) + 1) / 6;
                        case 3:
                        default:
                            return t * t * t / 6;
                    }
                }

                // x is within [xmin, xmax].
                static bool contained(
                        const point &xmin, const point &xmax, const point &x)
                {
                    for(size_t d = 0; d < NDIM; ++d) {
                        static const real eps = 1e-12;

                        if (x[d] - eps <  xmin[d]) return false;
                        if (x[d] + eps >= xmax[d]) return false;
                    }

                    return true;
                }

                // Get value of phi at index (i + d).
                template <class Shift>
                inline real get(const index &i, const Shift &d) const {
                    size_t idx = 0;

                    for(size_t k = 0; k < NDIM; ++k) {
                        size_t j = i[k] + d[k];

                        if (j >= n[k]) return 0;
                        idx += j * stride[k];
                    }

                    return phi[idx];
                }
        };
};

namespace traits {

template <>
struct is_vector_expr_terminal< mba_terminal > : std::true_type {};

template <>
struct proto_terminal_is_value< mba_terminal > : std::true_type {};

template <class MBA, class ExprTuple>
struct terminal_preamble< mba_interp<MBA, ExprTuple> > {
    static void get(backend::source_generator &src,
            const mba_interp<MBA, ExprTuple>&,
            const backend::command_queue&, const std::string &prm_name,
            detail::kernel_generator_state_ptr)
    {
        typedef typename MBA::value_type real;

        std::string B = prm_name + "_B";

        src.begin_function<real>(B + "0");
        src.begin_function_parameters();
        src.template parameter<real>("t");
        src.end_function_parameters();
        src.new_line() << "return (t * (t * (-t + 3) - 3) + 1) / 6;";
        src.end_function();

        src.begin_function<real>(B + "1");
        src.begin_function_parameters();
        src.template parameter<real>("t");
        src.end_function_parameters();
        src.new_line() << "return (t * t * (3 * t - 6) + 4) / 6;";
        src.end_function();

        src.begin_function<real>(B + "2");
        src.begin_function_parameters();
        src.template parameter<real>("t");
        src.end_function_parameters();
        src.new_line() << "return (t * (t * (-3 * t + 3) + 3) + 1) / 6;";
        src.end_function();

        src.begin_function<real>(B + "3");
        src.begin_function_parameters();
        src.template parameter<real>("t");
        src.end_function_parameters();
        src.new_line() << "return t * t * t / 6;";
        src.end_function();

        src.begin_function<real>(prm_name + "_mba");
        src.begin_function_parameters();

        for(size_t k = 0; k < MBA::ndim; ++k)
            src.template parameter<real>("x") << k;

        for(size_t k = 0; k < MBA::ndim; ++k) {
            src.template parameter<real>("c" + std::to_string(k));
            src.template parameter<real>("h" + std::to_string(k));
            src.parameter<size_t>("n" + std::to_string(k));
            src.parameter<size_t>("m" + std::to_string(k));
        }

        src.template parameter< global_ptr<const real> >("phi");
        src.end_function_parameters();
        src.new_line() << type_name<real>() << " u;";
        for(size_t k = 0; k < MBA::ndim; ++k) {
            src.new_line() << "u = (x" << k << " - c" << k << ") * h" << k << ";";
            src.new_line() << type_name<size_t>() << " i" << k << " = floor(u) - 1;";
            src.new_line() << type_name<real>() << " s" << k << " = u - floor(u);";
        }
        src.new_line() << type_name<real>() << " f = 0;";
        src.new_line() << type_name<size_t>() << " j, idx;";

        for(detail::scounter<4,MBA::ndim> d; d.valid(); ++d) {
            src.new_line() << "idx = 0;";
            for(size_t k = 0; k < MBA::ndim; ++k) {
                src.new_line() << "j = i" << k << " + " << d[k] << ";";
                src.new_line() << "if (j < n" << k << ")";
                src.open("{").new_line() << "idx += j * m" << k << ";";
            }

            src.new_line() << "f += ";
            for(size_t k = 0; k < MBA::ndim; ++k) {
                if (k) src << " * ";
                src << B << d[k] << "(s" << k << ")";
            }

            src << " * phi[idx];";

            for(size_t k = 0; k < MBA::ndim; ++k)
                src.close("}");
        }
        src.new_line() << "return f;";
        src.end_function();
    }
};

template <class MBA, class ExprTuple>
struct kernel_param_declaration< mba_interp<MBA, ExprTuple> > {
    static void get(backend::source_generator &src,
            const mba_interp<MBA, ExprTuple> &term,
            const backend::command_queue &queue, const std::string &prm_name,
            detail::kernel_generator_state_ptr state)
    {
        typedef typename MBA::value_type real;
        boost::fusion::for_each(term.coord, prmdecl(src, queue, prm_name, state));

        for(size_t k = 0; k < MBA::ndim; ++k) {
            src.template parameter<real>(prm_name + "_c" + std::to_string(k));
            src.template parameter<real>(prm_name + "_h" + std::to_string(k));
            src.parameter<size_t>(prm_name + "_n" + std::to_string(k));
            src.parameter<size_t>(prm_name + "_m" + std::to_string(k));
        }

        src.parameter< global_ptr<const real> >(prm_name + "_phi");
    }

    struct prmdecl {
        backend::source_generator &s;
        const backend::command_queue &queue;
        const std::string &prm_name;
        detail::kernel_generator_state_ptr state;
        mutable int pos;

        prmdecl(backend::source_generator &s,
                const backend::command_queue &queue, const std::string &prm_name,
                detail::kernel_generator_state_ptr state
            ) : s(s), queue(queue), prm_name(prm_name), state(state), pos(0)
        {}

        template <class Expr>
        void operator()(const Expr &expr) const {
            std::ostringstream prefix;
            prefix << prm_name << "_x" << pos;
            detail::declare_expression_parameter ctx(s, queue, prefix.str(), state);
            detail::extract_terminals()(boost::proto::as_child(expr), ctx);

            pos++;
        }
    };
};

template <class MBA, class ExprTuple>
struct local_terminal_init< mba_interp<MBA, ExprTuple> > {
    static void get(backend::source_generator &src,
            const mba_interp<MBA, ExprTuple> &term,
            const backend::command_queue &queue, const std::string &prm_name,
            detail::kernel_generator_state_ptr state)
    {
        boost::fusion::for_each(term.coord, local_init(src, queue, prm_name, state));
    }

    struct local_init {
        backend::source_generator &s;
        const backend::command_queue &queue;
        const std::string &prm_name;
        detail::kernel_generator_state_ptr state;
        mutable int pos;

        local_init(backend::source_generator &s,
                const backend::command_queue &queue, const std::string &prm_name,
                detail::kernel_generator_state_ptr state
            ) : s(s), queue(queue), prm_name(prm_name), state(state), pos(0)
        {}

        template <class Expr>
        void operator()(const Expr &expr) const {
            std::ostringstream prefix;
            prefix << prm_name << "_x" << pos;

            detail::output_local_preamble init_ctx(s, queue, prefix.str(), state);
            boost::proto::eval(boost::proto::as_child(expr), init_ctx);

            pos++;
        }
    };
};

template <class MBA, class ExprTuple>
struct partial_vector_expr< mba_interp<MBA, ExprTuple> > {
    static void get(backend::source_generator &src,
            const mba_interp<MBA, ExprTuple> &term,
            const backend::command_queue &queue, const std::string &prm_name,
            detail::kernel_generator_state_ptr state)
    {
        src << prm_name << "_mba(";

        boost::fusion::for_each(term.coord, buildexpr(src, queue, prm_name, state));

        for(size_t k = 0; k < MBA::ndim; ++k) {
            src << ", " << prm_name << "_c" << k
                << ", " << prm_name << "_h" << k
                << ", " << prm_name << "_n" << k
                << ", " << prm_name << "_m" << k;
        }

        src << ", " << prm_name << "_phi)";
    }

    struct buildexpr {
        backend::source_generator &s;
        const backend::command_queue &queue;
        const std::string &prm_name;
        detail::kernel_generator_state_ptr state;
        mutable int pos;

        buildexpr(backend::source_generator &s,
                const backend::command_queue &queue, const std::string &prm_name,
                detail::kernel_generator_state_ptr state
            ) : s(s), queue(queue), prm_name(prm_name), state(state), pos(0)
        {}

        template <class Expr>
        void operator()(const Expr &expr) const {
            if(pos) s << ", ";

            std::ostringstream prefix;
            prefix << prm_name << "_x" << pos;

            detail::vector_expr_context ctx(s, queue, prefix.str(), state);
            boost::proto::eval(boost::proto::as_child(expr), ctx);

            pos++;
        }
    };
};

template <class MBA, class ExprTuple>
struct kernel_arg_setter< mba_interp<MBA, ExprTuple> > {
    static void set(const mba_interp<MBA, ExprTuple> &term,
            backend::kernel &kernel, unsigned part, size_t index_offset,
            detail::kernel_generator_state_ptr state)
    {

        boost::fusion::for_each(term.coord,
                setargs(kernel, part, index_offset, state));

        for(size_t k = 0; k < MBA::ndim; ++k) {
            kernel.push_arg(term.cloud.xmin[k]);
            kernel.push_arg(term.cloud.hinv[k]);
            kernel.push_arg(term.cloud.n[k]);
            kernel.push_arg(term.cloud.stride[k]);
        }
        kernel.push_arg(term.cloud.phi[part]);
    }

    struct setargs {
        backend::kernel &kernel;
        unsigned part;
        size_t index_offset;
        detail::kernel_generator_state_ptr state;

        setargs(
                backend::kernel &kernel, unsigned part, size_t index_offset,
                detail::kernel_generator_state_ptr state
               )
            : kernel(kernel), part(part), index_offset(index_offset), state(state)
        {}

        template <class Expr>
        void operator()(const Expr &expr) const {
            detail::set_expression_argument ctx(kernel, part, index_offset, state);
            detail::extract_terminals()( boost::proto::as_child(expr), ctx);
        }
    };
};

template <class MBA, class ExprTuple>
struct expression_properties< mba_interp<MBA, ExprTuple> > {
    static void get(const mba_interp<MBA, ExprTuple> &term,
            std::vector<backend::command_queue> &queue_list,
            std::vector<size_t> &partition,
            size_t &size
            )
    {
        boost::fusion::for_each(term.coord, extrprop(queue_list, partition, size));
    }

    struct extrprop {
        std::vector<backend::command_queue> &queue_list;
        std::vector<size_t> &partition;
        size_t &size;

        extrprop(std::vector<backend::command_queue> &queue_list,
            std::vector<size_t> &partition, size_t &size
            ) : queue_list(queue_list), partition(partition), size(size)
        {}

        template <class Expr>
        void operator()(const Expr &expr) const {
            if (queue_list.empty()) {
                detail::get_expression_properties prop;
                detail::extract_terminals()(boost::proto::as_child(expr), prop);

                queue_list = prop.queue;
                partition  = prop.part;
                size       = prop.size;
            }
        }
    };
};

} //namespace traits

} // namespace vex


#endif
