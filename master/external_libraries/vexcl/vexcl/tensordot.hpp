#ifndef VEXCL_TENSORDOT_HPP
#define VEXCL_TENSORDOT_HPP

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
 * \file   vexcl/tensordot.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Tensor dot product along specified axes for multidimensional arrays.
 */

#include <vector>
#include <array>
#include <numeric>

#include <boost/config.hpp>
#ifdef BOOST_NO_VARIADIC_TEMPLATES
#  include <boost/preprocessor/arithmetic/mul.hpp>
#endif

#include <vexcl/operations.hpp>
#include <vexcl/vector_view.hpp>

namespace vex {

struct tensordot_terminal {};

typedef vector_expression<
    typename boost::proto::terminal< tensordot_terminal >::type
    > tensordot_terminal_expression;

template <class LHS, size_t LDIM, class RHS, size_t RDIM, size_t CDIM>
struct tensordot_expr : public tensordot_terminal_expression
{
    typedef typename detail::return_type<
        typename boost::proto::multiplies<LHS, RHS>::type
        >::type value_type;

    const LHS lhs;
    const RHS rhs;

    gslice<LDIM> lslice;
    gslice<RDIM> rslice;

    tensordot_expr(
            const LHS &lhs, const gslice<LDIM> &lsl,
            const RHS &rhs, const gslice<RDIM> &rsl,
            const std::array<std::array<size_t, 2>, CDIM> &common_axes
            )
        : lhs(lhs), rhs(rhs), lslice(lsl), rslice(rsl)
    {
        for(size_t k = 0; k < CDIM; ++k)
            precondition(
                lslice.length[common_axes[k][0]] == rslice.length[common_axes[k][1]],
                "Incompatible common dimensions in tensordot"
                );

        rearrange<0>(lslice, common_axes);
        rearrange<1>(rslice, common_axes);
    }

    // Rearrange strides/lengths so that non-common axes are in the beginning.
    template <size_t SIDE, size_t DIM>
    void rearrange(
            gslice<DIM> &slice,
            const std::array<std::array<size_t, 2>, CDIM> &common_axes
            )
    {
        bool common[DIM];

        for(size_t i = 0; i < DIM; ++i)  common[i] = false;
        for(size_t i = 0; i < CDIM; ++i) common[common_axes[i][SIDE]] = true;

        size_t length[DIM];
        size_t stride[DIM];

        for(size_t i = 0, j = 0, k = DIM - CDIM; i < DIM; ++i) {
            if (common[i]) {
                length[k] = slice.length[i];
                stride[k] = slice.stride[i];
                ++k;
            } else {
                length[j] = slice.length[i];
                stride[j] = slice.stride[i];
                ++j;
            }
        }

        for(size_t i = 0; i < DIM; ++i) {
            slice.length[i] = length[i];
            slice.stride[i] = stride[i];
        }
    }

};

namespace traits {

template <>
struct is_vector_expr_terminal<tensordot_terminal> : std::true_type {};

template <>
struct proto_terminal_is_value<tensordot_terminal> : std::true_type {};

template <class LHS, size_t LDIM, class RHS, size_t RDIM, size_t CDIM>
struct terminal_preamble< tensordot_expr<LHS, LDIM, RHS, RDIM, CDIM> > {
    static void get(backend::source_generator &src,
            const tensordot_expr<LHS, LDIM, RHS, RDIM, CDIM> &term,
            const backend::command_queue &queue, const std::string &prm_name,
            detail::kernel_generator_state_ptr state)
    {
        {
            detail::output_terminal_preamble termpream(src, queue, prm_name + "_lhs", state);
            boost::proto::eval(boost::proto::as_child(term.lhs), termpream);
        }

        {
            detail::output_terminal_preamble termpream(src, queue, prm_name + "_rhs", state);
            boost::proto::eval(boost::proto::as_child(term.rhs), termpream);
        }
    }
};

template <class LHS, size_t LDIM, class RHS, size_t RDIM, size_t CDIM>
struct local_terminal_init< tensordot_expr<LHS, LDIM, RHS, RDIM, CDIM> > {
    template <typename T>
    static typename std::enable_if<is_cl_vector<T>::value, std::string>::type
    init(int v) {
        std::ostringstream buf;
        buf << "{";
        for (size_t i = 0; i < cl_vector_length<T>::value; ++i) {
            if (i) buf << ", ";
            buf << v;
        }
        buf << "}";
        return buf.str();
    }

    template <typename T>
    static typename std::enable_if<!is_cl_vector<T>::value, std::string>::type
    init(int v) {
        return std::to_string(v);
    }

    static void get(backend::source_generator &src,
            const tensordot_expr<LHS, LDIM, RHS, RDIM, CDIM> &term,
            const backend::command_queue &queue, const std::string &prm_name,
            detail::kernel_generator_state_ptr state)
    {
        typedef typename tensordot_expr<LHS, LDIM, RHS, RDIM, CDIM>::value_type T;

        src.new_line() << type_name<T>() << " " << prm_name << "_sum = "
            << init<T>(0) << ";";
        src.open("{");

        src.new_line() << type_name<size_t>() << " lptr0 = " << prm_name << "_lhs_start;";
        src.new_line() << type_name<size_t>() << " rptr0 = " << prm_name << "_rhs_start;";

        src.open("{");
        src.new_line() << type_name<size_t>() << " pos = idx;";

        for(size_t i = RDIM - CDIM; i-->0; ) {
            src.open("{");
            src.new_line() << type_name<size_t>() << " i = pos % " << prm_name << "_rhs_length" << i << ";";
            src.new_line() << "pos /= " << prm_name << "_rhs_length" << i << ";";
            src.new_line() << "rptr0 += i * " << prm_name << "_rhs_stride" << i << ";";
            src.close("}");
        }
        for(size_t i = LDIM - CDIM; i-->0; ) {
            src.open("{");
            src.new_line() << type_name<size_t>() << " i = pos % " << prm_name << "_lhs_length" << i << ";";
            src.new_line() << "pos /= " << prm_name << "_lhs_length" << i << ";";
            src.new_line() << "lptr0 += i * " << prm_name << "_lhs_stride" << i << ";";
            src.close("}");
        }
        src.close("}");

        for(size_t i = 1, il = LDIM - CDIM, ir = RDIM - CDIM; i <= CDIM; ++i, ++il, ++ir) {
            src.new_line() << "for(" << type_name<size_t>() << " i" << i << " = 0,"
                " lptr" << i << " = lptr" << i - 1 << ", rptr" << i << " = rptr" << i - 1 << ";"
                " i" << i << " < " << prm_name << "_lhs_length" << il << "; ++i" << i << ","
                " lptr" << i << " += " << prm_name << "_lhs_stride" << il << ","
                " rptr" << i << " += " << prm_name << "_rhs_stride" << ir << ")";
            src.open("{");
        }

        src.new_line() << type_name<T>() << " prod = " << init<T>(1) << ";";

        {
            src.open("{");
            src.new_line() << type_name<size_t>() << " idx = lptr" << CDIM << ";";

            detail::output_local_preamble init_ctx(src, queue, prm_name + "_lhs", state);
            boost::proto::eval(boost::proto::as_child(term.lhs), init_ctx);

            src.new_line() << "prod *= ";

            detail::vector_expr_context expr_ctx(src, queue, prm_name + "_lhs", state);
            boost::proto::eval(boost::proto::as_child(term.lhs), expr_ctx);

            src << ";";
            src.close("}");
        }

        {
            src.open("{");
            src.new_line() << type_name<size_t>() << " idx = rptr" << CDIM << ";";

            detail::output_local_preamble init_ctx(src, queue, prm_name + "_rhs", state);
            boost::proto::eval(boost::proto::as_child(term.rhs), init_ctx);

            src.new_line() << "prod *= ";

            detail::vector_expr_context expr_ctx(src, queue, prm_name + "_rhs", state);
            boost::proto::eval(boost::proto::as_child(term.rhs), expr_ctx);

            src << ";";
            src.close("}");
        }

        src.new_line() << prm_name << "_sum += prod;";

        for(size_t i = 1; i <= CDIM; ++i) src.close("}");

        src.close("}");
    }
};

template <class LHS, size_t LDIM, class RHS, size_t RDIM, size_t CDIM>
struct kernel_param_declaration< tensordot_expr<LHS, LDIM, RHS, RDIM, CDIM> > {
    static void get(backend::source_generator &src,
            const tensordot_expr<LHS, LDIM, RHS, RDIM, CDIM> &term,
            const backend::command_queue &queue, const std::string &prm_name,
            detail::kernel_generator_state_ptr state)
    {
        {
            detail::declare_expression_parameter declare(src, queue, prm_name + "_lhs", state);
            detail::extract_terminals()(boost::proto::as_child(term.lhs), declare);
            term.lslice.parameter_declaration(src, prm_name + "_lhs", queue, state);
        }

        {
            detail::declare_expression_parameter declare(src, queue, prm_name + "_rhs", state);
            detail::extract_terminals()(boost::proto::as_child(term.rhs), declare);
            term.rslice.parameter_declaration(src, prm_name + "_rhs", queue, state);
        }
    }
};

template <class LHS, size_t LDIM, class RHS, size_t RDIM, size_t CDIM>
struct partial_vector_expr< tensordot_expr<LHS, LDIM, RHS, RDIM, CDIM> > {
    static void get(backend::source_generator &src,
            const tensordot_expr<LHS, LDIM, RHS, RDIM, CDIM>&,
            const backend::command_queue&, const std::string &prm_name,
            detail::kernel_generator_state_ptr)
    {
        src << prm_name << "_sum";
    }
};

template <class LHS, size_t LDIM, class RHS, size_t RDIM, size_t CDIM>
struct kernel_arg_setter< tensordot_expr<LHS, LDIM, RHS, RDIM, CDIM> > {
    static void set(const tensordot_expr<LHS, LDIM, RHS, RDIM, CDIM> &term,
            backend::kernel &kernel, unsigned part, size_t index_offset,
            detail::kernel_generator_state_ptr state)
    {
        detail::set_expression_argument setarg(kernel, part, index_offset, state);

        detail::extract_terminals()( boost::proto::as_child(term.lhs), setarg);
        term.lslice.setArgs(kernel, part, index_offset, state);

        detail::extract_terminals()( boost::proto::as_child(term.rhs), setarg);
        term.rslice.setArgs(kernel, part, index_offset, state);
    }
};

template <class LHS, size_t LDIM, class RHS, size_t RDIM, size_t CDIM>
struct expression_properties< tensordot_expr<LHS, LDIM, RHS, RDIM, CDIM> > {
    static void get(const tensordot_expr<LHS, LDIM, RHS, RDIM, CDIM> &term,
            std::vector<backend::command_queue> &queue_list,
            std::vector<size_t> &partition,
            size_t &size
            )
    {
        detail::get_expression_properties prop;
        detail::extract_terminals()(boost::proto::as_child(term.lhs), prop);
        // Sizes of lhs and rhs expressions most probably differ.
        // Avoid triggering size check exception:
        prop.size = 0;
        detail::extract_terminals()(boost::proto::as_child(term.rhs), prop);

        queue_list = prop.queue;
        partition  = std::vector<size_t>(2, 0);
        size       =
            std::accumulate(
                    term.lslice.length,
                    term.lslice.length + LDIM - CDIM,
                    static_cast<size_t>(1),
                    std::multiplies<size_t>()
                    ) *
            std::accumulate(
                    term.rslice.length,
                    term.rslice.length + RDIM - CDIM,
                    static_cast<size_t>(1),
                    std::multiplies<size_t>()
                    );

        partition.back() = size;
    }
};

} // namespace traits


#ifndef BOOST_NO_VARIADIC_TEMPLATES
/// Helper function for creating axes pairs.
/**
Example:
\code
auto axes = axes_pairs(a0, b0, a1, b1);
assert(axes[0][0] == a0 && axes[0][1] == b0);
assert(axes[1][0] == a1 && axes[1][1] == b1);
\endcode
 */
template <class... Args>
std::array< std::array<size_t, 2>, sizeof...(Args) / 2>
axes_pairs(Args... args) {
    static_assert(
            sizeof...(Args) % 2 == 0,
            "Odd number of arguments in axes_pairs"
            );

    return std::array<std::array<size_t, 2>, sizeof...(Args) / 2>{{static_cast<size_t>(args)...}};
}
#else

#define VEXCL_INIT_ARRAY(z, n, data) static_cast<size_t>(t ## n)
#define VEXCL_AXES_PAIRS(z, n, data)                                           \
  template <BOOST_PP_ENUM_PARAMS(BOOST_PP_MUL(n, 2), class T)>                 \
  std::array< std::array<size_t, 2>, n>                                        \
  axes_pairs(BOOST_PP_ENUM_BINARY_PARAMS(BOOST_PP_MUL(n, 2), T, t)) {          \
    std::array< std::array<size_t, 2>, n> a = { {                              \
        BOOST_PP_ENUM(BOOST_PP_MUL(n, 2), VEXCL_INIT_ARRAY, ~) } };            \
    return a;                                                                  \
  }

BOOST_PP_REPEAT_FROM_TO(1, VEXCL_MAX_ARITY, VEXCL_AXES_PAIRS, ~)

#undef VEXCL_AXES_PAIRS
#undef VEXCL_INIT_ARRAY

#endif

/// Tensor dot product along specified axes for multidimensional arrays.
#ifdef DOXYGEN
template <class SlicedExpr1, size_t SlicedExpr2, size_t CDIM>
auto tensordot(
        const SlicedExpr1 &lhs,
        const SlicedExpr2 &rhs,
        const std::array<std::array<size_t, 2>, CDIM> &common_axes
        );
#else
template <class LHS, size_t LDIM, class RHS, size_t RDIM, size_t CDIM>
tensordot_expr<LHS, LDIM, RHS, RDIM, CDIM> tensordot(
        const vector_view<LHS, gslice<LDIM>> &lhs,
        const vector_view<RHS, gslice<RDIM>> &rhs,
        const std::array<std::array<size_t, 2>, CDIM> &common_axes
        )
{
    return tensordot_expr<LHS, LDIM, RHS, RDIM, CDIM>(
            lhs.expr, lhs.slice, rhs.expr, rhs.slice, common_axes
            );
}
#endif

} // namespace vex

#endif
