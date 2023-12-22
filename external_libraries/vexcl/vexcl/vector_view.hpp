#ifndef VEXCL_VECTOR_VIEW_HPP
#define VEXCL_VECTOR_VIEW_HPP

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
 * \file   vexcl/vector_view.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Provides sub-view for an existing vex::vector.
 */

#include <vector>
#include <array>
#include <string>
#include <sstream>
#include <numeric>
#include <algorithm>

#include <vexcl/vector.hpp>
#include <vexcl/element_index.hpp>
#include <vexcl/tagged_terminal.hpp>
#if (VEXCL_CHECK_SIZES > 1)
#  include <vexcl/reductor.hpp>
#endif

#include <boost/fusion/container/vector.hpp>
#include <boost/fusion/container/vector/convert.hpp>
#include <boost/fusion/algorithm/transformation/push_back.hpp>
#include <boost/fusion/include/value_at.hpp>

namespace vex {

struct vector_view_terminal {};

typedef vector_expression<
    typename boost::proto::terminal< vector_view_terminal >::type
    > vector_view_terminal_expression;

template <class Expr, class Slice>
struct vector_view : public vector_view_terminal_expression
{
    typedef typename detail::return_type<Expr>::type value_type;
    const Expr   expr;
    const Slice  slice;

    vector_view(const Expr &expr, const Slice &slice)
        : expr(expr), slice(slice)
    { }

    const vector_view& operator=(const vector_view & other) {
        detail::assign_expression<assign::SET>(*this, other);
        return *this;
    }

#define VEXCL_ASSIGNMENT(cop, op)                                              \
  template <class RHS>                                                         \
  typename std::enable_if<                                                     \
      boost::proto::matches<                                                   \
          typename boost::proto::result_of::as_expr<RHS>::type,                \
          vector_expr_grammar>::value,                                         \
      const vector_view &>::type operator cop(const RHS & rhs) {               \
    detail::assign_expression<op>(*this, rhs);                                 \
    return *this;                                                              \
  }

    VEXCL_ASSIGNMENTS(VEXCL_ASSIGNMENT)

#undef VEXCL_ASSIGNMENT
};

// Allow vector_view to participate in vector expressions:
namespace traits {

template <>
struct is_vector_expr_terminal< vector_view_terminal > : std::true_type {};

template <>
struct proto_terminal_is_value< vector_view_terminal > : std::true_type {};

template <typename Expr, class Slice>
struct terminal_preamble< vector_view<Expr, Slice> > {
    static void get(backend::source_generator &src,
            const vector_view<Expr, Slice> &term,
            const backend::command_queue &queue, const std::string &prm_name,
            detail::kernel_generator_state_ptr state)
    {
        detail::output_terminal_preamble termpream(src, queue, prm_name + "_expr", state);
        boost::proto::eval(boost::proto::as_child(term.expr), termpream);

        term.slice.preamble(src, prm_name + "_slice", queue, state);
    }
};

template <typename Expr, class Slice>
struct kernel_param_declaration< vector_view<Expr, Slice> > {
    static void get(backend::source_generator &src,
            const vector_view<Expr, Slice> &term,
            const backend::command_queue &queue, const std::string &prm_name,
            detail::kernel_generator_state_ptr state)
    {
        detail::declare_expression_parameter declare(src, queue, prm_name + "_expr", state);
        detail::extract_terminals()(boost::proto::as_child(term.expr),  declare);

        term.slice.parameter_declaration(src, prm_name + "_slice", queue, state);
    }
};

template <typename Expr, class Slice>
struct local_terminal_init< vector_view<Expr, Slice> > {
    static void get(backend::source_generator &src,
            const vector_view<Expr, Slice> &term,
            const backend::command_queue &queue, const std::string &prm_name,
            detail::kernel_generator_state_ptr state)
    {
        term.slice.local_preamble(src, prm_name + "_slice", queue, state);

        src.new_line()
            << type_name<typename vector_view<Expr, Slice>::value_type>()
            << " " << prm_name << "_val;";

        src.open("{").new_line()
          << "size_t pos = ";

        term.slice.index(src, prm_name + "_slice", queue, state);
        src << ";";
        src.new_line() << "size_t idx = pos;";

        detail::output_local_preamble init_ctx(src, queue, prm_name + "_expr", state);
        boost::proto::eval(boost::proto::as_child(term.expr), init_ctx);

        src.new_line() << prm_name << "_val = ";
        detail::vector_expr_context ctx(src, queue, prm_name + "_expr", state);
        boost::proto::eval(boost::proto::as_child(term.expr), ctx);
        src << ";";
        src.close("}");
    }
};

template <typename T, class Slice>
struct local_terminal_init< vector_view<const vector<T>&, Slice> > {
    static void get(backend::source_generator &src,
            const vector_view<const vector<T>&, Slice> &term,
            const backend::command_queue &queue, const std::string &prm_name,
            detail::kernel_generator_state_ptr state)
    {
        term.slice.local_preamble(src, prm_name + "_slice", queue, state);
    }
};

template <typename T, size_t Tag, class Slice>
struct local_terminal_init< vector_view<tagged_terminal<Tag, const vector<T>&>, Slice> > {
    static void get(backend::source_generator &src,
            const vector_view<tagged_terminal<Tag, const vector<T>&>, Slice> &term,
            const backend::command_queue &queue, const std::string &prm_name,
            detail::kernel_generator_state_ptr state)
    {
        term.slice.local_preamble(src, prm_name + "_slice", queue, state);
    }
};

template <typename Expr, class Slice>
struct partial_vector_expr< vector_view<Expr, Slice> > {
    static void get(backend::source_generator &src,
            const vector_view<Expr, Slice>&,
            const backend::command_queue&, const std::string &prm_name,
            detail::kernel_generator_state_ptr)
    {
        src << prm_name << "_val";
    }
};

template <typename T, class Slice>
struct partial_vector_expr< vector_view<const vector<T>&, Slice> > {
    static void get(backend::source_generator &src,
            const vector_view<const vector<T>&, Slice> &term,
            const backend::command_queue &queue, const std::string &prm_name,
            detail::kernel_generator_state_ptr state)
    {
        src << prm_name << "_expr_1[";
        term.slice.index(src, prm_name + "_slice", queue, state);
        src << "]";
    }
};

template <typename T, size_t Tag, class Slice>
struct partial_vector_expr< vector_view<tagged_terminal<Tag, const vector<T>&>, Slice> > {
    static void get(backend::source_generator &src,
            const vector_view<tagged_terminal<Tag, const vector<T>&>, Slice> &term,
            const backend::command_queue &queue, const std::string &prm_name,
            detail::kernel_generator_state_ptr state)
    {
        src << "prm_tag_" << Tag << "_1[";
        term.slice.index(src, prm_name + "_slice", queue, state);
        src << "]";
    }
};

template <typename Expr, class Slice>
struct kernel_arg_setter< vector_view<Expr, Slice> > {
    static void set(const vector_view<Expr, Slice> &term,
            backend::kernel &kernel, unsigned part, size_t index_offset,
            detail::kernel_generator_state_ptr state)
    {
        assert(part == 0);

        detail::set_expression_argument setarg(kernel, part, index_offset, state);
        detail::extract_terminals()( boost::proto::as_child(term.expr),  setarg);
        term.slice.setArgs(kernel, part, index_offset, state);
    }
};

template <typename T, class Slice>
struct expression_properties< vector_view<T, Slice> > {
    static void get(const vector_view<T, Slice> &term,
            std::vector<backend::command_queue> &queue_list,
            std::vector<size_t> &partition,
            size_t &size
            )
    {
        detail::get_expression_properties prop;
        detail::extract_terminals()(boost::proto::as_child(term.expr), prop);

        queue_list = prop.queue;
        partition  = prop.part;
        size       = term.slice.size();

        assert(partition.size() == 2);
        partition.back() = size;
    }
};

} // namespace traits

/// Generalized slice selector.
/**
 * This is very similar to std::gslice.
 *
 * Index to base vector is obtained as start + sum(i_k * stride[k]), where i_k
 * is coordinate along each dimension of gslice.
 */
template <size_t NDIM>
struct gslice {
    static_assert(NDIM > 0, "Incorrect dimension for gslice");

    size_t    start;
    size_t    length[NDIM];
    ptrdiff_t stride[NDIM]; // Signed type allows reverse slicing.

#ifndef BOOST_NO_INITIALIZER_LISTS
    template <typename T1, typename T2>
    gslice(size_t start,
           const std::initializer_list<T1> &p_length,
           const std::initializer_list<T2> &p_stride
          ) : start(start)
    {
        assert(p_length.size() == NDIM);
        assert(p_stride.size() == NDIM);

        std::copy(p_length.begin(), p_length.end(), length);
        std::copy(p_stride.begin(), p_stride.end(), stride);
    }
#endif

    template <typename T1, typename T2>
    gslice(size_t start,
           const std::array<T1, NDIM> &p_length,
           const std::array<T2, NDIM> &p_stride
          ) : start(start)
    {
        std::copy(p_length.begin(), p_length.end(), length);
        std::copy(p_stride.begin(), p_stride.end(), stride);
    }

    template <typename T1, typename T2>
    gslice(size_t start,
           const T1 *p_length,
           const T2 *p_stride
          ) : start(start)
    {
        std::copy(p_length, p_length + NDIM, length);
        std::copy(p_stride, p_stride + NDIM, stride);
    }

    size_t size() const {
        return std::accumulate(length, length + NDIM,
            static_cast<size_t>(1), std::multiplies<size_t>());
    }

    void preamble(backend::source_generator &src,
            const std::string &prm_name,
            const backend::command_queue&, detail::kernel_generator_state_ptr) const
    {
        src.begin_function<size_t>(prm_name + "_func");
        src.begin_function_parameters();
        src.parameter<size_t>("start");

        for(size_t k = 0; k < NDIM; ++k) {
            std::ostringstream plen, pstr;
            plen << "length" << k;
            pstr << "stride" << k;
            src.parameter<size_t>(plen.str());
            src.parameter<ptrdiff_t>(pstr.str());
        }
        src.parameter<size_t>("idx");
        src.end_function_parameters();

        if (NDIM == 1) {
            src.new_line() << "return start + idx * stride0;";
        } else {
            src.new_line()
                << "size_t ptr = start + (idx % length" << NDIM - 1
                << ") * stride" << NDIM - 1 << ";";

            for(size_t k = NDIM - 1; k-- > 0;) {
                src.new_line() << "idx /= length" << k + 1 << ";";
                src.new_line() << "ptr += (idx % length" << k <<  ") * stride"
                    << k <<  ";";
            }
            src.new_line() << "return ptr;";
        }
        src.end_function();
    }

    void parameter_declaration(backend::source_generator &src,
            const std::string &prm_name,
            const backend::command_queue&, detail::kernel_generator_state_ptr) const
    {
        src.parameter<size_t>(prm_name + "_start");

        for(size_t k = 0; k < NDIM; ++k) {
            std::ostringstream plen, pstr;
            plen << prm_name << "_length" << k;
            pstr << prm_name << "_stride" << k;
            src.parameter<size_t>(plen.str());
            src.parameter<ptrdiff_t>(pstr.str());
        }
    }

    void local_preamble(backend::source_generator&,
            const std::string&/*prm_name*/,
            const backend::command_queue&, detail::kernel_generator_state_ptr) const
    {
    }

    void index(backend::source_generator &src,
            const std::string &prm_name,
            const backend::command_queue&, detail::kernel_generator_state_ptr) const
    {
        src << prm_name << "_func" << "("
            << prm_name << "_start";

        for(size_t k = 0; k < NDIM; ++k)
            src << ", " << prm_name << "_length" << k
                << ", " << prm_name << "_stride" << k;

        src << ", idx)";
    }

    void setArgs(backend::kernel &kernel, unsigned/*part*/, size_t/*index_offset*/,
            detail::kernel_generator_state_ptr) const
    {
        kernel.push_arg(start);
        for(size_t k = 0; k < NDIM; ++k) {
            kernel.push_arg(length[k]);
            kernel.push_arg(stride[k]);
        }
    }

    /// Returns sliced vector.
    template <class Expr>
    vector_view<
        typename boost::proto::result_of::as_child<const Expr, vector_domain>::type,
        gslice
        >
    operator()(const Expr &expr) const {
#if (VEXCL_CHECK_SIZES > 0)
        {
            size_t max_idx = start;
            for(size_t i = 0; i < NDIM; ++i)
                max_idx += (length[i] - 1) * stride[i];

            detail::get_expression_properties expr_prop;
            detail::extract_terminals()(boost::proto::as_child(expr), expr_prop);

            precondition(max_idx < expr_prop.size,
                    "Slice will result in array overrun");
        }
#endif
        return vector_view<
                    typename boost::proto::result_of::as_child<const Expr, vector_domain>::type,
                    gslice
                >(boost::proto::as_child<vector_domain>(expr), *this);
    }
};


/// An index range for use with slicer class.
struct range {
    ptrdiff_t start;
    ptrdiff_t stride;
    ptrdiff_t stop;

    /// Unbounded range (all elements along the current dimension).
    range () : start(0), stride(0), stop(0) {}

    /// Range with a single element.
    range(ptrdiff_t i)
        : start(i), stride(1), stop(i + 1) {}

    /// Elements from open interval with given stride.
    range(ptrdiff_t start, ptrdiff_t stride, ptrdiff_t stop)
        : start(start), stride(stride), stop(stop) {}

    /// Every element from open interval.
    range(ptrdiff_t start, ptrdiff_t stop)
        : start(start), stride(1), stop(stop) {}

    bool empty() const {
        return !(start || stride || stop);
    }
};

/// Placeholder for an unbounded range.
const range _;

template <size_t NDIM>
struct extent_gen {
    std::array<size_t, NDIM> dim;

    extent_gen() {}

    extent_gen<NDIM + 1> operator[](size_t new_dim) const {
        extent_gen<NDIM + 1> new_extent;
        std::copy(dim.begin(), dim.end(), new_extent.dim.begin());
        new_extent.dim.back() = new_dim;
        return new_extent;
    }

    size_t size() const {
        return std::accumulate(dim.begin(), dim.end(),
                static_cast<size_t>(1), std::multiplies<size_t>());
    }
};

/// Helper object for specifying slicer dimensions.
const extent_gen<0> extents;

template <size_t NR, class Dimensions = boost::fusion::vector<> >
struct index_gen {
    std::array<range, NR> ranges;

    index_gen() {}

    typedef
        typename boost::fusion::result_of::as_vector<
            typename boost::fusion::result_of::push_back<
                Dimensions,
                boost::mpl::size_t<NR>
            >::type
        >::type next_dim;

    index_gen<NR+1, next_dim> operator[](const range &r) const {
        return append_range<next_dim>(r);
    }

    index_gen<NR+1, Dimensions> operator[](size_t i) const {
        return append_range<Dimensions>(i);
    }

    private:
        template <class Dim>
        index_gen<NR+1, Dim> append_range(const range &r) const {
            index_gen<NR+1, Dim> idx;
            std::copy(ranges.begin(), ranges.end(), idx.ranges.begin());
            idx.ranges.back() = r;
            return idx;
        }
};

/// Helper object for specifying slicer shape.
const index_gen<0> indices;

/// Slicing operator.
/**
 * Provides information about shape of multidimensional vector expressions,
 * allows to slice the expressions.
 */
template <size_t NR>
struct slicer {
    std::array<size_t, NR> dim;
    std::array<size_t, NR> stride;

    /// Creates slicer with the given dimensions.
    template <typename T>
    slicer(const std::array<T, NR> &target_dimensions) {
        init(target_dimensions.data());
    }

    /// Creates slicer with the given dimensions.
    template <typename T>
    slicer(const T *target_dimensions) {
        init(target_dimensions);
    }

    /// Creates slicer with the given dimensions.
    slicer(const extent_gen<NR> &ext) {
        init(ext.dim.data());
    }

    template <class Dimensions>
    gslice<NR> operator()(const index_gen<NR, Dimensions> &idx) const {
        size_t start = 0;
        std::array<size_t, NR> len;
        std::array<size_t, NR> str;

        for(size_t i = 0; i < NR; ++i) {
            range r = idx.ranges[i].empty() ? range(0, dim[i]) : idx.ranges[i];

            start += r.start * stride[i];
            len[i] = (std::abs(r.stop - r.start) + std::abs(r.stride) - 1) / std::abs(r.stride);
            str[i] = r.stride * stride[i];
        }

        return gslice<NR>(start, len, str);
    }

    template <size_t C>
    struct slice : public gslice<NR> {
        const slicer &parent;

        slice(const slicer &parent, const range &r)
            : gslice<NR>(r.start * parent.stride[0], parent.dim, parent.stride),
              parent(parent)
        {
            static_assert(C == 0, "Wrong slice constructor!");

            this->length[0] = (std::abs(r.stop - r.start) + std::abs(r.stride) - 1) / std::abs(r.stride);
            this->stride[0] *= r.stride;
        }

        slice(const slice<C-1> &parent, const range &r)
            : gslice<NR>(parent),
              parent(parent.parent)
        {
            static_assert(C > 0, "Wrong slice constructor!");

            this->start += r.start * this->stride[C];
            this->length[C] = (std::abs(r.stop - r.start) + std::abs(r.stride) - 1) / std::abs(r.stride);
            this->stride[C] *= r.stride;
        }

#ifdef _MSC_VER
#  pragma warning(push)
#  pragma warning(disable: 4307)
#endif
        slice<C+1> operator[](const range &r) const {
            return slice<C+1>(*this, r.empty() ? range(0, parent.dim[C + 1]) : r);
        }
#ifdef _MSC_VER
#  pragma warning(pop)
#endif
    };

    slice<0> operator[](const range &r) const {
        return slice<0>(*this, r.empty() ? range(0, dim[0]) : r);
    }

    private:
        template <typename T>
        void init(const T *target_dim) {
            std::copy(target_dim, target_dim + NR, dim.begin());

            stride.back() = 1;
            for(size_t i = NR - 1; i-- > 0;)
                stride[i] = stride[i + 1] * dim[i + 1];
        }
};

// Expression-based permutation operator.
template <class Expr>
struct expr_permutation {
    const Expr expr;

    expr_permutation(const Expr &expr) : expr(expr) {}

    size_t size() const {
        detail::get_expression_properties prop;
        detail::extract_terminals()(expr, prop);
        return prop.size;
    }

    void preamble(backend::source_generator &src,
            const std::string &prm_name,
            const backend::command_queue &dev, detail::kernel_generator_state_ptr state) const
    {
        detail::output_terminal_preamble ctx(src, dev, prm_name, state);
        boost::proto::eval(boost::proto::as_child(expr), ctx);
    }

    void parameter_declaration(backend::source_generator &src,
            const std::string &prm_name,
            const backend::command_queue &dev, detail::kernel_generator_state_ptr state) const
    {
        detail::declare_expression_parameter ctx(src, dev, prm_name, state);
        detail::extract_terminals()(boost::proto::as_child(expr), ctx);
    }

    void local_preamble(backend::source_generator &src,
            const std::string &prm_name,
            const backend::command_queue &dev, detail::kernel_generator_state_ptr state) const
    {
        detail::output_local_preamble init_ctx(src, dev, prm_name, state);
        boost::proto::eval(boost::proto::as_child(expr), init_ctx);
    }

    void index(backend::source_generator &src,
            const std::string &prm_name,
            const backend::command_queue &dev, detail::kernel_generator_state_ptr state) const
    {
        detail::vector_expr_context ctx(src, dev, prm_name, state);
        boost::proto::eval(boost::proto::as_child(expr), ctx);
    }

    void setArgs(backend::kernel &kernel, unsigned part, size_t index_offset,
            detail::kernel_generator_state_ptr state) const
    {
        detail::extract_terminals()( boost::proto::as_child(expr),
                detail::set_expression_argument(kernel, part, index_offset, state));
    }

    template <class Base>
    vector_view<
        typename boost::proto::result_of::as_child<const Base, vector_domain>::type,
        expr_permutation
        >
    operator()(const Base &base) const {
#if (VEXCL_CHECK_SIZES > 1)
        {
            detail::get_expression_properties base_prop;
            detail::extract_terminals()(boost::proto::as_child(base), base_prop);

            precondition(!base_prop.queue.empty(),
                    "Can not permute stateless expression");

            detail::get_expression_properties expr_prop;
            detail::extract_terminals()(boost::proto::as_child(expr), expr_prop);

            vex::Reductor<size_t, vex::MAX> max(base_prop.queue);

            precondition(
                    (expr_prop.size == 0 && base_prop.size == 0) || max(expr) < base_prop.size,
                    "Permutation will result in array overrun");
        }
#endif
        return vector_view<
                    typename boost::proto::result_of::as_child<const Base, vector_domain>::type,
                    expr_permutation
                >(boost::proto::as_child<vector_domain>(base), *this);
    }
};

/// Returns permutation functor which is based on an integral expression.
template <class Expr>
auto permutation(const Expr &expr) ->
    typename std::enable_if<
        std::is_integral<typename detail::return_type<Expr>::type>::value,
        expr_permutation<
            typename boost::proto::result_of::as_child<const Expr, vector_domain>::type
        >
    >::type
{
    return expr_permutation<
        typename boost::proto::result_of::as_child<const Expr, vector_domain>::type
        >(boost::proto::as_child(expr));
}

//---------------------------------------------------------------------------
// Slice reduction
//---------------------------------------------------------------------------
struct reduced_vector_view_terminal {};

typedef vector_expression<
    typename boost::proto::terminal< reduced_vector_view_terminal >::type
    > reduced_vector_view_terminal_expression;

template <class Expr, size_t NDIM, size_t NR, class RDC>
struct reduced_vector_view : public reduced_vector_view_terminal_expression
{
    typedef typename detail::return_type<Expr>::type value_type;
    const Expr   expr;
    gslice<NDIM> slice;
    std::array<size_t, NR> reduce_dims;

    reduced_vector_view(
            const Expr &expr, const gslice<NDIM> slice, std::array<size_t, NR> dims
            ) : expr(expr), slice(slice)
    {
        std::copy(dims.begin(), dims.end(), reduce_dims.begin());
        std::sort(reduce_dims.begin(), reduce_dims.end());
    }
};

namespace traits {

template <>
struct is_vector_expr_terminal< reduced_vector_view_terminal >
    : std::true_type
{ };

template <>
struct proto_terminal_is_value< reduced_vector_view_terminal >
    : std::true_type
{ };

template <class Expr, size_t NDIM, size_t NR, class RDC>
struct terminal_preamble< reduced_vector_view<Expr, NDIM, NR, RDC> > {
    static void get(backend::source_generator &src,
            const reduced_vector_view<Expr, NDIM, NR, RDC> &term,
            const backend::command_queue &queue, const std::string &prm_name,
            detail::kernel_generator_state_ptr state)
    {
        detail::output_terminal_preamble termpream(src, queue, prm_name, state);
        boost::proto::eval(boost::proto::as_child(term.expr), termpream);

        typedef typename detail::return_type<Expr>::type T;
        typedef typename RDC::template impl<T>::result_type T_out;
        typedef typename RDC::template impl<T>::device_in fun_in;

        boost::proto::eval(boost::proto::as_child( fun_in() (T_out(), T()) ), termpream);
    }
};

template <typename Expr, size_t NDIM, size_t NR, class RDC>
struct local_terminal_init< reduced_vector_view<Expr, NDIM, NR, RDC> > {
    static void get(backend::source_generator &src,
            const reduced_vector_view<Expr, NDIM, NR, RDC> &term,
            const backend::command_queue &queue, const std::string &prm_name,
            detail::kernel_generator_state_ptr state)
    {
        typedef typename detail::return_type<Expr>::type T;
        typedef typename RDC::template impl<T>::result_type T_out;
        typedef typename RDC::template impl<T>::device_in fun;

        src.new_line() << type_name<T_out>() << " " << prm_name << "_sum = " << RDC::template impl<T>::initial() << ";";
        src.open("{");

        src.new_line()
            << "size_t pos = idx;";

        src.new_line()
            << "size_t ptr" << NDIM - NR << " = "
            << prm_name << "_start";
        if (NDIM > NR) src
            << " + (pos % " << prm_name << "_length"
            << NDIM - NR - 1 << ") * " << prm_name << "_stride"
            << NDIM - NR - 1;
        src << ";";

        for(size_t k = NDIM - NR; k-- > 1;) {
            src.new_line()
                << "pos /= " << prm_name << "_length" << k << ";";
            src.new_line()
                << "ptr" << NDIM - NR << " += (pos % " << prm_name
                << "_length" << k - 1 << ") * " << prm_name << "_stride" << k - 1 << ";";
        }

        for(size_t k = NDIM - NR; k < NDIM; ++k) {
            src.new_line()
                << "for(size_t i" << k << " = 0, ptr" << k + 1 << " = ptr"
                << k << "; i" << k << " < " << prm_name << "_length" << k
                << "; ++i" << k << ", ptr" << k + 1 << " += " << prm_name
                << "_stride" << k << ")";
            src.open("{");
        }

        src.new_line() << "size_t idx = ptr" << NDIM << ";";

        detail::output_local_preamble init_ctx(src, queue, prm_name, state);
        boost::proto::eval(boost::proto::as_child(term.expr), init_ctx);

        src.new_line()
            << prm_name << "_sum = " << fun::name() << "("
            << prm_name << "_sum, ";

        detail::vector_expr_context expr_ctx(src, queue, prm_name, state);
        boost::proto::eval(boost::proto::as_child(term.expr), expr_ctx);

        src << ");";

        for(size_t k = NDIM - NR; k < NDIM; ++k) src.close("}");
        src.close("}");
    }
};

template <typename Expr, size_t NDIM, size_t NR, class RDC>
struct kernel_param_declaration< reduced_vector_view<Expr, NDIM, NR, RDC> > {
    static void get(backend::source_generator &src,
            const reduced_vector_view<Expr, NDIM, NR, RDC> &term,
            const backend::command_queue &queue, const std::string &prm_name,
            detail::kernel_generator_state_ptr state)
    {
        detail::declare_expression_parameter declare(src, queue, prm_name, state);
        detail::extract_terminals()(boost::proto::as_child(term.expr), declare);
        term.slice.parameter_declaration(src, prm_name, queue, state);
    }
};

template <typename Expr, size_t NDIM, size_t NR, class RDC>
struct partial_vector_expr< reduced_vector_view<Expr, NDIM, NR, RDC> > {
    static void get(backend::source_generator &src,
            const reduced_vector_view<Expr, NDIM, NR, RDC>&,
            const backend::command_queue&, const std::string &prm_name,
            detail::kernel_generator_state_ptr)
    {
        src << prm_name << "_sum";
    }
};

template <typename Expr, size_t NDIM, size_t NR, class RDC>
struct kernel_arg_setter< reduced_vector_view<Expr, NDIM, NR, RDC> > {
    static void set(const reduced_vector_view<Expr, NDIM, NR, RDC> &term,
            backend::kernel &kernel, unsigned part, size_t index_offset,
            detail::kernel_generator_state_ptr state)
    {
        detail::set_expression_argument setarg(kernel, part, index_offset, state);
        detail::extract_terminals()( boost::proto::as_child(term.expr), setarg);

        kernel.push_arg(term.slice.start);

        for(size_t k = 0; k < NDIM; ++k) {
            if (!std::binary_search(term.reduce_dims.begin(), term.reduce_dims.end(), k)) {
                kernel.push_arg(term.slice.length[k]);
                kernel.push_arg(term.slice.stride[k]);
            }
        }

        for(size_t k = 0; k < NR; ++k) {
            kernel.push_arg(term.slice.length[term.reduce_dims[k]]);
            kernel.push_arg(term.slice.stride[term.reduce_dims[k]]);
        }
    }
};

template <typename Expr, size_t NDIM, size_t NR, class RDC>
struct expression_properties< reduced_vector_view<Expr, NDIM, NR, RDC> > {
    static void get(const reduced_vector_view<Expr, NDIM, NR, RDC> &term,
            std::vector<backend::command_queue> &queue_list,
            std::vector<size_t> &partition,
            size_t &size
            )
    {
        detail::get_expression_properties prop;
        detail::extract_terminals()(boost::proto::as_child(term.expr), prop);

        queue_list = prop.queue;
        partition  = std::vector<size_t>(2, 0);
        size       = 1;

        for(size_t k = 0; k < NDIM; ++k)
            if (!std::binary_search(term.reduce_dims.begin(), term.reduce_dims.end(), k))
                size *= term.slice.length[k];

        partition.back() = size;
    }
};

} // namespace traits

#ifdef DOXYGEN
/// Reduce multidimensional vector expression along specified dimensions.
template <class Reducer, class SlicedExpr, class ReduceDims>
auto reduce(const SlicedExpr &expr, const ReduceDims &reduce_dims);

/// Reduce multidimensional vector expression along specified dimensions.
template <class Reducer, class ExprShape, class Expr, class ReduceDims>
auto reduce(const ExprShape &shape, const Expr &expr, const ReduceDims &reduce_dims);
#else
/// Reduce vector_view along specified dimensions.
template <class RDC, typename Expr, size_t NDIM, size_t NR>
reduced_vector_view<Expr, NDIM, NR, RDC> reduce(
        const vector_view<Expr, gslice<NDIM> > &view,
        const std::array<size_t, NR> &reduce_dims
        )
{
    return reduced_vector_view<Expr, NDIM, NR, RDC>(view.expr, view.slice, reduce_dims);
}

/// Reduce vector_view along specified dimension.
template <class RDC, typename Expr, size_t NDIM>
reduced_vector_view<Expr, NDIM, 1, RDC> reduce(
        const vector_view<Expr, gslice<NDIM> > &view,
        size_t reduce_dim
        )
{
    std::array<size_t, 1> dim = {{reduce_dim}};
    return reduced_vector_view<Expr, NDIM, 1, RDC>(view.expr, view.slice, dim);
}

/// Reduce sliced expression along specified dimensions.
template <class RDC, typename Expr, size_t NDIM, size_t NR>
reduced_vector_view<
    typename boost::proto::result_of::as_child<const Expr, vector_domain>::type,
    NDIM, NR, RDC
> reduce(
        const gslice<NDIM> &slice,
        const Expr &expr,
        const std::array<size_t, NR> &reduce_dims
        )
{
    return reduced_vector_view<
        typename boost::proto::result_of::as_child<const Expr, vector_domain>::type,
        NDIM, NR, RDC>(boost::proto::as_child(expr), slice, reduce_dims);
}

/// Reduce sliced expression along specified dimension.
template <class RDC, typename Expr, size_t NDIM>
reduced_vector_view<
    typename boost::proto::result_of::as_child<const Expr, vector_domain>::type,
    NDIM, 1, RDC
> reduce(
        const gslice<NDIM> &slice,
        const Expr &expr,
        size_t reduce_dim
        )
{
    std::array<size_t, 1> dim = {{reduce_dim}};
    return reduced_vector_view<
        typename boost::proto::result_of::as_child<const Expr, vector_domain>::type,
        NDIM, 1, RDC>(boost::proto::as_child(expr), slice, dim);
}

/// Reduce sliced expression along specified dimensions.
template <class RDC, typename Expr, size_t NDIM, size_t NR>
reduced_vector_view<
    typename boost::proto::result_of::as_child<const Expr, vector_domain>::type,
    NDIM, NR, RDC
> reduce(
        const extent_gen<NDIM> &ext,
        const Expr &expr,
        const std::array<size_t, NR> &reduce_dims
        )
{
    return reduce<RDC>(slicer<NDIM>(ext)[_], expr, reduce_dims);
}

/// Reduce sliced expression along specified dimensions.
template <class RDC, typename Expr, size_t NDIM, size_t NR>
reduced_vector_view<
    typename boost::proto::result_of::as_child<const Expr, vector_domain>::type,
    NDIM, NR, RDC
> reduce(
        const extent_gen<NDIM> &ext,
        const Expr &expr,
        const extent_gen<NR> &reduce_dims
        )
{
    return reduce<RDC>(slicer<NDIM>(ext)[_], expr, reduce_dims.dim);
}

/// Reduce sliced expression along specified dimension.
template <class RDC, typename Expr, size_t NDIM>
reduced_vector_view<
    typename boost::proto::result_of::as_child<const Expr, vector_domain>::type,
    NDIM, 1, RDC
> reduce(
        const extent_gen<NDIM> &ext,
        const Expr &expr,
        size_t reduce_dim
        )
{
    return reduce<RDC>(slicer<NDIM>(ext)[_], expr, reduce_dim);
}
#endif

//---------------------------------------------------------------------------
// Matrix reshaper
//---------------------------------------------------------------------------

namespace detail {

template <size_t Nout, size_t Nin>
struct reshape_helper {
    std::array<size_t, Nout> dst_dims;
    std::array<size_t, Nin>  src_dims;
    std::array<size_t, Nout> dst_stride;
    std::array<size_t, Nin>  src_stride;

    decltype( vex::element_index() ) idx;

    reshape_helper(
            std::array<size_t, Nout> dst_dims,
            std::array<size_t, Nin>  src_dims
            )
        : dst_dims(dst_dims), src_dims(src_dims),
        idx(vex::element_index(0,
                    std::accumulate(dst_dims.begin(), dst_dims.end(),
                        static_cast<size_t>(1), std::multiplies<size_t>())
                    )
           )
    {
        src_stride[Nin - 1] = 1;
        for(size_t k = Nin - 1; k-- > 0;)
            src_stride[k] = src_stride[k + 1] * dst_dims[src_dims[k + 1]];

        dst_stride[Nout - 1] = 1;
        for(size_t k = Nout - 1; k-- > 0;)
            dst_stride[k] = dst_stride[k + 1] * dst_dims[k + 1];
    }

    template <class T>
    static T fake_instance();

    template <size_t I, class Enable = void>
    struct return_type {};

    template <size_t I>
    struct return_type<I, typename std::enable_if<I == 0>::type> {
        typedef
            decltype(
                    size_t(1) * ( ( vex::element_index() / size_t(1) ) % size_t(1) )
                    )
            type;
    };

    template <size_t I>
    struct return_type<I, typename std::enable_if<(I > 0)>::type> {
        typedef
            decltype(
                    size_t(1) * ( ( vex::element_index() / size_t(1) ) % size_t(1) )
                    + fake_instance<typename return_type<I - 1>::type>()
                    )
            type;
    };

    template <size_t I>
    auto tail_sum() const ->
        typename std::enable_if<
            (I == 0),
            typename return_type<I>::type
        >::type
    {
        size_t j = src_dims[I];
        return src_stride[I] * ( ( idx / dst_stride[j] ) % dst_dims[j] );
    }

    template <size_t I>
        typename std::enable_if<
            (I > 0),
            typename return_type<I>::type
        >::type
    tail_sum() const
    {
        size_t j = src_dims[I];
        return src_stride[I] * ( ( idx / dst_stride[j] ) % dst_dims[j] )
            + tail_sum<I - 1>();
    }

    typename return_type<Nin - 1>::type
    operator()() const
    {
        return tail_sum<Nin - 1>();
    }

};

} //namespace detail


#ifdef DOXYGEN
/// Reshapes the expression.
/**
 * Makes a multidimensional expression shaped as `dst_dims` from an input
 * expression shaped as `dst_dims[src_dims]`. `scr_dims` are specified as
 * indices into `dst_dims`.
 */
template <class Expr, class DstDims, class SrcDims>
auto reshape(const Expr &expr, const DstDims &dst_dims, const SrcDims &src_dims);
#else
/// Reshapes the expression.
template <class Expr, size_t Nout, size_t Nin>
auto reshape(
        const Expr &expr,
        const std::array<size_t, Nout> &dst_dims,
        const std::array<size_t, Nin>  &src_dims
        ) ->
    decltype(vex::permutation(detail::reshape_helper<Nout, Nin>(dst_dims, src_dims)())(expr))
{
    return vex::permutation(detail::reshape_helper<Nout, Nin>(dst_dims, src_dims)())(expr);
}

/// Reshapes the expression.
template <class Expr, size_t Nout, size_t Nin>
auto reshape(
        const Expr &expr,
        const extent_gen<Nout> &dst_dims,
        const extent_gen<Nin>  &src_dims
        ) ->
    decltype(vex::permutation(detail::reshape_helper<Nout, Nin>(dst_dims.dim, src_dims.dim)())(expr))
{
    return vex::permutation(detail::reshape_helper<Nout, Nin>(dst_dims.dim, src_dims.dim)())(expr);
}
#endif

} // namespace vex

#endif
