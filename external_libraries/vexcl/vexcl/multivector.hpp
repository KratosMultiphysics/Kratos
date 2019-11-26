#ifndef VEXCL_MULTIVECTOR_HPP
#define VEXCL_MULTIVECTOR_HPP

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
 * \file   vexcl/multivector.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  OpenCL device multi-vector.
 */

#include <array>
#include <vector>
#include <map>
#include <iostream>
#include <sstream>
#include <string>
#include <type_traits>
#include <boost/proto/proto.hpp>
#include <boost/iterator/iterator_facade.hpp>

#include <vexcl/util.hpp>
#include <vexcl/operations.hpp>
#include <vexcl/vector.hpp>

/// Vector expression template library for OpenCL.
namespace vex {

struct multivector_terminal {};

template <typename T, size_t N> class multivector;

namespace traits {

// Extract component directly from terminal rather than from value(terminal):
template <>
struct proto_terminal_is_value< multivector_terminal >
    : std::true_type
{ };

template <>
struct is_multivector_expr_terminal< multivector_terminal >
    : std::true_type
{ };

// Hold multivector terminals by reference:
template <class T>
struct hold_terminal_by_reference< T,
        typename std::enable_if<
            boost::proto::matches<
                typename boost::proto::result_of::as_expr< T >::type,
                boost::proto::terminal< multivector_terminal >
            >::value
        >::type
    >
    : std::true_type
{ };

template <typename T, size_t N>
struct number_of_components< multivector<T, N> > : boost::mpl::size_t<N> {};

template <size_t I, typename T, size_t N>
struct component< I, multivector<T, N> > {
    typedef const vector<T>& type;
};

} // namespace traits

template <size_t I, typename T, size_t N>
const vector<T>& get(const multivector<T, N> &mv) {
    static_assert(I < N, "Component number out of bounds");

    return mv(I);
}

template <size_t I, typename T, size_t N>
vector<T>& get(multivector<T, N> &mv) {
    static_assert(I < N, "Component number out of bounds");

    return mv(I);
}

namespace detail {

template <class OP, class LHS, class RHS>
void assign_multiexpression(LHS &lhs, const RHS &rhs,
        const std::vector<backend::command_queue> &queue,
        const std::vector<size_t> &part
        );
}


typedef multivector_expression<
    typename boost::proto::terminal< multivector_terminal >::type
    > multivector_terminal_expression;

/// Container for several equally sized instances of vex::vector<T>.
template <typename T, size_t N>
class multivector : public multivector_terminal_expression {
    public:
        typedef vex::vector<T>  subtype;
        typedef std::array<T,N> value_type;
        typedef T               sub_value_type;

        const static size_t NDIM = N;

        // Proxy class.
        class element {
            public:
                operator const value_type () const {
                    value_type val;
                    for(unsigned i = 0; i < N; i++) val[i] = vec(i)[index];
                    return val;
                }

                const value_type operator=(value_type val) {
                    for(unsigned i = 0; i < N; i++) vec(i)[index] = val[i];
                    return val;
                }
            private:
                element(multivector &vec, size_t index)
                    : vec(vec), index(index) {}

                multivector &vec;
                const size_t      index;

                friend class multivector;
        };

        // Proxy class.
        class const_element {
            public:
                operator const value_type () const {
                    value_type val;
                    for(unsigned i = 0; i < N; i++) val[i] = vec(i)[index];
                    return val;
                }
            private:
                const_element(const multivector &vec, size_t index)
                    : vec(vec), index(index) {}

                const multivector &vec;
                const size_t      index;

                friend class multivector;
        };

        template <class V, class E>
        class iterator_type
            : public boost::iterator_facade<
                        iterator_type<V, E>,
                        sub_value_type,
                        std::random_access_iterator_tag,
                        E
                     >
        {
            public:
                typedef boost::iterator_facade<
                            iterator_type<V, E>,
                            sub_value_type,
                            std::random_access_iterator_tag,
                            E
                         > super_type;
                typedef typename super_type::reference       reference;
                typedef typename super_type::difference_type difference_type;

                reference dereference() const {
                    return E(vec, pos);
                }

                bool equal(const iterator_type &it) const {
                    return pos == it.pos;
                }

                void increment() {
                    ++pos;
                }

                void decrement() {
                    --pos;
                }

                void advance(difference_type n) {
                    pos += n;
                }

                difference_type distance_to(const iterator_type &it) const {
                    return static_cast<difference_type>(it.pos - pos);
                }
            private:
                iterator_type(V &vec, size_t pos) : vec(vec), pos(pos) {}

                V      &vec;
                size_t pos;

                friend class multivector;
        };

        typedef iterator_type<multivector, element> iterator;
        typedef iterator_type<const multivector, const_element> const_iterator;

        multivector() {};

        /// Constructor.
        /**
         * The host vector data is divided equally between the created
         * multivector components. Each component gets continuous chunk of the
         * source vector.
         */
        multivector(const std::vector<backend::command_queue> &queue,
                const std::vector<T> &host,
                backend::mem_flags flags = backend::MEM_READ_WRITE)
        {
            static_assert(N > 0, "What's the point?");

            size_t size = host.size() / N;
            assert(N * size == host.size());

            for(size_t i = 0; i < N; ++i)
                vec[i].resize(queue, size, host.data() + i * size, flags);
        }

        /// Constructor.
        /**
         * If host pointer is not NULL, it is copied to the underlying vector
         * components of the multivector.  Each component gets continuous chunk
         * of the source vector.
         */
        multivector(const std::vector<backend::command_queue> &queue, size_t size,
                const T *host = 0, backend::mem_flags flags = backend::MEM_READ_WRITE)
        {
            static_assert(N > 0, "What's the point?");

            for(size_t i = 0; i < N; ++i)
                vec[i].resize(queue, size, host ? host + i * size : 0, flags);
        }

        /// Constructor.
        /**
         * Uses the most recently created VexCL context.
         */
        multivector(size_t size) {
            static_assert(N > 0, "What's the point?");

            for(size_t i = 0; i < N; ++i) vec[i].resize(size);
        }

#ifdef VEXCL_NO_COPY_CONSTRUCTORS
    private:
#endif
        /// Copy constructor.
        multivector(const multivector &mv) {
#ifdef VEXCL_SHOW_COPIES
            std::cout << "Copying vex::multivector<" << type_name<T>()
                      << ", " << N << "> of size " << size() << std::endl;
#endif
            for(size_t i = 0; i < N; ++i) vec[i].resize(mv(i));
        }
#ifdef VEXCL_NO_COPY_CONSTRUCTORS
    public:
#endif

        /// Move constructor.
        multivector(multivector &&mv) noexcept {
            for(size_t i = 0; i < N; ++i) vec[i].swap(mv.vec[i]);
        }

        /// Resizes the multivector.
        /**
         * This is equivalent to reconstructing the vector with the given
         * parameters.  Any data contained in the resized vector will be lost
         * as a result.
         */
        void resize(const std::vector<backend::command_queue> &queue, size_t size) {
            for(unsigned i = 0; i < N; i++) vec[i].resize(queue, size);
        }

        /// Resizes the multivector.
        /**
         * Uses the most recently created VexCL context.
         * This is equivalent to reconstructing the vector with the given
         * parameters.  Any data contained in the resized vector will be lost
         * as a result.
         */
        void resize(size_t size) {
            for(unsigned i = 0; i < N; i++) vec[i].resize(size);
        }

        /// Fills the multivector with zeros.
        void clear() {
            *this = static_cast<T>(0);
        }

        /// Returns size of the multivector (equals size of individual components).
        size_t size() const {
            return vec[0].size();
        }

        /// Returns i-th multivector component.
        const vex::vector<T>& operator()(size_t i) const {
            return vec[i];
        }

        /// Returns i-th multivector component.
        vex::vector<T>& operator()(size_t i) {
            return vec[i];
        }

        /// Returns const iterator to the first element of the multivector.
        const_iterator begin() const {
            return const_iterator(*this, 0);
        }

        /// Returns const iterator to the first element of the multivector.
        iterator begin() {
            return iterator(*this, 0);
        }

        /// Returns const iterator referring to the past-the-end element in the multivector.
        const_iterator end() const {
            return const_iterator(*this, size());
        }

        /// Returns iterator referring to the past-the-end element in the multivector.
        iterator end() {
            return iterator(*this, size());
        }

        /// Returns i-th elements of all components packed in a std::array<T,N>.
        const_element operator[](size_t i) const {
            return const_element(*this, i);
        }

        /// Assigns values from std::array<T,N> to i-th elements of all components.
        element operator[](size_t i) {
            return element(*this, i);
        }

        /// Returns reference to the multivector's queue list.
        const std::vector<backend::command_queue>& queue_list() const {
            return vec[0].queue_list();
        }

        /// Assignment operator
        const multivector& operator=(const multivector &mv) {
            if (this != &mv)
                detail::assign_multiexpression<assign::SET>(
                        *this, mv, vec[0].queue_list(), vec[0].partition());
            return *this;
        }

#define VEXCL_ASSIGNMENT(op, op_type)                                          \
  /** Assignment operator */                                                   \
  template <class Expr>                                                        \
  auto operator op(const Expr &expr) ->                                        \
      typename std::enable_if<                                                 \
          boost::proto::matches<                                               \
              typename boost::proto::result_of::as_expr<Expr>::type,           \
              multivector_expr_grammar>::value ||                              \
              is_tuple<Expr>::value,                                           \
          const multivector &>::type                                           \
  {                                                                            \
    detail::assign_multiexpression<op_type>(*this, expr, vec[0].queue_list(),  \
                                       vec[0].partition());                    \
    return *this;                                                              \
  }

        VEXCL_ASSIGNMENTS(VEXCL_ASSIGNMENT)

#undef VEXCL_ASSIGNMENT

#ifndef DOXYGEN
        template <class Expr>
        typename std::enable_if<
            boost::proto::matches<
                typename boost::proto::result_of::as_expr<Expr>::type,
                additive_multivector_transform_grammar
            >::value,
            const multivector&
        >::type
        operator=(const Expr &expr) {
            detail::apply_additive_transform</*append=*/false>(*this,
                    detail::simplify_additive_transform()( expr ));
            return *this;
        }

        template <class Expr>
        typename std::enable_if<
            boost::proto::matches<
                typename boost::proto::result_of::as_expr<Expr>::type,
                additive_multivector_transform_grammar
            >::value,
            const multivector&
        >::type
        operator+=(const Expr &expr) {
            detail::apply_additive_transform</*append=*/true>(*this,
                    detail::simplify_additive_transform()( expr ));
            return *this;
        }

        template <class Expr>
        typename std::enable_if<
            boost::proto::matches<
                typename boost::proto::result_of::as_expr<Expr>::type,
                additive_multivector_transform_grammar
            >::value,
            const multivector&
        >::type
        operator-=(const Expr &expr) {
            detail::apply_additive_transform</*append=*/true>(*this,
                    detail::simplify_additive_transform()( -expr ));
            return *this;
        }

        template <class Expr>
        typename std::enable_if<
            boost::proto::matches<
                typename boost::proto::result_of::as_expr<Expr>::type,
                multivector_full_grammar
            >::value &&
            !boost::proto::matches<
                typename boost::proto::result_of::as_expr<Expr>::type,
                multivector_expr_grammar
            >::value &&
            !boost::proto::matches<
                typename boost::proto::result_of::as_expr<Expr>::type,
                additive_multivector_transform_grammar
            >::value,
            const multivector&
        >::type
        operator=(const Expr &expr) {
            *this  = detail::extract_multivector_expressions()( expr );
            *this += detail::extract_additive_multivector_transforms()( expr );

            return *this;
        }

        template <class Expr>
        typename std::enable_if<
            boost::proto::matches<
                typename boost::proto::result_of::as_expr<Expr>::type,
                multivector_full_grammar
            >::value &&
            !boost::proto::matches<
                typename boost::proto::result_of::as_expr<Expr>::type,
                multivector_expr_grammar
            >::value &&
            !boost::proto::matches<
                typename boost::proto::result_of::as_expr<Expr>::type,
                additive_multivector_transform_grammar
            >::value,
            const multivector&
        >::type
        operator+=(const Expr &expr) {
            *this += detail::extract_multivector_expressions()( expr );
            *this += detail::extract_additive_multivector_transforms()( expr );

            return *this;
        }

        template <class Expr>
        typename std::enable_if<
            boost::proto::matches<
                typename boost::proto::result_of::as_expr<Expr>::type,
                multivector_full_grammar
            >::value &&
            !boost::proto::matches<
                typename boost::proto::result_of::as_expr<Expr>::type,
                multivector_expr_grammar
            >::value &&
            !boost::proto::matches<
                typename boost::proto::result_of::as_expr<Expr>::type,
                additive_multivector_transform_grammar
            >::value,
            const multivector&
        >::type
        operator-=(const Expr &expr) {
            *this -= detail::extract_multivector_expressions()( expr );
            *this -= detail::extract_additive_multivector_transforms()( expr );

            return *this;
        }
#endif

    private:
        std::array<vector<T>,N> vec;
};

/// Copy multivector to host vector.
template <class T, size_t N>
void copy(const multivector<T,N> &mv, std::vector<T> &hv) {
    for(size_t i = 0; i < N; ++i)
        vex::copy(mv(i).begin(), mv(i).end(), hv.begin() + i * mv.size());
}

/// Copy host vector to multivector.
template <class T, size_t N>
void copy(const std::vector<T> &hv, multivector<T,N> &mv) {
    for(size_t i = 0; i < N; ++i)
        vex::copy(hv.begin() + i * mv.size(), hv.begin() + (i + 1) * mv.size(),
                mv(i).begin());
}

/// Download and print the vector elements.
template<class T, size_t N>
std::ostream &operator<<(std::ostream &o, const vex::multivector<T, N> &t) {
    boost::io::ios_all_saver stream_state(o);
    const size_t n = t.size();
    const size_t chunk = std::max<int>(1, (std::is_integral<T>::value ? 12 : 6) / N);

    std::vector<T> data(n * N);
    copy(t, data);

    o << "{" << std::setprecision(6);
    for(size_t i = 0 ; i < n; i++) {
        if (i % chunk == 0) o << "\n" << std::setw(6) << i << ":";

        std::cout << " (";
        for(size_t j = 0; j < N; ++j) {
            if (std::is_integral<T>::value)
                o << " " << std::setw(6) << data[j * n + i];
            else if (std::is_arithmetic<T>::value)
                o << std::scientific << std::setw(14) << data[j * n + i];
            else
                o << " " << data[j * n + i];
        }
        std::cout << ")";
    }
    return o << "\n}\n";
}

} // namespace vex

namespace boost { namespace fusion { namespace traits {

template <class T, size_t N>
struct is_sequence< vex::multivector<T, N> > : std::false_type
{};

} } }


#endif
