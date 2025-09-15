#ifndef VEXCL_ENQUEUE_HPP
#define VEXCL_ENQUEUE_HPP

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
 * \file   vexcl/enqueue.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Assignment expression wrapper allowing to explicitly set the command queue to use.
 */

#include <vexcl/operations.hpp>

namespace vex {

template <class LHS>
struct enqueue_vector_impl {
    LHS &lhs;

    std::vector<backend::command_queue> q;
    std::vector<size_t> p;

    enqueue_vector_impl(LHS &lhs, const std::vector<backend::command_queue> &q)
        : lhs(lhs), q(q)
    {
        detail::get_expression_properties prop;
        detail::extract_terminals()(boost::proto::as_child(lhs), prop);

        p = prop.part;
    }

#ifdef DOXYGEN
#define VEXCL_ASSIGNMENT(op, op_type)                                          \
    /** Expression assignment operator. */                                     \
    template <class RHS> const LHS& operator op(const RHS &rhs);
#else
#define VEXCL_ASSIGNMENT(op, op_type)                                          \
    template <class RHS>                                                       \
    auto operator op(const RHS &rhs) ->                                        \
        typename std::enable_if<                                               \
            boost::proto::matches<                                             \
                typename boost::proto::result_of::as_expr<RHS>::type,          \
                vector_expr_grammar>::value,                                   \
            const LHS&>::type                                                  \
    {                                                                          \
        detail::assign_expression<op_type>(lhs, rhs, q, p);                    \
        return lhs;                                                            \
    }
#endif

    VEXCL_ASSIGNMENTS(VEXCL_ASSIGNMENT)

#undef VEXCL_ASSIGNMENT
};


template <class LHS>
struct enqueue_multiex_impl {
    LHS &lhs;

    std::vector<backend::command_queue> q;
    std::vector<size_t> p;

    enqueue_multiex_impl(LHS &lhs, const std::vector<backend::command_queue> &q)
        : lhs(lhs), q(q)
    {
        detail::get_expression_properties prop;
        detail::extract_terminals()(detail::subexpression<0>::get(lhs), prop);

        p = prop.part;
    }

#ifdef DOXYGEN
#define VEXCL_ASSIGNMENT(op, op_type)                                          \
    /** Expression assignment operator. */                                     \
    template <class RHS> const LHS& operator op(const RHS &rhs);
#else
#define VEXCL_ASSIGNMENT(op, op_type)                                          \
    template <class RHS>                                                       \
    auto operator op(const RHS &rhs) ->                                        \
        typename std::enable_if<                                               \
            boost::proto::matches<                                             \
                typename boost::proto::result_of::as_expr<RHS>::type,          \
                multivector_expr_grammar>::value || is_tuple<LHS>::value,      \
            const LHS&>::type                                                  \
    {                                                                          \
        detail::assign_multiexpression<op_type>(lhs, rhs, q, p);                    \
        return lhs;                                                            \
    }
#endif

    VEXCL_ASSIGNMENTS(VEXCL_ASSIGNMENT)

#undef VEXCL_ASSIGNMENT
};



/// Assignment operation proxy.
/**
 * Allows to explicitly specify the command queue to use for the assignment
 */
template <class LHS>
auto enqueue(const std::vector<backend::command_queue> &q, LHS &lhs) ->
    typename std::enable_if<
        boost::proto::matches<
            typename boost::proto::result_of::as_expr<LHS>::type,
            vector_expr_grammar>::value,
        enqueue_vector_impl<LHS> >::type
{
    return enqueue_vector_impl<LHS>(lhs, q);
}

template <class LHS>
auto enqueue(const std::vector<backend::command_queue> &q, LHS &lhs) ->
    typename std::enable_if<
        boost::proto::matches<
            typename boost::proto::result_of::as_expr<LHS>::type,
            multivector_expr_grammar>::value || is_tuple<LHS>::value,
        enqueue_multiex_impl<LHS> >::type
{
    return enqueue_multiex_impl<LHS>(lhs, q);
}

} // namespace vex

#endif
