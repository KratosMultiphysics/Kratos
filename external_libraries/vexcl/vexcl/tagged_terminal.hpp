#ifndef VEXCL_TAGGED_TERMINAL_HPP
#define VEXCL_TAGGED_TERMINAL_HPP

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
 * \file   vexcl/tagged_terminal.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Tagged terminal wrapper.
 */

#include <vector>
#include <string>
#include <sstream>
#include <map>
#include <set>

#include <vexcl/operations.hpp>

namespace vex {

struct tagged_terminal_terminal {};

typedef vector_expression<
    typename boost::proto::terminal< tagged_terminal_terminal >::type
    > tagged_terminal_expression;

template <size_t Tag, class Term>
struct tagged_terminal : tagged_terminal_expression
{
    typedef typename detail::return_type<Term>::type value_type;

    Term term;

    tagged_terminal(const Term &term) : term(term) {}

    // Expression assignments.
#define VEXCL_ASSIGNMENT(cop, op)                                              \
  template <class Expr>                                                        \
  typename std::enable_if<                                                     \
      boost::proto::matches<                                                   \
          typename boost::proto::result_of::as_expr<Expr>::type,               \
          vector_expr_grammar>::value,                                         \
      const tagged_terminal &>::type operator cop(const Expr & expr) const {   \
    detail::assign_expression<op>(*this, expr);                                \
    return *this;                                                              \
  }

    VEXCL_ASSIGNMENTS(VEXCL_ASSIGNMENT)

#undef VEXCL_ASSIGNMENT
};

namespace traits {

template <>
struct is_vector_expr_terminal< tagged_terminal_terminal > : std::true_type {};

template <>
struct proto_terminal_is_value< tagged_terminal_terminal > : std::true_type {};

template <size_t Tag, class Term>
struct terminal_preamble< tagged_terminal<Tag, Term> > {
    static void get(backend::source_generator &src,
            const tagged_terminal<Tag, Term> &term,
            const backend::command_queue &queue, const std::string&/*prm_name*/,
            detail::kernel_generator_state_ptr state)
    {
        auto s = state->find("tag_pream");

        if (s == state->end()) {
            s = state->insert(std::make_pair(
                        std::string("tag_pream"),
                        boost::any(std::set<size_t>())
                        )).first;
        }

        auto &pos = boost::any_cast< std::set<size_t>& >(s->second);
        auto p = pos.find(Tag);

        if (p == pos.end()) {
            pos.insert(Tag);

            std::ostringstream prm_name;
            prm_name << "prm_tag_" << Tag;

            detail::output_terminal_preamble termpream(src, queue, prm_name.str(), state);
            boost::proto::eval(boost::proto::as_child(term.term), termpream);
        }
    }
};

template <size_t Tag, class Term>
struct kernel_param_declaration< tagged_terminal<Tag, Term> > {
    static void get(backend::source_generator &src,
            const tagged_terminal<Tag, Term> &term,
            const backend::command_queue &queue, const std::string&/*prm_name*/,
            detail::kernel_generator_state_ptr state)
    {
        auto s = state->find("tag_param");

        if (s == state->end()) {
            s = state->insert(std::make_pair(
                        std::string("tag_param"),
                        boost::any(std::set<size_t>())
                        )).first;
        }

        auto &pos = boost::any_cast< std::set<size_t>& >(s->second);
        auto p = pos.find(Tag);

        if (p == pos.end()) {
            pos.insert(Tag);

            std::ostringstream prm_name;
            prm_name << "prm_tag_" << Tag;

            detail::declare_expression_parameter declare(src, queue, prm_name.str(), state);
            detail::extract_terminals()(boost::proto::as_child(term.term),  declare);
        }
    }
};

template <size_t Tag, class Term>
struct local_terminal_init< tagged_terminal<Tag, Term> > {
    static void get(backend::source_generator &src,
            const tagged_terminal<Tag, Term> &term,
            const backend::command_queue &queue,
            const std::string&/*prm_name*/,
            detail::kernel_generator_state_ptr state)
    {
        auto s = state->find("tag_locinit");

        if (s == state->end()) {
            s = state->insert(std::make_pair(
                        std::string("tag_locinit"),
                        boost::any(std::set<size_t>())
                        )).first;
        }

        auto &pos = boost::any_cast< std::set<size_t>& >(s->second);
        auto p = pos.find(Tag);

        if (p == pos.end()) {
            pos.insert(Tag);

            std::ostringstream prm_name;
            prm_name << "prm_tag_" << Tag;

            detail::output_local_preamble init_ctx(src, queue, prm_name.str(), state);
            boost::proto::eval(boost::proto::as_child(term.term), init_ctx);
        }
    }
};

template <size_t Tag, class Term>
struct partial_vector_expr< tagged_terminal<Tag, Term> > {
    static void get(backend::source_generator &src,
            const tagged_terminal<Tag, Term> &term,
            const backend::command_queue &queue,
            const std::string&/*prm_name*/,
            detail::kernel_generator_state_ptr state)
    {
        std::ostringstream prm_name;
        prm_name << "prm_tag_" << Tag;

        detail::vector_expr_context expr_ctx(src, queue, prm_name.str(), state);
        boost::proto::eval(boost::proto::as_child(term.term), expr_ctx);
    }
};

template <size_t Tag, class Term>
struct kernel_arg_setter< tagged_terminal<Tag, Term> > {
    static void set(const tagged_terminal<Tag, Term> &term,
            backend::kernel &kernel, unsigned part, size_t index_offset,
            detail::kernel_generator_state_ptr state)
    {
        auto s = state->find("tag_args");

        if (s == state->end()) {
            s = state->insert(std::make_pair(
                        std::string("tag_args"),
                        boost::any(std::set<size_t>())
                        )).first;
        }

        auto &pos = boost::any_cast< std::set<size_t>& >(s->second);
        auto p = pos.find(Tag);

        if (p == pos.end()) {
            pos.insert(Tag);

            detail::set_expression_argument setarg(kernel, part, index_offset, state);
            detail::extract_terminals()( boost::proto::as_child(term.term),  setarg);
        }
    }
};

template <size_t Tag, class Term>
struct expression_properties< tagged_terminal<Tag, Term> > {
    static void get(const tagged_terminal<Tag, Term> &term,
            std::vector<backend::command_queue> &queue_list,
            std::vector<size_t> &partition,
            size_t &size
            )
    {
        detail::get_expression_properties prop;
        detail::extract_terminals()(boost::proto::as_child(term.term), prop);

        queue_list = prop.queue;
        partition  = prop.part;
        size       = prop.size;
    }
};

} // namespace traits

/// Tags terminal with a unique (in a single expression) tag.
/**
 * By tagging terminals user guarantees that the terminals with same tags
 * actually refer to the same data. VexCL is able to use this information in
 * order to reduce number of kernel parameters and unnecessary global memory
 * I/O operations.
 */
template <size_t Tag, class Expr>
auto tag(const Expr& expr)
        -> const tagged_terminal<
                        Tag,
                        decltype(boost::proto::as_child<vector_domain>(expr))
                >
{
        static_assert(
                boost::proto::matches<
                        typename boost::proto::result_of::as_expr<Expr>::type,
                        boost::proto::terminal<boost::proto::_>
                >::value,
                "Tagging non-terminals is not allowed"
                );

    return tagged_terminal<
                Tag,
                decltype(boost::proto::as_child<vector_domain>(expr))
                >(boost::proto::as_child<vector_domain>(expr));
}

} //namespace vex
#endif
