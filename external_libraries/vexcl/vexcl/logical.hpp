#ifndef VEXCL_LOGICAL_HPP
#define VEXCL_LOGICAL_HPP

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
 * \file   vexcl/logical.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Implementation of any_of and all_of primitives.
 */

#include <vector>
#include <vexcl/operations.hpp>

namespace vex {

/// Functor that checks if any element in a vector expression is true.
/**
 * Example:
\code
vex::any_of any_of(ctx);
bool have_zeros = any_of(x == 0);
\endcode
 */
class any_of {
    public:
        any_of(const std::vector<backend::command_queue> &queue
#ifndef VEXCL_NO_STATIC_CONTEXT_CONSTRUCTORS
                = current_context().queue()
#endif
              ) : queue(queue)
        {
            result.reserve(queue.size());
            for(auto q = queue.begin(); q != queue.end(); ++q)
                result.push_back( backend::device_vector<char>(*q, 1) );
        }

        template <class Expr>
        bool operator()(const Expr &expr) const {
            using namespace detail;

            get_expression_properties prop;
            extract_terminals()(expr, prop);

            for(unsigned d = 0; d < queue.size(); ++d) {
                if (size_t psize = prop.part_size(d)) {
                    backend::kernel& k = make_kernel(queue[d], expr);
                    k.push_arg(psize);
                    extract_terminals()(expr,
                            set_expression_argument(
                                k, d, prop.part_start(d), empty_state()
                                )
                            );
                    k.push_arg(result[d]);
                    k(queue[d]);
                }
            }

            for(unsigned d = 0; d < queue.size(); ++d) {
                if (prop.part_size(d)) {
                    char r;
                    result[d].read(queue[d], 0, 1, &r, true);
                    if (r) return true;
                }
            }
            return false;
        }

    private:
        std::vector<backend::command_queue> queue;
        std::vector<backend::device_vector<char>> result;

        template <class Expr>
        static backend::kernel& make_kernel(
                const backend::command_queue &q,
                const Expr &expr
                )
        {
            using namespace detail;
            static kernel_cache cache;

            auto kernel = cache.find(q);

            if (kernel == cache.end()) {
                backend::source_generator src(q);

                output_terminal_preamble gpre(src, q, "prm", empty_state());
                boost::proto::eval(boost::proto::as_child(expr), gpre);

                src.begin_kernel("vexcl_any_of_kernel");
                src.begin_kernel_parameters();
                src.parameter<size_t>("n");

                extract_terminals()(expr,
                        declare_expression_parameter(
                            src, q, "prm", empty_state()
                            )
                        );

                src.template parameter< global_ptr<char> >("result");

                src.end_kernel_parameters();
                src.new_line() << "for(ulong idx = 0; idx < n; ++idx)";
                src.open("{");

                output_local_preamble lpre(src, q, "prm", empty_state());
                boost::proto::eval(boost::proto::as_child(expr), lpre);

                src.new_line() << "if (";

                vector_expr_context expr_ctx(src, q, "prm", empty_state());
                boost::proto::eval(boost::proto::as_child(expr), expr_ctx);

                src << ")";
                src.open("{");
                src.new_line() << "result[0] = 1;";
                src.new_line() << "return;";
                src.close("}");

                src.close("}");
                src.new_line() << "result[0] = 0;";
                src.end_kernel();

                kernel = cache.insert(q, backend::kernel(
                            q, src.str(), "vexcl_any_of_kernel"
                            ));
                kernel->second.config(1, 1);
            }

            return kernel->second;
        }
};

/// Functor that checks if all elements in a vector expression are true.
/**
 * Example:
\code
vex::all_of all_of(ctx);
bool odds_only = all_of(x % 2 == 1);
\endcode
 */
class all_of : public any_of {
    public:
        all_of(const std::vector<backend::command_queue> &queue
#ifndef VEXCL_NO_STATIC_CONTEXT_CONSTRUCTORS
                = current_context().queue()
#endif
              ) : any_of(queue) {}

        template <class Expr>
        bool operator()(const Expr &expr) const {
            return !any_of::operator()(!expr);
        }
};

} // namespace vexcl


#endif
