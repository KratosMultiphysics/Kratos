#ifndef VEXCL_EVAL_HPP
#define VEXCL_EVAL_HPP

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
 * \file   vexcl/eval.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Evaluate vector expression.
 */

#include <vexcl/operations.hpp>

namespace vex {

template <class Expr>
void eval(const Expr &expr,
        const std::vector<backend::command_queue> &queue,
        const std::vector<size_t> &part
        )
{
    using namespace vex::detail;

#if (VEXCL_CHECK_SIZES > 0)
    {
        get_expression_properties prop;
        extract_terminals()(boost::proto::as_child(expr), prop);

        precondition(
                prop.queue.empty() || prop.queue.size() == queue.size(),
                "Incompatible queue lists"
                );

        precondition(
                prop.size == 0 || prop.size == part.back(),
                "Incompatible expression sizes"
                );
    }
#endif
    static kernel_cache cache;

    for(unsigned d = 0; d < queue.size(); d++) {
        auto kernel = cache.find(queue[d]);

        backend::select_context(queue[d]);

        if (kernel == cache.end()) {
            backend::source_generator source(queue[d]);

            output_terminal_preamble termpream(source, queue[d], "prm", empty_state());

            boost::proto::eval(boost::proto::as_child(expr), termpream);

            source.begin_kernel("vexcl_eval_kernel");
            source.begin_kernel_parameters();
            source.parameter<size_t>("n");

            declare_expression_parameter declare(source, queue[d], "prm", empty_state());
            extract_terminals()(boost::proto::as_child(expr), declare);

            source.end_kernel_parameters();
            source.grid_stride_loop().open("{");

            output_local_preamble loc_init(source, queue[d], "prm", empty_state());
            boost::proto::eval(boost::proto::as_child(expr), loc_init);

            source.new_line();
            vector_expr_context expr_ctx(source, queue[d], "prm", empty_state());
            boost::proto::eval(boost::proto::as_child(expr), expr_ctx);
            source << ";";
            source.close("}");
            source.end_kernel();

            kernel = cache.insert(queue[d], backend::kernel(
                        queue[d], source.str(), "vexcl_eval_kernel"));
        }

        if (size_t psize = part[d + 1] - part[d]) {
            auto &K = kernel->second;
            K.push_arg(psize);
            set_expression_argument setarg(K, d, part[d], empty_state());
            extract_terminals()( boost::proto::as_child(expr), setarg);
            K(queue[d]);
        }
    }
}

template <class Expr>
void eval(const Expr &expr) {
    using namespace vex::detail;

    get_expression_properties prop;
    extract_terminals()(boost::proto::as_child(expr), prop);

    precondition(!prop.queue.empty() && !prop.part.empty(),
            "Can not determine expression size and queue list"
            );

    eval(expr, prop.queue, prop.part);
}

}

#endif
