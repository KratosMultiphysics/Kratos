#ifndef VEXCL_SPARSE_SPMV_OPS_HPP
#define VEXCL_SPARSE_SPMV_OPS_HPP

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
 * \file   vexcl/sparse/spmv_ops.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Default low-level operations for sparse matrix-vector product.
 */

#include <vexcl/operations.hpp>

namespace vex {
namespace sparse {

template <class mat_type, class vec_type, class enable = void>
struct spmv_ops_impl {
    static void decl_accum_var(backend::source_generator &src, const std::string &name)
    {
        typedef decltype(std::declval<mat_type>() * std::declval<vec_type>()) res_type;
        src.new_line() << type_name<res_type>() << " " << name << " = " << res_type() << ";";
    }

    static void append(backend::source_generator &src,
            const std::string &sum, const std::string &val)
    {
        src.new_line() << sum << " += " << val << ";";
    }

    static void append_product(backend::source_generator &src,
            const std::string &sum, const std::string &mat_val, const std::string &vec_val)
    {
        src.new_line() << sum << " += " << mat_val << " * " << vec_val << ";";
    }
};

} // namespace sparse
} // namespace vex


#endif
