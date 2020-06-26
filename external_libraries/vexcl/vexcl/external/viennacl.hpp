#ifndef VEXCL_VIENNACL_HPP
#define VEXCL_VIENNACL_HPP

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
 * \file   external/viennacl.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Enables use of ViennaCL (http://viennacl.sourceforge.net) iterative solvers.
 */

namespace viennacl {
    namespace linalg {
        /// Inner product of two vex::vector instances.
        template< typename T1, typename T2 >
        decltype(T1() * T2())
        inner_prod(const vex::vector<T1> &v1, const vex::vector<T2> &v2) {
            vex::Reductor<decltype(T1() * T2()), vex::SUM> sum(v1.queue_list());
            return sum(v1 * v2);
        }

        /// Matrix-vector product.
        template <typename T, typename C, typename I>
        auto prod(const vex::SpMat<T, C, I> &A, const vex::vector<T> &x)
            -> decltype(A * x)
        {
            return A * x;
        }

        /// L2-norm.
        template <typename T>
        T norm_2(const vex::vector<T> &x) {
            return sqrt(inner_prod(x, x));
        }
    }

    namespace traits {
        /// Zeros-out the given vector.
        template <class T>
        void clear(vex::vector<T> &vec) {
            vec = 0;
        }
    }
}

#endif
