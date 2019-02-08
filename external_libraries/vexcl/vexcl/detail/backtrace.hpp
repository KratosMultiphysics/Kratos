#ifndef VEXCL_DETAIL_BACKTRACE_HPP
#define VEXCL_DETAIL_BACKTRACE_HPP

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
 * \file   vexcl/detail/backtrace.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Print backtrace if possible.
 */

#include <iostream>

#ifdef __linux__
#  include <execinfo.h>
#  include <stdlib.h>
#endif

namespace vex {
namespace detail {

#ifdef __linux__
inline void print_backtrace() {
    const size_t nbuf = 100;
    void *buffer[nbuf];

    int nptrs = backtrace(buffer, nbuf);

    if (char **strings = backtrace_symbols(buffer, nptrs)) {
        for(int i = 0; i < nptrs; ++i)
            std::cerr << strings[i] << "\n";

        std::cerr << std::endl;

        free(strings);
    }
}
#else
inline void print_backtrace() {}
#endif

} // namespace detail
} // namespace vex

#endif
