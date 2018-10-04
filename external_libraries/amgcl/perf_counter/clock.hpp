#ifndef AMGCL_PERF_COUNTER_CLOCK_HPP
#define AMGCL_PERF_COUNTER_CLOCK_HPP

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
 * \file   amgcl/perf_counter/clock.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Clock class.
 *
 * A minimal wrapper around either omp_get_wtime() or std::clock().
 */

#ifdef _OPENMP
#  include <omp.h>
#else
#  include <ctime>
#endif

namespace amgcl {

/// Performance counters for use with amgcl::profiler
namespace perf_counter {

/// Clock class.
/**
 * Designed to interchangeable (in context of amgcl::profiler) with either
 * std::chrono or boost::chrono clocks.
 *
 * Uses omp_get_wtime() when available, std::clock() otherwise.
 */
struct clock {
    typedef double value_type;

    static const char* units() {
        return "s";
    }

    /// Current time point.
    static double current() {
#ifdef _OPENMP
        return omp_get_wtime();
#else
        return std::clock() / static_cast<double>(CLOCKS_PER_SEC);
#endif
    }
};

} // namespace perf_counter
} // namespace amgcl


#endif
