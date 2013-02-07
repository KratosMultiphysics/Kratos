#ifndef AMGCL_TICTOC_HPP
#define AMGCL_TICTOC_HPP

/*
The MIT License

Copyright (c) 2012-2013 Denis Demidov <ddemidov@ksu.ru>

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
 * \file   tictoc.hpp
 * \author Denis Demidov <ddemidov@ksu.ru>
 * \brief  TIC and TOC macros for internal profiling.
 *
 * If AMGCL_PROFILING macro is defined at compilation, then TIC(name) and
 * TOC(name) macros correspond to prof.tic(name) and prof.toc(name).
 * amgcl::prof should be an instance of amgcl::profiler<> defined in a user
 * code similar to:
 * \code
 * namespace amgcl { profiler<> prof; }
 * \endcode
 * If AMGCL_PROFILING is undefined, then TIC and TOC are noop macros.
 */

#ifdef AMGCL_PROFILING
#  include <amgcl/profiler.hpp>
#  define TIC(name) prof.tic(name);
#  define TOC(name) prof.toc(name);
namespace amgcl {
extern profiler<> prof;
}
#else
#  define TIC(name)
#  define TOC(name)
#endif

#endif
