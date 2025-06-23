#ifndef VEXCL_BACKEND_JIT_FILTER_HPP
#define VEXCL_BACKEND_JIT_FILTER_HPP

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
 * \file   vexcl/backend/jit/filter.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Device filters for the JIT backend.
 */

#include <string>
#include <vector>
#include <functional>
#include <cstdlib>

namespace vex {
namespace Filter {

struct DummyFilter {
    bool v;
    DummyFilter(bool v) : v(v) {}
    bool operator()(const backend::device &d) const { return v; }
};

const DummyFilter GPU(false);
const DummyFilter CPU(true);
const DummyFilter Accelerator(false);
const DummyFilter DoublePrecision(true);

struct Name {
    explicit Name(std::string name) : devname(std::move(name)) {}

    bool operator()(const backend::device &d) const {
        return d.name().find(devname) != std::string::npos;
    }

    private:
        std::string devname;
};

// The JIT backend does not support the exclusive filter functionality,
// but provides the filter for compatibility.
template <class Filter>
Filter Exclusive(Filter&& filter) {
    return std::forward<Filter>(filter);
}

inline std::vector< std::function<bool(const backend::device&)> >
backend_env_filters()
{
    std::vector< std::function<bool(const backend::device&)> > filter;

    if (const char *name = getenv("OCL_DEVICE"))
        filter.push_back(Name(name));

    return filter;
}

} // namespace Filter
} // namespace vex


#endif
