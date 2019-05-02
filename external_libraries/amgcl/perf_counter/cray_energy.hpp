#ifndef AMGCL_PERF_COUNTER_CRAY_ENERGY_HPP
#define AMGCL_PERF_COUNTER_CRAY_ENERGY_HPP

/*
The MIT License

Copyright (c) 2012-2019 Denis Demidov <dennis.demidov@gmail.com>
Copyright (c) 2016 Mohammad Siahatgar <siahatgar@luis.uni-hannover.de>

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
 * \file   amgcl/perf_counter/cray_energy.hpp
 * \author Mohammad Siahatgar <siahatgar@luis.uni-hannover.de>
 * \brief  Cray energy counter.
 */

#include <fstream>
#include <amgcl/util.hpp>

namespace amgcl {
namespace perf_counter {

class cray_energy {
    public:
        typedef long long value_type;

        cray_energy() : f(attribute_path()) {
            precondition(f, std::string("Failed to open ") + attribute_path());
        }

        static const char* units() {
            return "J";
        }

        value_type current() {
            f.clear();
            f.seekg(0, std::ios::beg);

            value_type v;
            f >> v;
            return v;
        }
    private:
        static const char* attribute_path() {
            return "/sys/cray/pm_counters/energy";
        }

        std::ifstream f;
};

} // namespace perf_counter
} // namespace amgcl

#endif
