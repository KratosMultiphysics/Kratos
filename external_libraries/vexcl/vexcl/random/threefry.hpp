#ifndef VEXCL_RANDOM_THREEFRY_HPP
#define VEXCL_RANDOM_THREEFRY_HPP

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
 * \file   vexcl/random/threefry.hpp
 * \author Pascal Germroth <pascal@ensieve.org>
 * \brief  Threefry RNG.

Threefry, based on the Threefish cipher, is a non cryptographic algorithm
for pseudorandom number generation from the Random123 suite,
see <http://www.deshawresearch.com/resources_random123.html>

The original code came with the following copyright notice:

\verbatim
Copyright 2010-2011, D. E. Shaw Research.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

* Redistributions of source code must retain the above copyright
  notice, this list of conditions, and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions, and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

* Neither the name of D. E. Shaw Research nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
\endverbatim
*/

namespace vex {
namespace random {

namespace detail {
    template <size_t bits, size_t w>
    struct rotation_table;

    template<>
    struct rotation_table<32, 2> {
        static std::array<unsigned, 8> get() {
            static const std::array<unsigned, 8> R = {{
                13, 15, 26,  6, 17, 29, 16, 24
            }};
            return R;
        }
    };

    template <>
    struct rotation_table<32, 4> {
        static std::array<unsigned, 16> get() {
            static const std::array<unsigned, 16> R = {{
                10, 26, 11, 21, 13, 27, 23,  5,
                 6, 20, 17, 11, 25, 10, 18, 20
            }};
            return R;
        }
    };

    template <>
    struct rotation_table<64, 2> {
        static std::array<unsigned, 8> get() {
            static const std::array<unsigned, 8> R = {{
                16, 42, 12, 31, 16, 32, 24, 21
            }};
            return R;
        }
    };

    template <>
    struct rotation_table<64, 4> {
        static std::array<unsigned, 16> get() {
            static const std::array<unsigned, 16> R = {{
                14, 16, 52, 57, 23, 40,  5, 37,
                25, 33, 46, 12, 58, 22, 32, 32
            }};
            return R;
        }
    };
}

/// Threefry random number generator.
/**
 * Threefry, based on the Threefish cipher, is a non cryptographic algorithm
 * for pseudorandom number generation from the Random123 suite.
 * \see http://www.deshawresearch.com/resources_random123.html
 * \sa vex::Random
 * \sa vex::RandomNormal
 */
struct threefry {
    static std::string name() {
        return "threefry";
    }

    // Generates function threefry(ctr, key).
    // ctr will be modified, containing the random output.
    // key will be preserved.
    template <class T, size_t N, size_t R = 20>
    struct function {
        static const size_t K = N;

        static_assert(
                N == 2 || N == 4,
                "Only supports vectors with 2 or 4 components."
                );

        static_assert(
                std::is_same<T, cl_uint>::value ||
                std::is_same<T, cl_ulong>::value,
                "Only supports 32 or 64 bit integers."
                );

        static std::string name() {
            std::ostringstream s;
            s << "threefry_" << type_name<T>() << "_" << N << "_" << R;
            return s.str();
        }

        static void define(backend::source_generator &src) {
            const size_t bits = sizeof(T) * 8;
            auto rot = detail::rotation_table<bits, N>::get();

            src.begin_function<void>( name() );
            src.begin_function_parameters();
            src.template parameter< regstr_ptr<T> >("ctr");
            src.template parameter< regstr_ptr<T> >("key");
            src.end_function_parameters();

#if defined(VEXCL_BACKEND_CUDA) || defined(VEXCL_BACKEND_JIT)
            src.new_line() << "#define rotate(x, b) "
                "(((x) << (b)) | ((x) >> (sizeof(x)*8 - (b))))";
#endif

            src.new_line() << "const " << type_name<T>() << " p = "
                << (bits == 32 ? "0x1BD11BDA" : "0x1BD11BDAA9FC1A22");
            for(size_t i = 0; i < N; ++i)
                src << " ^ key[" << i << "]";
            src << ";";

        // Insert initial key before round 0
        for(size_t i = 0; i < N; ++i)
            src.new_line() << "ctr[" << i << "] += key[" << i << "];";

        for(size_t round = 0; round < R; ++round) {
            if(N == 2) {
                src.new_line()
                    << "ctr[0] += ctr[1]; "
                    << "ctr[1] = rotate(ctr[1], " << rot[round % 8] << "u); "
                    << "ctr[1] ^= ctr[0];";
            } else {
                const size_t r = 2 * (round % 8),
                    r0 = r + (round % 2),
                    r1 = r + ((round + 1) % 2);
                src.new_line()
                    << "ctr[0] += ctr[1]; "
                    << "ctr[1] = rotate(ctr[1], " << rot[r0] << "u); "
                    << "ctr[1] ^= ctr[0];";
                src.new_line()
                    << "ctr[2] += ctr[3]; "
                    << "ctr[3] = rotate(ctr[3], " << rot[r1] << "u); "
                    << "ctr[3] ^= ctr[2];";
            }

            // inject key
            if((round + 1) % 4 == 0) {
                const size_t j = round / 4 + 1;
                for(size_t i = 0; i < N; ++i) {
                    const size_t ii = ((j + i) % (N + 1));
                    src.new_line() << "ctr[" << i << "] += ";
                    if(ii == N) src << "p; ";
                    else src << "key[" << ii << "]; ";
                }
                src << "ctr[" << (N - 1) << "] += " << j << ";";
            }
        }

#ifdef VEXCL_BACKEND_CUDA
            src.new_line() << "#undef rotate";
#endif

            src.end_function();
        }
    };
};


} // namespace random
} // namespace vex

#endif
