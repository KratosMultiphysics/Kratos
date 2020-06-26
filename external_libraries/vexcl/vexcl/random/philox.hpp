#ifndef VEXCL_RANDOM_PHILOX_HPP
#define VEXCL_RANDOM_PHILOX_HPP

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
 * \file   vexcl/random/philox.hpp
 * \author Pascal Germroth <pascal@ensieve.org>
 * \brief  Philox RNG.

Philox is an integer multiplication based, non cryptographic algorithm
for pseudorandom number generation from the Random123 suite,
see <http://www.deshawresearch.com/resources_random123.html>.

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

/// Random generators.
namespace random {


/// Philox random number generator.
/**
 * Philox is an integer multiplication based, non cryptographic algorithm
 * for pseudorandom number generation from the Random123 suite.
 * \see http://www.deshawresearch.com/resources_random123.html
 * \sa vex::Random
 * \sa vex::RandomNormal
 */
struct philox {

    static std::string name() {
        return "philox";
    }

    // Generates function philox(ctr, key);
    // modifies both inputs, uses the components of ctr for randomness.
    template <class T, size_t N, size_t R = 10>
    struct function {
        static const size_t K = N / 2;

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
            s << "philox_" << type_name<T>() << "_" << N << "_" << R;
            return s.str();
        }

        static void define(backend::source_generator &src) {
            std::string M[2], W[2];
            if(std::is_same<T, cl_uint>::value) { // 32
                W[0] = "(" + type_name<T>() + ")0x9E3779B9";
                if(N == 2) {
                    M[0] = "(" + type_name<T>() + ")0xD256D193";
                } else {
                    W[1] = "(" + type_name<T>() + ")0xBB67AE85";
                    M[0] = "(" + type_name<T>() + ")0xD2511F53";
                    M[1] = "(" + type_name<T>() + ")0xCD9E8D57";
                }
            } else { // 64
                W[0] = "(" + type_name<T>() + ")0x9E3779B97F4A7C15"; // golden ratio
                if(N == 2) {
                    M[0] = "(" + type_name<T>() + ")0xD2B74407B1CE6E93";
                } else {
                    W[1] = "(" + type_name<T>() + ")0xBB67AE8584CAA73B"; // sqrt(3)-1
                    M[0] = "(" + type_name<T>() + ")0xD2E7470EE14C6C93";
                    M[1] = "(" + type_name<T>() + ")0xCA5A826395121157";
                }
            }

            src.begin_function<void>( name() );
            src.begin_function_parameters();
            src.template parameter< regstr_ptr<T> >("ctr");
            src.template parameter< regstr_ptr<T> >("key");
            src.end_function_parameters();

            src.new_line() << type_name<T>() << " m[" << N << "];";

#if defined(VEXCL_BACKEND_CUDA)
            src.new_line() << "#define mul_hi __umulhi";
#elif defined(VEXCL_BACKEND_JIT)
            src.new_line() << "#define mul_hi mulhi";
#endif

            for(size_t round = 0; round < R; ++round) {
                if(round > 0) { // bump key
                    for(size_t i = 0; i < K; ++i)
                        src.new_line() << "key[" << i <<"] += " << W[i] << ";";
                }
                // next round
                if(N == 2) {
                    src.new_line() << "m[0] = mul_hi(" << M[0] << ", ctr[0]);";
                    src.new_line() << "m[1] = " << M[0] << " * ctr[0];";
                    src.new_line() << "ctr[0] = m[0] ^ key[0] ^ ctr[1];";
                    src.new_line() << "ctr[1] = m[1];";
                } else {
                    src.new_line() << "m[0] = mul_hi(" << M[0] << ", ctr[0]);";
                    src.new_line() << "m[1] = " << M[0] << " * ctr[0];";
                    src.new_line() << "m[2] = mul_hi(" << M[1] << ", ctr[2]);";
                    src.new_line() << "m[3] = " << M[1] << " * ctr[2];";
                    src.new_line() << "ctr[0] = m[2] ^ ctr[1] ^ key[0];";
                    src.new_line() << "ctr[1] = m[3];";
                    src.new_line() << "ctr[2] = m[0] ^ ctr[3] ^ key[1];";
                    src.new_line() << "ctr[3] = m[1];";
                }
            }

#ifdef VEXCL_BACKEND_CUDA
            src.new_line() << "#undef mul_hi";
#endif

            src.end_function();
        }
    };
};


} // namespace random
} // namespace vex

#endif
