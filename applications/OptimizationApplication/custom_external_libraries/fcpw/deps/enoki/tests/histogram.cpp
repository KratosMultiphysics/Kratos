/*
    test/histogram.cpp -- Test which uses transform_<> to build a histogram
    of a set of normally distributed pseudorandom samples

    Enoki is a C++ template library that enables transparent vectorization
    of numerical kernels using SIMD instruction sets available on current
    processor architectures.

    Copyright (c) 2019 Wenzel Jakob <wenzel.jakob@epfl.ch>

    All rights reserved. Use of this source code is governed by a BSD-style
    license that can be found in the LICENSE file.
*/

#if defined(NDEBUG)
#  undef NDEBUG
#endif

#include <enoki/random.h>
#include <enoki/special.h>

using namespace enoki;

int main(int /* argc */, char * /* argv */[]) {
    using UInt32      = Packet<uint32_t>;
    using UInt32Mask  = mask_t<UInt32>;
    using RNG         = PCG32<UInt32>;
    using Float32     = RNG::Float32;
    using UInt64      = RNG::UInt64;

    /* Bin configuration */
    const float min_value = -4;
    const float max_value =  4;
    const uint32_t bin_count = 31;
    uint32_t bins[bin_count] { };

    for (size_t j = 0; j < 16 / UInt32::Size; ++j) {
        RNG rng(PCG32_DEFAULT_STATE, arange<UInt64>() + (j * UInt32::Size));

        for (size_t i = 0; i < 1024 * 1024; ++i) {
            /* Generate a uniform variate */
            Float32 x = rng.next_float32();

            /* Importance sample a normal distribution */
            Float32 y = float(M_SQRT2) * erfinv(2.f*x - 1.f);

            /* Compute bin index */
            UInt32 idx((y - min_value) * float(bin_count) / (max_value - min_value));

            /* Discard samples that are out of bounds */
            UInt32Mask mask = idx >= zero<UInt32>() && idx < bin_count;

            /* Increment the bin indices */
            scatter_add(bins, UInt32(1), idx, mask);
        }
    }

    uint32_t sum = 0;
    for (uint32_t i = 0; i < bin_count; ++i) {
        std::cout << "bin[" << i << "] = ";
        for (uint32_t j = 0; j < bins[i] / 50000; ++j)
            std::cout << "*";
        std::cout << " " << bins[i] << std::endl;
        sum += bins[i];
    }

#if defined(__aarch64__)
    assert(std::abs(int(16 * 1024 * 1024 - sum) - 743) <= 200);
#else
    assert(std::abs(int(16 * 1024 * 1024 - sum) - 743) <= 3);
#endif
    assert(bins[1] == 2558);
    assert(bins[2] == 6380);

    return 0;
}
