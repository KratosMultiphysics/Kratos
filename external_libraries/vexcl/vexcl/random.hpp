#ifndef VEXCL_RANDOM_HPP
#define VEXCL_RANDOM_HPP

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
 * \file   vexcl/random.hpp
 * \author Pascal Germroth <pascal@ensieve.org>
 * \brief  Random number generators.
 */

#include <vexcl/operations.hpp>
#include <boost/math/constants/constants.hpp>
#include <vexcl/random/philox.hpp>
#include <vexcl/random/threefry.hpp>


namespace vex {

/// Returns uniformly distributed random numbers.
/**
 * For integral types, generated values span the complete range.
 *
 * For floating point types, generated values are in \f$[0,1]\f$.
 *
 * Uses Random123 generators which provide 64(2x32), 128(4x32, 2x64)
 * and 256(4x64) random bits, this limits the supported output types,
 * which means `cl_double8` (512bit) is not supported, but `cl_uchar2` is.
 *
 * Supported generator families are ``random::philox`` (based on integer
 * multiplication, default) and ``random::threefry`` (based on the Threefish
 * encryption function). Both satisfy rigorous statistical testing (passing
 * BigCrush in TestU01), vectorize and parallelize well (each generator can
 * produce at least \f$2^{64}\f$ independent streams), have long periods (the
 * period of each stream is at least \f$2^{128}\f$), require little or no
 * memory or state, and have excellent performance (a few clock cycles per byte
 * of random output).
 */
template <class T, class Generator = random::philox>
struct Random : UserFunction<Random<T, Generator>, T(cl_ulong, cl_ulong)> {
    static_assert(
            sizeof(T) == 1  ||
            sizeof(T) == 2  ||
            sizeof(T) == 4  ||
            sizeof(T) == 8  ||
            sizeof(T) == 16 ||
            sizeof(T) == 32,
            "Unsupported random output type."
            );

    // TODO: parameter should be same size as ctr_t
    // to allow using full range of the generator.
    typedef typename cl_scalar_of<T>::type Ts;

    static std::string name() {
        return "random_" + type_name<T>() + "_" + Generator::name();
    }

    static void define(backend::source_generator &src) {
        define( src, name() );
    }

    static void define(backend::source_generator &src, const std::string &fname)
    {
        const size_t N = cl_vector_length<T>::value;

        typedef typename std::conditional<
                    sizeof(T) < 32, cl_uint, cl_ulong
                >::type ctr_t;

        const size_t ctr_n = sizeof(T) <= 8 ? 2 : 4;

        typedef typename Generator::template function<ctr_t, ctr_n> generator;

        const size_t key_n = generator::K;

        generator::define(src);

        src.begin_function<T>(fname);
        src.begin_function_parameters();
        src.template parameter<cl_ulong>("prm1");
        src.template parameter<cl_ulong>("prm2");
        src.end_function_parameters();

        src.new_line() << "union ";
        src.open("{");
        src.new_line() << type_name<ctr_t>() << " ctr[" << ctr_n << "];";
        if (std::is_same<Ts, cl_float>::value) {
            src.new_line()
                << type_name<cl_uint>() << " res_i[" << N << "];";
            src.new_line()
                << type_name<cl_float>() << " res_f[" << N << "];";
        } else if (std::is_same<Ts, cl_double>::value) {
            src.new_line()
                << type_name<cl_ulong>() << " res_i[" << N << "];";
            src.new_line()
                << type_name<cl_double>() << " res_f[" << N << "];";
        }
        src.new_line() << type_name<T>() << " res;";
        src.close("} ctr;");

        src.new_line() << type_name<ctr_t>() << " key[" << key_n << "];";

        for(size_t i = 0; i < ctr_n; i += 2)
            src.new_line()
                << "ctr.ctr[" << i     << "] = prm1; "
                << "ctr.ctr[" << i + 1 << "] = prm2;";

        for(size_t i = 0; i < key_n; ++i)
            src.new_line() << "key[" << i << "] = 0x12345678;";

        src.new_line() << generator::name() << "(ctr.ctr, key);";

        if(std::is_same<Ts, cl_float>::value) {
            for(size_t i = 0; i < N; ++i)
                src.new_line()
                    << "ctr.res_f[" << i << "] = ctr.res_i[" << i
                    << "] / " << std::numeric_limits<cl_uint>::max()
                    << ".0f;";
        } else if (std::is_same<Ts, cl_double>::value) {
            for(size_t i = 0; i < N; ++i)
                src.new_line()
                    << "ctr.res_f[" << i << "] = ctr.res_i[" << i
                    << "] / " << std::numeric_limits<cl_ulong>::max()
                    << ".0;";
        }
        src.new_line() << "return ctr.res;";

        src.end_function();
    }
};


/// Returns normally distributed random numbers.
/** Uses Box-Muller transform. */
template <class T, class Generator = random::philox>
struct RandomNormal : UserFunction<RandomNormal<T,Generator>, T(cl_ulong, cl_ulong)> {
    typedef typename cl_scalar_of<T>::type Ts;
    static_assert(
            std::is_same<Ts, cl_float>::value ||
            std::is_same<Ts, cl_double>::value,
            "Must use float or double vector or scalar."
            );

    static std::string name() {
        return "random_normal_" + type_name<T>() + "_" + Generator::name();
    }

    static void define(backend::source_generator &src) {
        define( src, name() );
    }

    static void define(backend::source_generator &src, const std::string &fname)
    {
        const size_t N        = cl_vector_length<T>::value;
        const bool   is_float = std::is_same<Ts, cl_float>::value;
        const size_t ctr_n    = is_float ? 2 : 4;

        typedef typename Generator::template function<cl_uint, ctr_n> generator;

        const size_t key_n = generator::K;

        generator::define(src);

        src.begin_function<T>(fname);
        src.begin_function_parameters();
        src.template parameter<cl_ulong>("prm1");
        src.template parameter<cl_ulong>("prm2");
        src.end_function_parameters();

#if defined(VEXCL_BACKEND_JIT)
        src.new_line() << "#define cospi(x) cos(M_PI * (x))";
#endif

        src.new_line() << "union ";
        src.open("{");
        src.new_line() << type_name<cl_uint>() << " ctr[" << ctr_n << "];";
        if (is_float) {
            src.new_line() << type_name<cl_uint>()  << " res_i[2];";
        } else {
            src.new_line() << type_name<cl_ulong>()  << " res_i[2];";
        }
        src.close("} ctr;");
        src.new_line() << type_name<Ts>() << " u[2];";

        src.new_line() << type_name<cl_uint>() << " key[" << key_n << "];";

        for(size_t i = 0; i < ctr_n; i += 2)
            src.new_line()
                << "ctr.ctr[" << i     << "] = prm1; "
                << "ctr.ctr[" << i + 1 << "] = prm2;";

        for(size_t i = 0; i < key_n; ++i)
            src.new_line() << "key[" << i << "] = 0x12345678;";

        if (N > 1) {
            src.new_line() << "union ";
            src.open("{");
            src.new_line() << type_name<Ts>() << " z[" << N << "];";
            src.new_line() << type_name<T>() << " v;";
            src.close("} res;");
        }

        for(size_t i = 0 ; i < N ; i += 2) {
            src.new_line() << generator::name() << "(ctr.ctr, key);";

            if(is_float) {
                for(size_t i = 0; i < 2; ++i)
                    src.new_line()
                        << "u[" << i << "] = ctr.res_i[" << i
                        << "] / " << std::numeric_limits<cl_uint>::max()
                        << ".0f;";
            } else {
                for(size_t i = 0; i < 2; ++i)
                    src.new_line()
                        << "u[" << i << "] = ctr.res_i[" << i
                        << "] / " << std::numeric_limits<cl_ulong>::max()
                        << ".0;";
            }

            if(N == 1) {
                src.new_line()
                    << "return sqrt(-2 * log(u[0])) * cospi(2 * u[1]);\n";
            } else {
                src.open("{");

                src.new_line() << type_name<Ts>()
                    << " l = sqrt(-2 * log(u[0])), cs, sn;";

#if defined(VEXCL_BACKEND_CUDA)
                src.new_line() << "sincospi(2 * u[1], &sn, &cs);";
#else
                src.new_line() << "sn = sincos("
                    << std::setprecision(16)
                    << boost::math::constants::two_pi<double>()
                    << " * u[1], &cs);";
#endif
                src.new_line() << "res.z[" << i     << "] = l * cs;";
                src.new_line() << "res.z[" << i + 1 << "] = l * sn;";

                src.close("}");
            }
        }

        if (N > 1)
            src.new_line() << "return res.v;";

        src.end_function();
    }
};




} // namespace vex



#endif
