#ifndef VEXCL_FFT_PLAN_HPP
#define VEXCL_FFT_PLAN_HPP

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
 * \file   vexcl/fft/plan.hpp
 * \author Pascal Germroth <pascal@ensieve.org>
 * \brief  FFT plan, stores kernels and buffers for one configuration.
 */

#include <cmath>
#include <queue>
#include <numeric>

#include <vexcl/profiler.hpp>
#include <vexcl/vector.hpp>
#include <vexcl/function.hpp>
#include <vexcl/fft/unrolled_dft.hpp>
#include <vexcl/fft/kernels.hpp>

namespace vex {

/// Fast Fourier Transform
namespace fft {

/// FFT direction.
enum direction {
    forward, ///< Forward transform.
    inverse, ///< Inverse transform.
    none     ///< Specifies dimension(s) to do batch transform.
};

// Returns successive prime numbers on each call.
struct prime_generator {
    typedef std::pair<size_t, size_t> P;

    std::priority_queue<P, std::vector<P>, std::greater<P>> cross;
    size_t x;
    prime_generator() : x(2) {}

    size_t operator()() {
        for(;; x++) {
            if(cross.empty() || cross.top().first != x) { // is a prime.
                cross.push(P(x * x, x));
                return x++;
            } else { // has prime factor. cross out.
                while(cross.top().first == x) {
                    size_t prime = cross.top().second;
                    cross.pop();
                    cross.push(P(x + prime, prime));
                }
            }
        }
    }
};

// Returns the prime factors of a number.
inline std::vector<pow> prime_factors(size_t n) {
    std::vector<pow> fs;
    if(n != 0) {
        prime_generator next;
        while(n != 1) {
            const auto prime = next();
            size_t exp = 0;
            while(n % prime == 0) {
                n /= prime;
                exp++;
            }
            if(exp != 0)
                fs.push_back(pow(prime, exp));
        }
    }
    return fs;
}

// for arbitrary x, return smallest y=a^b c^d e^f...>=x where a,c,e are iterated over (assumed to be prime)
template<class Iterator>
inline size_t next_prime_power(Iterator begin, Iterator end, size_t target, size_t n = 1, size_t best = -1) {
    const size_t prime = *begin++;
    for(; n < target ; n *= prime) {
        if(begin != end)
            best = next_prime_power(begin, end, target, n, best);
    }
    return std::min(n, best);
}



struct planner {
    const size_t max_size;
    std::vector<size_t> primes;

    planner(size_t s = 25) : max_size(std::min(s, supported_kernel_sizes().back())) {
        auto ps = supported_primes();
        for(auto i = ps.begin() ; i != ps.end() ; i++)
            if(*i <= s) primes.push_back(*i);
    }

    // returns the size the data must be padded to.
    size_t best_size(size_t n) const {
        return next_prime_power(primes.begin(), primes.end(), n);
    }

    // splits n into a list of powers 2^a 2^b 2^c 3^d 5^e...
    // exponents are limited by available kernels,
    // if no kernel for prime is available, exponent will be 0, interpret as 1.
    std::vector<pow> factor(size_t n) const {
        std::vector<pow> out, factors = prime_factors(n);
        for(auto f = factors.begin() ; f != factors.end() ; f++) {
            if(std::find(primes.begin(), primes.end(), f->base) != primes.end()) {
                // split exponent into reasonable parts.
                auto qs = stages(*f);
                // use smallest radix first
                std::copy(qs.rbegin(), qs.rend(), std::back_inserter(out));
            } else {
                // unsupported prime.
                for(size_t i = 0 ; i != f->exponent ; i++)
                    out.push_back(pow(f->base, 0));
            }
        }
        return out;
    }

  private:
    std::vector<pow> stages(pow p) const {
        size_t t = static_cast<size_t>(std::log(max_size + 1.0) / std::log(static_cast<double>(p.base)));
        std::vector<pow> fs;
#ifdef FFT_SIMPLE_PLANNER
        // use largest radixes, i.e. 2^4 2^4 2^1
        for(size_t e = p.exponent ; ; ) {
            if(e > t) {
                fs.push_back(pow(p.base, t));
                e -= t;
            } else if(e <= t) {
                fs.push_back(pow(p.base, e));
                break;
            }
        }
#else
        // avoid very small radixes, i.e. 2^3 2^3 2^3
        // number of parts
        size_t m = (p.exponent + t - 1) / t;
        // two levels.
        size_t r = t * m - p.exponent;
        size_t u = m * (r / m);
        size_t v = t - (r / m);
        for(size_t i = 0 ; i < m - r + u ; i++)
            fs.push_back(pow(p.base, v));
        for(size_t i = 0 ; i < r - u ; i++)
            fs.push_back(pow(p.base, v - 1));
#endif
        return fs;
    }
};


template <class Tv, class Planner = planner>
struct plan {
    typedef typename cl_scalar_of<Tv>::type Ts;
    static_assert(
            std::is_same<Ts, cl_float >::value ||
            std::is_same<Ts, cl_double>::value,
            "Only float and double data supported."
            );

    typedef typename cl_vector_of<Ts, 2>::type T2;

    VEX_FUNCTION_S(T2, r2c, (Ts, v), type_name<T2>() + " r = {v, 0}; return r;");
    VEX_FUNCTION_S(Ts, c2r, (T2, v), "return v.x;");
    VEX_FUNCTION_S(T2, scl, (T2, v)(Ts, s), "v.x *= s; v.y *= s; return v;");

    const std::vector<backend::command_queue> &queues;
    Planner planner;
    Ts scale;
    const std::vector<size_t> sizes;

    std::vector<kernel_call> kernels;

    size_t input, output;
    std::vector< vex::vector<T2> > bufs;

    profiler<> *profile;

    // \param sizes
    //  1D case: {n}.
    //  2D case: {h, w} in row-major format: x + y * w. (like FFTw)
    //  etc.
    plan(const std::vector<backend::command_queue> &_queues, const std::vector<size_t> sizes,
        const std::vector<direction> dirs, const Planner &planner = Planner())
        : queues(_queues), planner(planner), sizes(sizes), profile(NULL)
    {
        assert(sizes.size() >= 1);
        assert(sizes.size() == dirs.size());

        precondition(
                queues.size() == 1,
                "FFT is only supported for single-device contexts."
                );

        auto queue   = queues[0];

        size_t total_n = std::accumulate(sizes.begin(), sizes.end(),
            static_cast<size_t>(1), std::multiplies<size_t>());
        size_t current = bufs.size(); bufs.push_back(vex::vector<T2>(queues, total_n));
        size_t other   = bufs.size(); bufs.push_back(vex::vector<T2>(queues, total_n));

        size_t inv_n = 1;
        for(size_t i = 0 ; i < sizes.size() ; i++)
            if(dirs[i] == inverse)
                inv_n *= sizes[i];
        scale = (Ts)1 / inv_n;

        // Build the list of kernels.
        input = current;
        for(size_t i = 1 ; i <= sizes.size() ; i++) {
            const size_t j = sizes.size() - i;
            const size_t w = sizes[j], h = total_n / w;
            if(w > 1) {
                // 1D, each row.
                if(dirs[j] != none)
                    plan_cooley_tukey(dirs[j] == inverse, w, h, current, other, false);

                if(h > 1 && !(dirs.size() == 2 && dirs[0] == none)) {
                    kernels.push_back(transpose_kernel<Ts>(queue, w, h, bufs[current](), bufs[other]()));
                    std::swap(current, other);
                }
            }
        }
        output = current;
    }

    void plan_cooley_tukey(bool inverse, size_t n, size_t batch, size_t &current, size_t &other, bool once) {
        size_t p = 1;
        auto rs = planner.factor(n);
        for(auto r = rs.begin() ; r != rs.end() ; r++) {
            if(r->exponent == 0) {
                plan_bluestein(n, batch, inverse, r->base, p, current, other);
                p *= r->base;
            } else {
                kernels.push_back(radix_kernel<Ts>(once, queues[0], n, batch,
                    inverse, *r, p, bufs[current](), bufs[other]()));
                std::swap(current, other);
                p *= r->value;
            }
        }
    }

    void plan_bluestein(size_t width, size_t batch, bool inverse, size_t n, size_t p, size_t &current, size_t &other) {
        size_t conv_n = planner.best_size(2 * n);
        size_t threads = width / n;

        size_t b_twiddle = bufs.size(); bufs.push_back(vex::vector<T2>(queues, n));
        size_t b_other   = bufs.size(); bufs.push_back(vex::vector<T2>(queues, conv_n));
        size_t b_current = bufs.size(); bufs.push_back(vex::vector<T2>(queues, conv_n));
        size_t a_current = bufs.size(); bufs.push_back(vex::vector<T2>(queues, conv_n * batch * threads));
        size_t a_other   = bufs.size(); bufs.push_back(vex::vector<T2>(queues, conv_n * batch * threads));

        // calculate twiddle factors
        kernels.push_back(bluestein_twiddle<Ts>(queues[0], n, inverse,
            bufs[b_twiddle]())); // once

        // first part of the convolution
        kernels.push_back(bluestein_pad_kernel<Ts>(queues[0], n, conv_n,
            bufs[b_twiddle](), bufs[b_current]())); // once

        plan_cooley_tukey(false, conv_n, 1, b_current, b_other, /*once*/true);

        // other part of convolution
        kernels.push_back(bluestein_mul_in<Ts>(queues[0], inverse, batch, n, p, threads, conv_n,
            bufs[current](), bufs[b_twiddle](), bufs[a_current]()));

        plan_cooley_tukey(false, conv_n, threads * batch, a_current, a_other, false);

        // calculate convolution
        kernels.push_back(bluestein_mul<Ts>(queues[0], conv_n, threads * batch,
            bufs[a_current](), bufs[b_current](), bufs[a_other]()));
        std::swap(a_current, a_other);

        plan_cooley_tukey(true, conv_n, threads * batch, a_current, a_other, false);

        // twiddle again
        kernels.push_back(bluestein_mul_out<Ts>(queues[0], batch, p, n, threads, conv_n,
            bufs[a_current](), bufs[b_twiddle](), bufs[other]()));
        std::swap(current, other);
    }

    // Execute the complete transformation.
    // Converts real-valued input and output, supports multiply-adding to output.
    template<class Expr>
    void transform(const Expr &in) {
        if(profile) {
            std::ostringstream prof_name;
            prof_name << "fft(n={";
            for(size_t i = 0 ; i < sizes.size() ; i++) {
                if(i != 0) prof_name << ", ";
                prof_name << sizes[i];
            }
            prof_name << "}, scale=" << scale << ")";
            profile->tic_cl(prof_name.str());
            profile->tic_cl("in");
        }
        vector<T2> &in_c = bufs[input];
        if(cl_vector_length<Tv>::value == 1) in_c = r2c(in);
        else in_c = in;
        if(profile) profile->toc("in");
        for(auto run = kernels.begin(); run != kernels.end(); ++run) {
            if(!run->once || run->count == 0) {
#ifdef FFT_DUMP_ARRAYS
                for(auto b = bufs.begin(); b != bufs.end(); ++b)
                    std::cerr << "   " << std::setprecision(2) << (*b)()() << " = " << (*b) << std::endl;
                std::cerr << "run " << run->desc << std::endl;
#endif
                if(profile) {
                    std::ostringstream s;
                    s << "(k" << (run - kernels.begin()) << ")";
                    if(run->once) s << " ONCE";
                    s << " " << run->desc;
                    profile->tic_cl(s.str());
                }
                run->kernel(queues[0]);
                run->count++;
                if(profile) profile->toc("");
            }
        }
#ifdef FFT_DUMP_ARRAYS
        for(auto b = bufs.begin(); b != bufs.end(); ++b)
            std::cerr << "   " << (*b)()() << " = " << (*b) << std::endl;
#endif
        if (profile) profile->toc("");
    }

    template <typename Tout, class Expr>
    auto apply(const Expr &expr) ->
        typename std::enable_if<
            cl_vector_length<Tout>::value == 1,
            decltype( scale * c2r(bufs[output]) )
        >::type
    {
        transform(expr);
        return scale * c2r(bufs[output]);
    }

    template <typename Tout, class Expr>
    auto apply(const Expr &expr) ->
        typename std::enable_if<
            cl_vector_length<Tout>::value == 2,
            decltype( scl(bufs[output], scale) )
        >::type
    {
        transform(expr);
        return scl(bufs[output], scale);
    }

    std::string desc() const {
        std::ostringstream o;
        o << "FFT(";
        // sizes
        for(auto n = sizes.begin() ; n != sizes.end() ; n++) {
            if(n != sizes.begin()) o << " x ";
            o << *n;
            auto fs = prime_factors(*n);
            if(fs.size() > 1) {
                o << '=';
                for(auto f = fs.begin() ; f != fs.end() ; f++) {
                    if(f != fs.begin()) o << '*';
                    o << *f;
                }
            }
        }
        o << ")";
        return o.str();
    }
};


template <class T, class P>
inline std::ostream &operator<<(std::ostream &o, const plan<T,P> &p) {
    o << p.desc() << "{\n";
    for(auto k = p.kernels.begin() ; k != p.kernels.end() ; k++) {
        o << "  ";
        if(k->once) o << "once: ";
        o << k->desc << "\n";
    }
    return o << "}";
}

} // namespace fft
} // namespace vex
#endif
