#ifndef VEXCL_FFT_KERNELS_HPP
#define VEXCL_FFT_KERNELS_HPP

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
 * \file   vexcl/fft/kernels.hpp
 * \author Pascal Germroth <pascal@ensieve.org>
 * \brief  Kernel generator for FFT.
 */

#include <boost/math/constants/constants.hpp>

#ifndef VEX_FAST_MATH_OPTS
#  if defined(VEXCL_BACKEND_OPENCL) || defined(VEXCL_BACKEND_COMPUTE)
#    define VEX_FAST_MATH_OPTS "-cl-mad-enable -cl-fast-relaxed-math"
#  elif defined(VEXCL_BACKEND_CUDA)
#    define VEX_FAST_MATH_OPTS "--use_fast_math"
#  elif defined(VEXCL_BACKEND_JIT)
#    define VEX_FAST_MATH_OPTS "-ffast-math"
#  endif
#endif

namespace vex {
namespace fft {

struct pow {
    size_t base, exponent, value;
    pow(size_t b, size_t e) : base(b), exponent(e),
        value(static_cast<size_t>(std::pow(static_cast<double>(b), static_cast<double>(e)))) {}
};

inline std::ostream &operator<<(std::ostream &o, const pow &p) {
    o << p.base;
    if(p.exponent != 1) o << '^' << p.exponent;
    return o;
}

struct kernel_call {
    bool once;
    size_t count;
    std::string desc;
    backend::kernel kernel;
    kernel_call(bool o, std::string d, backend::kernel k)
        : once(o), count(0), desc(d), kernel(k)
    {}
};



// generates "(prefix vfrom,vfrom+1,...,vto)"
inline void param_list(backend::source_generator &o,
        std::string prefix, size_t from, size_t to, size_t step = 1)
{
    o << '(';
    for(size_t i = from ; i != to ; i += step) {
        if(i != from) o << ", ";
        o << prefix << 'v' << i;
    } o << ')';
}

template <class T, class T2>
inline void kernel_radix(backend::source_generator &o, pow radix, bool invert) {
    o << in_place_dft(radix.value, invert);

    // kernel.
    o.begin_kernel("radix");
    o.begin_kernel_parameters();
    o.template parameter< global_ptr<const T2> >("x");
    o.template parameter< global_ptr<      T2> >("y");
    o.template parameter< cl_uint              >("p");
    o.template parameter< cl_uint              >("threads");
    o.end_kernel_parameters();

    o.new_line() << "const size_t i = " << o.global_id(0) << ";";
    o.new_line() << "if(i >= threads) return;";

    // index in input sequence, in 0..P-1
    o.new_line() << "const size_t k = i % p;";
    o.new_line() << "const size_t batch_offset = " << o.global_id(1) << " * threads * " << radix.value << ";";

    // read
    o.new_line() << "x += i + batch_offset;";
    for(size_t i = 0; i < radix.value; ++i)
        o.new_line() << type_name<T2>() << " v" << i << " = x[" << i << " * threads];";

    // twiddle
    o.new_line() << "if(p != 1)";
    o.open("{");
    for(size_t i = 1; i < radix.value; ++i) {
        const T alpha = -boost::math::constants::two_pi<T>() * i / radix.value;
        o.new_line() << "v" << i << " = mul(v" << i << ", twiddle("
          << "(" << type_name<T>() << ")" << std::setprecision(16) << alpha << " * k / p));";
    }
    o.close("}");

    // inplace DFT
    o.new_line() << "dft" << radix.value;
    param_list(o, "&", 0, radix.value);
    o << ";";

    // write back
    o.new_line() << "const size_t j = k + (i - k) * " << radix.value << ";";
    o.new_line() << "y += j + batch_offset;";
    for(size_t i = 0; i < radix.value; i++)
        o.new_line() << "y[" << i << " * p] = v" << i << ";";
    o.end_kernel();
}


template <class T>
inline std::string fft_kernel_header() {
    std::ostringstream src;
    src
#ifndef VEXCL_BACKEND_CUDA
      << "#define DEVICE\n"
#else
      << "#define DEVICE __device__\n"
#endif
      << "typedef " << type_name<T>() << " real_t;\n"
      << "typedef " << type_name<T>() << "2 real2_t;\n";

    return src.str();
}

// Return A*B (complex multiplication)
template <class T2>
inline void mul_code(backend::source_generator &o, bool invert) {
    o.begin_function<T2>("mul");
    o.begin_function_parameters();
    o.template parameter<T2>("a");
    o.template parameter<T2>("b");
    o.end_function_parameters();

    if(invert) { // conjugate b
        o.new_line() << type_name<T2>() << " r = {"
            "a.x * b.x + a.y * b.y, "
            "a.y * b.x - a.x * b.y};";
    } else {
        o.new_line() << type_name<T2>() << " r = {"
            "a.x * b.x - a.y * b.y, "
            "a.y * b.x + a.x * b.y};";
    }

    o.new_line() << "return r;";
    o.end_function();
}

// A * exp(alpha * I) == A  * (cos(alpha) + I * sin(alpha))
// native_cos(), native_sin() is a *lot* faster than sincos, on nVidia.
template <class T, class T2>
inline void twiddle_code(backend::source_generator &o) {
    o.begin_function<T2>("twiddle");
    o.begin_function_parameters();
    o.template parameter<T>("alpha");
    o.end_function_parameters();

    if(std::is_same<T, cl_double>::value) {
        // use sincos with double since we probably want higher precision
#if defined(VEXCL_BACKEND_OPENCL) || defined(VEXCL_BACKEND_COMPUTE)
        o.new_line() << type_name<T>() << " cs, sn = sincos(alpha, &cs);";
#else
        o.new_line() << type_name<T>() << " sn, cs;";
        o.new_line() << "sincos(alpha, &sn, &cs);";
#endif
        o.new_line() << type_name<T2>() << " r = {cs, sn};";
    } else {
        // use native with float since we probably want higher performance
#if defined(VEXCL_BACKEND_OPENCL) || defined(VEXCL_BACKEND_COMPUTE)
        o.new_line() << type_name<T2>() << " r = {"
            "native_cos(alpha), native_sin(alpha)};";
#elif defined(VEXCL_BACKEND_CUDA)
        o.new_line() << type_name<T>() << " sn, cs;";
        o.new_line() << "__sincosf(alpha, &sn, &cs);";
        o.new_line() << type_name<T2>() << " r = {cs, sn};";
#elif defined(VEXCL_BACKEND_JIT)
        o.new_line() << type_name<T>() << " sn, cs;";
        o.new_line() << "sincosf(alpha, &sn, &cs);";
        o.new_line() << type_name<T2>() << " r = {cs, sn};";
#else
#  error Unsupported backend!
#endif
    }

    o.new_line() << "return r;";
    o.end_function();
}


template <class T, class T2>
inline kernel_call radix_kernel(
        bool once, const backend::command_queue &queue, size_t n, size_t batch,
        bool invert, pow radix, size_t p,
        const backend::device_vector<T2> &in,
        const backend::device_vector<T2> &out
        )
{
    scoped_program_header header(queue, fft_kernel_header<T>());

    backend::source_generator o(queue);
    o << std::setprecision(25);

    mul_code<T2>(o, invert);
    twiddle_code<T, T2>(o);

    const size_t m = n / radix.value;
    kernel_radix<T, T2>(o, radix, invert);

    backend::kernel kernel(queue, o.str(), "radix", 0, VEX_FAST_MATH_OPTS);

    kernel.push_arg(in);
    kernel.push_arg(out);
    kernel.push_arg(static_cast<cl_uint>(p));
    kernel.push_arg(static_cast<cl_uint>(m));

    const size_t wg_mul = kernel.preferred_work_group_size_multiple(queue);
    //const size_t max_cu = device.getInfo<CL_DEVICE_MAX_COMPUTE_UNITS>();
    //const size_t max_wg = device.getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>();
    size_t wg = wg_mul;
    //while(wg * max_cu < max_wg) wg += wg_mul;
    //wg -= wg_mul;
    const size_t threads = (m + wg - 1) / wg;

    kernel.config(backend::ndrange(threads, batch), backend::ndrange(wg, 1));

    std::ostringstream desc;
    desc << "dft{r=" << radix << ", p=" << p << ", n=" << n << ", batch=" << batch << ", threads=" << m << "(" << threads << "), wg=" << wg << "}";

    return kernel_call(once, desc.str(), kernel);
}


template <class T, class T2>
inline kernel_call transpose_kernel(
        const backend::command_queue &queue, size_t width, size_t height,
        const backend::device_vector<T2> &in,
        const backend::device_vector<T2> &out
        )
{
    scoped_program_header header(queue, fft_kernel_header<T>());

    backend::source_generator o(queue);
    o << std::setprecision(25);

    // determine max block size to fit into local memory/workgroup
    size_t block_size = is_cpu(queue) ? 1 : 128;
    {
#if defined(VEXCL_BACKEND_OPENCL) || defined(VEXCL_BACKEND_COMPUTE)
        cl_device_id dev = backend::get_device_id(queue);
        cl_ulong local_size;
        size_t workgroup;
        clGetDeviceInfo(dev, CL_DEVICE_LOCAL_MEM_SIZE, sizeof(cl_ulong), &local_size, NULL);
        clGetDeviceInfo(dev, CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(size_t), &workgroup, NULL);
#else
        const auto local_size = queue.device().max_shared_memory_per_block();
        const auto workgroup = queue.device().max_threads_per_block();
#endif
        while(block_size * block_size * sizeof(T) * 2 > local_size) block_size /= 2;
        while(block_size * block_size > workgroup) block_size /= 2;
    }

    // from NVIDIA SDK.
    o.begin_kernel("transpose");
    o.begin_kernel_parameters();
    o.template parameter< global_ptr<const T2> >("input");
    o.template parameter< global_ptr<      T2> >("output");
    o.template parameter< cl_uint              >("width");
    o.template parameter< cl_uint              >("height");
    o.end_kernel_parameters();

    o.new_line() << "const size_t global_x = " << o.global_id(0) << ";";
    o.new_line() << "const size_t global_y = " << o.global_id(1) << ";";
    o.new_line() << "const size_t local_x  = " << o.local_id(0)  << ";";
    o.new_line() << "const size_t local_y  = " << o.local_id(1)  << ";";
    o.new_line() << "const bool range = global_x < width && global_y < height;";

    // local memory
    {
        std::ostringstream s;
        s << "block[" << block_size * block_size << "]";
        o.smem_static_var(type_name<T2>(), s.str());
    }

    // copy from input to local memory
    o.new_line() << "if(range) "
        << "block[local_x + local_y * " << block_size << "] = input[global_x + global_y * width];";

    // wait until the whole block is filled
    o.new_line().barrier();

    // transpose local block to target
    o.new_line() << "if(range) "
      << "output[global_x * height + global_y] = block[local_x + local_y * " << block_size << "];";

    o.end_kernel();

    backend::kernel kernel(queue, o.str(), "transpose");

    kernel.push_arg(in);
    kernel.push_arg(out);
    kernel.push_arg(static_cast<cl_uint>(width));
    kernel.push_arg(static_cast<cl_uint>(height));

    // range multiple of wg size, last block maybe not completely filled.
    size_t r_w = (width  + block_size - 1) / block_size;
    size_t r_h = (height + block_size - 1) / block_size;

    kernel.config(backend::ndrange(r_w, r_h), backend::ndrange(block_size, block_size));

    std::ostringstream desc;
    desc << "transpose{"
         << "w=" << width << "(" << r_w << "), "
         << "h=" << height << "(" << r_h << "), "
         << "bs=" << block_size << "}";

    return kernel_call(false, desc.str(), kernel);
}



template <class T, class T2>
inline kernel_call bluestein_twiddle(
        const backend::command_queue &queue, size_t n, bool inverse,
        const backend::device_vector<T2> &out
        )
{
    scoped_program_header header(queue, fft_kernel_header<T>());

    backend::source_generator o(queue);
    o << std::setprecision(25);

    twiddle_code<T, T2>(o);

    o.begin_kernel("bluestein_twiddle");
    o.begin_kernel_parameters();
    o.template parameter< size_t         >("n");
    o.template parameter< global_ptr<T2> >("output");
    o.end_kernel_parameters();

    o.new_line() << "const size_t x = " << o.global_id(0) << ";";

    o.new_line() << "const size_t xx = ((ulong)x * x) % (2 * n);";
    o.new_line() << "if (x < n) output[x] = twiddle("
        << std::setprecision(16)
        << (inverse ? 1 : -1) * boost::math::constants::pi<T>()
        << " * xx / n);";

    o.end_kernel();

    backend::kernel kernel(queue, o.str(), "bluestein_twiddle");
    kernel.push_arg(n);
    kernel.push_arg(out);

    size_t ws = kernel.preferred_work_group_size_multiple(queue);
    size_t gs = (n + ws - 1) / ws;

    kernel.config(gs, ws);

    std::ostringstream desc;
    desc << "bluestein_twiddle{n=" << n << ", inverse=" << inverse << "}";
    return kernel_call(true, desc.str(), kernel);
}

template <class T, class T2>
inline kernel_call bluestein_pad_kernel(
        const backend::command_queue &queue, size_t n, size_t m,
        const backend::device_vector<T2> &in,
        const backend::device_vector<T2> &out
        )
{
    scoped_program_header header(queue, fft_kernel_header<T>());

    backend::source_generator o(queue);
    o << std::setprecision(25);

    o.begin_function<T2>("conj");
    o.begin_function_parameters();
    o.template parameter<T2>("v");
    o.end_function_parameters();
    o.new_line() << type_name<T2>() << " r = {v.x, -v.y};";
    o.new_line() << "return r;";
    o.end_function();

    o.begin_kernel("bluestein_pad_kernel");
    o.begin_kernel_parameters();
    o.template parameter< global_ptr<const T2> >("input");
    o.template parameter< global_ptr<      T2> >("output");
    o.template parameter< cl_uint              >("n");
    o.template parameter< cl_uint              >("m");
    o.end_kernel_parameters();
    o.new_line() << "const uint x = " << o.global_id(0) << ";";
    o.new_line() << "if (x < m)";
    o.open("{");
    o.new_line() << "if(x < n || m - x < n)";
    o.open("{");
    o.new_line() << "output[x] = conj(input[min(x, m - x)]);";
    o.close("}");
    o.new_line() << "else";
    o.open("{");
    o.new_line() << type_name<T2>() << " r = {0,0};";
    o.new_line() << "output[x] = r;";
    o.close("}");
    o.close("}");
    o.end_kernel();

    backend::kernel kernel(queue, o.str(), "bluestein_pad_kernel");
    kernel.push_arg(in);
    kernel.push_arg(out);
    kernel.push_arg(static_cast<cl_uint>(n));
    kernel.push_arg(static_cast<cl_uint>(m));

    size_t ws = kernel.preferred_work_group_size_multiple(queue);
    size_t gs = (m + ws - 1) / ws;

    kernel.config(gs, ws);

    std::ostringstream desc;
    desc << "bluestein_pad_kernel{n=" << n << ", m=" << m << "}";
    return kernel_call(true, desc.str(), kernel);
}

template <class T, class T2>
inline kernel_call bluestein_mul_in(
        const backend::command_queue &queue, bool inverse, size_t batch,
        size_t radix, size_t p, size_t threads, size_t stride,
        const backend::device_vector<T2> &data,
        const backend::device_vector<T2> &exp,
        const backend::device_vector<T2> &out
        )
{
    scoped_program_header header(queue, fft_kernel_header<T>());

    backend::source_generator o(queue);
    o << std::setprecision(25);

    mul_code<T2>(o, false);
    twiddle_code<T, T2>(o);

    o.begin_kernel("bluestein_mul_in");
    o.begin_kernel_parameters();
    o.template parameter< global_ptr<const T2> >("data");
    o.template parameter< global_ptr<const T2> >("exp");
    o.template parameter< global_ptr<      T2> >("output");
    o.template parameter< cl_uint              >("radix");
    o.template parameter< cl_uint              >("p");
    o.template parameter< cl_uint              >("out_stride");
    o.end_kernel_parameters();

    o.new_line() << "const size_t thread  = " << o.global_id(0)   << ";";
    o.new_line() << "const size_t threads = " << o.global_size(0) << ";";
    o.new_line() << "const size_t batch   = " << o.global_id(1)   << ";";
    o.new_line() << "const size_t element = " << o.global_id(2)   << ";";

    o.new_line() << "if(element < out_stride)";
    o.open("{");

    o.new_line() << "const size_t in_off  = thread + batch * radix * threads + element * threads;";
    o.new_line() << "const size_t out_off = thread * out_stride + batch * out_stride * threads + element;";

    o.new_line() << "if(element < radix)";
    o.open("{");

    o.new_line() << type_name<T2>() << " w = exp[element];";

    o.new_line() << "if(p != 1)";
    o.open("{");

    o.new_line() << "ulong a = (ulong)element * (thread % p);";
    o.new_line() << "ulong b = (ulong)radix * p;";
    o.new_line() << type_name<T2>() << " t = twiddle(" << std::setprecision(16)
        << (inverse ? 1 : -1) * boost::math::constants::two_pi<T>()
        << " * (a % (2 * b)) / b);";
    o.new_line() << "w = mul(w, t);";
    o.close("}");

    o.new_line() << "output[out_off] = mul(data[in_off], w);";

    o.close("}");
    o.new_line() << "else";
    o.open("{");

    o.new_line() << type_name<T2>() << " r = {0,0};";
    o.new_line() << "output[out_off] = r;";

    o.close("}");
    o.close("}");
    o.end_kernel();

    backend::kernel kernel(queue, o.str(), "bluestein_mul_in");
    kernel.push_arg(data);
    kernel.push_arg(exp);
    kernel.push_arg(out);
    kernel.push_arg(static_cast<cl_uint>(radix));
    kernel.push_arg(static_cast<cl_uint>(p));
    kernel.push_arg(static_cast<cl_uint>(stride));

    const size_t wg = kernel.preferred_work_group_size_multiple(queue);
    const size_t stride_pad = (stride + wg - 1) / wg;

    kernel.config(
            backend::ndrange(threads, batch, stride_pad),
            backend::ndrange(      1,     1,         wg)
            );

    std::ostringstream desc;
    desc << "bluestein_mul_in{batch=" << batch << ", radix=" << radix << ", p=" << p << ", threads=" << threads << ", stride=" << stride << "(" << stride_pad << "), wg=" << wg << "}";
    return kernel_call(false, desc.str(), kernel);
}

template <class T, class T2>
inline kernel_call bluestein_mul_out(
        const backend::command_queue &queue, size_t batch, size_t p,
        size_t radix, size_t threads, size_t stride,
        const backend::device_vector<T2> &data,
        const backend::device_vector<T2> &exp,
        const backend::device_vector<T2> &out
        )
{
    scoped_program_header header(queue, fft_kernel_header<T>());

    backend::source_generator o(queue);
    o << std::setprecision(25);

    mul_code<T2>(o, false);

    o.begin_function<T2>("scale");
    o.begin_function_parameters();
    o.template parameter<T2>("x");
    o.template parameter<T >("a");
    o.end_function_parameters();

    o.new_line() << type_name<T2>() << " r = {x.x * a, x.y * a};";
    o.new_line() << "return r;";
    o.end_function();

    o.begin_kernel("bluestein_mul_out");
    o.begin_kernel_parameters();
    o.template parameter< global_ptr<const T2> >("data");
    o.template parameter< global_ptr<const T2> >("exp");
    o.template parameter< global_ptr<      T2> >("output");
    o.template parameter< T                    >("div");
    o.template parameter< cl_uint              >("p");
    o.template parameter< cl_uint              >("in_stride");
    o.template parameter< cl_uint              >("radix");
    o.end_kernel_parameters();

    o.new_line() << "const size_t i = " << o.global_id(0) << ";";
    o.new_line() << "const size_t threads = " << o.global_size(0) << ";";
    o.new_line() << "const size_t b = " << o.global_id(1) << ";";
    o.new_line() << "const size_t l = " << o.global_id(2) << ";";

    o.new_line() << "if(l < radix)";
    o.open("{");

    o.new_line() << "const size_t k = i % p;";
    o.new_line() << "const size_t j = k + (i - k) * radix;";
    o.new_line() << "const size_t in_off = i * in_stride + b * in_stride * threads + l;";
    o.new_line() << "const size_t out_off = j + b * threads * radix + l * p;";

    o.new_line() << "output[out_off] = mul(scale(data[in_off], div), exp[l]);";

    o.close("}");
    o.end_kernel();

    backend::kernel kernel(queue, o.str(), "bluestein_mul_out");
    kernel.push_arg(data);
    kernel.push_arg(exp);
    kernel.push_arg(out);
    kernel.push_arg(static_cast<T>(1.0 / stride));
    kernel.push_arg(static_cast<cl_uint>(p));
    kernel.push_arg(static_cast<cl_uint>(stride));
    kernel.push_arg(static_cast<cl_uint>(radix));

    const size_t wg = kernel.preferred_work_group_size_multiple(queue);
    const size_t radix_pad = (radix + wg - 1) / wg;

    kernel.config(
            backend::ndrange(threads, batch, radix_pad),
            backend::ndrange(      1,     1,        wg)
            );

    std::ostringstream desc;
    desc << "bluestein_mul_out{r=" << radix << "(" << radix_pad << "), wg=" << wg << ", batch=" << batch << ", p=" << p << ", thr=" << threads << ", stride=" << stride << "}";
    return kernel_call(false, desc.str(), kernel);
}

template <class T, class T2>
inline kernel_call bluestein_mul(
        const backend::command_queue &queue, size_t n, size_t batch,
        const backend::device_vector<T2> &data,
        const backend::device_vector<T2> &exp,
        const backend::device_vector<T2> &out
        )
{
    scoped_program_header header(queue, fft_kernel_header<T>());

    backend::source_generator o(queue);
    o << std::setprecision(25);

    mul_code<T2>(o, false);

    o.begin_kernel("bluestein_mul");
    o.begin_kernel_parameters();
    o.template parameter< global_ptr<const T2> >("data");
    o.template parameter< global_ptr<const T2> >("exp");
    o.template parameter< global_ptr<      T2> >("output");
    o.template parameter< cl_uint              >("stride");
    o.end_kernel_parameters();

    o.new_line() << "const size_t x = " << o.global_id(0) << ";";
    o.new_line() << "const size_t y = " << o.global_id(1) << ";";

    o.new_line() << "if(x < stride)";
    o.open("{");

    o.new_line() << "const size_t off = x + stride * y;";
    o.new_line() << "output[off] = mul(data[off], exp[x]);";

    o.close("}");
    o.end_kernel();

    backend::kernel kernel(queue, o.str(), "bluestein_mul");
    kernel.push_arg(data);
    kernel.push_arg(exp);
    kernel.push_arg(out);
    kernel.push_arg(static_cast<cl_uint>(n));

    const size_t wg = kernel.preferred_work_group_size_multiple(queue);
    const size_t threads = (n + wg - 1) / wg;

    kernel.config(backend::ndrange(threads, batch), backend::ndrange(wg, 1));

    std::ostringstream desc;
    desc << "bluestein_mul{n=" << n << "(" << threads << "), wg=" << wg << ", batch=" << batch << "}";
    return kernel_call(false, desc.str(), kernel);
}

} // namespace fft
} // namespace vex


#endif
