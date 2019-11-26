#ifndef VEXCL_STENCIL_HPP
#define VEXCL_STENCIL_HPP

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
 * \file   stencil.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Stencil convolution.
 */

#include <vector>
#include <map>
#include <sstream>
#include <cassert>
#include <vexcl/vector.hpp>

namespace vex {

template <typename T>
class stencil_base {
    protected:
        template <class Iterator>
        stencil_base(
                const std::vector<backend::command_queue> &queue,
                unsigned width, unsigned center, Iterator begin, Iterator end
                );

        void exchange_halos(const vex::vector<T> &x) const;

        mutable std::vector<backend::command_queue> queue;

        mutable std::vector<T>  hbuf;
        std::vector< backend::device_vector<T> > dbuf;
        std::vector< backend::device_vector<T> > s;

        int lhalo;
        int rhalo;
};

template <typename T> template <class Iterator>
stencil_base<T>::stencil_base(
        const std::vector<backend::command_queue> &q,
        unsigned width, unsigned center, Iterator begin, Iterator end
        )
    : queue(q), hbuf(q.size() * (width - 1)),
      dbuf(q.size()), s(q.size()),
      lhalo(center), rhalo(width - center - 1)
{
    assert(queue.size());
    assert(lhalo >= 0);
    assert(rhalo >= 0);
    assert(width);
    assert(center < width);

    for(unsigned d = 0; d < queue.size(); d++) {
        if (begin != end)
            s[d] = backend::device_vector<T>(queue[d], end - begin, &begin[0], backend::MEM_READ_ONLY);

        // Allocate one element more than needed, to be sure size is nonzero.
        dbuf[d] = backend::device_vector<T>(queue[d], width);
    }

    for(unsigned d = 0; d < queue.size(); d++) queue[d].finish();
}

template <typename T>
void stencil_base<T>::exchange_halos(const vex::vector<T> &x) const {
    int width = lhalo + rhalo;

    if ((queue.size() <= 1) || (width <= 0)) return;

    // Get halos from neighbours.
    for(unsigned d = 0; d < queue.size(); d++) {
        if (!x.part_size(d)) continue;

        // Get halo from left neighbour.
        if (d > 0 && lhalo > 0) {
            size_t end   = x.part_start(d);
            size_t begin = end >= static_cast<unsigned>(lhalo) ?  end - lhalo : 0;
            size_t size  = end - begin;
            x.read_data(begin, size, &hbuf[d * width + lhalo - size], false);
        }

        // Get halo from right neighbour.
        if (d + 1 < queue.size() && rhalo > 0) {
            size_t begin = x.part_start(d + 1);
            size_t end   = std::min(begin + rhalo, x.size());
            size_t size  = end - begin;
            x.read_data(begin, size, &hbuf[d * width + lhalo], false);
        }
    }

    // Wait for the end of transfer.
    for(unsigned d = 0; d < queue.size(); d++) queue[d].finish();

    // Write halos to a local buffer.
    for(unsigned d = 0; d < queue.size(); d++) {
        if (!x.part_size(d)) continue;

        if (d > 0 && lhalo > 0) {
            size_t end   = x.part_start(d);
            size_t begin = end >= static_cast<unsigned>(lhalo) ?  end - lhalo : 0;
            size_t size  = end - begin;
            if (size)
                std::fill(&hbuf[d * width], &hbuf[d * width + lhalo - size], hbuf[d * width + lhalo - size]);
            else
                std::fill(&hbuf[d * width], &hbuf[d * width + lhalo - size], static_cast<T>(x[0]));
        }

        if (d + 1 < queue.size() && rhalo > 0) {
            size_t begin = x.part_start(d + 1);
            size_t end   = std::min(begin + rhalo, x.size());
            size_t size  = end - begin;
            if (size)
                std::fill(&hbuf[d * width + lhalo + size], &hbuf[(d + 1) * width], hbuf[d * width + lhalo + size - 1]);
            else
                std::fill(&hbuf[d * width + lhalo + size], &hbuf[(d + 1) * width], static_cast<T>(x[x.size()-1]));

        }

        if ((d > 0 && lhalo > 0) || (d + 1 < queue.size() && rhalo > 0))
            dbuf[d].write(queue[d], 0, width, &hbuf[d * width]);
    }

    // Wait for the end of transfer.
    for(unsigned d = 0; d < queue.size(); d++) queue[d].finish();
}

/// Stencil.
/**
 * Should be used for stencil convolutions with vex::vectors as in
 \code
 void convolve(
          const vex::stencil<double> &s,
          const vex::vector<double>  &x,
          vex::vector<double> &y)
 {
     y = x * s;
 }
 \endcode
 * Stencil should be small enough to fit into local memory of all compute
 * devices it resides on.
 */
template <typename T>
class stencil : private stencil_base<T> {
    public:
        typedef T value_type;

        /// Costructor.
        /**
         * \param queue  vector of queues. Each queue represents one
         *               compute device.
         * \param st     vector holding stencil values.
         * \param center center of the stencil.
         */
        stencil(const std::vector<backend::command_queue> &queue,
                const std::vector<T> &st, unsigned center
                )
            : stencil_base<T>(queue, static_cast<unsigned>(st.size()), center, st.begin(), st.end()),
              conv(queue.size()), smem(queue.size())
        {
            init(static_cast<unsigned>(st.size()));
        }

        /// Costructor.
        /**
         * \param queue  vector of queues. Each queue represents one
         *               compute device.
         * \param begin  iterator to begin of sequence holding stencil data.
         * \param end    iterator to end of sequence holding stencil data.
         * \param center center of the stencil.
         */
        template <class Iterator>
        stencil(const std::vector<backend::command_queue> &queue,
                Iterator begin, Iterator end, unsigned center
                )
            : stencil_base<T>(queue, static_cast<unsigned>(end - begin), center, begin, end),
              conv(queue.size()), smem(queue.size())
        {
            init(static_cast<unsigned>(end - begin));
        }

#ifndef BOOST_NO_INITIALIZER_LISTS
        /// Costructor.
        /**
         * \param queue  vector of queues. Each queue represents one
         *               compute device.
         * \param list   intializer list holding stencil values.
         * \param center center of the stencil.
         */
        stencil(const std::vector<backend::command_queue> &queue,
                std::initializer_list<T> list, unsigned center
                )
            : stencil_base<T>(queue, list.size(), center, list.begin(), list.end()),
              conv(queue.size()), smem(queue.size())
        {
            init(list.size());
        }
#endif

        /// Convolve stencil with a vector.
        /**
         * y = alpha * conv(x) + y;
         * \param x input vector.
         * \param y output vector.
         * \param alpha Scaling coefficient in front of y.
         * \param append whether to append the result to the output vector
         *               (alternative is to replace the output vector).
         */
        void apply(const vex::vector<T> &x, vex::vector<T> &y,
                T alpha = 1, bool append = false) const;
    private:
        typedef stencil_base<T> Base;

        using Base::queue;
        using Base::hbuf;
        using Base::dbuf;
        using Base::s;
        using Base::lhalo;
        using Base::rhalo;

        mutable std::vector<backend::kernel> conv;
        std::vector<size_t>  smem;

        void init(unsigned width);

        static const backend::kernel& slow_conv(const backend::command_queue &queue);
        static const backend::kernel& fast_conv(const backend::command_queue &queue);
};

namespace detail {

template <typename T>
inline void define_read_x(backend::source_generator &source) {
    source.begin_function<T>("read_x");
    source.begin_function_parameters();
    source.template parameter<ptrdiff_t>("g_id");
    source.template parameter<size_t>("n");
    source.template parameter<char>("has_left");
    source.template parameter<char>("has_right");
    source.template parameter<int>("lhalo");
    source.template parameter<int>("rhalo");
    source.template parameter< global_ptr<const T> >("xloc");
    source.template parameter< global_ptr<const T> >("xrem");
    source.end_function_parameters();

    source.new_line() << "if (g_id >= 0 && g_id < n)";
    source.open("{");
    source.new_line() << "return xloc[g_id];";
    source.close("}");
    source.new_line() << "else if (g_id < 0)";
    source.open("{");
    source.new_line() << "if (has_left) "
        "return (lhalo + g_id >= 0) ? xrem[lhalo + g_id] : 0;";
    source.new_line() << "else return xloc[0];";
    source.close("}");
    source.new_line() << "else";
    source.open("{");
    source.new_line() << "if (has_right) "
        "return (g_id < n + rhalo) ? xrem[lhalo + g_id - n] : 0;";
    source.new_line() << "else return xloc[n - 1];";
    source.close("}");
    source.end_function();
}

}

template <typename T>
const backend::kernel& stencil<T>::slow_conv(const backend::command_queue &queue) {
    using namespace detail;

    static kernel_cache cache;

    auto kernel = cache.find(queue);

    backend::select_context(queue);

    if (kernel == cache.end()) {
        backend::source_generator source(queue);

        define_read_x<T>(source);

        source.begin_kernel("slow_conv");
        source.begin_kernel_parameters();
        source.template parameter<size_t>("n");
        source.template parameter<char>("has_left");
        source.template parameter<char>("has_right");
        source.template parameter<int>("lhalo");
        source.template parameter<int>("rhalo");
        source.template parameter< global_ptr<const T> >("s");
        source.template parameter< global_ptr<const T> >("xloc");
        source.template parameter< global_ptr<const T> >("xrem");
        source.template parameter< global_ptr<T> >("y");
        source.template parameter<T>("alpha");
        source.template parameter<T>("beta");
        source.end_kernel_parameters();

        source.grid_stride_loop().open("{");

        source.new_line() << type_name<T>() << " sum = 0;";
        source.new_line() << "for(int j = -lhalo; j <= rhalo; j++)";
        source.open("{");
        source.new_line() << "sum += s[lhalo + j] * read_x(("
            << type_name<ptrdiff_t>()
            << ")idx + j, n, has_left, has_right, lhalo, rhalo, xloc, xrem);";
        source.close("}");
        source.new_line() << "if (alpha) y[idx] = alpha * y[idx] + beta * sum;";
        source.new_line() << "else y[idx] = beta * sum;";
        source.close("}");
        source.end_kernel();

        kernel = cache.insert(queue, backend::kernel(
                    queue, source.str(), "slow_conv"));
    }

    return kernel->second;
}

template <typename T>
const backend::kernel& stencil<T>::fast_conv(const backend::command_queue &queue) {
    using namespace detail;

    static kernel_cache cache;

    auto kernel = cache.find(queue);

    backend::select_context(queue);

    if (kernel == cache.end()) {
        backend::source_generator source(queue);

        define_read_x<T>(source);

        source.begin_kernel("fast_conv");
        source.begin_kernel_parameters();
        source.template parameter<size_t>("n");
        source.template parameter<char>("has_left");
        source.template parameter<char>("has_right");
        source.template parameter<int>("lhalo");
        source.template parameter<int>("rhalo");
        source.template parameter< global_ptr<const T> >("s");
        source.template parameter< global_ptr<const T> >("xloc");
        source.template parameter< global_ptr<const T> >("xrem");
        source.template parameter< global_ptr<T> >("y");
        source.template parameter<T>("alpha");
        source.template parameter<T>("beta");
        source.template smem_parameter< T >();
        source.end_kernel_parameters();

        source.smem_declaration<T>();
        source.new_line() << type_name< shared_ptr<T> >() << " S = smem;";
        source.new_line() << type_name< shared_ptr<T> >() << " X = smem + lhalo + rhalo + 1;";

        source.new_line() << "size_t grid_size = " << source.global_size(0) << ";";
        source.new_line() << "int l_id = " << source.local_id(0) << ";";
        source.new_line() << "int block_size = " << source.local_size(0) << ";";
        source.new_line() << "for(int i = l_id; i < rhalo + lhalo + 1; i += block_size) S[i] = s[i];";
        source.new_line() << "for(long g_id = " << source.global_id(0) << ", pos = 0; pos < n; g_id += grid_size, pos += grid_size)";
        source.open("{");
        source.new_line() << "for(int i = l_id, j = g_id - lhalo; i < block_size + lhalo + rhalo; i += block_size, j += block_size)";
        source.open("{");
        source.new_line() << "X[i] = read_x(j, n, has_left, has_right, lhalo, rhalo, xloc, xrem);";
        source.close("}");
        source.new_line().barrier();
        source.new_line() << "if (g_id < n)";
        source.open("{");
        source.new_line() << type_name<T>() << " sum = 0;";
        source.new_line() << "for(int j = -lhalo; j <= rhalo; j++)";
        source.open("{");
        source.new_line() << "sum += S[lhalo + j] * X[lhalo + l_id + j];";
        source.close("}");
        source.new_line() << "if (alpha) "
            "y[g_id] = alpha * y[g_id] + beta * sum;";
        source.new_line() << "else y[g_id] = beta * sum;";
        source.close("}");
        source.new_line().barrier();
        source.close("}");
        source.end_kernel();

        kernel = cache.insert(queue, backend::kernel(
                    queue, source.str(), "fast_conv"));
    }

    return kernel->second;
}

template <typename T>
void stencil<T>::init(unsigned width) {
    for (unsigned d = 0; d < queue.size(); d++) {
        auto   fast_krn = fast_conv(queue[d]);
        size_t max_smem = fast_krn.max_shared_memory_per_block(queue[d]);

        auto smem_required = [width](size_t wgs) { return sizeof(T) * (2 * width + wgs - 1); };

        if (backend::is_cpu(queue[d]) || max_smem < smem_required(64)) {
            conv[d]  = slow_conv(queue[d]);
            conv[d].config(queue[d], [](size_t){ return 0; });
            smem[d] = 0;
        } else {
            conv[d] = fast_krn;
            conv[d].config(queue[d], smem_required);
            smem[d] = smem_required(conv[d].workgroup_size());
        }
    }
}

template <typename T>
void stencil<T>::apply(const vex::vector<T> &x, vex::vector<T> &y,
        T alpha, bool append) const
{
    Base::exchange_halos(x);

    T beta = static_cast<T>(append ? 1 : 0);

    for(unsigned d = 0; d < queue.size(); d++) {
        if (size_t psize = x.part_size(d)) {
            char has_left  = d > 0;
            char has_right = d + 1 < queue.size();

            conv[d].push_arg(psize);
            conv[d].push_arg(has_left);
            conv[d].push_arg(has_right);
            conv[d].push_arg(lhalo);
            conv[d].push_arg(rhalo);
            conv[d].push_arg(s[d]);
            conv[d].push_arg(x(d));
            conv[d].push_arg(dbuf[d]);
            conv[d].push_arg(y(d));
            conv[d].push_arg(beta);
            conv[d].push_arg(alpha);

            if (smem[d]) conv[d].set_smem([&](size_t){ return smem[d]; });

            conv[d](queue[d]);
        }
    }
}

/// Convolve the stencil with the vector.
template <typename T>
additive_operator< stencil<T>, vector<T> >
operator*( const stencil<T> &s, const vector<T> &x ) {
    return additive_operator< stencil<T>, vector<T> >(s, x);
}

/// Convolve the stencil with the vector.
template <typename T>
additive_operator< stencil<T>, vector<T> >
operator*(const vector<T> &x, const stencil<T> &s) {
    return additive_operator< stencil<T>, vector<T> >(s, x);
}

#ifdef VEXCL_MULTIVECTOR_HPP

/// Convolve the stencil with the multivector.
template <typename T, size_t N>
multiadditive_operator< stencil<T>, multivector<T, N> >
operator*( const stencil<T> &s, const multivector<T, N> &x ) {
    return multiadditive_operator< stencil<T>, multivector<T, N> >(s, x);
}

/// Convolve the stencil with the multivector.
template <typename T, size_t N>
multiadditive_operator< stencil<T>, multivector<T, N> >
operator*( const multivector<T, N> &x, const stencil<T> &s ) {
    return multiadditive_operator< stencil<T>, multivector<T, N> >(s, x);
}

#endif

/// User-defined stencil operator
/**
 * Is used to define custom stencil operator. For example, to implement the
 * following nonlinear operator:
 \code
 y[i] = x[i] + pow3(x[i-1] + x[i+1]);
 \endcode
 * one has to write:
 \code
 extern const char pow3_oper_body[] = "return X[0] + pow(X[-1] + X[1], 3);";
 StencilOperator<double, 3, 1, pow3_oper_body> pow3_oper(ctx);

 y = pow3_oper(x);
 \endcode
 */
template <typename T, unsigned width, unsigned center, class Impl>
class StencilOperator : private stencil_base<T> {
    public:
        typedef T value_type;

        StencilOperator(const std::vector<backend::command_queue> &queue);

        additive_operator< StencilOperator, vector<T> >
        operator()(const vector<T> &x) const {
            return additive_operator< StencilOperator, vector<T> >(*this, x);
        }

#ifdef VEXCL_MULTIVECTOR_HPP
        template <size_t N>
        multiadditive_operator< StencilOperator, multivector<T, N> >
        operator()(const multivector<T, N> &x) const {
            return multiadditive_operator< StencilOperator, multivector<T, N> >(*this, x);
        }
#endif

        void apply(const vex::vector<T> &x, vex::vector<T> &y,
                T alpha = 1, bool append = true) const;
    private:
        typedef stencil_base<T> Base;

        using Base::queue;
        using Base::hbuf;
        using Base::dbuf;
        using Base::lhalo;
        using Base::rhalo;
};

template <typename T, unsigned width, unsigned center, class Impl>
StencilOperator<T, width, center, Impl>::StencilOperator(
        const std::vector<backend::command_queue> &queue)
    : Base(queue, width, center, static_cast<T*>(0), static_cast<T*>(0))
{ }

template <typename T, unsigned width, unsigned center, class Impl>
void StencilOperator<T, width, center, Impl>::apply(
        const vex::vector<T> &x, vex::vector<T> &y, T alpha, bool append) const
{
    using namespace detail;

    T beta = append ? 1 : 0;

    static kernel_cache cache;
    static std::map<backend::context_id, size_t> lmem;

    Base::exchange_halos(x);

    for(unsigned d = 0; d < queue.size(); d++) {
        backend::select_context(queue[d]);

        auto key    = backend::get_context_id(queue[d]);
        auto kernel = cache.find(queue[d]);

        if (kernel == cache.end()) {
            backend::source_generator source(queue[d]);

            define_read_x<T>(source);

            source.begin_function<T>("stencil_oper");
            source.begin_function_parameters();
            source.template parameter< shared_ptr<const T> >("X");
            source.end_function_parameters();
            source.new_line() << Impl::body();
            source.end_function();

            source.begin_kernel("convolve");
            source.begin_kernel_parameters();
            source.template parameter<size_t>("n");
            source.template parameter<char>("has_left");
            source.template parameter<char>("has_right");
            source.template parameter<int>("lhalo");
            source.template parameter<int>("rhalo");
            source.template parameter< global_ptr<const T> >("xloc");
            source.template parameter< global_ptr<const T> >("xrem");
            source.template parameter< global_ptr<T> >("y");
            source.template parameter<T>("alpha");
            source.template parameter<T>("beta");
            source.template smem_parameter<T>();
            source.end_kernel_parameters();

            source.smem_declaration<T>();
            source.new_line() << type_name< shared_ptr<T> >() << " X = smem;";

            source.new_line() << "size_t grid_size = " << source.global_size(0) << ";";
            source.new_line() << "int l_id = " << source.local_id(0) << ";";
            source.new_line() << "int block_size = " << source.local_size(0) << ";";
            source.new_line() << "for(long g_id = " << source.global_id(0)
                << ", pos = 0; pos < n; g_id += grid_size, pos += grid_size)";
            source.open("{");
            source.new_line() << "for(int i = l_id, j = g_id - lhalo; i < block_size + lhalo + rhalo; i += block_size, j += block_size)";
            source.open("{");
            source.new_line() << "X[i] = read_x(j, n, has_left, has_right, lhalo, rhalo, xloc, xrem);";
            source.close("}");
            source.new_line().barrier();
            source.new_line() << "if (g_id < n)";
            source.open("{");
            source.new_line() << type_name<T>() << " sum = stencil_oper(X + lhalo + l_id);";
            source.new_line() << "if (alpha) y[g_id] = alpha * y[g_id] + beta * sum;";
            source.new_line() << "else y[g_id] = beta * sum;";
            source.close("}");
            source.new_line().barrier();
            source.close("}");
            source.end_kernel();

            kernel = cache.insert(queue[d], backend::kernel(
                        queue[d], source.str(), "convolve",
                        [](size_t wgs) { return (width + wgs - 1) * sizeof(T); }
                        ));

            lmem[key] = sizeof(T) * (kernel->second.workgroup_size() + width - 1);
        }

        if (size_t psize = x.part_size(d)) {
            char has_left  = d > 0;
            char has_right = d + 1 < queue.size();

            kernel->second.push_arg(psize);
            kernel->second.push_arg(has_left);
            kernel->second.push_arg(has_right);
            kernel->second.push_arg(lhalo);
            kernel->second.push_arg(rhalo);
            kernel->second.push_arg(x(d));
            kernel->second.push_arg(dbuf[d]);
            kernel->second.push_arg(y(d));
            kernel->second.push_arg(beta);
            kernel->second.push_arg(alpha);

            size_t smem_bytes = lmem[key];
            kernel->second.set_smem([smem_bytes](size_t){ return smem_bytes; });

            kernel->second(queue[d]);
        }
    }
}

/// Macro to declare a user-defined stencil operator type.
/**
 \code
 VEX_STENCIL_OPERATOR_TYPE(pow3_oper_t, double, 3, 1, "return X[0] + pow(X[-1] + X[1], 3.0);");
 pow3_oper_t pow3_oper(ctx);
 output = pow3_oper(input);
 \endcode
 *
 * \note Should be used in case same operator is used in several places (to
 * save on OpenCL kernel recompilations). Otherwise VEX_STENCIL_OPERATOR should
 * be used locally.
 */
#define VEX_STENCIL_OPERATOR_TYPE(name, type, width, center, body_str) \
    struct name : vex::StencilOperator<type, width, center, name> { \
        name(const std::vector<vex::backend::command_queue> &q) : vex::StencilOperator<type, width, center, name>(q) {} \
        static std::string body() { return body_str; } \
    }

/// Macro to declare a user-defined stencil operator.
/**
 \code
 VEX_STENCIL_OPERATOR(pow3_oper, double, 3, 1, "return X[0] + pow(X[-1] + X[1], 3.0);", queue);
 output = pow3_oper(input);
 \endcode
 */
#define VEX_STENCIL_OPERATOR(name, type, width, center, body, queue)           \
  VEX_STENCIL_OPERATOR_TYPE(                                                   \
          stencil_operator_##name##_t, type, width, center,  body)             \
    const name(queue)

} // namespace vex

#endif
