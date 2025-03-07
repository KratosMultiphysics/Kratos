/*
    enoki/autodiff.h -- Reverse mode automatic differentiation

    Enoki is a C++ template library that enables transparent vectorization
    of numerical kernels using SIMD instruction sets available on current
    processor architectures.

    Copyrighe (c) 2019 Wenzel Jakob <wenzel.jakob@epfl.ch>

    All rights reserved. Use of this source code is governed by a BSD-style
    license that can be found in the LICENSE file.
*/

#pragma once

#include <enoki/array.h>
#include <vector>

#define ENOKI_AUTODIFF_H 1

NAMESPACE_BEGIN(enoki)

template <typename Type> struct Tape {
private:
    template <typename T> friend struct DiffArray;

    struct Detail;
    struct Node;
    struct Edge;
    struct Special;
    struct SimplificationLock;

    using Index = uint32_t;
    using Mask = mask_t<Type>;
    using Int64 = int64_array_t<Type>;

    Tape();

    // -----------------------------------------------------------------------
    //! @{ \name Append unary/binary/ternary operations to the tape
    // -----------------------------------------------------------------------

    Index append(const char *label, size_t size, Index i1, const Type &w1);

    Index append(const char *label, size_t size, Index i1, Index i2,
                 const Type &w1, const Type &w2);

    Index append(const char *label, size_t size, Index i1, Index i2, Index i3,
                 const Type &w1, const Type &w2, const Type &w3);

    Index append_psum(Index i);
    Index append_reverse(Index i);

    Index append_gather(const Int64 &offset, const Mask &mask);

    void append_scatter(Index index, const Int64 &offset, const Mask &mask,
                        bool scatter_add);

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Append nodes and edges to the tape
    // -----------------------------------------------------------------------

    Index append_node(size_t size, const char *label);
    Index append_leaf(size_t size);
    void append_edge(Index src, Index dst, const Type &weight);
    void append_edge_prod(Index src, Index dst, const Type &weight1,
                          const Type &weight2);

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Reference counting
    // -----------------------------------------------------------------------

    void dec_ref_ext(Index index);
    void inc_ref_ext(Index index);
    void dec_ref_int(Index index, Index from);
    void inc_ref_int(Index index, Index from);
    void free_node(Index index);

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Other operations
    // -----------------------------------------------------------------------

    void set_scatter_gather_operand(Index *index, size_t size, bool permute);
    void push_prefix(const char *);
    void pop_prefix();
    void backward(bool free_graph);
    void forward(bool free_graph);
    void backward(Index index, bool free_graph);
    void forward(Index index, bool free_graph);
    void set_gradient(Index index, const Type &value,
                      bool backward = true);
    void set_label(Index index, const char *name);
    const Type &gradient(Index index);
    std::string graphviz(const std::vector<Index> &indices);
    /// Current log level (0 == none, 1 == minimal, 2 == moderate, 3 == high, 4 == everything)
    void set_log_level(uint32_t);
    uint32_t log_level() const;
    void set_graph_simplification(bool);
    void simplify_graph();
    std::string whos() const;
    static void cuda_callback(void*);

    //! @}
    // -----------------------------------------------------------------------

    static Tape* get() ENOKI_PURE;

public:
    ~Tape();

private:

    static std::unique_ptr<Tape> s_tape;
    Detail *d;
};

template <typename Type>
struct DiffArray : ArrayBase<value_t<Type>, DiffArray<Type>> {
public:
    using Base = enoki::ArrayBase<value_t<Type>, DiffArray<Type>>;
    using typename Base::Scalar;
    using Tape = enoki::Tape<Type>;
    using Index = uint32_t;

    using UnderlyingType = Type;
    using ArrayType = DiffArray;
    using MaskType = DiffArray<mask_t<Type>>;

    static constexpr size_t Size = is_scalar_v<Type> ? 1 : array_size_v<Type>;
    static constexpr size_t Depth = is_scalar_v<Type> ? 1 : array_depth_v<Type>;
    static constexpr bool IsMask = is_mask_v<Type>;
    static constexpr bool IsCUDA = is_cuda_array_v<Type>;
    static constexpr bool IsDiff = true;
    static constexpr bool Enabled =
        std::is_floating_point_v<scalar_t<Type>> && !is_mask_v<Type>;

    template <typename T>
    using ReplaceValue = DiffArray<replace_scalar_t<Type, T, false>>;

    static_assert(array_depth_v<Type> <= 1,
                  "DiffArray requires a scalar or (non-nested) static or "
                  "dynamic Enoki array as template parameter.");

    // -----------------------------------------------------------------------
    //! @{ \name Constructors / destructors
    // -----------------------------------------------------------------------

    DiffArray() = default;

    ~DiffArray() {
        if constexpr (Enabled)
            tape()->dec_ref_ext(m_index);
    }

    DiffArray(const DiffArray &a) : m_value(a.m_value), m_index(a.m_index) {
        if constexpr (Enabled)
            tape()->inc_ref_ext(m_index);
    }

    DiffArray(DiffArray &&a) : m_value(std::move(a.m_value)) {
        if constexpr (Enabled) {
            m_index = a.m_index;
            a.m_index = 0;
        }
    }

    template <typename T>
    DiffArray(const DiffArray<T> &v, detail::reinterpret_flag) :
        m_value(v.value_(), detail::reinterpret_flag()) { /* no derivatives */ }

    template <typename Type2, enable_if_t<!std::is_same_v<Type, Type2>> = 0>
    DiffArray(const DiffArray<Type2> &a) : m_value(a.value_()) { }

    template <typename Type2, enable_if_t<!std::is_same_v<Type, Type2>> = 0>
    DiffArray(DiffArray<Type2> &&a) : m_value(std::move(a.value_())) { }

    DiffArray(Type &&value) : m_value(std::move(value)) { }

    template <typename... Args,
             enable_if_t<sizeof...(Args) != 0 && std::conjunction_v<
                  std::negation<is_diff_array<Args>>...>> = 0>
    DiffArray(Args&&... args) : m_value(std::forward<Args>(args)...) { }

    DiffArray &operator=(const DiffArray &a) {
        m_value = a.m_value;
        if constexpr (Enabled) {
            auto t = tape();
            t->inc_ref_ext(a.m_index);
            t->dec_ref_ext(m_index);
            m_index = a.m_index;
        }
        return *this;
    }

    DiffArray &operator=(DiffArray &&a) {
        m_value = std::move(a.m_value);
        if constexpr (Enabled)
            std::swap(m_index, a.m_index);
        return *this;
    }


    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Vertical operations
    // -----------------------------------------------------------------------

    DiffArray add_(const DiffArray &a) const {
        if constexpr (is_mask_v<Type> || std::is_pointer_v<Scalar>) {
            fail_unsupported("add_");
        } else {
            Index index_new = 0;
            Type result = m_value + a.m_value;
            if constexpr (Enabled)
                index_new = tape()->append("add", slices(result), m_index,
                                           a.m_index, 1.f, 1.f);
            return DiffArray::create(index_new, std::move(result));
        }
    }

    DiffArray sub_(const DiffArray &a) const {
        if constexpr (is_mask_v<Type> || std::is_pointer_v<Scalar>) {
            fail_unsupported("sub_");
        } else {
            Index index_new = 0;
            Type result = m_value - a.m_value;
            if constexpr (Enabled)
                index_new = tape()->append("sub", slices(result), m_index,
                                           a.m_index, 1.f, -1.f);
            return DiffArray::create(index_new, std::move(result));
        }
    }

    DiffArray mul_(const DiffArray &a) const {
        if constexpr (is_mask_v<Type> || std::is_pointer_v<Scalar>) {
            fail_unsupported("mul_");
        } else {
            Index index_new = 0;
            Type result = m_value * a.m_value;
            if constexpr (Enabled)
                index_new = tape()->append("mul", slices(result), m_index,
                                           a.m_index, a.m_value, m_value);
            return DiffArray::create(index_new, std::move(result));
        }
    }

    DiffArray div_(const DiffArray &a) const {
        if constexpr (is_mask_v<Type> || std::is_pointer_v<Scalar>) {
            fail_unsupported("div_");
        } else {
            Index index_new = 0;
            Type result = m_value / a.m_value;
            if constexpr (Enabled) {
                Type rcp_a = rcp(a.m_value);
                index_new = tape()->append("div", slices(result),
                                           m_index, a.m_index, rcp_a,
                                           -m_value * sqr(rcp_a));
            }
            return DiffArray::create(index_new, std::move(result));
        }
    }

    DiffArray fmadd_(const DiffArray &a, const DiffArray &b) const {
        if constexpr (is_mask_v<Type>) {
            fail_unsupported("fmadd_");
        } else {
            Index index_new = 0;
            Type result = fmadd(m_value, a.m_value, b.m_value);
            if constexpr (Enabled)
                index_new = tape()->append("fmadd", slices(result),
                                           m_index, a.m_index, b.m_index,
                                           a.m_value, m_value, 1);
            return DiffArray::create(index_new, std::move(result));
        }
    }

    DiffArray fmsub_(const DiffArray &a, const DiffArray &b) const {
        if constexpr (is_mask_v<Type>) {
            fail_unsupported("fmsub_");
        } else {
            Type result = fmsub(m_value, a.m_value, b.m_value);
            Index index_new = 0;
            if constexpr (Enabled)
                index_new = tape()->append("fmsub", slices(result),
                                           m_index, a.m_index, b.m_index,
                                           a.m_value, m_value, -1);
            return DiffArray::create(index_new, std::move(result));
        }
    }

    DiffArray fnmadd_(const DiffArray &a, const DiffArray &b) const {
        if constexpr (is_mask_v<Type>) {
            fail_unsupported("fnmadd_");
        } else {
            Type result = fnmadd(m_value, a.m_value, b.m_value);
            Index index_new = 0;
            if constexpr (Enabled)
                index_new = tape()->append("fnmadd", slices(result),
                                           m_index, a.m_index, b.m_index,
                                           -a.m_value, -m_value, 1);
            return DiffArray::create(index_new, std::move(result));
        }
    }

    DiffArray fnmsub_(const DiffArray &a, const DiffArray &b) const {
        if constexpr (is_mask_v<Type>) {
            fail_unsupported("fnmsub_");
        } else {
            Index index_new = 0;
            Type result = fnmsub(m_value, a.m_value, b.m_value);
            if constexpr (Enabled)
                index_new = tape()->append("fnmsub", slices(result),
                                           m_index, a.m_index, b.m_index,
                                           -a.m_value, -m_value, -1);
            return DiffArray::create(index_new, std::move(result));
        }
    }

    DiffArray neg_() const {
        if constexpr (is_mask_v<Type> || std::is_pointer_v<Scalar>) {
            fail_unsupported("neg_");
        } else {
            Index index_new = 0;
            if constexpr (Enabled)
                index_new = tape()->append("neg", slices(m_value), m_index, -1.f);
            return DiffArray::create(index_new, -m_value);
        }
    }

    DiffArray abs_() const {
        if constexpr (is_mask_v<Type> || std::is_pointer_v<Scalar>) {
            fail_unsupported("abs_");
        } else {
            Index index_new = 0;
            if constexpr (Enabled)
                index_new = tape()->append("abs", slices(m_value), m_index,
                                           sign(m_value));
            return DiffArray::create(index_new, abs(m_value));
        }
    }

    DiffArray sqrt_() const {
        if constexpr (is_mask_v<Type> || !std::is_floating_point_v<Scalar>) {
            fail_unsupported("sqrt_");
        } else {
            Index index_new = 0;
            Type result = sqrt(m_value);
            if constexpr (Enabled)
                index_new = tape()->append("sqrt", slices(result), m_index,
                                           .5f / result);
            return DiffArray::create(index_new, std::move(result));
        }
    }

    DiffArray cbrt_() const {
        if constexpr (is_mask_v<Type> || !std::is_floating_point_v<Scalar>) {
            fail_unsupported("cbrt_");
        } else {
            Index index_new = 0;
            Type result = cbrt(m_value);
            if constexpr (Enabled)
                index_new = tape()->append("cbrt", slices(result), m_index,
                                           1.f / (3 * sqr(result)));
            return DiffArray::create(index_new, std::move(result));
        }
    }

    DiffArray rcp_() const {
        if constexpr (is_mask_v<Type> || !std::is_floating_point_v<Scalar>) {
            fail_unsupported("rcp_");
        } else {
            Index index_new = 0;
            Type result = rcp(m_value);
            if constexpr (Enabled)
                index_new = tape()->append("rcp", slices(result), m_index,
                                           -sqr(result));
            return DiffArray::create(index_new, std::move(result));
        }
    }

    DiffArray rsqrt_() const {
        if constexpr (is_mask_v<Type> || !std::is_floating_point_v<Scalar>) {
            fail_unsupported("rsqrt_");
        } else {
            Index index_new = 0;
            Type result = rsqrt(m_value);
            if constexpr (Enabled) {
                Type rsqrt_2 = sqr(result), rsqrt_3 = result * rsqrt_2;
                index_new = tape()->append("rsqrt", slices(result), m_index,
                                           -.5f * rsqrt_3);
            }
            return DiffArray::create(index_new, std::move(result));
        }
    }

    DiffArray min_(const DiffArray &a) const {
        if constexpr (is_mask_v<Type>) {
            fail_unsupported("min_");
        } else {
            Index index_new = 0;
            Type result = min(m_value, a.m_value);
            if constexpr (Enabled) {
                mask_t<Type> m = m_value < a.m_value;
                index_new = tape()->append("min", slices(result),
                                           m_index, a.m_index,
                                           select(m, Type(1), Type(0)),
                                           select(m, Type(0), Type(1)));
            }
            return DiffArray::create(index_new, std::move(result));
        }
    }

    DiffArray max_(const DiffArray &a) const {
        if constexpr (is_mask_v<Type>) {
            fail_unsupported("max_");
        } else {
            Index index_new = 0;
            Type result = max(m_value, a.m_value);
            if constexpr (Enabled) {
                mask_t<Type> m = m_value > a.m_value;
                index_new = tape()->append("max", slices(result),
                                           m_index, a.m_index,
                                           select(m, Type(1), Type(0)),
                                           select(m, Type(0), Type(1)));
            }
            return DiffArray::create(index_new, std::move(result));
        }
    }

    static DiffArray select_(const DiffArray<mask_t<Type>> &m,
                             const DiffArray &t,
                             const DiffArray &f) {
        Index index_new = 0;
        Type result = select(m.value_(), t.m_value, f.m_value);
        if constexpr (Enabled) {
            index_new =
                tape()->append("select", slices(result), t.m_index, f.m_index,
                               select(m.value_(), Type(1), Type(0)),
                               select(m.value_(), Type(0), Type(1)));
        }
        return DiffArray::create(index_new, std::move(result));
    }

    DiffArray floor_() const {
        if constexpr (is_mask_v<Type> || !std::is_floating_point_v<Scalar>)
            fail_unsupported("floor_");
        else
            return DiffArray::create(0, floor(m_value));
    }

    DiffArray ceil_() const {
        if constexpr (is_mask_v<Type> || !std::is_floating_point_v<Scalar>)
            fail_unsupported("ceil_");
        else
            return DiffArray::create(0, ceil(m_value));
    }

    DiffArray trunc_() const {
        if constexpr (is_mask_v<Type> || !std::is_floating_point_v<Scalar>)
            fail_unsupported("trunc_");
        else
            return DiffArray::create(0, trunc(m_value));
    }

    DiffArray round_() const {
        if constexpr (is_mask_v<Type> || !std::is_floating_point_v<Scalar>)
            fail_unsupported("round_");
        else
            return DiffArray::create(0, round(m_value));
    }

    template <typename T> T ceil2int_() const {
        return T(ceil2int<typename T::UnderlyingType>(m_value));
    }

    template <typename T> T floor2int_() const {
        return T(floor2int<typename T::UnderlyingType>(m_value));
    }

    DiffArray sin_() const {
        if constexpr (is_mask_v<Type> || !std::is_floating_point_v<Scalar>) {
            fail_unsupported("sin_");
        } else {
            Index index_new = 0;
            auto [s, c] = sincos(m_value);
            if constexpr (Enabled)
                index_new = tape()->append("sin", slices(m_value), m_index, c);
            return DiffArray::create(index_new, std::move(s));
        }
    }

    DiffArray cos_() const {
        if constexpr (is_mask_v<Type> || !std::is_floating_point_v<Scalar>) {
            fail_unsupported("cos_");
        } else {
            Index index_new = 0;
            auto [s, c] = sincos(m_value);
            if constexpr (Enabled)
                index_new = tape()->append("cos", slices(m_value), m_index, -s);
            return DiffArray::create(index_new, std::move(c));
        }
    }

    std::pair<DiffArray, DiffArray> sincos_() const {
        if constexpr (is_mask_v<Type> || !std::is_floating_point_v<Scalar>) {
            fail_unsupported("sincos_");
        } else {
            Index index_new_s = 0, index_new_c = 0;
            auto [s, c] = sincos(m_value);
            if constexpr (Enabled) {
                index_new_s = tape()->append("sin", slices(m_value), m_index,  c);
                index_new_c = tape()->append("cos", slices(m_value), m_index, -s);
            }
            return {
                DiffArray::create(index_new_s, std::move(s)),
                DiffArray::create(index_new_c, std::move(c))
            };
        }
    }

    DiffArray tan_() const {
        if constexpr (is_mask_v<Type> || !std::is_floating_point_v<Scalar>) {
            fail_unsupported("tan_");
        } else {
            Index index_new = 0;
            if constexpr (Enabled)
                index_new = tape()->append("tan", slices(m_value), m_index,
                                           sqr(sec(m_value)));
            return DiffArray::create(index_new, tan(m_value));
        }
    }

    DiffArray csc_() const {
        if constexpr (is_mask_v<Type> || !std::is_floating_point_v<Scalar>) {
            fail_unsupported("csc_");
        } else {
            Index index_new = 0;
            Type csc_value = csc(m_value);
            if constexpr (Enabled)
                index_new = tape()->append("csc", slices(m_value), m_index,
                                           -csc_value * cot(m_value));
            return DiffArray::create(index_new, std::move(csc_value));
        }
    }

    DiffArray sec_() const {
        if constexpr (is_mask_v<Type> || !std::is_floating_point_v<Scalar>) {
            fail_unsupported("sec_");
        } else {
            Index index_new = 0;
            Type sec_value = sec(m_value);
            if constexpr (Enabled)
                index_new = tape()->append("sec", slices(m_value), m_index,
                                           sec_value * tan(m_value));
            return DiffArray::create(index_new, std::move(sec_value));
        }
    }

    DiffArray cot_() const {
        if constexpr (is_mask_v<Type> || !std::is_floating_point_v<Scalar>) {
            fail_unsupported("cot_");
        } else {
            Index index_new = 0;
            if constexpr (Enabled)
                index_new = tape()->append("cot", slices(m_value), m_index,
                                           -sqr(csc(m_value)));
            return DiffArray::create(index_new, cot(m_value));
        }
    }

    DiffArray asin_() const {
        if constexpr (is_mask_v<Type> || !std::is_floating_point_v<Scalar>) {
            fail_unsupported("asin_");
        } else {
            Index index_new = 0;
            if constexpr (Enabled)
                index_new = tape()->append("asin", slices(m_value), m_index,
                                           rsqrt(1 - sqr(m_value)));
            return DiffArray::create(index_new, asin(m_value));
        }
    }

    DiffArray acos_() const {
        if constexpr (is_mask_v<Type> || !std::is_floating_point_v<Scalar>) {
            fail_unsupported("acos_");
        } else {
            Index index_new = 0;
            if constexpr (Enabled)
                index_new = tape()->append("acos", slices(m_value), m_index,
                                           -rsqrt(1 - sqr(m_value)));
            return DiffArray::create(index_new, acos(m_value));
        }
    }

    DiffArray atan_() const {
        if constexpr (is_mask_v<Type> || !std::is_floating_point_v<Scalar>) {
            fail_unsupported("atan_");
        } else {
            Index index_new = 0;
            if constexpr (Enabled)
                index_new = tape()->append("atan", slices(m_value), m_index,
                                           rcp(1 + sqr(m_value)));
            return DiffArray::create(index_new, atan(m_value));
        }
    }

    DiffArray atan2_(const DiffArray &x) const {
        if constexpr (is_mask_v<Type> || !std::is_floating_point_v<Scalar>) {
            fail_unsupported("atan2_");
        } else {
            Index index_new = 0;

            if constexpr (Enabled) {
                Type il2 = rcp(sqr(m_value) + sqr(x.m_value));
                index_new = tape()->append("atan2", slices(il2),
                                           m_index, x.m_index,
                                           il2 * x.m_value, -il2 * m_value);
            }

            return DiffArray::create(index_new, atan2(m_value, x.m_value));
        }
    }

    DiffArray sinh_() const {
        if constexpr (is_mask_v<Type> || !std::is_floating_point_v<Scalar>) {
            fail_unsupported("sinh_");
        } else {
            Index index_new = 0;
            auto [s, c] = sincosh(m_value);
            if constexpr (Enabled)
                index_new = tape()->append("sinh", slices(m_value), m_index, c);
            return DiffArray::create(index_new, std::move(s));
        }
    }

    DiffArray cosh_() const {
        if constexpr (is_mask_v<Type> || !std::is_floating_point_v<Scalar>) {
            fail_unsupported("cosh_");
        } else {
            Index index_new = 0;
            auto [s, c] = sincosh(m_value);
            if constexpr (Enabled)
                index_new = tape()->append("cosh", slices(m_value), m_index, s);
            return DiffArray::create(index_new, std::move(c));
        }
    }

    DiffArray csch_() const {
        if constexpr (is_mask_v<Type> || !std::is_floating_point_v<Scalar>) {
            fail_unsupported("csch_");
        } else {
            Index index_new = 0;
            Type result = csch(m_value);
            if constexpr (Enabled)
                index_new = tape()->append("csch", slices(m_value), m_index,
                                           -result * coth(m_value));
            return DiffArray::create(index_new, std::move(result));
        }
    }

    DiffArray sech_() const {
        if constexpr (is_mask_v<Type> || !std::is_floating_point_v<Scalar>) {
            fail_unsupported("sech_");
        } else {
            Index index_new = 0;
            Type result = sech(m_value);
            if constexpr (Enabled)
                index_new = tape()->append("sech", slices(m_value), m_index,
                                           -result * tanh(m_value));
            return DiffArray::create(index_new, std::move(result));
        }
    }

    DiffArray tanh_() const {
        if constexpr (is_mask_v<Type> || !std::is_floating_point_v<Scalar>) {
            fail_unsupported("tanh_");
        } else {
            Index index_new = 0;
            Type result = tanh(m_value);
            if constexpr (Enabled)
                index_new = tape()->append("index", slices(m_value), m_index,
                                           sqr(sech(m_value)));
            return DiffArray::create(index_new, std::move(result));
        }
    }

    DiffArray asinh_() const {
        if constexpr (is_mask_v<Type> || !std::is_floating_point_v<Scalar>) {
            fail_unsupported("asinh_");
        } else {
            Index index_new = 0;
            if constexpr (Enabled)
                index_new = tape()->append("asinh", slices(m_value), m_index,
                                           rsqrt((Scalar) 1 + sqr(m_value)));
            return DiffArray::create(index_new, asinh(m_value));
        }
    }

    DiffArray acosh_() const {
        if constexpr (is_mask_v<Type> || !std::is_floating_point_v<Scalar>) {
            fail_unsupported("acosh_");
        } else {
            Index index_new = 0;
            if constexpr (Enabled)
                index_new = tape()->append("acosh", slices(m_value), m_index,
                                           rsqrt(sqr(m_value) - (Scalar) 1));
            return DiffArray::create(index_new, acosh(m_value));
        }
    }

    DiffArray atanh_() const {
        if constexpr (is_mask_v<Type> || !std::is_floating_point_v<Scalar>) {
            fail_unsupported("atanh_");
        } else {
            Index index_new = 0;
            if constexpr (Enabled)
                index_new = tape()->append("atanh", slices(m_value), m_index,
                                           rcp((Scalar) 1 - sqr(m_value)));
            return DiffArray::create(index_new, atanh(m_value));
        }
    }

    DiffArray exp_() const {
        if constexpr (is_mask_v<Type> || !std::is_floating_point_v<Scalar>) {
            fail_unsupported("exp_");
        } else {
            Index index_new = 0;
            Type result = exp(m_value);
            if constexpr (Enabled)
                index_new = tape()->append("exp", slices(m_value),
                                           m_index, result);
            return DiffArray::create(index_new, std::move(result));
        }
    }

    DiffArray log_() const {
        if constexpr (is_mask_v<Type> || !std::is_floating_point_v<Scalar>) {
            fail_unsupported("log_");
        } else {
            Index index_new = 0;
            if constexpr (Enabled)
                index_new = tape()->append("log", slices(m_value), m_index,
                                           rcp(m_value));
            return DiffArray::create(index_new, log(m_value));
        }
    }

    DiffArray or_(const DiffArray &m) const {
        if constexpr (!is_mask_v<Type> && !std::is_integral_v<Scalar>)
            fail_unsupported("or_");
        else
            return DiffArray::create(0, m_value | m.value_());
    }

    template <typename Mask> DiffArray or_(const Mask &m) const {
        Index index_new = 0;
        if constexpr (Enabled && is_mask_v<Mask>)
            index_new = tape()->append("or", slices(m_value), m_index, 1);
        return DiffArray::create(index_new, m_value | m.value_());
    }

    DiffArray and_(const DiffArray &m) const {
        if constexpr (!is_mask_v<Type> && !std::is_integral_v<Scalar>)
            fail_unsupported("and_");
        else
            return DiffArray::create(0, m_value & m.value_());
    }

    template <typename Mask>
    DiffArray and_(const Mask &m) const {
        Index index_new = 0;
        if constexpr (Enabled && is_mask_v<Mask>)
            index_new = tape()->append("and", slices(m_value), m_index,
                                       select(m.value_(), Type(1), Type(0)));
        return DiffArray::create(index_new, m_value & m.value_());
    }

    DiffArray xor_(const DiffArray &m) const {
        if constexpr (!is_mask_v<Type> && !std::is_integral_v<Scalar>)
            fail_unsupported("xor_");
        else
            return DiffArray::create(0, m_value ^ m.value_());
    }

    template <typename Mask>
    DiffArray xor_(const Mask &m) const {
        if (Enabled && m_index != 0)
            fail_unsupported("xor_ -- gradients are not implemented.");
        return DiffArray(m_value ^ m.value_());
    }

    DiffArray andnot_(const DiffArray &m) const {
        if constexpr (!is_mask_v<Type> && !std::is_integral_v<Scalar>)
            fail_unsupported("andnot_");
        else
            return DiffArray::create(0, andnot(m_value, m.value_()));
    }

    template <typename Mask>
    DiffArray andnot_(const Mask &m) const {
        if (Enabled && m_index != 0)
            fail_unsupported("andnot_ -- gradients are not implemented.");
        return DiffArray(andnot(m_value, m.value_()));
    }

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Operations that don't require derivatives
    // -----------------------------------------------------------------------

    DiffArray mod_(const DiffArray &a) const {
        if constexpr (!std::is_integral_v<Scalar>)
            fail_unsupported("mod_");
        else
            return m_value % a.m_value;
    }

    DiffArray mulhi_(const DiffArray &a) const {
        if constexpr (!std::is_integral_v<Scalar>)
            fail_unsupported("mulhi_");
        else
            return mulhi(m_value, a.m_value);
    }

    DiffArray not_() const {
        if constexpr ((!is_mask_v<Type> && !std::is_integral_v<Scalar>) ||
                      std::is_pointer_v<Scalar>)
            fail_unsupported("not_");
        else
            return DiffArray::create(0, ~m_value);
    }

    template <typename Mask>
    ENOKI_INLINE value_t<Type> extract_(const Mask &mask) const {
        if constexpr (is_mask_v<Type> || Enabled)
            fail_unsupported("extract_");
        else
            return extract(m_value, mask.value_());
    }

    DiffArray lzcnt_() const {
        if constexpr ((!is_mask_v<Type> && !std::is_integral_v<Scalar>) ||
                      std::is_pointer_v<Scalar>)
            fail_unsupported("lzcnt_");
        else
            return DiffArray::create(0, lzcnt(m_value));
    }

    DiffArray tzcnt_() const {
        if constexpr ((!is_mask_v<Type> && !std::is_integral_v<Scalar>) ||
                      std::is_pointer_v<Scalar>)
            fail_unsupported("tzcnt_");
        else
            return DiffArray::create(0, tzcnt(m_value));
    }

    DiffArray popcnt_() const {
        if constexpr ((!is_mask_v<Type> && !std::is_integral_v<Scalar>) ||
                      std::is_pointer_v<Scalar>)
            fail_unsupported("popcnt_");
        else
            return DiffArray::create(0, popcnt(m_value));
    }

    template <size_t Imm> DiffArray sl_() const {
        if constexpr (is_mask_v<Type> || !std::is_integral_v<Scalar>)
            fail_unsupported("sl_");
        else
            return DiffArray::create(0, sl<Imm>(m_value));
    }

    template <size_t Imm> DiffArray sr_() const {
        if constexpr (is_mask_v<Type> || !std::is_integral_v<Scalar>)
            fail_unsupported("sr_");
        else
            return DiffArray::create(0, sr<Imm>(m_value));
    }

    DiffArray sl_(const DiffArray &a) const {
        if constexpr (is_mask_v<Type> || !std::is_integral_v<Scalar>)
            fail_unsupported("sl_");
        else
            return DiffArray::create(0, m_value << a.m_value);
    }

    DiffArray sr_(const DiffArray &a) const {
        if constexpr (is_mask_v<Type> || !std::is_integral_v<Scalar>)
            fail_unsupported("sr_");
        else
            return DiffArray::create(0, m_value >> a.m_value);
    }

    DiffArray sl_(size_t size) const {
        if constexpr (is_mask_v<Type> || !std::is_integral_v<Scalar>)
            fail_unsupported("sl_");
        else
            return DiffArray::create(0, m_value << size);
    }

    DiffArray sr_(size_t size) const {
        if constexpr (is_mask_v<Type> || !std::is_integral_v<Scalar>)
            fail_unsupported("sr_");
        else
            return DiffArray::create(0, m_value >> size);
    }

    template <size_t Imm> DiffArray rol_() const {
        if constexpr (is_mask_v<Type> || !std::is_integral_v<Scalar>)
            fail_unsupported("rol_");
        else
            return DiffArray::create(0, rol<Imm>(m_value));
    }

    template <size_t Imm> DiffArray ror_() const {
        if constexpr (is_mask_v<Type> || !std::is_integral_v<Scalar>)
            fail_unsupported("ror_");
        else
            return DiffArray::create(0, ror<Imm>(m_value));
    }

    DiffArray rol_(const DiffArray &a) const {
        if constexpr (is_mask_v<Type> || !std::is_integral_v<Scalar>)
            fail_unsupported("rol_");
        else
            return DiffArray::create(0, rol(m_value, a.m_value));
    }

    DiffArray ror_(const DiffArray &a) const {
        if constexpr (is_mask_v<Type> || !std::is_integral_v<Scalar>)
            fail_unsupported("ror_");
        else
            return DiffArray::create(0, ror(m_value, a.m_value));
    }

    auto eq_ (const DiffArray &d) const { return MaskType(eq(m_value, d.m_value)); }
    auto neq_(const DiffArray &d) const { return MaskType(neq(m_value, d.m_value)); }
    auto lt_ (const DiffArray &d) const { return MaskType(m_value < d.m_value); }
    auto le_ (const DiffArray &d) const { return MaskType(m_value <= d.m_value); }
    auto gt_ (const DiffArray &d) const { return MaskType(m_value > d.m_value); }
    auto ge_ (const DiffArray &d) const { return MaskType(m_value >= d.m_value); }

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Scatter/gather operations
    // -----------------------------------------------------------------------

    template <size_t Stride, typename Offset, typename Mask>
    static DiffArray gather_(const void *ptr, const Offset &offset,
                             const Mask &mask) {
        static_assert(!Enabled || Stride == sizeof(Scalar),
                      "Differentiable gather: unsupported stride!");


        Type result = gather<Type, Stride>(ptr, offset.value_(), mask.value_());

        Index index_new = 0;
        if constexpr (Enabled)
            index_new = tape()->append_gather(offset.value_(), mask.value_());

        return DiffArray::create(index_new, std::move(result));
    }

    template <size_t Stride, typename Offset, typename Mask>
    void scatter_(void *ptr, const Offset &offset, const Mask &mask) const {
        static_assert(!Enabled || Stride == sizeof(Scalar),
                      "Differentiable scatter: unsupported stride!");

        scatter<Stride>(ptr, m_value, offset.value_(), mask.value_());

        if constexpr (Enabled)
            tape()->append_scatter(m_index, offset.value_(), mask.value_(), false);
    }

    template <size_t Stride, typename Offset, typename Mask>
    void scatter_add_(void *ptr, const Offset &offset, const Mask &mask) const {
        static_assert(!Enabled || Stride == sizeof(Scalar),
                      "Differentiable scatter_add: unsupported stride!");

        scatter_add<Stride>(ptr, m_value, offset.value_(), mask.value_());

        if constexpr (Enabled)
            tape()->append_scatter(m_index, offset.value_(), mask.value_(), true);
    }

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Horizontal operations
    // -----------------------------------------------------------------------

    auto all_() const {
        if constexpr (!is_mask_v<Type>)
            fail_unsupported("all_");
        else
            return all(m_value);
    }

    auto any_() const {
        if constexpr (!is_mask_v<Type>)
            fail_unsupported("any_");
        else
            return any(m_value);
    }

    auto count_() const {
        if constexpr (!is_mask_v<Type>)
            fail_unsupported("count_");
        else
            return count(m_value);
    }

    DiffArray reverse_() const {
        if constexpr (is_mask_v<Type> || std::is_pointer_v<Scalar>) {
            fail_unsupported("reverse_");
        } else {
            Index index_new = 0;
            if constexpr (Enabled)
                index_new = tape()->append_reverse(m_index);

            return DiffArray::create(index_new, reverse(m_value));
        }
    }

    DiffArray psum_() const {
        if constexpr (is_mask_v<Type> || std::is_pointer_v<Scalar>) {
            fail_unsupported("psum_");
        } else {
            Index index_new = 0;
            if constexpr (Enabled)
                index_new = tape()->append_psum(m_index);

            return DiffArray::create(index_new, psum(m_value));
        }
    }

    DiffArray hsum_() const {
        if constexpr (is_mask_v<Type> || std::is_pointer_v<Scalar>) {
            fail_unsupported("hsum_");
        } else {
            Index index_new = 0;
            if constexpr (Enabled)
                index_new = tape()->append("hsum", 1, m_index, 1.f);

            return DiffArray::create(index_new, hsum(m_value));
        }
    }

    DiffArray hprod_() const {
        if constexpr (is_mask_v<Type> || std::is_pointer_v<Scalar>) {
            fail_unsupported("hprod_");
        } else {
            Index index_new = 0;
            Type result = hprod(m_value);
            if constexpr (Enabled)
                index_new = tape()->append(
                    "hprod", 1, m_index,
                    select(eq(m_value, (Scalar) 0), (Scalar) 0, result / m_value));
            return DiffArray::create(index_new, std::move(result));
        }
    }

    DiffArray hmax_() const {
        if constexpr (is_mask_v<Type> || std::is_pointer_v<Scalar>) {
            fail_unsupported("hmax_");
        } else {
            if (Enabled && m_index != 0)
                fail_unsupported("hmax_: gradients not yet implemented!");
            return DiffArray::create(0, hmax(m_value));
        }
    }

    DiffArray hmin_() const {
        if constexpr (is_mask_v<Type> || std::is_pointer_v<Scalar>) {
            fail_unsupported("hmin_");
        } else {
            if (Enabled && m_index != 0)
                fail_unsupported("hmin_: gradients not yet implemented!");
            return DiffArray::create(0, hmin(m_value));
        }
    }

    template <typename T = Scalar, enable_if_t<std::is_pointer_v<T>> = 0>
    auto partition_() const {
        std::vector<std::pair<T, uint32_array_t<DiffArray, false>>> result;

        auto p = partition(m_value);
        result.reserve(p.size());

        for (auto &kv : p)
            result.emplace_back(kv.first, std::move(kv.second));

        return result;
    }

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Access to internals
    // -----------------------------------------------------------------------

    void set_index_(Index index) {
        if constexpr (Enabled) {
            auto t = tape();
            t->inc_ref_ext(index);
            t->dec_ref_ext(m_index);
        }
        m_index = index;
    }
    Index index_() const { return m_index; }
    Type &value_() { return m_value; }
    const Type &value_() const { return m_value; }

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Coefficient access
    // -----------------------------------------------------------------------

    ENOKI_INLINE size_t size() const {
        if constexpr (is_scalar_v<Type>)
            return 1;
        else
            return slices(m_value);
    }

    ENOKI_INLINE bool empty() const {
        if constexpr (is_scalar_v<Type>)
            return false;
        else
            return slices(m_value) == 0;
    }

    ENOKI_NOINLINE void resize(size_t size) {
        ENOKI_MARK_USED(size);
        if constexpr (!is_scalar_v<Type>)
            m_value.resize(size);
    }

    ENOKI_INLINE Scalar *data() {
        if constexpr (is_scalar_v<Type>)
            return &m_value;
        else
            return m_value.data();
    }

    ENOKI_INLINE const Scalar *data() const {
        if constexpr (is_scalar_v<Type>)
            return &m_value;
        else
            return m_value.data();
    }

    template <typename... Args>
    ENOKI_INLINE decltype(auto) coeff(Args... args) {
        static_assert(sizeof...(Args) == Depth, "coeff(): Invalid number of arguments!");
        if constexpr (is_scalar_v<Type>)
            return m_value;
        else
            return m_value.coeff((size_t) args...);
    }

    template <typename... Args>
    ENOKI_INLINE decltype(auto) coeff(Args... args) const {
        static_assert(sizeof...(Args) == Depth, "coeff(): Invalid number of arguments!");
        if constexpr (is_scalar_v<Type>)
            return m_value;
        else
            return m_value.coeff((size_t) args...);
    }

    const Type &gradient_() const {
        if constexpr (!Enabled)
            fail_unsupported("gradient_");
        else
            return tape()->gradient(m_index);
    }

    static const Type &gradient_static_(Index index) {
        if constexpr (!Enabled)
            fail_unsupported("gradient_static_");
        else
            return tape()->gradient(index);
    }

    void set_gradient_(const Type &value, bool backward = true) {
        if constexpr (!Enabled)
            fail_unsupported("set_gradient_");
        else
            return tape()->set_gradient(m_index, value, backward);
    }

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Standard initializers
    // -----------------------------------------------------------------------

    template <typename... Args>
    static DiffArray empty_(Args... args) { return enoki::empty<Type>(args...); }
    template <typename... Args>
    static DiffArray zero_(Args... args) { return zero<Type>(args...); }
    template <typename... Args>
    static DiffArray arange_(Args... args) { return arange<Type>(args...); }
    template <typename... Args>
    static DiffArray linspace_(Args... args) { return linspace<Type>(args...); }
    template <typename... Args>
    static DiffArray full_(Args... args) { return full<Type>(args...); }

    //! @}
    // -----------------------------------------------------------------------

    void set_requires_gradient_(bool value) {
        if constexpr (!Enabled) {
            fail_unsupported("set_requires_gradient_");
        } else {
            if (value && m_index == 0) {
                m_index = tape()->append_leaf(slices(m_value));
            } else if (!value && m_index != 0) {
                tape()->dec_ref_ext(m_index);
                m_index = 0;
            }
        }
    }

    bool requires_gradient_() const {
        return Enabled && m_index != 0;
    }

    void set_label_(const char *label) const {
        ENOKI_MARK_USED(label);
        if constexpr (Enabled)
            tape()->set_label(m_index, label);
        set_label(m_value, label);
    }

    void backward_(bool free_graph) const {
        if constexpr (!Enabled) {
            fail_unsupported("backward_");
        } else {
            tape()->backward(m_index, free_graph);
        }
    }

    void forward_(bool free_graph) const {
        if constexpr (!Enabled) {
            fail_unsupported("forward_");
        } else {
            tape()->forward(m_index, free_graph);
        }
    }

    static void backward_static_(bool free_graph) {
        tape()->backward(free_graph);
    }

    static void forward_static_(bool free_graph) {
        tape()->forward(free_graph);
    }

    static std::string graphviz_(const std::vector<Index> &indices) {
        if constexpr (!Enabled)
            fail_unsupported("graphviz_");
        else
            return tape()->graphviz(indices);
    }

    static void push_prefix_(const char *label) {
        if constexpr (Enabled)
            tape()->push_prefix(label);
    }

    static void pop_prefix_() {
        if constexpr (Enabled)
            tape()->pop_prefix();
    }

    static void inc_ref_ext_(Index index) {
        if constexpr (Enabled)
            tape()->inc_ref_ext(index);
    }

    static void dec_ref_ext_(Index index) {
        if constexpr (Enabled)
            tape()->dec_ref_ext(index);
    }

    static void set_scatter_gather_operand_(const DiffArray &v, bool permute) {
        ENOKI_MARK_USED(v);
        ENOKI_MARK_USED(permute);
        if constexpr (Enabled)
            tape()->set_scatter_gather_operand(const_cast<Index *>(&v.m_index),
                                               v.size(), permute);
    }

    static void clear_scatter_gather_operand_() {
        if constexpr (Enabled)
            tape()->set_scatter_gather_operand(nullptr, 0, false);
    }

    static void set_log_level_(uint32_t level) {
        if constexpr (Enabled)
            tape()->set_log_level(level);
    }

    static uint32_t log_level_() {
        if constexpr (Enabled)
            return tape()->log_level();
        else
            return 0;
    }


    static void set_graph_simplification_(uint32_t level) {
        if constexpr (Enabled)
            tape()->set_graph_simplification(level);
    }

    static void simplify_graph_() {
        if constexpr (Enabled)
            tape()->simplify_graph();
    }

    static std::string whos_() {
        if constexpr (!Enabled)
            fail_unsupported("whos");
        else
            return tape()->whos();
    }

    static DiffArray map(void *ptr, size_t size, bool dealloc = false) {
        if constexpr (!is_dynamic_array_v<Type>)
            fail_unsupported("map");
        else
            return DiffArray::create(0, Type::map(ptr, size, dealloc));
    }

    static DiffArray copy(const void *ptr, size_t size) {
        if constexpr (!is_dynamic_array_v<Type>)
            fail_unsupported("copy");
        else
            return DiffArray::create(0, Type::copy(ptr, size));
    }

    DiffArray &managed() {
        if constexpr (is_cuda_array_v<Type>)
            m_value.managed();
        return *this;
    }

    const DiffArray &managed() const {
        if constexpr (is_cuda_array_v<Type>)
            m_value.managed();
        return *this;
    }


    DiffArray &eval() {
        if constexpr (is_cuda_array_v<Type>)
            m_value.eval();
        return *this;
    }

    const DiffArray &eval() const {
        if constexpr (is_cuda_array_v<Type>)
            m_value.eval();
        return *this;
    }

    auto operator->() const {
        using BaseType = std::decay_t<std::remove_pointer_t<Scalar>>;
        return call_support<BaseType, DiffArray>(*this);
    }

private:
    ENOKI_INLINE static Tape* tape() { return Tape::get(); }

    using Arg = std::conditional_t<std::is_scalar_v<Type>, Type, Type&&>;

    ENOKI_INLINE static DiffArray create(Index index, Arg value) {
        DiffArray result(std::move(value));
        result.m_index = index;
        return result;
    }

    [[noreturn]]
    ENOKI_NOINLINE static void fail_unsupported(const char *msg) {
        fprintf(stderr, "DiffArray: unsupported operation for type %s", msg);
        exit(EXIT_FAILURE);
    }

    Type m_value;
    Index m_index = 0;
};

template <typename T, enable_if_t<is_diff_array_v<T>> = 0>
ENOKI_INLINE void set_label(const T& a, const char *label) {
    if constexpr (array_depth_v<T> >= 2) {
        for (size_t i = 0; i < T::Size; ++i)
            set_label(a.coeff(i), (std::string(label) + "." + std::to_string(i)).c_str());
    } else {
        a.set_label_(label);
    }
}

template <typename T> ENOKI_INLINE bool requires_gradient(T& a) {
    if constexpr (is_diff_array_v<T>) {
        if constexpr (array_depth_v<T> >= 2) {
            for (size_t i = 0; i < a.size(); ++i) {
                if (requires_gradient(a.coeff(i)))
                    return true;
            }
            return false;
        } else {
            return a.requires_gradient_();
        }
    }
    return false;
}

template <typename T> ENOKI_INLINE void set_requires_gradient(T& a, bool value = true) {
    if constexpr (is_diff_array_v<T>) {
        if constexpr (array_depth_v<T> >= 2) {
            for (size_t i = 0; i < a.size(); ++i)
                set_requires_gradient(a.coeff(i), value);
        } else {
            a.set_requires_gradient_(value);
        }
    }
}

template <typename T> auto gradient_index(const T &a) {
    if constexpr (array_depth_v<T> >= 2) {
        using Result = std::array<decltype(gradient_index(a.coeff(0))), T::Size>;
        Result result;
        for (size_t i = 0; i < T::Size; ++i)
            result[i] = gradient_index(a.coeff(i));
        return result;
    } else if constexpr (is_diff_array_v<T>) {
        return a.index_();
    } else {
        static_assert(detail::false_v<T>, "The given array does not support derivatives.");
    }
}

template <typename T1, typename T2> void set_gradient(T1 &a, const T2 &b, bool backward = true) {
    if constexpr (array_depth_v<T1> >= 2) {
        for (size_t i = 0; i < array_size_v<T1>; ++i)
            set_gradient(a[i], b[i], backward);
    } else if constexpr (is_diff_array_v<T1>) {
        a.set_gradient_(b, backward);
    } else {
        static_assert(detail::false_v<T1, T2>, "The given array does not support derivatives.");
    }
}

template <typename T1> void reattach(T1 &a, const T1 &b) {
    if constexpr (array_depth_v<T1> >= 2) {
        for (size_t i = 0; i < array_size_v<T1>; ++i)
            reattach(a[i], b[i]);
    } else if constexpr (is_diff_array_v<T1>) {
        a.set_index_(b.index_());
    } else {
        static_assert(detail::false_v<T1>, "The given array does not support derivatives.");
    }
}

template <typename T> void forward(const T& a, bool free_graph = true) {
    a.forward_(free_graph);
}

template <typename T> void backward(const T& a, bool free_graph = true) {
    a.backward_(free_graph);
}

template <typename T> void backward(bool free_graph = true) {
    T::backward_static_(free_graph);
}

template <typename T> void forward(bool free_graph = true) {
    T::forward_static_(free_graph);
}

namespace detail {
    template <typename T>
    void collect_indices(const T &value, std::vector<uint32_t> &indices) {
        if constexpr (is_diff_array_v<T>) {
            if constexpr (array_depth_v<T> == 1) {
                if (value.index_() != 0)
                    indices.push_back(value.index_());
            } else {
                for (size_t i = 0; i < T::Size; ++i)
                    collect_indices(value.coeff(i), indices);
            }
        }
    }
};

namespace detail {
    template <typename T, typename = int> struct diff_type {
        using type = T;
    };
    template <typename T> using diff_type_t = typename diff_type<T>::type;
    template <typename T> struct diff_type<T, enable_if_t<is_diff_array_v<value_t<T>>>> {
        using type = diff_type_t<value_t<T>>;
    };
}

template <typename T> std::string graphviz(const T &value) {
    std::vector<uint32_t> indices;
    detail::collect_indices(value, indices);
    return detail::diff_type_t<T>::graphviz_(indices);
}

#if defined(ENOKI_AUTODIFF_BUILD)
#  define ENOKI_AUTODIFF_EXTERN extern
#  define ENOKI_AUTODIFF_EXPORT ENOKI_EXPORT
#else
#  define ENOKI_AUTODIFF_EXPORT ENOKI_IMPORT
#  if defined(_MSC_VER)
#    define ENOKI_AUTODIFF_EXTERN
#else
#    define ENOKI_AUTODIFF_EXTERN extern
#  endif
#endif

#if !defined(ENOKI_BUILD)
    ENOKI_AUTODIFF_EXTERN template struct ENOKI_AUTODIFF_EXPORT Tape<float>;
    ENOKI_AUTODIFF_EXTERN template struct ENOKI_AUTODIFF_EXPORT DiffArray<float>;

    ENOKI_AUTODIFF_EXTERN template struct ENOKI_AUTODIFF_EXPORT Tape<double>;
    ENOKI_AUTODIFF_EXTERN template struct ENOKI_AUTODIFF_EXPORT DiffArray<double>;

#  if defined(ENOKI_DYNAMIC_H)
        ENOKI_AUTODIFF_EXTERN template struct ENOKI_AUTODIFF_EXPORT Tape<DynamicArray<Packet<float>>>;
        ENOKI_AUTODIFF_EXTERN template struct ENOKI_AUTODIFF_EXPORT DiffArray<DynamicArray<Packet<float>>>;

        ENOKI_AUTODIFF_EXTERN template struct ENOKI_AUTODIFF_EXPORT Tape<DynamicArray<Packet<double>>>;
        ENOKI_AUTODIFF_EXTERN template struct ENOKI_AUTODIFF_EXPORT DiffArray<DynamicArray<Packet<double>>>;
#  endif

#  if defined(ENOKI_CUDA_H)
        ENOKI_AUTODIFF_EXTERN template struct ENOKI_AUTODIFF_EXPORT Tape<CUDAArray<float>>;
        ENOKI_AUTODIFF_EXTERN template struct ENOKI_AUTODIFF_EXPORT DiffArray<CUDAArray<float>>;

        ENOKI_AUTODIFF_EXTERN template struct ENOKI_AUTODIFF_EXPORT Tape<CUDAArray<double>>;
        ENOKI_AUTODIFF_EXTERN template struct ENOKI_AUTODIFF_EXPORT DiffArray<CUDAArray<double>>;
#  endif
#endif

NAMESPACE_END(enoki)
