/*
    enoki/quaternion.h -- Matrix data structure

    Enoki is a C++ template library that enables transparent vectorization
    of numerical kernels using SIMD instruction sets available on current
    processor architectures.

    Copyright (c) 2019 Wenzel Jakob <wenzel.jakob@epfl.ch>

    All rights reserved. Use of this source code is governed by a BSD-style
    license that can be found in the LICENSE file.
*/

#pragma once

#include <enoki/array.h>

NAMESPACE_BEGIN(enoki)

/// Value trait to access the column type of a matrix
template <typename T> using column_t = typename std::decay_t<T>::Column;

/// Value trait to access the entry type of a matrix
template <typename T> using entry_t = value_t<column_t<T>>;

/// SFINAE helper for matrixs
template <typename T> using is_matrix_helper = enable_if_t<std::decay_t<T>::IsMatrix>;
template <typename T> constexpr bool is_matrix_v = is_detected_v<is_matrix_helper, T>;
template <typename T> using enable_if_matrix_t = enable_if_t<is_matrix_v<T>>;
template <typename T> using enable_if_not_matrix_t = enable_if_t<!is_matrix_v<T>>;

template <typename Value_, size_t Size_>
struct Matrix : StaticArrayImpl<Array<Value_, Size_>, Size_, false, Matrix<Value_, Size_>> {

    using Entry = Value_;
    using Column = Array<Entry, Size_>;

    using Base = StaticArrayImpl<Column, Size_, false, Matrix<Value_, Size_>>;
    using Base::coeff;

    ENOKI_ARRAY_IMPORT_BASIC(Base, Matrix);
    using Base::operator=;

    static constexpr bool IsMatrix = true;
    static constexpr bool IsVector = false;

    using ArrayType = Matrix;
    using MaskType = Mask<mask_t<Column>, Size_>;

    template <typename T> using ReplaceValue = Matrix<value_t<T>, Size_>;

    Matrix() = default;

    /// Initialize from a incompatible matrix
    template <typename Value2, size_t Size2, enable_if_t<Size2 == Size_> = 0>
    ENOKI_INLINE Matrix(const Matrix<Value2, Size2> &m)
     : Base(m) { }

    /// Initialize from an incompatible matrix
    template <size_t Size2, enable_if_t<Size2 != Size_> = 0>
    ENOKI_INLINE Matrix(const Matrix<Value_, Size2> &m) {
        if constexpr (Size2 > Size) {
            /// Other matrix is bigger -- retain the top left part
            for (size_t i = 0; i < Size; ++i)
                coeff(i) = head<Size>(m.coeff(i));
        } else {
            /// Other matrix is smaller -- copy the top left part and set remainder to identity
            using Remainder = Array<Value_, Size - Size2>;
            for (size_t i = 0; i < Size2; ++i)
                coeff(i) = concat(m.coeff(i), zero<Remainder>());
            for (size_t i = Size2; i < Size; ++i) {
                auto col = zero<Column>();
                col.coeff(i) = 1;
                coeff(i) = col;
            }
        }
    }

    template <typename T, enable_if_t<(array_depth_v<T> <= Base::Depth - 2)> = 0,
                          enable_if_not_matrix_t<T> = 0>
    ENOKI_INLINE Matrix(T&& v) {
        for (size_t i = 0; i < Size; ++i) {
            coeff(i) = zero<Column>();
            coeff(i, i) = v;
        }
    }

    template <typename T, enable_if_t<(array_depth_v<T> == Base::Depth)> = 0,
                          enable_if_not_matrix_t<T> = 0>
    ENOKI_INLINE Matrix(T&& v) : Base(std::forward<T>(v)) { }

    /// Initialize the matrix from a list of columns
    template <typename... Args, enable_if_t<sizeof...(Args) == Size_ &&
              std::conjunction_v<std::is_constructible<Column, Args>...>> = 0>
    ENOKI_INLINE Matrix(const Args&... args) : Base(args...) { }

    /// Initialize the matrix from a list of entries in row-major order
    template <typename... Args, enable_if_t<sizeof...(Args) == Size_ * Size_ &&
              std::conjunction_v<std::is_constructible<Entry, Args>...>> = 0>
    ENOKI_INLINE Matrix(const Args&... args) {
        alignas(alignof(Column)) Entry values[sizeof...(Args)] = { Entry(args)... };
        for (size_t j = 0; j < Size; ++j)
            for (size_t i = 0; i < Size; ++i)
                coeff(j, i) = values[i * Size + j];
    }

    template <typename... Column>
    ENOKI_INLINE static Matrix from_cols(const Column&... cols) {
        return Matrix(cols...);
    }

    template <typename... Row>
    ENOKI_INLINE static Matrix from_rows(const Row&... rows) {
        return transpose(Matrix(rows...));
    }

    ENOKI_INLINE Column& col(size_t index) { return coeff(index); }
    ENOKI_INLINE const Column& col(size_t index) const { return coeff(index); }

    ENOKI_INLINE Column row(size_t index) const {
        using Index = Array<uint32_t, Size>;
        return gather<Column>(coeff(0).data() + index,
                              arange<Index>() * uint32_t(Size));
    }

    /// Return a reference to the (i, j) element
    ENOKI_INLINE decltype(auto) operator()(size_t i, size_t j) { return coeff(j, i); }

    /// Return a reference to the (i, j) element (const)
    ENOKI_INLINE decltype(auto) operator()(size_t i, size_t j) const { return coeff(j, i); }

    static ENOKI_INLINE Derived zero_(size_t size) {
        Derived result;
        for (size_t i = 0; i < Size; ++i)
            result.coeff(i) = zero<Column>(size);
        return result;
    }

    static ENOKI_INLINE Derived empty_(size_t size) {
        Derived result;
        for (size_t i = 0; i < Size; ++i)
            result.coeff(i) = empty<Column>(size);
        return result;
    }

    template <typename T> ENOKI_INLINE static Matrix full_(const T &value, size_t size) {
        return Array<Column, Size>::full_(value, size);
    }
};

template <typename T0, typename T1, size_t Size,
          typename Result = Matrix<expr_t<T0, T1>, Size>,
          typename Column = column_t<Result>>
ENOKI_INLINE Result operator*(const Matrix<T0, Size> &m0,
                              const Matrix<T1, Size> &m1) {
    Result result;
    /* 4x4 case reduced to 4 multiplications, 12 fused multiply-adds,
       and 16 broadcasts (also fused on AVX512VL) */
    for (size_t j = 0; j < Size; ++j) {
        Column sum = m0.coeff(0) * Column::full_(m1(0, j), 1);
        for (size_t i = 1; i < Size; ++i)
            sum = fmadd(m0.coeff(i), Column::full_(m1(i, j), 1), sum);
        result.coeff(j) = sum;
    }

    return result;
}

template <typename T0, typename T1, size_t Size, enable_if_t<!T1::IsMatrix> = 0>
ENOKI_INLINE auto operator*(const Matrix<T0, Size> &m, const T1 &s) {
    if constexpr (array_size_v<T1> == Size && T1::IsVector) {
        using EValue  = expr_t<T0, value_t<T1>>;
        using EVector = Array<EValue, Size>;
        EVector sum = m.coeff(0) * EVector::full_(s.coeff(0), 1);
        for (size_t i = 1; i < Size; ++i)
            sum = fmadd(m.coeff(i), EVector::full_(s.coeff(i), 1), sum);
        return sum;
    } else {
        using EValue  = expr_t<T0, T1>;
        using EArray  = Array<Array<EValue, Size>, Size>;
        using EMatrix = Matrix<EValue, Size>;

        return EMatrix(EArray(m) * EArray::full_(EValue(s), 1));
    }
}

template <typename T0, typename T1, size_t Size, enable_if_t<!T0::IsMatrix> = 0>
ENOKI_INLINE auto operator*(const T0 &s, const Matrix<T1, Size> &m) {
    using EValue  = expr_t<T0, T1>;
    using EArray  = Array<Array<EValue, Size>, Size>;
    using EMatrix = Matrix<EValue, Size>;

    return EMatrix(EArray::full_(EValue(s), 1) * EArray(m));
}

template <typename T0, typename T1, size_t Size, enable_if_t<!T1::IsMatrix> = 0>
ENOKI_INLINE auto operator/(const Matrix<T0, Size> &m, const T1 &s) {
    using EValue  = expr_t<T0, T1>;
    using EArray  = Array<Array<EValue, Size>, Size>;
    using EMatrix = Matrix<EValue, Size>;

    return EMatrix(EArray(m) * EArray::full_(rcp(EValue(s)), 1));
}

template <typename Value, size_t Size>
ENOKI_INLINE expr_t<Value> trace(const Matrix<Value, Size> &m) {
    expr_t<Value> result = m.coeff(0, 0);
    for (size_t i = 1; i < Size; ++i)
        result += m(i, i);
    return result;
}

template <typename Value, size_t Size>
ENOKI_INLINE expr_t<Value> frob(const Matrix<Value, Size> &matrix) {
    expr_t<column_t<Matrix<Value, Size>>> result = sqr(matrix.coeff(0));
    for (size_t i = 1; i < Size; ++i)
        result = fmadd(matrix.coeff(i), matrix.coeff(i), result);
    return hsum(result);
}

template <typename T, enable_if_matrix_t<T> = 0>
ENOKI_INLINE T identity(size_t size = 1) {
    T result = zero<T>(size);
    for (size_t i = 0; i < T::Size; ++i)
        result(i, i) = full<typename T::Entry>(scalar_t<T>(1.f), size);
    return result;
}


template <typename Matrix, enable_if_matrix_t<Matrix> = 0>
ENOKI_INLINE Matrix diag(const column_t<Matrix> &value) {
    Matrix result = zero<Matrix>();
    for (size_t i = 0; i < Matrix::Size; ++i)
        result(i, i) = value.coeff(i);
    return result;
}

template <typename Matrix, enable_if_matrix_t<Matrix> = 0>
ENOKI_INLINE column_t<expr_t<Matrix>> diag(const Matrix &value) {
    column_t<expr_t<Matrix>> result;
    for (size_t i = 0; i < Matrix::Size; ++i)
        result.coeff(i) = value(i, i);
    return result;
}

template <typename T, typename E = expr_t<T>>
ENOKI_INLINE Matrix<E, 1> inverse(const Matrix<T, 1> &m) {
    return rcp(m(0, 0));
}

template <typename T, typename E = expr_t<T>>
ENOKI_INLINE Matrix<E, 1>
inverse_transpose(const Matrix<T, 1> &m) {
    return rcp(m(0, 0));
}

template <typename T, typename E = expr_t<T>>
ENOKI_INLINE E det(const Matrix<T, 1> &m) {
    return m(0, 0);
}

template <typename T, typename E = expr_t<T>>
ENOKI_INLINE Matrix<E, 2> inverse(const Matrix<T, 2> &m) {
    E inv_det = rcp(fmsub(m(0, 0), m(1, 1), m(0, 1) * m(1, 0)));
    return Matrix<E, 2>(
        m(1, 1) * inv_det, -m(0, 1) * inv_det,
       -m(1, 0) * inv_det,  m(0, 0) * inv_det
    );
}

template <typename T, typename E = expr_t<T>>
ENOKI_INLINE E det(const Matrix<T, 2> &m) {
    return fmsub(m(0, 0), m(1, 1), m(0, 1) * m(1, 0));
}

template <typename T, typename E = expr_t<T>>
ENOKI_INLINE Matrix<E, 2>
inverse_transpose(const Matrix<T, 2> &m) {
    E inv_det = rcp(fmsub(m(0, 0), m(1, 1), m(0, 1) * m(1, 0)));
    return Matrix<E, 2>(
        m(1, 1) * inv_det, -m(1, 0) * inv_det,
       -m(0, 1) * inv_det,  m(0, 0) * inv_det
    );
}

template <typename T, typename E = expr_t<T>>
ENOKI_INLINE Matrix<E, 3>
inverse_transpose(const Matrix<T, 3> &m) {
    using Vector = Array<E, 3>;

    Vector col0 = m.coeff(0),
           col1 = m.coeff(1),
           col2 = m.coeff(2);

    Vector row0 = cross(col1, col2),
           row1 = cross(col2, col0),
           row2 = cross(col0, col1);

    Vector inv_det = Vector(rcp(dot(col0, row0)));

    return Matrix<E, 3>(
        row0 * inv_det,
        row1 * inv_det,
        row2 * inv_det
    );
}

template <typename T, typename E = expr_t<T>>
ENOKI_INLINE Matrix<E, 3> inverse(const Matrix<T, 3> &m) {
    return transpose(inverse_transpose(m));
}

template <typename T, typename E = expr_t<T>>
ENOKI_INLINE E det(const Matrix<T, 3> &m) {
    return dot(m.coeff(0), cross(m.coeff(1), m.coeff(2)));
}

template <typename T, typename E = expr_t<T>>
ENOKI_INLINE Matrix<E, 4>
inverse_transpose(const Matrix<T, 4> &m) {
    using Vector = Array<E, 4>;

    Vector col0 = m.coeff(0), col1 = m.coeff(1),
           col2 = m.coeff(2), col3 = m.coeff(3);

    col1 = shuffle<2, 3, 0, 1>(col1);
    col3 = shuffle<2, 3, 0, 1>(col3);

    Vector tmp, row0, row1, row2, row3;

    tmp = shuffle<1, 0, 3, 2>(col2 * col3);
    row0 = col1 * tmp;
    row1 = col0 * tmp;
    tmp = shuffle<2, 3, 0, 1>(tmp);
    row0 = fmsub(col1, tmp, row0);
    row1 = shuffle<2, 3, 0, 1>(fmsub(col0, tmp, row1));

    tmp = shuffle<1, 0, 3, 2>(col1 * col2);
    row0 = fmadd(col3, tmp, row0);
    row3 = col0 * tmp;
    tmp = shuffle<2, 3, 0, 1>(tmp);
    row0 = fnmadd(col3, tmp, row0);
    row3 = shuffle<2, 3, 0, 1>(fmsub(col0, tmp, row3));

    tmp = shuffle<1, 0, 3, 2>(shuffle<2, 3, 0, 1>(col1) * col3);
    col2 = shuffle<2, 3, 0, 1>(col2);
    row0 = fmadd(col2, tmp, row0);
    row2 = col0 * tmp;
    tmp = shuffle<2, 3, 0, 1>(tmp);
    row0 = fnmadd(col2, tmp, row0);
    row2 = shuffle<2, 3, 0, 1>(fmsub(col0, tmp, row2));

    tmp = shuffle<1, 0, 3, 2>(col0 * col1);
    row2 = fmadd(col3, tmp, row2);
    row3 = fmsub(col2, tmp, row3);
    tmp = shuffle<2, 3, 0, 1>(tmp);
    row2 = fmsub(col3, tmp, row2);
    row3 = fnmadd(col2, tmp, row3);

    tmp = shuffle<1, 0, 3, 2>(col0 * col3);
    row1 = fnmadd(col2, tmp, row1);
    row2 = fmadd(col1, tmp, row2);
    tmp = shuffle<2, 3, 0, 1>(tmp);
    row1 = fmadd(col2, tmp, row1);
    row2 = fnmadd(col1, tmp, row2);

    tmp = shuffle<1, 0, 3, 2>(col0 * col2);
    row1 = fmadd(col3, tmp, row1);
    row3 = fnmadd(col1, tmp, row3);
    tmp = shuffle<2, 3, 0, 1>(tmp);
    row1 = fnmadd(col3, tmp, row1);
    row3 = fmadd(col1, tmp, row3);

    Vector inv_det = Vector(rcp(dot(col0, row0)));

    return Matrix<E, 4>(
        row0 * inv_det, row1 * inv_det,
        row2 * inv_det, row3 * inv_det
    );
}

template <typename T, typename E = expr_t<T>>
ENOKI_INLINE Matrix<E, 4> inverse(const Matrix<T, 4> &m) {
    return transpose(inverse_transpose(m));
}

template <typename T, typename E = expr_t<T>>
ENOKI_INLINE E det(const Matrix<T, 4> &m) {
    using Vector = Array<E, 4>;

    Vector col0 = m.coeff(0), col1 = m.coeff(1),
           col2 = m.coeff(2), col3 = m.coeff(3);

    col1 = shuffle<2, 3, 0, 1>(col1);
    col3 = shuffle<2, 3, 0, 1>(col3);

    Vector tmp, row0;

    tmp = shuffle<1, 0, 3, 2>(col2 * col3);
    row0 = col1 * tmp;
    tmp = shuffle<2, 3, 0, 1>(tmp);
    row0 = fmsub(col1, tmp, row0);

    tmp = shuffle<1, 0, 3, 2>(col1 * col2);
    row0 = fmadd(col3, tmp, row0);
    tmp = shuffle<2, 3, 0, 1>(tmp);
    row0 = fnmadd(col3, tmp, row0);

    col1 = shuffle<2, 3, 0, 1>(col1);
    col2 = shuffle<2, 3, 0, 1>(col2);
    tmp = shuffle<1, 0, 3, 2>(col1 * col3);
    row0 = fmadd(col2, tmp, row0);
    tmp = shuffle<2, 3, 0, 1>(tmp);
    row0 = fnmadd(col2, tmp, row0);

    return dot(col0, row0);
}

template <typename Value, size_t Size, bool IsMask_, typename Derived>
ENOKI_INLINE auto transpose(const StaticArrayBase<Value, Size, IsMask_, Derived> &a) {
    static_assert(Value::Size == Size && array_depth<Derived>::value >= 2,
                  "Array must be a square matrix!");
    using Column = value_t<Derived>;

    if constexpr (Column::IsNative) {
        #if defined(ENOKI_X86_SSE42)
            if constexpr (std::is_same_v<value_t<Column>, float> && Size == 3) {
                __m128 c0 = a.derived().coeff(0).m,
                       c1 = a.derived().coeff(1).m,
                       c2 = a.derived().coeff(2).m;

                __m128 t0 = _mm_unpacklo_ps(c0, c1);
                __m128 t1 = _mm_unpacklo_ps(c2, c2);
                __m128 t2 = _mm_unpackhi_ps(c0, c1);
                __m128 t3 = _mm_unpackhi_ps(c2, c2);

                return Derived(
                    _mm_movelh_ps(t0, t1),
                    _mm_movehl_ps(t1, t0),
                    _mm_movelh_ps(t2, t3)
                );
            } else if constexpr (std::is_same_v<value_t<Column>, float> && Size == 4) {
                __m128 c0 = a.derived().coeff(0).m, c1 = a.derived().coeff(1).m,
                       c2 = a.derived().coeff(2).m, c3 = a.derived().coeff(3).m;

                __m128 t0 = _mm_unpacklo_ps(c0, c1);
                __m128 t1 = _mm_unpacklo_ps(c2, c3);
                __m128 t2 = _mm_unpackhi_ps(c0, c1);
                __m128 t3 = _mm_unpackhi_ps(c2, c3);

                return Derived(
                    _mm_movelh_ps(t0, t1),
                    _mm_movehl_ps(t1, t0),
                    _mm_movelh_ps(t2, t3),
                    _mm_movehl_ps(t3, t2)
                );
            }
        #endif

        #if defined(ENOKI_X86_AVX)
            if constexpr (std::is_same_v<value_t<Column>, double> && Size == 3) {
                __m256d c0 = a.derived().coeff(0).m,
                        c1 = a.derived().coeff(1).m,
                        c2 = a.derived().coeff(2).m;

                __m256d t3 = _mm256_shuffle_pd(c2, c2, 0b0000),
                        t2 = _mm256_shuffle_pd(c2, c2, 0b1111),
                        t1 = _mm256_shuffle_pd(c0, c1, 0b0000),
                        t0 = _mm256_shuffle_pd(c0, c1, 0b1111);

                return Derived(
                    _mm256_permute2f128_pd(t1, t3, 0b0010'0000),
                    _mm256_permute2f128_pd(t0, t2, 0b0010'0000),
                    _mm256_permute2f128_pd(t1, t3, 0b0011'0001)
                );
            } else if constexpr (std::is_same_v<value_t<Column>, double> && Size == 4) {
                __m256d c0 = a.derived().coeff(0).m, c1 = a.derived().coeff(1).m,
                        c2 = a.derived().coeff(2).m, c3 = a.derived().coeff(3).m;

                __m256d t3 = _mm256_shuffle_pd(c2, c3, 0b0000),
                        t2 = _mm256_shuffle_pd(c2, c3, 0b1111),
                        t1 = _mm256_shuffle_pd(c0, c1, 0b0000),
                        t0 = _mm256_shuffle_pd(c0, c1, 0b1111);

                return Derived(
                    _mm256_permute2f128_pd(t1, t3, 0b0010'0000),
                    _mm256_permute2f128_pd(t0, t2, 0b0010'0000),
                    _mm256_permute2f128_pd(t1, t3, 0b0011'0001),
                    _mm256_permute2f128_pd(t0, t2, 0b0011'0001)
                );
            }
        #endif

        #if defined(ENOKI_ARM_NEON)
            if constexpr (std::is_same_v<value_t<Column>, float> && Size == 3) {
                float32x4x2_t v01 = vtrnq_f32(a.derived().coeff(0).m, a.derived().coeff(1).m);
                float32x4x2_t v23 = vtrnq_f32(a.derived().coeff(2).m, a.derived().coeff(2).m);

                return Derived(
                    vcombine_f32(vget_low_f32 (v01.val[0]), vget_low_f32 (v23.val[0])),
                    vcombine_f32(vget_low_f32 (v01.val[1]), vget_low_f32 (v23.val[1])),
                    vcombine_f32(vget_high_f32(v01.val[0]), vget_high_f32(v23.val[0]))
                );
            } else if constexpr (std::is_same_v<value_t<Column>, float> && Size == 4) {
                float32x4x2_t v01 = vtrnq_f32(a.derived().coeff(0).m, a.derived().coeff(1).m);
                float32x4x2_t v23 = vtrnq_f32(a.derived().coeff(2).m, a.derived().coeff(3).m);

                return Derived(
                    vcombine_f32(vget_low_f32 (v01.val[0]), vget_low_f32 (v23.val[0])),
                    vcombine_f32(vget_low_f32 (v01.val[1]), vget_low_f32 (v23.val[1])),
                    vcombine_f32(vget_high_f32(v01.val[0]), vget_high_f32(v23.val[0])),
                    vcombine_f32(vget_high_f32(v01.val[1]), vget_high_f32(v23.val[1]))
                );
            }
        #endif
    }

    ENOKI_CHKSCALAR("transpose");

    Derived result;
    for (size_t i = 0; i < Size; ++i)
        for (size_t j = 0; j < Size; ++j)
            result.coeff(i, j) = a.derived().coeff(j, i);
    return result;
}

template <typename T, size_t Size, typename Expr = expr_t<T>,
          typename Matrix = Matrix<Expr, Size>>
std::pair<Matrix, Matrix> ENOKI_INLINE
polar_decomp(const enoki::Matrix<T, Size> &A, size_t it = 10) {
    using Arr = Array<Array<Expr, Size>, Size>;
    Matrix Q = A;
    for (size_t i = 0; i < it; ++i) {
        Matrix Qi = inverse_transpose(Q);
        Expr gamma = sqrt(frob(Qi) / frob(Q));
        Q = fmadd(Arr(Q), gamma * .5f, Arr(Qi) * (rcp(gamma) * 0.5f));
    }
    return std::make_pair(Q, transpose(Q) * A);
}

// =======================================================================
//! @{ \name Enoki accessors for static & dynamic vectorization
// =======================================================================

template <typename T, size_t Size>
struct struct_support<Matrix<T, Size>,
                      enable_if_static_array_t<Matrix<T, Size>>> {
    static constexpr bool IsDynamic = enoki::is_dynamic_v<T>;
    using Dynamic = Matrix<enoki::make_dynamic_t<T>, Size>;
    using Value = Matrix<T, Size>;
    using Column = column_t<Value>;

    static ENOKI_INLINE size_t slices(const Value &value) {
        return enoki::slices(value.coeff(0, 0));
    }

    static ENOKI_INLINE size_t packets(const Value &value) {
        return enoki::packets(value.coeff(0, 0));
    }

    static ENOKI_INLINE void set_slices(Value &value, size_t size) {
        for (size_t i = 0; i < Size; ++i)
            enoki::set_slices(value.coeff(i), size);
    }

    template <typename T2>
    static ENOKI_INLINE auto packet(T2&& value, size_t i) {
        return packet(value, i, std::make_index_sequence<Size>());
    }

    template <typename T2>
    static ENOKI_INLINE auto slice(T2&& value, size_t i) {
        return slice(value, i, std::make_index_sequence<Size>());
    }

    template <typename T2>
    static ENOKI_INLINE auto slice_ptr(T2&& value, size_t i) {
        return slice_ptr(value, i, std::make_index_sequence<Size>());
    }

    template <typename T2>
    static ENOKI_INLINE auto ref_wrap(T2&& value) {
        return ref_wrap(value, std::make_index_sequence<Size>());
    }

    template <typename T2>
    static ENOKI_INLINE auto detach(T2&& value) {
        return detach(value, std::make_index_sequence<Size>());
    }

    template <typename T2>
    static ENOKI_INLINE auto gradient(T2&& value) {
        return gradient(value, std::make_index_sequence<Size>());
    }

    static ENOKI_INLINE Value zero(size_t size) {
        return Value::zero_(size);
    }

    static ENOKI_INLINE Value empty(size_t size) {
        return Value::empty_(size);
    }

    template <typename T2, typename Mask,
              enable_if_t<array_size<T2>::value == array_size<Mask>::value> = 0>
    static ENOKI_INLINE auto masked(T2 &value, const Mask &mask) {
        return detail::MaskedArray<T2>{ value, mask_t<T2>(mask) };
    }

    template <typename T2, typename Mask,
              enable_if_t<array_size<T2>::value != array_size<Mask>::value> = 0>
    static ENOKI_INLINE auto masked(T2 &value, const Mask &mask) {
        using Arr = Array<Array<T, Size>, Size>;
        return enoki::masked((Arr&) value, mask_t<Arr>(mask));
    }

private:
    template <typename T2, size_t... Index>
    static ENOKI_INLINE auto packet(T2&& value, size_t i, std::index_sequence<Index...>) {
        return Matrix<decltype(enoki::packet(value.coeff(0, 0), i)), Size>(
            enoki::packet(value.coeff(Index), i)...);
    }

    template <typename T2, size_t... Index>
    static ENOKI_INLINE auto slice(T2&& value, size_t i, std::index_sequence<Index...>) {
        return Matrix<decltype(enoki::slice(value.coeff(0, 0), i)), Size>(
            enoki::slice(value.coeff(Index), i)...);
    }

    template <typename T2, size_t... Index>
    static ENOKI_INLINE auto slice_ptr(T2&& value, size_t i, std::index_sequence<Index...>) {
        return Matrix<decltype(enoki::slice_ptr(value.coeff(0, 0), i)), Size>(
            enoki::slice_ptr(value.coeff(Index), i)...);
    }

    template <typename T2, size_t... Index>
    static ENOKI_INLINE auto ref_wrap(T2&& value, std::index_sequence<Index...>) {
        return Matrix<decltype(enoki::ref_wrap(value.coeff(0, 0))), Size>(
            enoki::ref_wrap(value.coeff(Index))...);
    }

    template <typename T2, size_t... Index>
    static ENOKI_INLINE auto detach(T2&& value, std::index_sequence<Index...>) {
        return Matrix<decltype(enoki::detach(value.coeff(0, 0))), Size>(
            enoki::detach(value.coeff(Index))...);
    }

    template <typename T2, size_t... Index>
    static ENOKI_INLINE auto gradient(T2&& value, std::index_sequence<Index...>) {
        return Matrix<decltype(enoki::gradient(value.coeff(0, 0))), Size>(
            enoki::gradient(value.coeff(Index))...);
    }
};

//! @}
// =======================================================================

NAMESPACE_END(enoki)
