/*
    enoki/transform.h -- 3D homogeneous coordinate transformations

    Enoki is a C++ template library that enables transparent vectorization
    of numerical kernels using SIMD instruction sets available on current
    processor architectures.

    Copyright (c) 2019 Wenzel Jakob <wenzel.jakob@epfl.ch>

    All rights reserved. Use of this source code is governed by a BSD-style
    license that can be found in the LICENSE file.
*/

#pragma once

#include <enoki/quaternion.h>

NAMESPACE_BEGIN(enoki)

template <typename Matrix, typename Vector> ENOKI_INLINE Matrix translate(const Vector &v) {
    Matrix trafo = identity<Matrix>();
    trafo.coeff(Matrix::Size - 1) = concat(v, scalar_t<Matrix>(1));
    return trafo;
}

template <typename Matrix, typename Vector> ENOKI_INLINE Matrix scale(const Vector &v) {
    return diag<Matrix>(concat(v, scalar_t<Matrix>(1)));
}

template <typename Matrix, enable_if_t<Matrix::IsMatrix && Matrix::Size == 3> = 0>
ENOKI_INLINE Matrix rotate(const entry_t<Matrix> &angle) {
    entry_t<Matrix> z(0.f), o(1.f);
    auto [s, c] = sincos(angle);
    return Matrix(c, -s, z, s, c, z, z, z, o);
}

template <typename Matrix, typename Vector3, enable_if_t<Matrix::IsMatrix && Matrix::Size == 4> = 0>
ENOKI_INLINE Matrix rotate(const Vector3 &axis, const entry_t<Matrix> &angle) {
    using Value = entry_t<Matrix>;
    using Vector4 = column_t<Matrix>;

    auto [sin_theta, cos_theta] = sincos(angle);
    Value cos_theta_m = 1.f - cos_theta;

    auto shuf1 = shuffle<1, 2, 0>(axis),
         shuf2 = shuffle<2, 0, 1>(axis),
         tmp0  = fmadd(axis * axis, cos_theta_m, cos_theta),
         tmp1  = fmadd(axis * shuf1, cos_theta_m, shuf2 * sin_theta),
         tmp2  = fmsub(axis * shuf2, cos_theta_m, shuf1 * sin_theta);

    return Matrix(
        Vector4(tmp0.x(), tmp1.x(), tmp2.x(), 0.f),
        Vector4(tmp2.y(), tmp0.y(), tmp1.y(), 0.f),
        Vector4(tmp1.z(), tmp2.z(), tmp0.z(), 0.f),
        Vector4(0.f, 0.f, 0.f, 1.f)
    );
}

template <typename Matrix>
ENOKI_INLINE Matrix perspective(const entry_t<Matrix> &fov,
                                const entry_t<Matrix> &near_,
                                const entry_t<Matrix> &far_,
                                const entry_t<Matrix> &aspect = 1.f) {
    static_assert(Matrix::Size == 4, "Matrix::perspective(): implementation assumes 4x4 matrix output");

    auto recip = rcp(near_ - far_);
    auto c = cot(.5f * fov);

    Matrix trafo = diag<Matrix>(
        column_t<Matrix>(c / aspect, c, (near_ + far_) * recip, 0.f));

    trafo(2, 3) = 2.f * near_ * far_ * recip;
    trafo(3, 2) = -1.f;

    return trafo;
}

template <typename Matrix>
ENOKI_INLINE Matrix frustum(const entry_t<Matrix> &left,
                            const entry_t<Matrix> &right,
                            const entry_t<Matrix> &bottom,
                            const entry_t<Matrix> &top,
                            const entry_t<Matrix> &near_,
                            const entry_t<Matrix> &far_) {
    static_assert(Matrix::Size == 4, "Matrix::frustum(): implementation assumes 4x4 matrix output");

    auto rl = rcp(right - left),
         tb = rcp(top - bottom),
         fn = rcp(far_ - near_);

    Matrix trafo = zero<Matrix>();
    trafo(0, 0) = (2.f * near_) * rl;
    trafo(1, 1) = (2.f * near_) * tb;
    trafo(0, 2) = (right + left) * rl;
    trafo(1, 2) = (top + bottom) * tb;
    trafo(2, 2) = -(far_ + near_) * fn;
    trafo(3, 2) = -1.f;
    trafo(2, 3) = -2.f * far_ * near_ * fn;

    return trafo;
}

template <typename Matrix>
ENOKI_INLINE Matrix ortho(const entry_t<Matrix> &left,
                          const entry_t<Matrix> &right,
                          const entry_t<Matrix> &bottom,
                          const entry_t<Matrix> &top,
                          const entry_t<Matrix> &near_,
                          const entry_t<Matrix> &far_) {
    static_assert(Matrix::Size == 4, "Matrix::ortho(): implementation assumes 4x4 matrix output");

    auto rl = rcp(right - left),
         tb = rcp(top - bottom),
         fn = rcp(far_ - near_);

    Matrix trafo = zero<Matrix>();

    trafo(0, 0) = 2.f * rl;
    trafo(1, 1) = 2.f * tb;
    trafo(2, 2) = -2.f * fn;
    trafo(3, 3) = 1.f;
    trafo(0, 3) = -(right + left) * rl;
    trafo(1, 3) = -(top + bottom) * tb;
    trafo(2, 3) = -(far_ + near_) * fn;

    return trafo;
}

template <typename Matrix, typename Point, typename Vector>
Matrix look_at(const Point &origin, const Point &target, const Vector &up) {
    static_assert(Matrix::Size == 4, "Matrix::look_at(): implementation "
                                     "assumes 4x4 matrix output");

    auto dir = normalize(target - origin);
    auto left = normalize(cross(dir, up));
    auto new_up = cross(left, dir);
    using Scalar = scalar_t<Matrix>;

    return Matrix(
        concat(left, Scalar(0)),
        concat(new_up, Scalar(0)),
        concat(-dir, Scalar(0)),
        column_t<Matrix>(
            -dot(left, origin),
            -dot(new_up, origin),
             dot(dir, origin),
             1.f
        )
    );
}

template <typename T,
          typename E       = expr_t<T>,
          typename Matrix3 = Matrix<E, 3>,
          typename Vector3 = Array<E, 3>,
          typename Quat    = Quaternion<E>>
std::tuple<Matrix3, Quat, Vector3> transform_decompose(const Matrix<T, 4> &A, size_t it = 10) {
    Matrix3 A_sub(A), Q, P;
    std::tie(Q, P) = polar_decomp(A_sub, it);

    if (ENOKI_UNLIKELY(any(enoki::isnan(Q(0, 0)))))
        Q = identity<Matrix3>();

    auto sign_q = det(Q);
    Q = mulsign(Array<Vector3, 3>(Q), sign_q);
    P = mulsign(Array<Vector3, 3>(P), sign_q);

    return std::make_tuple(
        P,
        matrix_to_quat(Q),
        head<3>(A.col(3))
    );
}

template <typename T,
          typename E = expr_t<T>,
          typename Matrix3 = Matrix<E, 3>,
          typename Matrix4 = Matrix<E, 4>,
          typename Vector3>
Matrix4 transform_compose(const Matrix<T, 3> &S,
                          const Quaternion<T> &q,
                          const Vector3 &t) {
    Matrix4 result = Matrix4(quat_to_matrix<Matrix3>(q) * S);
    result.coeff(3) = concat(t, scalar_t<Matrix4>(1));
    return result;
}

template <typename T,
          typename E = expr_t<T>,
          typename Matrix3 = Matrix<E, 3>,
          typename Matrix4 = Matrix<E, 4>,
          typename Vector3>
Matrix4 transform_compose_inverse(const Matrix<T, 3> &S,
                                  const Quaternion<T> &q,
                                  const Vector3 &t) {
    auto inv_m = inverse(quat_to_matrix<Matrix3>(q) * S);
    Matrix4 result = Matrix4(inv_m);
    result.coeff(3) = concat(inv_m * -t, scalar_t<Matrix4>(1));
    return result;
}

NAMESPACE_END(enoki)
