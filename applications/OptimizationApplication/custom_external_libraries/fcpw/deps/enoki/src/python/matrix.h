#pragma once
#include "common.h"
#include <enoki/matrix.h>

using Matrix2f  = Matrix<Float32 , 2>;
using Matrix2d  = Matrix<Float64 , 2>;
using Matrix3f  = Matrix<Float32 , 3>;
using Matrix3d  = Matrix<Float64 , 3>;
using Matrix4f  = Matrix<Float32 , 4>;
using Matrix4d  = Matrix<Float64 , 4>;
using Matrix44f = Matrix<Vector4f, 4>;
using Matrix41f = Matrix<Vector1f, 4>;
using Matrix44d = Matrix<Vector4d, 4>;
using Matrix41d = Matrix<Vector1d, 4>;
using Matrix2m  = mask_t<Matrix2f>;
using Matrix3m  = mask_t<Matrix3f>;
using Matrix4m  = mask_t<Matrix4f>;
using Matrix44m = mask_t<Matrix44f>;
using Matrix41m = mask_t<Matrix41f>;

using Matrix2fX  = Matrix<Float32X, 2>;
using Matrix2dX  = Matrix<Float64X, 2>;
using Matrix3fX  = Matrix<Float32X, 3>;
using Matrix3dX  = Matrix<Float64X, 3>;
using Matrix4fX  = Matrix<Float32X, 4>;
using Matrix4dX  = Matrix<Float64X, 4>;
using Matrix44fX = Matrix<Vector4fX, 4>;
using Matrix44dX = Matrix<Vector4dX, 4>;
using Matrix2mX  = mask_t<Matrix2fX>;
using Matrix3mX  = mask_t<Matrix3fX>;
using Matrix4mX  = mask_t<Matrix4fX>;
using Matrix44mX = mask_t<Matrix44fX>;

#if defined(ENOKI_CUDA)
using Matrix2fC  = Matrix<Float32C, 2>;
using Matrix2dC  = Matrix<Float64C, 2>;
using Matrix3fC  = Matrix<Float32C, 3>;
using Matrix3dC  = Matrix<Float64C, 3>;
using Matrix4fC  = Matrix<Float32C, 4>;
using Matrix4dC  = Matrix<Float64C, 4>;
using Matrix44fC = Matrix<Vector4fC , 4>;
using Matrix44dC = Matrix<Vector4dC , 4>;
using Matrix2mC  = mask_t<Matrix2fC>;
using Matrix3mC  = mask_t<Matrix3fC>;
using Matrix4mC  = mask_t<Matrix4fC>;
using Matrix44mC = mask_t<Matrix44fC>;
#endif

#if defined(ENOKI_AUTODIFF)
using Matrix2fD  = Matrix<Float32D, 2>;
using Matrix2dD  = Matrix<Float64D, 2>;
using Matrix3fD  = Matrix<Float32D, 3>;
using Matrix3dD  = Matrix<Float64D, 3>;
using Matrix4fD  = Matrix<Float32D, 4>;
using Matrix4dD  = Matrix<Float64D, 4>;
using Matrix44fD = Matrix<Vector4fD , 4>;
using Matrix44dD = Matrix<Vector4dD , 4>;
using Matrix2mD  = mask_t<Matrix2fD>;
using Matrix3mD  = mask_t<Matrix3fD>;
using Matrix4mD  = mask_t<Matrix4fD>;
using Matrix44mD = mask_t<Matrix44fD>;
#endif

template <typename Matrix>
py::class_<Matrix> bind_matrix(py::module &m, py::module &s, const char *name) {
    using Value  = typename Matrix::Entry;
    using Vector = typename Matrix::Column;
    using Array  = Array<Vector, Matrix::Size>;
    using Scalar = scalar_t<Matrix>;
    using Mask   = mask_t<Matrix>;

    static constexpr bool IsDiff    = is_diff_array_v<Matrix>;
    static constexpr bool IsDynamic = is_dynamic_v<Value>;

    auto cl = py::class_<Matrix>(s, name)
        .def(py::init<>())
        .def(py::init<const Value &>())
        .def(py::init([](const py::list &list) -> Matrix {
            size_t size = list.size();
            if (size != Matrix::Size)
                throw py::reference_cast_error();
            Matrix result;
            for (size_t i = 0; i < size; ++i)
                result[i] = py::cast<Vector>(list[i]);
            return transpose(result);
        }))
        .def(py::self == py::self)
        .def(py::self != py::self)
        .def("__repr__", [](const Matrix &a) -> std::string {
            if (*implicit_conversion)
                return "";
            std::ostringstream oss;
            oss << a;
            return oss.str();
        })
        .def("__getitem__", [](const Matrix &a, std::pair<size_t, size_t> index) {
            if (index.first >= Matrix::Size || index.second >= Matrix::Size)
                throw py::index_error();
            return a.coeff(index.second, index.first);
        })
        .def("__getitem__", [](const Matrix &a, size_t index) {
            if (index >= Matrix::Size)
                throw py::index_error();
            return a[index];
        })
        .def("__setitem__", [](Matrix &a, std::pair<size_t, size_t> index, const Value &value) {
            if (index.first >= Matrix::Size || index.second >= Matrix::Size)
                throw py::index_error();
            a.coeff(index.second, index.first) = value;
        })
        .def("__setitem__", [](Matrix &a, const mask_t<Value> &m, const Matrix &b) {
            a[m] = b;
        })
        .def(py::self - py::self)
        .def(Value() - py::self)
        .def(py::self + py::self)
        .def(Value() + py::self)
        .def(py::self * Value())
        .def(py::self / Value())
        .def("__lt__",
            [](const Matrix &a, const Matrix &b) -> Mask {
                return a < b;
            }, py::is_operator())
        .def("__lt__",
            [](const Value &a, const Matrix &b) -> Mask {
                return a < b;
            }, py::is_operator())
        .def("__gt__",
            [](const Matrix &a, const Matrix &b) -> Mask {
                return a > b;
            }, py::is_operator())
        .def("__gt__",
            [](const Value &a, const Matrix &b) -> Mask {
                return a > b;
            }, py::is_operator())
        .def("__le__",
            [](const Matrix &a, const Matrix &b) -> Mask {
                return a <= b;
            }, py::is_operator())
        .def("__le__",
            [](const Value &a, const Matrix &b) -> Mask {
                return a <= b;
            }, py::is_operator())
        .def("__ge__",
            [](const Matrix &a, const Matrix &b) -> Mask {
                return a >= b;
            }, py::is_operator())
        .def("__ge__",
            [](const Value &a, const Matrix &b) -> Mask {
                return a >= b;
            }, py::is_operator())
        .def(-py::self)
        .def("__mul__", [](const Matrix &a, const Matrix &b) {
            return Matrix(Array(a) * Array(b));
        });

    cl.def("__matmul__",
            [](const Matrix &a, const Matrix &b) { return a * b; });

    if constexpr (array_depth_v<Value> == (IsDynamic ? 1 : 0))
        cl.def("__matmul__",
                [](const Matrix &a, const Vector &b) { return a * b; });

    cl.def_static("identity", [](size_t size) { return identity<Matrix>(size); }, "size"_a = 1);
    cl.def_static("zero", [](size_t size) { return zero<Matrix>(size); }, "size"_a = 1);
    cl.def_static("full", [](Scalar value, size_t size) { return full<Matrix>(value, size); },
                  "value"_a, "size"_a = 1);

    cl.def("__len__", [](const Matrix &m) { return Matrix::Size; });

    cl.def(py::init([](const py::ndarray &obj) -> Matrix {
            using T = expr_t<decltype(detach(std::declval<Matrix &>()))>;
            return numpy_to_enoki<T>(obj);
        }))
    .def("numpy", [](const Matrix &a, bool eval) { return enoki_to_numpy(detach(a), eval); },
        "eval"_a = true)
    .def_property_readonly("__array_interface__", [](const py::object &o) {
        py::object np_array = o.attr("numpy")();
        py::dict result;
        result["data"]    = np_array;
        result["shape"]   = np_array.attr("shape");
        result["version"] = 3;
        char typestr[4]   = { '<', 0, '0' + sizeof(Scalar), 0 };
        if (std::is_floating_point_v<Scalar>)
            typestr[1] = 'f';
        else if (std::is_same_v<Scalar, bool>)
            typestr[1] = 'b';
        else if (std::is_unsigned_v<Scalar>)
            typestr[1] = 'u';
        else
            typestr[1] = 'i';
        result["typestr"] = typestr;
        return result;
    });

    cl.def(py::init([](const py::torch_tensor &obj) -> Matrix {
        using T = expr_t<decltype(detach(std::declval<Matrix&>()))>;
        return torch_to_enoki<T>(obj);
    }))
    .def("torch", [](const Matrix &a, bool eval) {
        return enoki_to_torch(detach(a), eval); },
        "eval"_a = true
    );

    m.def("shape", [](const Matrix &a) { return shape(a); });
    m.def("slices", [](const Matrix &a) { return slices(a); });
    m.def("set_slices", [](Matrix &a, size_t size) { set_slices(a, size); }, "array"_a, "slices"_a);

    if constexpr (Matrix::Size == 2) {
        cl.def(py::init<const Value &, const Value &,
                        const Value &, const Value &>());
        if constexpr (array_depth_v<Value> == (IsDynamic ? 1 : 0))
            cl.def(py::init<const Vector &, const Vector &>());
    } else if constexpr (Matrix::Size == 3) {
        cl.def(py::init<const Value &, const Value &, const Value &,
                        const Value &, const Value &, const Value &,
                        const Value &, const Value &, const Value &>());
        if constexpr (array_depth_v<Value> == (IsDynamic ? 1 : 0))
            cl.def(py::init<const Vector &, const Vector &, const Vector &>());
    } else if constexpr (Matrix::Size == 4) {
        using Vector3 = enoki::Array<Value, 3>;

        cl.def(py::init<const Value &, const Value &, const Value &, const Value &,
                        const Value &, const Value &, const Value &, const Value &,
                        const Value &, const Value &, const Value &, const Value &,
                        const Value &, const Value &, const Value &, const Value &>());

        if constexpr (array_depth_v<Value> == (IsDynamic ? 1 : 0)) {
            cl.def(py::init<const Vector &, const Vector &, const Vector &, const Vector &>())
            .def_static("translate", [](const Vector3 &v) { return translate<Matrix>(v); })
            .def_static("scale", [](const Vector3 &v) { return scale<Matrix>(v); })
            .def_static("rotate", [](const Vector3 &axis, const Value &angle) {
                    return rotate<Matrix>(axis, angle);
                },
                "axis"_a, "angle"_a
            )
            .def_static("look_at", [](const Vector3 &origin,
                                        const Vector3 &target,
                                        const Vector3 &up) {
                    return look_at<Matrix>(origin, target, up);
                },
                "origin"_a, "target"_a, "up"_a
            );
        }
    }

    auto transpose_m = [](const Matrix &m) { return transpose(m); };
    cl.def_property_readonly("T", transpose_m);

    m.def("transpose", transpose_m);
    m.def("det", [](const Matrix &m) { return det(m); });
    m.def("inverse", [](const Matrix &m) { return inverse(m); });
    m.def("inverse_transpose", [](const Matrix &m) { return inverse_transpose(m); });

    m.def("abs", [](const Matrix &m) { return enoki::abs(m); });

    m.def("isfinite", [](const Matrix &m) -> Mask { return enoki::isfinite(m); });
    m.def("isnan",    [](const Matrix &m) -> Mask { return enoki::isnan(m); });
    m.def("isinf",    [](const Matrix &m) -> Mask { return enoki::isinf(m); });

    m.def("eq",  [](const Matrix &a, const Matrix &b) -> Mask { return eq(a, b); });
    m.def("neq", [](const Matrix &a, const Matrix &b) -> Mask { return neq(a, b); });


    if constexpr (IsDiff) {
        using Detached = expr_t<decltype(detach(std::declval<Matrix&>()))>;

        m.def("detach", [](const Matrix &a) -> Detached { return detach(a); });
        m.def("requires_gradient",
              [](const Matrix &a) { return requires_gradient(a); },
              "array"_a);

        m.def("set_requires_gradient",
              [](Matrix &a, bool value) { set_requires_gradient(a, value); },
              "array"_a, "value"_a = true);

        m.def("gradient", [](Matrix &a) { return eval(gradient(a)); });
        m.def("set_gradient",
              [](Matrix &a, const Detached &g, bool b) { set_gradient(a, g, b); },
              "array"_a, "gradient"_a, "backward"_a = true);

        m.def("graphviz", [](const Matrix &a) { return graphviz(a); });

        m.def("set_label", [](const Matrix &a, const char *label) {
            set_label(a, label);
        });
    }

    register_implicit_casts<Matrix, Value>();

    return cl;
}


template <typename Mask>
py::class_<Mask> bind_matrix_mask(py::module &m, py::module &s, const char *name) {
    auto cl = py::class_<Mask>(s, name)
        .def(py::init<>())
        .def(py::self == py::self)
        .def(py::self != py::self)
        .def(py::self | py::self)
        .def(py::self & py::self)
        .def(py::self ^ py::self)
        .def(!py::self)
        .def(~py::self)
        .def("__repr__", [](const Mask &a) -> std::string {
            if (*implicit_conversion)
                return "";
            std::ostringstream oss;
            oss << a;
            return oss.str();
        });

    m.def("shape", [](const Mask &a) { return shape(a); });
    m.def("slices", [](const Mask &a) { return slices(a); });
    m.def("set_slices", [](Mask &a, size_t size) { set_slices(a, size); }, "array"_a, "slices"_a);

    m.def("any",  [](const Mask &m) { return enoki::any(m); });
    m.def("none", [](const Mask &m) { return enoki::none(m); });
    m.def("all",  [](const Mask &m) { return enoki::all(m); });

    m.def("any_nested",   [](const Mask &m) { return enoki::any_nested(m); });
    m.def("none_nested",  [](const Mask &m) { return enoki::none_nested(m); });
    m.def("all_nested",   [](const Mask &m) { return enoki::all_nested(m); });

    m.def("eq",  [](const Mask &a, const Mask &b) -> Mask { return eq(a, b); });
    m.def("neq", [](const Mask &a, const Mask &b) -> Mask { return neq(a, b); });

    return cl;
}
