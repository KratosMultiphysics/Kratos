#pragma once

#if defined(ENOKI_CUDA)
#include <enoki/cuda.h>
#endif
#if defined(ENOKI_AUTODIFF)
#include <enoki/autodiff.h>
#endif
#include <enoki/dynamic.h>
#include <enoki/transform.h>
#include <enoki/special.h>
#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <sstream>

using namespace enoki;
namespace py = pybind11;
using namespace py::literals;

using Float32 = float;
using Float64 = double;
using Int32   = int32_t;
using UInt32  = uint32_t;
using Int64   = int64_t;
using UInt64  = uint64_t;

using Vector0f  = Array<Float32 , 0>;
using Vector0d  = Array<Float64 , 0>;
using Vector0i  = Array<Int32 , 0>;
using Vector0u  = Array<UInt32 , 0>;
using Vector0m  = mask_t<Vector0f >;
using Vector1f  = Array<Float32 , 1>;
using Vector1d  = Array<Float64 , 1>;
using Vector1i  = Array<Int32 , 1>;
using Vector1u  = Array<UInt32 , 1>;
using Vector1m  = mask_t<Vector1f >;
using Vector2f  = Array<Float32 , 2>;
using Vector2d  = Array<Float64 , 2>;
using Vector2i  = Array<Int32 , 2>;
using Vector2u  = Array<UInt32 , 2>;
using Vector2m  = mask_t<Vector2f >;
using Vector3f  = Array<Float32 , 3>;
using Vector3d  = Array<Float64 , 3>;
using Vector3i  = Array<Int32 , 3>;
using Vector3u  = Array<UInt32 , 3>;
using Vector3m  = mask_t<Vector3f >;
using Vector4f  = Array<Float32 , 4>;
using Vector4d  = Array<Float64 , 4>;
using Vector4i  = Array<Int32 , 4>;
using Vector4u  = Array<UInt32 , 4>;
using Vector4m  = mask_t<Vector4f >;
using Complex2f   = Complex<Float32 >;
using Complex2d   = Complex<Float64 >;
using Complex24f  = Complex<Vector4f >;
using Complex24d  = Complex<Vector4d >;
using Quaternion4f = Quaternion<float>;
using Quaternion4d = Quaternion<double>;

using Float32X = DynamicArray<Packet<Float32>>;
using Float64X = DynamicArray<Packet<Float64>>;
using Int32X   = DynamicArray<Packet<Int32>>;
using Int64X   = DynamicArray<Packet<Int64>>;
using UInt32X  = DynamicArray<Packet<UInt32>>;
using UInt64X  = DynamicArray<Packet<UInt64>>;
using MaskX    = mask_t<Float32X>;
using Mask64X  = mask_t<Float64X>;

using Vector0fX = Array<Float32X, 0>;
using Vector0dX = Array<Float64X, 0>;
using Vector0iX = Array<Int32X, 0>;
using Vector0uX = Array<UInt32X, 0>;
using Vector0mX = mask_t<Vector0fX>;
using Vector1fX = Array<Float32X, 1>;
using Vector1dX = Array<Float64X, 1>;
using Vector1iX = Array<Int32X, 1>;
using Vector1uX = Array<UInt32X, 1>;
using Vector1mX = mask_t<Vector1fX>;
using Vector2fX = Array<Float32X, 2>;
using Vector2dX = Array<Float64X, 2>;
using Vector2iX = Array<Int32X, 2>;
using Vector2uX = Array<UInt32X, 2>;
using Vector2mX = mask_t<Vector2fX>;
using Vector3fX = Array<Float32X, 3>;
using Vector3dX = Array<Float64X, 3>;
using Vector3iX = Array<Int32X, 3>;
using Vector3uX = Array<UInt32X, 3>;
using Vector3mX = mask_t<Vector3fX>;
using Vector4fX = Array<Float32X, 4>;
using Vector4dX = Array<Float64X, 4>;
using Vector4iX = Array<Int32X, 4>;
using Vector4uX = Array<UInt32X, 4>;
using Vector4mX = mask_t<Vector4fX>;

using Complex2fX  = Complex<Float32X>;
using Complex2dX  = Complex<Float64X>;
using Complex24fX = Complex<Vector4fX>;
using Complex24dX = Complex<Vector4dX>;

#if defined(ENOKI_CUDA)
using Float32C = CUDAArray<Float32>;
using Float64C = CUDAArray<Float64>;
using Int32C   = CUDAArray<Int32>;
using Int64C   = CUDAArray<Int64>;
using UInt32C  = CUDAArray<UInt32>;
using UInt64C  = CUDAArray<UInt64>;
using MaskC    = mask_t<Float32C>;

using Vector0fC = Array<Float32C, 0>;
using Vector0dC = Array<Float64C, 0>;
using Vector0iC = Array<Int32C, 0>;
using Vector0uC = Array<UInt32C, 0>;
using Vector0mC = mask_t<Vector0fC>;
using Vector1fC = Array<Float32C, 1>;
using Vector1dC = Array<Float64C, 1>;
using Vector1iC = Array<Int32C, 1>;
using Vector1uC = Array<UInt32C, 1>;
using Vector1mC = mask_t<Vector1fC>;
using Vector2fC = Array<Float32C, 2>;
using Vector2dC = Array<Float64C, 2>;
using Vector2iC = Array<Int32C, 2>;
using Vector2uC = Array<UInt32C, 2>;
using Vector2mC = mask_t<Vector2fC>;
using Vector3fC = Array<Float32C, 3>;
using Vector3dC = Array<Float64C, 3>;
using Vector3iC = Array<Int32C, 3>;
using Vector3uC = Array<UInt32C, 3>;
using Vector3mC = mask_t<Vector3fC>;
using Vector4fC = Array<Float32C, 4>;
using Vector4dC = Array<Float64C, 4>;
using Vector4iC = Array<Int32C, 4>;
using Vector4uC = Array<UInt32C, 4>;
using Vector4mC = mask_t<Vector4fC>;

using Complex2fC  = Complex<Float32C>;
using Complex2dC  = Complex<Float64C>;
using Complex24fC = Complex<Vector4fC>;
using Complex24dC = Complex<Vector4dC>;
#endif

#if defined(ENOKI_AUTODIFF)
using Float32D = DiffArray<Float32C>;
using Float64D = DiffArray<Float64C>;
using Int32D   = DiffArray<Int32C>;
using Int64D   = DiffArray<Int64C>;
using UInt32D  = DiffArray<UInt32C>;
using UInt64D  = DiffArray<UInt64C>;
using MaskD    = mask_t<Float32D>;

using Vector0fD = Array<Float32D, 0>;
using Vector0dD = Array<Float64D, 0>;
using Vector0iD = Array<Int32D, 0>;
using Vector0uD = Array<UInt32D, 0>;
using Vector0mD = mask_t<Vector0fD>;
using Vector1fD = Array<Float32D, 1>;
using Vector1dD = Array<Float64D, 1>;
using Vector1iD = Array<Int32D, 1>;
using Vector1uD = Array<UInt32D, 1>;
using Vector1mD = mask_t<Vector1fD>;
using Vector2fD = Array<Float32D, 2>;
using Vector2dD = Array<Float64D, 2>;
using Vector2iD = Array<Int32D, 2>;
using Vector2uD = Array<UInt32D, 2>;
using Vector2mD = mask_t<Vector2fD>;
using Vector3fD = Array<Float32D, 3>;
using Vector3dD = Array<Float64D, 3>;
using Vector3iD = Array<Int32D, 3>;
using Vector3uD = Array<UInt32D, 3>;
using Vector3mD = mask_t<Vector3fD>;
using Vector4fD = Array<Float32D, 4>;
using Vector4dD = Array<Float64D, 4>;
using Vector4iD = Array<Int32D, 4>;
using Vector4uD = Array<UInt32D, 4>;
using Vector4mD = mask_t<Vector4fD>;

using Complex2fD  = Complex<Float32D>;
using Complex2dD  = Complex<Float64D>;
using Complex24fD = Complex<Vector4fD>;
using Complex24dD = Complex<Vector4dD>;
#endif

namespace pybind11 {
    inline bool NDArray_Check(PyObject *obj) {
        return strcmp(obj->ob_type->tp_name, "numpy.ndarray") == 0;
    }

    inline bool TorchTensor_Check(PyObject *obj) {
        return strcmp(obj->ob_type->tp_name, "Tensor") == 0;
    }

    class ndarray : public object {
    public:
        PYBIND11_OBJECT_DEFAULT(ndarray, object, NDArray_Check)
    };

    class torch_tensor : public object {
    public:
        PYBIND11_OBJECT_DEFAULT(torch_tensor, object, TorchTensor_Check)
    };
};

template <bool IsCUDA> struct Buffer {
    Buffer(size_t size) {
        if constexpr (IsCUDA) {
            #if defined(ENOKI_CUDA)
                ptr = cuda_managed_malloc(size);
            #endif
        } else {
            #if defined(__APPLE__)
                if (posix_memalign(&ptr, 64, size))
                    throw std::runtime_error("Buffer: allocation failure!");
            #else
                #if defined(_MSC_VER)
                    ptr = _aligned_malloc(size, 64);
                #else
                    ptr = std::aligned_alloc(64, size);
                #endif
            #endif
        }
    }

    ~Buffer() {
        if constexpr (IsCUDA) {
            #if defined(ENOKI_CUDA)
                cuda_free(ptr);
            #endif
        } else {
            #if defined(_MSC_VER)
                _aligned_free(ptr);
            #else
                std::free(ptr);
            #endif
        }
    }

    void *ptr = nullptr;
    bool cuda_managed = false;
};

template <typename Array> py::object enoki_to_torch(const Array &array, bool eval);
template <typename Array> py::object enoki_to_numpy(const Array &array, bool eval);
template <typename Array> Array torch_to_enoki(py::object src);
template <typename Array> Array numpy_to_enoki(py::array src);

extern bool *implicit_conversion;

struct set_flag {
    bool &flag, backup;
    set_flag(bool &flag_, bool value = true) : flag(flag_), backup(flag_) { flag = value; }
    ~set_flag() { flag = backup; }
};

/// Customized version of pybind11::implicitly_convertible() which disables
/// __repr__ during implicit casts (this can be triggered at implicit cast
/// failures and causes a totally unnecessary/undesired cuda_eval() invocation)
template <typename InputType, typename OutputType> void implicitly_convertible() {
    auto implicit_caster = [](PyObject *obj, PyTypeObject *type) -> PyObject * {
        if (*implicit_conversion) return nullptr;
        set_flag flag_helper(*implicit_conversion);
        if (!py::detail::make_caster<InputType>().load(obj, std::is_scalar_v<InputType>)) {
            return nullptr;
        }
        py::tuple args(1);
        args[0] = obj;
        PyObject *result = PyObject_CallObject((PyObject *) type, args.ptr());
        if (result == nullptr)
            PyErr_Clear();
        return result;
    };

    if (auto tinfo = py::detail::get_type_info(typeid(OutputType)))
        tinfo->implicit_conversions.push_back(implicit_caster);
    else
        py::pybind11_fail("implicitly_convertible: Unable to find type " +
                          pybind11::type_id<OutputType>());
}

template <typename Array, typename Value> void register_implicit_casts() {
    using Scalar = std::conditional_t<!is_mask_v<Array>, scalar_t<Array>, bool>;
    static constexpr size_t Size = array_size_v<Array>;

    auto implicit_caster = [](PyObject *obj, PyTypeObject *type) -> PyObject * {
        const char *tp_name_src = obj->ob_type->tp_name,
                    *tp_name_dst = type->tp_name;

        if (*implicit_conversion) // limit nesting of implicit conversions
            return nullptr;
        set_flag flag_helper(*implicit_conversion);

        bool pass = false;

        if (PyList_CheckExact(obj)) {
            pass = Size == Dynamic || Size == PyList_GET_SIZE(obj);
        } else if (PyTuple_CheckExact(obj)) {
            pass = Size == Dynamic || Size == PyTuple_GET_SIZE(obj);
        } else if (Size != 0 && PyNumber_Check(obj)) {
            pass = true;
        } else if (Size != 0 &&
                    (strcmp(tp_name_src, "numpy.ndarray") == 0 ||
                    strcmp(tp_name_src, "Tensor") == 0)) {
            pass = true;
        } else {
            // Convert from a different vector type (Vector4f -> Vector4fX)
            if (strncmp(tp_name_src, "enoki.", 6) == 0 &&
                strncmp(tp_name_dst, "enoki.", 6) == 0) {
                const char *dot_src = strchr(tp_name_src + 6, '.'),
                            *dot_dst = strchr(tp_name_dst + 6, '.');

                if (dot_src && dot_dst)
                    pass |= strcmp(dot_src, dot_dst) == 0;
            }

            if constexpr (!std::is_same_v<Scalar, Value>) {
                auto tinfo = py::detail::get_global_type_info(typeid(Value));
                if (tinfo)
                    pass |= strcmp(tp_name_src, tinfo->type->tp_name) == 0;
            }
        }

        if (!pass)
            return nullptr;

        PyObject *args = PyTuple_New(1);
        Py_INCREF(obj);
        PyTuple_SET_ITEM(args, 0, obj);
        PyObject *result = PyObject_CallObject((PyObject *) type, args);
        if (result == nullptr)
            PyErr_Clear();
        Py_DECREF(args);
        return result;
    };

    auto tinfo = py::detail::get_type_info(typeid(Array));
    tinfo->implicit_conversions.push_back(implicit_caster);
}

template <typename Array>
py::class_<Array> bind(py::module &m, py::module &s, const char *name) {
    using Scalar = std::conditional_t<
        !is_mask_v<Array>, scalar_t<Array>, bool>;

    using Value  = std::conditional_t<
        is_mask_v<Array> && array_depth_v<Array> == 1,
        bool, value_t<Array>>;

    using Mask = std::conditional_t<
        !is_mask_v<Array>,
        mask_t<float32_array_t<Array>>,
        Array
    >;

    static constexpr bool IsMask    = is_mask_v<Array>;
    static constexpr bool IsFloat   = std::is_floating_point_v<Scalar>;
    static constexpr bool IsCUDA    = is_cuda_array_v<Array>;
    static constexpr bool IsDiff    = is_diff_array_v<Array>;
    static constexpr bool IsDynamic = is_dynamic_v<Array>;
    static constexpr bool IsKMask   = IsMask && !is_cuda_array_v<Array>;
    static constexpr size_t Size    = array_size_v<Array>;

    py::class_<Array> cl(s, name);

    cl.def(py::init<>())
      .def(py::self == py::self)
      .def(py::self != py::self)
      .def("managed", py::overload_cast<>(&Array::managed), py::return_value_policy::reference)
      .def("eval", py::overload_cast<>(&Array::eval), py::return_value_policy::reference)
      .def("__repr__", [](const Array &a) -> std::string {
          if (*implicit_conversion)
              return "";
          std::ostringstream oss;
          oss << a;
          return oss.str();
      })
      .def_static("zero", [](size_t size) { return zero<Array>(size); },
                  "size"_a = 1)
      .def_static("empty", [](size_t size) { return empty<Array>(size); },
                  "size"_a = 1)
      .def_static("full", [](Scalar value, size_t size) { return full<Array>(value, size); },
                  "value"_a, "size"_a = 1);

    cl.def(py::init<const Scalar &>());
    if constexpr (!std::is_same_v<Value, Scalar>)
        cl.def(py::init<const Value &>());
    cl.def(py::init<const Array &>());

    if constexpr (array_depth_v<Array> == 1 && is_cuda_array_v<Array>) {
        cl.def(
            py::init([](Scalar scalar, bool literal) -> Array {
                if (literal)
                    return Array(scalar);
                else
                    return Array::copy(&scalar, 1);
            }),
            "scalar"_a, "literal"_a = true
        );
    }

    if constexpr (!IsKMask) {
        cl.def(py::init([](const py::ndarray &obj) -> Array {
                   using T = expr_t<decltype(detach(std::declval<Array &>()))>;
                   return numpy_to_enoki<T>(obj);
               }))
          .def("numpy", [](const Array &a, bool eval) { return enoki_to_numpy(detach(a), eval); },
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
        cl.def(py::init([](const py::torch_tensor &obj) -> Array {
            using T = expr_t<decltype(detach(std::declval<Array&>()))>;
            return torch_to_enoki<T>(obj);
        }))
        .def("torch", [](const Array &a, bool eval) {
            return enoki_to_torch(detach(a), eval); },
            "eval"_a = true
        );
    }

    cl.def(py::init([](const py::list &list) -> Array {
        size_t size = list.size();

        if constexpr (IsDynamic && IsKMask && array_depth_v<Array> == 1) {
            using IntArray = replace_scalar_t<typename Array::ArrayType, int>;
            std::unique_ptr<int[]> result(new int[size]);
            for (size_t i = 0; i < size; ++i)
                result[i] = py::cast<int>(list[i]);
            return eq(IntArray::copy(result.get(), size), 1);
        } else if constexpr (IsDynamic && array_depth_v<Array> == 1) {
            std::unique_ptr<Value[]> result(new Value[size]);
            for (size_t i = 0; i < size; ++i)
                result[i] = py::cast<Value>(list[i]);
            return Array::copy(result.get(), size);
        } else {
            if (size != Array::Size)
                throw py::reference_cast_error();

            // allow nested implicit conversions
            set_flag flag_helper(*implicit_conversion, false);
            Array result;
            for (size_t i = 0; i < size; ++i)
                result[i] = py::cast<Value>(list[i]);
            return result;
        }
    }));

    cl.def(py::init([](const py::tuple &tuple) -> Array {
        size_t size = tuple.size();

        if constexpr (IsDynamic && array_depth_v<Array> == 1) {
            std::unique_ptr<Value[]> result(new Value[size]);
            for (size_t i = 0; i < size; ++i)
                result[i] = py::cast<Value>(tuple[i]);
            return Array::copy(result.get(), size);
        } else {
            if (size != Array::Size)
                throw py::reference_cast_error();

            // allow nested implicit conversions
            set_flag flag_helper(*implicit_conversion, false);
            Array result;
            for (size_t i = 0; i < size; ++i)
                result[i] = py::cast<Value>(tuple[i]);
            return result;
        }
    }));

    if constexpr (!IsMask) {
        cl.def(py::self + py::self)
          .def(Value() + py::self)
          .def(py::self - py::self)
          .def(Value() - py::self)
          .def(py::self * py::self)
          .def(Value() * py::self)
          .def("__lt__",
               [](const Array &a, const Array &b) -> Mask {
                   return a < b;
               }, py::is_operator())
          .def("__lt__",
               [](const Value &a, const Array &b) -> Mask {
                   return a < b;
               }, py::is_operator())
          .def("__gt__",
               [](const Array &a, const Array &b) -> Mask {
                   return a > b;
               }, py::is_operator())
          .def("__gt__",
               [](const Value &a, const Array &b) -> Mask {
                   return a > b;
               }, py::is_operator())
          .def("__le__",
               [](const Array &a, const Array &b) -> Mask {
                   return a <= b;
               }, py::is_operator())
          .def("__le__",
               [](const Value &a, const Array &b) -> Mask {
                   return a <= b;
               }, py::is_operator())
          .def("__ge__",
               [](const Array &a, const Array &b) -> Mask {
                   return a >= b;
               }, py::is_operator())
          .def("__ge__",
               [](const Value &a, const Array &b) -> Mask {
                   return a >= b;
               }, py::is_operator())
          .def(-py::self);

        if constexpr (std::is_integral_v<Scalar>) {
            cl.def("__floordiv__", [](const Array &a1, Scalar a2) {
                  return (a2 == 1) ? a1 : (a1 / a2);
              })
              .def("__floordiv__", [](const Array &a1, const Array &a2) {
                  return a1 / a2;
              })
              .def("__mod__", [](const Array &a1, Scalar a2) {
                  return (a2 == 1) ? 1 : (a1 % a2);
              })
              .def("__mod__", [](const Array &a1, const Array &a2) {
                  return a1 % a2;
              });
        } else {
            cl.def("__truediv__", [](const Array &a1, Scalar a2) {
                  return a1 / a2;
               })
              .def("__truediv__", [](const Array &a1, const Array &a2) {
                  return a1 / a2;
               });
        }
    }

    if constexpr (!IsFloat) {
        cl.def(py::self | py::self)
          .def(py::self & py::self)
          .def(py::self ^ py::self)
          .def(!py::self)
          .def(~py::self);
    }

    if constexpr (!IsMask && array_depth_v<Array> == 1 && array_size_v<Array> == -1) {
        if (IsFloat)
            cl.def_static("linspace",
                  [](Scalar min, Scalar max, size_t size) {
                      return linspace<Array>(min, max, size);
                  },
                  "min"_a, "max"_a, "size"_a);

        cl.def_static("arange",
              [](size_t size) { return arange<Array>(size); }, "size"_a);

        cl.def_static("arange",
                      [](size_t start, size_t end, size_t step) {
                          return arange<Array>(start, end, step);
                      },
                      "start"_a, "end"_a, "step"_a = 1);
    }

    cl.def("__getitem__", [](const Array &a, size_t index) -> Value {
        if (index >= a.size())
            throw py::index_error();
        return a.coeff(index);
    }, "index"_a);

    if constexpr (array_depth_v<Array> == 1 && IsDynamic && !IsKMask) {
        cl.def("__getitem__", [](const Array &s, py::slice slice) {
            ssize_t start, stop, step, slicelength;
            if (!slice.compute(s.size(), &start, &stop, &step, &slicelength))
                throw py::error_already_set();

            if (slicelength == 0)
                return Array();
            else if (step == 1 && slicelength == (ssize_t) s.size())
                return s; // Fast path

            using Int32 = int32_array_t<Array>;
            Int32 indices =
                arange<Int32>((uint32_t) slicelength) * (uint32_t) step +
                (uint32_t) start;

            return gather<Array>(s, indices);
        }, "slice"_a);

        cl.def("__setitem__", [](Array &dst, py::slice slice,
                                 const Array &src) {
            ssize_t start, stop, step, slicelength;
            if (!slice.compute(dst.size(), &start, &stop, &step, &slicelength))
                throw py::error_already_set();

            if (step == 1 && slicelength == (ssize_t) dst.size()) {
                dst = src; // Fast path
            } else {
                if (slicelength != src.size() && src.size() != 1)
                    throw py::index_error(
                        "Size mismatch: tried to assign an array of size " +
                        std::to_string(src.size()) + " to a slice of size " +
                        std::to_string(slicelength) + "!");

                if (step == 0)
                    return;

                using Int32 = int32_array_t<Array>;
                Int32 indices =
                    arange<Int32>((int32_t) slicelength) * (int32_t) step +
                    (int32_t) start;

                scatter(dst, src, indices);
            }
        }, "src"_a, "value"_a);
    } else if constexpr (!IsKMask) {
        cl.def("__getitem__", [](const Array &src, py::slice slice) {
            ssize_t start, stop, step, slicelength;
            if (!slice.compute(src.size(), &start, &stop, &step, &slicelength))
                throw py::error_already_set();

            py::list result;
            for (ssize_t i = 0, j = start; j < stop; ++i, j += step)
                result.append(py::cast(src[j], py::return_value_policy::copy));
            return result;
        }, "slice"_a);

        cl.def("__setitem__", [](Array &dst, py::slice slice, py::object src) {
            ssize_t start, stop, step, slicelength, src_size = -1;
            if (!slice.compute(dst.size(), &start, &stop, &step, &slicelength))
                throw py::error_already_set();

            py::list list;
            if (py::isinstance<py::list>(src)) {
                list = py::list(src);
                src_size = (ssize_t) list.size();
            }

            if (slicelength != src_size && src_size != -1)
                throw py::index_error(
                    "Size mismatch: tried to assign an array of size " +
                    std::to_string(src_size) + " to a slice of size " +
                    std::to_string(slicelength) + "!");

            try {
                if (src_size == -1) {
                    Value v = py::cast<Value>(src);
                    for (ssize_t j = start; j < stop; j += step)
                        dst[j] = v;
                } else {
                    for (ssize_t i = 0, j = start; j < stop; ++i, j += step)
                        dst[j] = py::cast<Value>(list[i]);
                }
            } catch (const py::cast_error &) {
                throw py::reference_cast_error();
            }
        }, "slice"_a, "value"_a);
    }

    struct Iterator {
        Array &array;
        size_t pos = 0;

        Iterator(Array &array) : array(array), pos(0) { }
    };

    py::class_<Iterator>(cl, "Iterator")
        .def("__iter__",
            [](Iterator &it) -> Iterator & { return it; }
        )
        .def("__next__",
            [](Iterator &it) -> Value {
                if (it.pos >= it.array.size())
                    throw py::stop_iteration();
                return it.array[it.pos++];
            }
        );

    cl.def("__iter__",
           [](Array &array) { return Iterator(array); },
           py::keep_alive<0, 1>());

    cl.def("__len__", [](const Array &a) { return a.size(); });

    m.def("shape", [](const Array &a) { return shape(a); });
    m.def("slices", [](const Array &a) { return slices(a); });
    m.def("set_slices", [](Array &a, size_t size) { set_slices(a, size); }, "array"_a, "size"_a);

    if constexpr (IsDynamic)
        cl.def("resize", [](Array &a, size_t size) { a.resize(size); });

    cl.def("__getitem__", [](Array &a, const Mask &mask) {
        return select(mask, a, Value(0.f));
    }, "mask"_a);

    cl.def("__setitem__", [](Array &a, const Mask &mask, const Array &value) {
        a[mask] = value;
    }, "mask"_a, "value"_a);

    if constexpr (Size != Dynamic) {
        cl.def("__setitem__", [](Array &a, size_t index, const Value &value) {
            if (index >= Array::Size)
                throw py::index_error();
            a.coeff(index) = value;
        }, "index"_a, "value"_a);

        if constexpr (!IsMask || array_depth_v<Array> > 1) {
            if constexpr (Size == 2)
                cl.def(py::init<Value, Value>(), "x"_a, "y"_a);
            else if constexpr (Size == 3)
                cl.def(py::init<Value, Value, Value>(), "x"_a, "y"_a, "z"_a);
            else if constexpr (Size == 4)
                cl.def(py::init<Value, Value, Value, Value>(), "x"_a, "y"_a, "z"_a, "w"_a);
        }

        if constexpr (Size >= 1)
            cl.def_property("x", [](const Array &a) { return a.x(); },
                                 [](Array &a, const Value &v) { a.x() = v; });
        if constexpr (Size >= 2)
            cl.def_property("y", [](const Array &a) { return a.y(); },
                                 [](Array &a, const Value &v) { a.y() = v; });
        if constexpr (Size >= 3)
            cl.def_property("z", [](const Array &a) { return a.z(); },
                                 [](Array &a, const Value &v) { a.z() = v; });
        if constexpr (Size >= 4)
            cl.def_property("w", [](const Array &a) { return a.w(); },
                                 [](Array &a, const Value &v) { a.w() = v; });
    } else {
        cl.def_property_readonly("data", [](const Array &a) { return (uintptr_t) a.data(); });

        if constexpr (IsCUDA || IsDiff)
            cl.def_property_readonly("index", [](const Array &a) { return a.index_(); });
    }

    if constexpr (!IsMask) {
        cl.def("__matmul__", [](const Array &a, const Array &b) { return enoki::dot(a, b); });
        m.def("dot", [](const Array &a, const Array &b) { return enoki::dot(a, b); });
        m.def("abs_dot", [](const Array &a, const Array &b) { return enoki::abs_dot(a, b); });
        m.def("normalize", [](const Array &a) { return enoki::normalize(a); });
        m.def("squared_norm", [](const Array &a) { return enoki::squared_norm(a); });
        m.def("norm", [](const Array &a) { return enoki::norm(a); });

        if constexpr (Size == 3)
            m.def("cross", [](const Array &a, const Array &b) { return enoki::cross(a, b); });
    }

    using Index = array_t<uint32_array_t<
        std::conditional_t<array_depth_v<Array> == 1, array_t<Array>, value_t<array_t<Array>>>>>;

    using IndexMask = mask_t<float32_array_t<Index>>;

    // Scatter/gather currently not supported for dynamic CPU arrays containing masks
    if constexpr (IsDynamic && !IsKMask) {
        m.def("gather",
              [](const Array &source, const Index &index, const IndexMask &mask) {
                  return gather<Array>(source, index, mask);
              },
              "source"_a, "index"_a, "mask"_a = true);

        m.def("scatter",
              [](Array &target, const Array &source,
                 const Index &index,
                 const IndexMask &mask) { scatter(target, source, index, mask); },
              "target"_a, "source"_a, "index"_a, "mask"_a = true);

        m.def("scatter_add",
            [](Array &target, const Array &source,
               const Index &index,
               const IndexMask &mask) { scatter_add(target, source, index, mask); },
            "target"_a, "source"_a, "index"_a, "mask"_a = true);

        if constexpr (!IsDiff) {
            m.def("compress", [](const Array &array, const IndexMask &mask) {
                return compress(array, mask);
            }, "array"_a, "mask"_a);
        }
    }

    if constexpr (!IsMask) {
        m.def("sqr",        [](const Array &a) { return enoki::sqr(a); });
        m.def("clamp", [](const Array &value, const Array &min, const Array &max) {
            return enoki::clamp(value, min, max);
        });

        if (std::is_signed_v<Scalar>)
            m.def("abs",    [](const Array &a) { return enoki::abs(a); });
    }

    if constexpr (IsFloat) {
        m.def("sqrt",       [](const Array &a) { return enoki::sqrt(a); });
        m.def("safe_sqrt",  [](const Array &a) { return enoki::safe_sqrt(a); });
        m.def("cbrt",       [](const Array &a) { return enoki::cbrt(a); });
        m.def("rcp",        [](const Array &a) { return enoki::rcp(a); });
        m.def("rsqrt",      [](const Array &a) { return enoki::rsqrt(a); });
        m.def("safe_rsqrt", [](const Array &a) { return enoki::safe_rsqrt(a); });

        m.def("ceil",       [](const Array &a) { return enoki::ceil(a); });
        m.def("floor",      [](const Array &a) { return enoki::floor(a); });
        m.def("round",      [](const Array &a) { return enoki::round(a); });
        m.def("trunc",      [](const Array &a) { return enoki::trunc(a); });

        m.def("sign",         [](const Array &a) { return enoki::sign(a); });
        m.def("copysign",     [](const Array &a, const Array &b) { return enoki::copysign(a, b); });
        m.def("copysign_neg", [](const Array &a, const Array &b) {
            return enoki::copysign_neg(a, b);
        });
        m.def("mulsign",      [](const Array &a, const Array &b) { return enoki::mulsign(a, b); });
        m.def("mulsign_neg",  [](const Array &a, const Array &b) {
            return enoki::mulsign_neg(a, b);
        });

        m.def("lerp", [](const Array &a, const Array &b, const Array &t) {
            return enoki::lerp(a, b, t);
        });

        m.def("isfinite", [](const Array &a) -> Mask { return enoki::isfinite(a); });
        m.def("isnan", [](const Array &a) -> Mask { return enoki::isnan(a); });
        m.def("isinf", [](const Array &a) -> Mask { return enoki::isinf(a); });
    }

    if constexpr (!IsFloat && !IsMask && Size == -1) {
        m.def("popcnt", [](const Array &a) { return enoki::popcnt(a); });
        m.def("lzcnt", [](const Array &a) { return enoki::lzcnt(a); });
        m.def("tzcnt", [](const Array &a) { return enoki::tzcnt(a); });
        m.def("log2i", [](const Array &a) { return enoki::log2i(a); });
        m.def("mulhi", [](const Array &a, const Array &b) { return enoki::mulhi(a, b); });
    }

    if constexpr (IsFloat && Size == -1) {
        m.def("sin",        [](const Array &a) { return enoki::sin(a); });
        m.def("cos",        [](const Array &a) { return enoki::cos(a); });
        m.def("sincos",     [](const Array &a) { return enoki::sincos(a); });
        m.def("tan",        [](const Array &a) { return enoki::tan(a); });
        m.def("sec",        [](const Array &a) { return enoki::sec(a); });
        m.def("csc",        [](const Array &a) { return enoki::csc(a); });
        m.def("cot",        [](const Array &a) { return enoki::cot(a); });
        m.def("asin",       [](const Array &a) { return enoki::asin(a); });
        m.def("safe_asin",  [](const Array &a) { return enoki::safe_asin(a); });
        m.def("acos",       [](const Array &a) { return enoki::acos(a); });
        m.def("safe_acos",  [](const Array &a) { return enoki::safe_acos(a); });
        m.def("atan",       [](const Array &a) { return enoki::atan(a); });
        m.def("atan2",      [](const Array &x, const Array &y) {
            return enoki::atan2(x, y);
        }, "y"_a, "x"_a);

        m.def("sinh",    [](const Array &a) { return enoki::sinh(a); });
        m.def("cosh",    [](const Array &a) { return enoki::cosh(a); });
        m.def("sincosh", [](const Array &a) { return enoki::sincosh(a); });
        m.def("tanh",    [](const Array &a) { return enoki::tanh(a); });
        m.def("sech",    [](const Array &a) { return enoki::sech(a); });
        m.def("csch",    [](const Array &a) { return enoki::csch(a); });
        m.def("coth",    [](const Array &a) { return enoki::coth(a); });
        m.def("asinh",   [](const Array &a) { return enoki::asinh(a); });
        m.def("acosh",   [](const Array &a) { return enoki::acosh(a); });
        m.def("atanh",   [](const Array &a) { return enoki::atanh(a); });

        m.def("log",    [](const Array &a) { return enoki::log(a); });
        m.def("exp",    [](const Array &a) { return enoki::exp(a); });
        m.def("erfinv", [](const Array &a) { return enoki::erfinv(a); });
        m.def("erf",    [](const Array &a) { return enoki::erf(a); });
        m.def("lgamma", [](const Array &a) { return enoki::lgamma(a); });
        m.def("tgamma", [](const Array &a) { return enoki::tgamma(a); });
        m.def("pow",    [](const Array &a, const Array &b) {
            return enoki::pow(a, b);
        });

        m.def("pow",    [](const Array &a, int b) {
            return enoki::pow(a, b);
        });

        m.def("pow",    [](const Array &a, Scalar b) {
            return enoki::pow(a, b);
        });

        cl.def("__pow__", [](const Array &a, const Array &b) {
            return enoki::pow(a, b);
        });

        cl.def("__pow__", [](const Array &a, int b) {
            return enoki::pow(a, b);
        });

        cl.def("__pow__", [](const Array &a, Scalar b) {
            return enoki::pow(a, b);
        });
    }

    m.def("reverse", [](const Array &a) { return enoki::reverse(a); });

    if constexpr (!IsMask) {
        m.def("max", [](const Array &a, const Array &b) { return enoki::max(a, b); });
        m.def("min", [](const Array &a, const Array &b) { return enoki::min(a, b); });

        m.def("psum", [](const Array &a) { return enoki::psum(a); });

        m.def("hsum", [](const Array &a) { return enoki::hsum(a); });
        m.def("hprod", [](const Array &a) { return enoki::hprod(a); });
        m.def("hmin", [](const Array &a) { return enoki::hmin(a); });
        m.def("hmax", [](const Array &a) { return enoki::hmax(a); });
        m.def("hmean", [](const Array &a) { return enoki::hmean(a); });

        m.def("hsum_nested", [](const Array &a) { return enoki::hsum_nested(a); });
        m.def("hprod_nested", [](const Array &a) { return enoki::hprod_nested(a); });
        m.def("hmin_nested", [](const Array &a) { return enoki::hmin_nested(a); });
        m.def("hmax_nested", [](const Array &a) { return enoki::hmax_nested(a); });
        m.def("hmean_nested", [](const Array &a) { return enoki::hmean_nested(a); });

        m.def("fmadd", [](const Array &a, const Array &b, const Array &c) {
            return enoki::fmadd(a, b, c);
        });
        m.def("fmsub", [](const Array &a, const Array &b, const Array &c) {
            return enoki::fmsub(a, b, c);
        });
        m.def("fnmadd", [](const Array &a, const Array &b, const Array &c) {
            return enoki::fnmadd(a, b, c);
        });
        m.def("fnmsub", [](const Array &a, const Array &b, const Array &c) {
            return enoki::fnmsub(a, b, c);
        });
    } else {
        m.def("any", [](const Array &a) { return enoki::any(a); });
        m.def("none", [](const Array &a) { return enoki::none(a); });
        m.def("all", [](const Array &a) { return enoki::all(a); });

        using Count = value_t<uint64_array_t<array_t<Array>>>;
        m.def("count", [](const Array &a) -> Count { return enoki::count(a); });

        m.def("any_nested", [](const Array &a) { return enoki::any_nested(a); });
        m.def("none_nested", [](const Array &a) { return enoki::none_nested(a); });
        m.def("all_nested", [](const Array &a) { return enoki::all_nested(a); });
        m.def("count_nested", [](const Array &a) { return enoki::count_nested(a); });
    }

    m.def("eq", [](const Array &a, const Array &b) -> Mask { return eq(a, b); });
    m.def("neq", [](const Array &a, const Array &b) -> Mask { return neq(a, b); });

    if constexpr (!IsMask) {
        m.def("select", [](const Mask &a, const Array &b, const Array &c) {
            return enoki::select(a, b, c);
        });
    } else {
        m.def("select", [](const Array &a, const Array &b, const Array &c) -> Array {
            return enoki::select(a, b, c);
        });
    }

    if constexpr (!IsMask && IsDiff) {
        using Detached = decltype(eval(detach(std::declval<Array&>())));

        m.def("detach", [](const Array &a) -> Detached { return detach(a); });

        if constexpr (IsFloat) {
            m.def("reattach", [](Array &a, Array &b) { reattach(a, b); });
            m.def("requires_gradient",
                  [](const Array &a) { return requires_gradient(a); },
                  "array"_a);

            m.def("set_requires_gradient",
                  [](Array &a, bool value) { set_requires_gradient(a, value); },
                  "array"_a, "value"_a = true);

            m.def("gradient", [](Array &a) -> Detached { return gradient(a); });
            m.def("gradient_index", [](Array &a) { return gradient_index(a); });

            m.def("set_gradient",
                  [](Array &a, const Detached &g, bool backward) { set_gradient(a, g, backward); },
                  "array"_a, "gradient"_a, "backward"_a = true);

            m.def("graphviz", [](const Array &a) { return graphviz(a); });

            if constexpr (array_depth_v<Array> == 1) {
                m.def("backward",
                      [](Array &a, bool free_graph) { backward(a, free_graph); },
                      "array"_a, "free_graph"_a = true);
                m.def("forward",
                      [](Array &a, bool free_graph) { return forward(a, free_graph); },
                      "array"_a, "free_graph"_a = true);
            }
        }
    }

    if constexpr (is_cuda_array_v<Array>) {
        m.def("set_label", [](const Array &a, const char *label) {
            set_label(a, label);
        });
    }

    register_implicit_casts<Array, value_t<Array>>();

    return cl;
}

template <typename Scalar> py::object torch_dtype() {
    py::object torch = py::module::import("torch");
    const char *name = nullptr;

    if (std::is_same_v<Scalar, enoki::half>) {
        name = "float16";
    } else if (std::is_same_v<Scalar, float>) {
        name = "float32";
    } else if (std::is_same_v<Scalar, double>) {
        name = "float64";
    } else if (std::is_integral_v<Scalar>) {
        if (sizeof(Scalar) == 1)
            name = std::is_signed_v<Scalar> ? "int8" : "uint8";
        else if (sizeof(Scalar) == 2)
            name = "int16";
        else if (sizeof(Scalar) == 4)
            name = "int32";
        else if (sizeof(Scalar) == 8)
            name = "int64";
    }

    if (name == nullptr)
        throw std::runtime_error("pytorch_dtype(): Unsupported type");

    return torch.attr(name);
}

template <bool Scatter, size_t Index, size_t Dim, typename Source, typename Target>
static void copy_array(size_t offset,
                       const std::array<size_t, Dim> &shape,
                       const std::array<size_t, Dim> &strides,
                       const Source &source, Target &target) {
    using namespace enoki;

    size_t cur_shape = shape[Index],
           cur_stride = strides[Index];

    if constexpr (Index == Dim - 1) {
        if constexpr (is_cuda_array_v<Source>) {
            using UInt32 = uint32_array_t<Source>;
            UInt32 index = fmadd(arange<UInt32>((uint32_t) cur_shape),
                                 (uint32_t) cur_stride, (uint32_t) offset);
            if constexpr (Scatter)
                scatter(target, source, index);
            else
                target = gather<Target>(source, index);
        } else {
            for (size_t i = 0; i < cur_shape; ++i) {
                if constexpr (Scatter)
                    target[offset + cur_stride*i] = source[i];
                else
                    target[i] = source[offset + cur_stride*i];
            }
        }
    } else {
        for (size_t i = 0; i < cur_shape; ++i) {
            if constexpr (Scatter)
                copy_array<Scatter, Index + 1, Dim>(offset, shape, strides,
                                                    source.coeff(i), target);
            else
                copy_array<Scatter, Index + 1, Dim>(offset, shape, strides,
                                                    source, target.coeff(i));
            offset += cur_stride;
        }
    }
}

template <typename Array>
ENOKI_NOINLINE py::object enoki_to_torch(const Array &src, bool eval) {
    constexpr size_t Depth = array_depth_v<Array>;
    using Scalar = scalar_t<Array>;

    if (enoki::ragged(src))
        throw std::runtime_error("Enoki array is ragged -- cannot convert to PyTorch format!");

    std::array<size_t, Depth> shape = enoki::shape(src),
                              shape_rev = shape,
                              strides;
    std::reverse(shape_rev.begin(), shape_rev.end());

    py::object torch = py::module::import("torch");
    py::object dtype_obj = torch_dtype<Scalar>();

    py::object result = torch.attr("empty")(
        py::cast(shape_rev),
        "dtype"_a = dtype_obj,
        "device"_a = is_cuda_array_v<Array> ? "cuda" : "cpu");

    size_t size = 1;
    for (size_t i : shape)
        size *= i;

    strides = py::cast<std::array<size_t, Depth>>(result.attr("stride")());
    std::reverse(strides.begin(), strides.end());
    if (size > 0) {
        using T = std::conditional_t<
            is_cuda_array_v<Array>,
            CUDAArray<Scalar>,
            DynamicArray<Packet<Scalar>>
        >;
        T target = T::map(
            (Scalar *) py::cast<uintptr_t>(result.attr("data_ptr")()), size);
        copy_array</* Scatter = */ true, 0>(0, shape, strides, src, target);

#if defined(ENOKI_CUDA)
        if constexpr (is_cuda_array_v<Array>) {
            if (eval) {
                cuda_eval();
                cuda_sync();
            }
        }
#endif
    }
    return result;
}

template <typename Array>
ENOKI_NOINLINE Array torch_to_enoki(py::object src) {
    constexpr size_t Depth = array_depth_v<Array>;
    using Scalar = scalar_t<Array>;

    py::tuple shape_obj = src.attr("shape");
    py::object dtype_obj = src.attr("dtype");
    py::object target_dtype = torch_dtype<Scalar>();
    std::string device = (std::string) (py::str) src.attr("device");

    if (is_cuda_array_v<Array> && strncmp(device.c_str(), "cuda", 4) != 0)
        throw std::runtime_error("Attempted to cast a Torch CPU tensor to a Enoki GPU array!");
    if (!is_cuda_array_v<Array> && strncmp(device.c_str(), "cpu", 3) != 0)
        throw std::runtime_error("Attempted to cast a Torch GPU tensor to a Enoki CPU array!");

    if (shape_obj.size() != Depth || !dtype_obj.is(target_dtype))
        throw py::reference_cast_error();

    auto shape = py::cast<std::array<size_t, Depth>>(shape_obj);
    auto strides = py::cast<std::array<size_t, Depth>>(src.attr("stride")());
    std::reverse(shape.begin(), shape.end());
    std::reverse(strides.begin(), strides.end());

    size_t size = 1;
    for (size_t i : shape)
        size *= i;

    Array result;

    if constexpr (!is_cuda_array_v<Array>)
        set_shape(result, shape);

    if (size > 0) {
        using T = std::conditional_t<
            is_cuda_array_v<Array>,
            CUDAArray<Scalar>,
            DynamicArray<Packet<Scalar>>
        >;

        T source = T::map(
            (Scalar *) py::cast<uintptr_t>(src.attr("data_ptr")()), size);

        copy_array</* Scatter = */ false, 0>(0, shape, strides, source, result);
    }
    return result;
}

template <typename Array>
ENOKI_NOINLINE py::object enoki_to_numpy(const Array &src, bool eval) {
    constexpr size_t Depth = array_depth_v<Array>;
    using Scalar = scalar_t<Array>;

    if (enoki::ragged(src))
        throw std::runtime_error("Enoki array is ragged -- cannot convert to NumPy format!");

    std::array<size_t, Depth> shape = enoki::shape(src),
                              shape_rev = shape, strides;
    std::reverse(shape_rev.begin(), shape_rev.end());

    size_t size = 1, stride = sizeof(Scalar);
    for (ssize_t i = (ssize_t) Depth - 1; i >= 0; --i) {
        size *= shape_rev[i];
        strides[i] = stride;
        stride *= shape_rev[i];
    }

    using BufferType = Buffer<is_cuda_array_v<Array>>;
    BufferType *buf = new BufferType(stride);
    py::object buf_py = py::cast(buf, py::return_value_policy::take_ownership);

    py::array result(py::dtype::of<Scalar>(), shape_rev, strides, buf->ptr, buf_py);
    for (ssize_t i = (ssize_t) Depth - 1; i >= 0; --i)
        strides[i] /= sizeof(Scalar);
    std::reverse(strides.begin(), strides.end());

    if (size > 0) {
        using T = std::conditional_t<
            is_cuda_array_v<Array>,
            CUDAArray<Scalar>,
            DynamicArray<Packet<Scalar>>
        >;
        T target = T::map(buf->ptr, stride);
        copy_array</* Scatter = */ true, 0>(0, shape, strides, src, target);

#if defined(ENOKI_CUDA)
        if constexpr (is_cuda_array_v<Array>) {
            if (eval) {
                cuda_eval();
                cuda_sync();
            }
        }
#endif
    }

    return result;
}

template <typename Array>
ENOKI_NOINLINE Array numpy_to_enoki(py::array src) {
    constexpr size_t Depth = array_depth_v<Array>;
    using SizeArray = std::array<size_t, Depth>;
    using Scalar = scalar_t<Array>;

    py::tuple shape_obj = src.attr("shape");
    py::object dtype_obj = src.attr("dtype");
    py::object target_dtype = py::dtype::of<Scalar>();

    if (shape_obj.size() != Depth)
        throw py::reference_cast_error();

    if (!dtype_obj.is(target_dtype))
        src = py::array_t<Scalar, py::array::forcecast>::ensure(src);

    SizeArray shape   = py::cast<SizeArray>(shape_obj),
              strides = py::cast<SizeArray>(src.attr("strides"));

    std::reverse(shape.begin(), shape.end());
    std::reverse(strides.begin(), strides.end());

    size_t size = 1;
    for (size_t i : shape)
        size *= i;

    for (ssize_t i = 0; i < Depth; ++i)
        strides[i] /= sizeof(Scalar);

    Array result;

    if constexpr (!is_cuda_array_v<Array>)
        set_shape(result, shape);

    if (size > 0) {
        using T = std::conditional_t<
            is_cuda_array_v<Array>,
            CUDAArray<Scalar>,
            DynamicArray<Packet<Scalar>>
        >;
        const T source = T::copy(src.data(), size);

        copy_array</* Scatter = */ false, 0>(0, shape, strides, source, result);
    }

    return result;
}
