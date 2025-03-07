#include "common.h"
#include <pybind11/functional.h>

void bind_cuda_autodiff_1d(py::module& m, py::module& s) {
    auto mask_class = bind<mask_t<Float32D>>(m, s, "Mask");
    auto uint32_class = bind<UInt32D>(m, s, "UInt32");
    auto uint64_class = bind<UInt64D>(m, s, "UInt64");
    auto int32_class = bind<Int32D>(m, s, "Int32");
    auto int64_class = bind<Int64D>(m, s, "Int64");
    auto float32_class = bind<Float32D>(m, s, "Float32");
    auto float64_class = bind<Float64D>(m, s, "Float64");

    mask_class
        .def(py::init<const MaskC &>());

    float32_class
        .def(py::init<const Float32C &>())
        .def(py::init<const Float64D &>())
        .def(py::init<const Int64D &>())
        .def(py::init<const Int32D &>())
        .def(py::init<const UInt64D &>())
        .def(py::init<const UInt32D &>())
        .def("set_graph_simplification", [](bool value) { Float32D::set_graph_simplification_(value); })
        .def("whos", []() { py::print(Float32D::whos_()); })
        .def_static("set_log_level",
             [](int log_level) { Float32D::set_log_level_(log_level); },
             "Sets the current log level (0 == none, 1 == minimal, 2 == moderate, 3 == high, 4 == everything)")
        .def_static("log_level", []() { return Float32D::log_level_(); })
        .def_static("simplify_graph", []() { Float32D::simplify_graph_(); })
        .def_static("backward",
                    [](bool free_graph) { backward<Float32D>(free_graph); },
                    "free_graph"_a = true)
        .def_static("forward",
                    [](bool free_graph) { forward<Float32D>(free_graph); },
                    "free_graph"_a = true);

    float64_class
        .def(py::init<const Float64C &>())
        .def(py::init<const Float32D &>())
        .def(py::init<const Int64D &>())
        .def(py::init<const Int32D &>())
        .def(py::init<const UInt64D &>())
        .def(py::init<const UInt32D &>())
        .def("set_graph_simplification", [](bool value) { Float64D::set_graph_simplification_(value); })
        .def("whos", []() { py::print(Float64D::whos_()); })
        .def_static("set_log_level",
             [](int log_level) { Float64D::set_log_level_(log_level); },
             "Sets the current log level (0 == none, 1 == minimal, 2 == moderate, 3 == high, 4 == everything)")
        .def_static("log_level", []() { return Float64D::log_level_(); })
        .def_static("simplify_graph", []() { Float64D::simplify_graph_(); })
        .def_static("backward",
                    [](bool free_graph) { backward<Float64D>(free_graph); },
                    "free_graph"_a = true)
        .def_static("forward",
                    [](bool free_graph) { forward<Float64D>(free_graph); },
                    "free_graph"_a = true);

    int32_class
        .def(py::init<const Int32C &>())
        .def(py::init<const Int64D &>())
        .def(py::init<const UInt32D &>())
        .def(py::init<const UInt64D &>())
        .def(py::init<const Float32D &>())
        .def(py::init<const Float64D &>());

    int64_class
        .def(py::init<const Int32D &>())
        .def(py::init<const Int64C &>())
        .def(py::init<const UInt32D &>())
        .def(py::init<const UInt64D &>())
        .def(py::init<const Float32D &>())
        .def(py::init<const Float64D &>());

    uint32_class
        .def(py::init<const Int32D &>())
        .def(py::init<const Int64D &>())
        .def(py::init<const UInt32C &>())
        .def(py::init<const UInt64D &>())
        .def(py::init<const Float32D &>())
        .def(py::init<const Float64D &>());

    uint64_class
        .def(py::init<const Int32D &>())
        .def(py::init<const Int64D &>())
        .def(py::init<const UInt32D &>())
        .def(py::init<const UInt64C &>())
        .def(py::init<const Float32D &>())
        .def(py::init<const Float64D &>());

    auto vector1m_class = bind<Vector1mD>(m, s, "Vector1m");
    auto vector1i_class = bind<Vector1iD>(m, s, "Vector1i");
    auto vector1u_class = bind<Vector1uD>(m, s, "Vector1u");
    auto vector1f_class = bind<Vector1fD>(m, s, "Vector1f");
    auto vector1d_class = bind<Vector1dD>(m, s, "Vector1d");

    vector1f_class
        .def(py::init<const Vector1f  &>())
        .def(py::init<const Vector1fC &>())
        .def(py::init<const Vector1dD &>())
        .def(py::init<const Vector1uD &>())
        .def(py::init<const Vector1iD &>());

    vector1d_class
        .def(py::init<const Vector1d  &>())
        .def(py::init<const Vector1dC &>())
        .def(py::init<const Vector1fD &>())
        .def(py::init<const Vector1uD &>())
        .def(py::init<const Vector1iD &>());

    vector1i_class
        .def(py::init<const Vector1i  &>())
        .def(py::init<const Vector1iC &>())
        .def(py::init<const Vector1fD &>())
        .def(py::init<const Vector1dD &>())
        .def(py::init<const Vector1uD &>());

    vector1u_class
        .def(py::init<const Vector1u  &>())
        .def(py::init<const Vector1uC &>())
        .def(py::init<const Vector1fD &>())
        .def(py::init<const Vector1dD &>())
        .def(py::init<const Vector1iD &>());

    m.def(
        "binary_search",
        [](uint32_t start,
           uint32_t end,
           const std::function<MaskD(UInt32D)> &pred) {
            return enoki::binary_search(start, end, pred);
        },
        "start"_a, "end"_a, "pred"_a);

    m.def("meshgrid", [](const Float32D &x, const Float32D &y) {
        auto result = meshgrid(x, y);
        return std::make_pair(std::move(result.x()), std::move(result.y()));
    });

    m.def("meshgrid", [](const Float64D &x, const Float64D &y) {
        auto result = meshgrid(x, y);
        return std::make_pair(std::move(result.x()), std::move(result.y()));
    });

    struct Scope {
        Scope(const std::string &name) : name(name) { }

        void enter() { Float32D::push_prefix_(name.c_str()); }
        void exit(py::handle, py::handle, py::handle) { Float32D::pop_prefix_(); }

        std::string name;
    };

    py::class_<Scope>(float32_class, "Scope")
        .def(py::init<const std::string &>())
        .def("__enter__", &Scope::enter)
        .def("__exit__", &Scope::exit);
}
