//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya (https://github.com/sunethwarna)
//

// System includes

// External includes
#include <pybind11/pybind11.h>

// Project includes

// Application includes
#include "custom_methods/temporal_methods.h"

// Include base h
#include "custom_python/add_custom_temporal_non_historical_nodal_methods_to_python.h"

namespace Kratos
{
namespace Python
{
void AddCustomTemporalNonHistoricalNodalMethodsToPython(pybind11::module& m)
{
    namespace py = pybind11;

    using type_double = double;
    using type_array = array_1d<double, 3>;
    using type_vector = Vector;
    using type_matrix = Matrix;

    auto container_module = m.def_submodule("Nodes");

    auto value_module = container_module.def_submodule("ValueMethods");
    {
        auto method_module = value_module.def_submodule("Sum");
        using current_method = TemporalMethods::NodalNonHistoricalTemporalMethods::SumMethod;
        py::class_<current_method::ValueMethod<type_double>, current_method::ValueMethod<type_double>::Pointer,
                   TemporalMethods::TemporalMethod>(method_module, "Double")
            .def(py::init<ModelPart&, const std::string&, const Variable<type_double>&,
                          const Variable<type_double>&>());
        py::class_<current_method::ValueMethod<type_array>, current_method::ValueMethod<type_array>::Pointer,
                   TemporalMethods::TemporalMethod>(method_module, "Array")
            .def(py::init<ModelPart&, const std::string&, const Variable<type_array>&,
                          const Variable<type_array>&>());
        py::class_<current_method::ValueMethod<type_vector>, current_method::ValueMethod<type_vector>::Pointer,
                   TemporalMethods::TemporalMethod>(method_module, "Vector")
            .def(py::init<ModelPart&, const std::string&, const Variable<type_vector>&,
                          const Variable<type_vector>&>());
        py::class_<current_method::ValueMethod<type_matrix>, current_method::ValueMethod<type_matrix>::Pointer,
                   TemporalMethods::TemporalMethod>(method_module, "Matrix")
            .def(py::init<ModelPart&, const std::string&, const Variable<type_matrix>&,
                          const Variable<type_matrix>&>());
    }
    {
        auto method_module = value_module.def_submodule("Mean");
        using current_method = TemporalMethods::NodalNonHistoricalTemporalMethods::MeanMethod;
        py::class_<current_method::ValueMethod<type_double>, current_method::ValueMethod<type_double>::Pointer,
                   TemporalMethods::TemporalMethod>(method_module, "Double")
            .def(py::init<ModelPart&, const std::string&, const Variable<type_double>&,
                          const Variable<type_double>&>());
        py::class_<current_method::ValueMethod<type_array>, current_method::ValueMethod<type_array>::Pointer,
                   TemporalMethods::TemporalMethod>(method_module, "Array")
            .def(py::init<ModelPart&, const std::string&, const Variable<type_array>&,
                          const Variable<type_array>&>());
        py::class_<current_method::ValueMethod<type_vector>, current_method::ValueMethod<type_vector>::Pointer,
                   TemporalMethods::TemporalMethod>(method_module, "Vector")
            .def(py::init<ModelPart&, const std::string&, const Variable<type_vector>&,
                          const Variable<type_vector>&>());
        py::class_<current_method::ValueMethod<type_matrix>, current_method::ValueMethod<type_matrix>::Pointer,
                   TemporalMethods::TemporalMethod>(method_module, "Matrix")
            .def(py::init<ModelPart&, const std::string&, const Variable<type_matrix>&,
                          const Variable<type_matrix>&>());
    }
    {
        auto method_module = value_module.def_submodule("Variance");
        using current_method = TemporalMethods::NodalNonHistoricalTemporalMethods::VarianceMethod;
        py::class_<current_method::ValueMethod<type_double>, current_method::ValueMethod<type_double>::Pointer,
                   TemporalMethods::TemporalMethod>(method_module, "Double")
            .def(py::init<ModelPart&, const std::string&, const Variable<type_double>&,
                          const Variable<type_double>&, const Variable<type_double>&>());
        py::class_<current_method::ValueMethod<type_array>, current_method::ValueMethod<type_array>::Pointer,
                   TemporalMethods::TemporalMethod>(method_module, "Array")
            .def(py::init<ModelPart&, const std::string&, const Variable<type_array>&,
                          const Variable<type_array>&, const Variable<type_array>&>());
        py::class_<current_method::ValueMethod<type_vector>, current_method::ValueMethod<type_vector>::Pointer,
                   TemporalMethods::TemporalMethod>(method_module, "Vector")
            .def(py::init<ModelPart&, const std::string&, const Variable<type_vector>&,
                          const Variable<type_vector>&, const Variable<type_vector>&>());
        py::class_<current_method::ValueMethod<type_matrix>, current_method::ValueMethod<type_matrix>::Pointer,
                   TemporalMethods::TemporalMethod>(method_module, "Matrix")
            .def(py::init<ModelPart&, const std::string&, const Variable<type_matrix>&,
                          const Variable<type_matrix>&, const Variable<type_matrix>&>());
    }

    auto norm_module = container_module.def_submodule("NormMethods");
    {
        auto method_module = norm_module.def_submodule("Sum");
        using current_method = TemporalMethods::NodalNonHistoricalTemporalMethods::SumMethod;
        py::class_<current_method::NormMethod<type_double>, current_method::NormMethod<type_double>::Pointer,
                   TemporalMethods::TemporalMethod>(method_module, "Double")
            .def(py::init<ModelPart&, const std::string&, const Variable<type_double>&,
                          const Variable<double>&>());
        py::class_<current_method::NormMethod<type_array>, current_method::NormMethod<type_array>::Pointer,
                   TemporalMethods::TemporalMethod>(method_module, "Array")
            .def(py::init<ModelPart&, const std::string&, const Variable<type_array>&,
                          const Variable<double>&>());
        py::class_<current_method::NormMethod<type_vector>, current_method::NormMethod<type_vector>::Pointer,
                   TemporalMethods::TemporalMethod>(method_module, "Vector")
            .def(py::init<ModelPart&, const std::string&, const Variable<type_vector>&,
                          const Variable<double>&>());
        py::class_<current_method::NormMethod<type_matrix>, current_method::NormMethod<type_matrix>::Pointer,
                   TemporalMethods::TemporalMethod>(method_module, "Matrix")
            .def(py::init<ModelPart&, const std::string&, const Variable<type_matrix>&,
                          const Variable<double>&>());
    }
    {
        auto method_module = norm_module.def_submodule("Mean");
        using current_method = TemporalMethods::NodalNonHistoricalTemporalMethods::MeanMethod;
        py::class_<current_method::NormMethod<type_double>, current_method::NormMethod<type_double>::Pointer,
                   TemporalMethods::TemporalMethod>(method_module, "Double")
            .def(py::init<ModelPart&, const std::string&, const Variable<type_double>&,
                          const Variable<double>&>());
        py::class_<current_method::NormMethod<type_array>, current_method::NormMethod<type_array>::Pointer,
                   TemporalMethods::TemporalMethod>(method_module, "Array")
            .def(py::init<ModelPart&, const std::string&, const Variable<type_array>&,
                          const Variable<double>&>());
        py::class_<current_method::NormMethod<type_vector>, current_method::NormMethod<type_vector>::Pointer,
                   TemporalMethods::TemporalMethod>(method_module, "Vector")
            .def(py::init<ModelPart&, const std::string&, const Variable<type_vector>&,
                          const Variable<double>&>());
        py::class_<current_method::NormMethod<type_matrix>, current_method::NormMethod<type_matrix>::Pointer,
                   TemporalMethods::TemporalMethod>(method_module, "Matrix")
            .def(py::init<ModelPart&, const std::string&, const Variable<type_matrix>&,
                          const Variable<double>&>());
    }
    {
        auto method_module = norm_module.def_submodule("Variance");
        using current_method = TemporalMethods::NodalNonHistoricalTemporalMethods::VarianceMethod;
        py::class_<current_method::NormMethod<type_double>, current_method::NormMethod<type_double>::Pointer,
                   TemporalMethods::TemporalMethod>(method_module, "Double")
            .def(py::init<ModelPart&, const std::string&, const Variable<type_double>&,
                          const Variable<double>&, const Variable<double>&>());
        py::class_<current_method::NormMethod<type_array>, current_method::NormMethod<type_array>::Pointer,
                   TemporalMethods::TemporalMethod>(method_module, "Array")
            .def(py::init<ModelPart&, const std::string&, const Variable<type_array>&,
                          const Variable<double>&, const Variable<double>&>());
        py::class_<current_method::NormMethod<type_vector>, current_method::NormMethod<type_vector>::Pointer,
                   TemporalMethods::TemporalMethod>(method_module, "Vector")
            .def(py::init<ModelPart&, const std::string&, const Variable<type_vector>&,
                          const Variable<double>&, const Variable<double>&>());
        py::class_<current_method::NormMethod<type_matrix>, current_method::NormMethod<type_matrix>::Pointer,
                   TemporalMethods::TemporalMethod>(method_module, "Matrix")
            .def(py::init<ModelPart&, const std::string&, const Variable<type_matrix>&,
                          const Variable<double>&, const Variable<double>&>());
    }
    {
        auto method_module = norm_module.def_submodule("Max");
        using current_method = TemporalMethods::NodalNonHistoricalTemporalMethods::MaxMethod;
        py::class_<current_method::NormMethod<type_double>, current_method::NormMethod<type_double>::Pointer,
                   TemporalMethods::TemporalMethod>(method_module, "Double")
            .def(py::init<ModelPart&, const std::string&, const Variable<type_double>&,
                          const Variable<double>&, const Variable<double>&>());
        py::class_<current_method::NormMethod<type_array>, current_method::NormMethod<type_array>::Pointer,
                   TemporalMethods::TemporalMethod>(method_module, "Array")
            .def(py::init<ModelPart&, const std::string&, const Variable<type_array>&,
                          const Variable<double>&, const Variable<double>&>());
        py::class_<current_method::NormMethod<type_vector>, current_method::NormMethod<type_vector>::Pointer,
                   TemporalMethods::TemporalMethod>(method_module, "Vector")
            .def(py::init<ModelPart&, const std::string&, const Variable<type_vector>&,
                          const Variable<double>&, const Variable<double>&>());
        py::class_<current_method::NormMethod<type_matrix>, current_method::NormMethod<type_matrix>::Pointer,
                   TemporalMethods::TemporalMethod>(method_module, "Matrix")
            .def(py::init<ModelPart&, const std::string&, const Variable<type_matrix>&,
                          const Variable<double>&, const Variable<double>&>());
    }
    {
        auto method_module = norm_module.def_submodule("Min");
        using current_method = TemporalMethods::NodalNonHistoricalTemporalMethods::MinMethod;
        py::class_<current_method::NormMethod<type_double>, current_method::NormMethod<type_double>::Pointer,
                   TemporalMethods::TemporalMethod>(method_module, "Double")
            .def(py::init<ModelPart&, const std::string&, const Variable<type_double>&,
                          const Variable<double>&, const Variable<double>&>());
        py::class_<current_method::NormMethod<type_array>, current_method::NormMethod<type_array>::Pointer,
                   TemporalMethods::TemporalMethod>(method_module, "Array")
            .def(py::init<ModelPart&, const std::string&, const Variable<type_array>&,
                          const Variable<double>&, const Variable<double>&>());
        py::class_<current_method::NormMethod<type_vector>, current_method::NormMethod<type_vector>::Pointer,
                   TemporalMethods::TemporalMethod>(method_module, "Vector")
            .def(py::init<ModelPart&, const std::string&, const Variable<type_vector>&,
                          const Variable<double>&, const Variable<double>&>());
        py::class_<current_method::NormMethod<type_matrix>, current_method::NormMethod<type_matrix>::Pointer,
                   TemporalMethods::TemporalMethod>(method_module, "Matrix")
            .def(py::init<ModelPart&, const std::string&, const Variable<type_matrix>&,
                          const Variable<double>&, const Variable<double>&>());
    }
}

} // namespace Python.
} // Namespace Kratos
