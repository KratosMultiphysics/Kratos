//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes

// External includes

// Project includes
#include "includes/define_python.h"
#include "future/input_output/vtu_output.h"

// Include base h
#include "add_io_to_python.h"

namespace Kratos::Future::Python
{

void  AddIOToPython(pybind11::module& m)
{
    namespace py = pybind11;

    auto vtu_output = py::class_<VtuOutput, VtuOutput::Pointer>(m, "VtuOutput");

    py::enum_<VtuOutput::WriterFormat>(vtu_output, "WriterFormat")
        .value("ASCII", VtuOutput::WriterFormat::ASCII)
        .value("BINARY", VtuOutput::WriterFormat::BINARY)
        .value("RAW", VtuOutput::WriterFormat::RAW)
        .value("COMPRESSED_RAW", VtuOutput::WriterFormat::COMPRESSED_RAW)
        .export_values()
        ;

    vtu_output
        .def(py::init<ModelPart&, const Globals::Configuration, const VtuOutput::WriterFormat, const std::size_t, const bool, const bool, const IndexType>(),
            py::arg("model_part"),
            py::arg("configuration") = Globals::Configuration::Initial,
            py::arg("output_format") = VtuOutput::WriterFormat::COMPRESSED_RAW,
            py::arg("precision") = 9,
            py::arg("output_sub_model_parts") = false,
            py::arg("write_ids") = false,
            py::arg("echo_level") = 0)
        .def("AddFlag", &VtuOutput::AddFlag, py::arg("flag_name"), py::arg("flag"), py::arg("data_location"))
        .def("AddVariable", &VtuOutput::AddVariable<int>, py::arg("variable"), py::arg("data_location"))
        .def("AddVariable", &VtuOutput::AddVariable<double>, py::arg("variable"), py::arg("data_location"))
        .def("AddVariable", &VtuOutput::AddVariable<array_1d<double, 3>>, py::arg("variable"), py::arg("data_location"))
        .def("AddVariable", &VtuOutput::AddVariable<array_1d<double, 4>>, py::arg("variable"), py::arg("data_location"))
        .def("AddVariable", &VtuOutput::AddVariable<array_1d<double, 6>>, py::arg("variable"), py::arg("data_location"))
        .def("AddVariable", &VtuOutput::AddVariable<array_1d<double, 9>>, py::arg("variable"), py::arg("data_location"))
        .def("AddVariable", &VtuOutput::AddVariable<Vector>, py::arg("variable"), py::arg("data_location"))
        .def("AddVariable", &VtuOutput::AddVariable<Matrix>, py::arg("variable"), py::arg("data_location"))
        .def("AddIntegrationPointVariable", &VtuOutput::AddIntegrationPointVariable<int>, py::arg("variable"), py::arg("data_location"))
        .def("AddIntegrationPointVariable", &VtuOutput::AddIntegrationPointVariable<double>, py::arg("variable"), py::arg("data_location"))
        .def("AddIntegrationPointVariable", &VtuOutput::AddIntegrationPointVariable<array_1d<double, 3>>, py::arg("variable"), py::arg("data_location"))
        .def("AddIntegrationPointVariable", &VtuOutput::AddIntegrationPointVariable<array_1d<double, 4>>, py::arg("variable"), py::arg("data_location"))
        .def("AddIntegrationPointVariable", &VtuOutput::AddIntegrationPointVariable<array_1d<double, 6>>, py::arg("variable"), py::arg("data_location"))
        .def("AddIntegrationPointVariable", &VtuOutput::AddIntegrationPointVariable<array_1d<double, 9>>, py::arg("variable"), py::arg("data_location"))
        .def("AddIntegrationPointVariable", &VtuOutput::AddIntegrationPointVariable<Vector>, py::arg("variable"), py::arg("data_location"))
        .def("AddIntegrationPointVariable", &VtuOutput::AddIntegrationPointVariable<Matrix>, py::arg("variable"), py::arg("data_location"))
        .def("AddTensorAdaptor", &VtuOutput::AddTensorAdaptor<TensorAdaptor<bool>::Pointer>, py::arg("tensor_adaptor_name"), py::arg("tensor_adaptor"))
        .def("AddTensorAdaptor", &VtuOutput::AddTensorAdaptor<TensorAdaptor<int>::Pointer>, py::arg("tensor_adaptor_name"), py::arg("tensor_adaptor"))
        .def("AddTensorAdaptor", &VtuOutput::AddTensorAdaptor<TensorAdaptor<double>::Pointer>, py::arg("tensor_adaptor_name"), py::arg("tensor_adaptor"))
        .def("ReplaceTensorAdaptor", &VtuOutput::ReplaceTensorAdaptor<TensorAdaptor<bool>::Pointer>, py::arg("tensor_adaptor_name"), py::arg("tensor_adaptor"))
        .def("ReplaceTensorAdaptor", &VtuOutput::ReplaceTensorAdaptor<TensorAdaptor<int>::Pointer>, py::arg("tensor_adaptor_name"), py::arg("tensor_adaptor"))
        .def("ReplaceTensorAdaptor", &VtuOutput::ReplaceTensorAdaptor<TensorAdaptor<double>::Pointer>, py::arg("tensor_adaptor_name"), py::arg("tensor_adaptor"))
        .def("EmplaceTensorAdaptor", &VtuOutput::EmplaceTensorAdaptor<TensorAdaptor<bool>::Pointer>, py::arg("tensor_adaptor_name"), py::arg("tensor_adaptor"))
        .def("EmplaceTensorAdaptor", &VtuOutput::EmplaceTensorAdaptor<TensorAdaptor<int>::Pointer>, py::arg("tensor_adaptor_name"), py::arg("tensor_adaptor"))
        .def("EmplaceTensorAdaptor", &VtuOutput::EmplaceTensorAdaptor<TensorAdaptor<double>::Pointer>, py::arg("tensor_adaptor_name"), py::arg("tensor_adaptor"))
        .def("GetModelPart", &VtuOutput::GetModelPart, py::return_value_policy::reference)
        .def("GetOutputContainerList", &VtuOutput::GetOutputContainerList)
        .def("PrintOutput", &VtuOutput::PrintOutput, py::arg("output_file_name_prefix"))
        .def("__str__", PrintObject<VtuOutput>)
        ;
}
}  // namespace Kratos::Python.


