//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main author:     Suneth Warnakulasuriya
//

// System includes

// External includes
#include <pybind11/stl.h>

// Project includes

// Application includes
#include "custom_utilities/tensor_adaptors/properties_variable_tensor_adaptor.h"

// Include base h
#include "add_custom_tensor_adaptors_to_python.h"

namespace Kratos {
namespace Python {

void  AddCustomTensorAdaptorsToPython(pybind11::module& m)
{
    namespace py = pybind11;

    auto tensor_adaptor_modules = m.def_submodule("TensorAdaptors");

    py::class_<PropertiesVariableTensorAdaptor, PropertiesVariableTensorAdaptor::Pointer, PropertiesVariableTensorAdaptor::BaseType>(tensor_adaptor_modules, "PropertiesVariableTensorAdaptor")
        .def(py::init<const PropertiesVariableTensorAdaptor::BaseType&, PropertiesVariableTensorAdaptor::VariablePointerType, const bool>(),
            py::arg("tensor_adaptor"),
            py::arg("variable"),
            py::arg("copy") = true)
        .def(py::init<ModelPart::ConditionsContainerType::Pointer, PropertiesVariableTensorAdaptor::VariablePointerType>(),
            py::arg("container"),
            py::arg("variable"))
        .def(py::init<ModelPart::ConditionsContainerType::Pointer, PropertiesVariableTensorAdaptor::VariablePointerType, const std::vector<unsigned int>&>(),
            py::arg("container"),
            py::arg("variable"),
            py::arg("data_shape"))
        .def(py::init<ModelPart::ElementsContainerType::Pointer, PropertiesVariableTensorAdaptor::VariablePointerType>(),
            py::arg("container"),
            py::arg("variable"))
        .def(py::init<ModelPart::ElementsContainerType::Pointer, PropertiesVariableTensorAdaptor::VariablePointerType, const std::vector<unsigned int>&>(),
            py::arg("container"),
            py::arg("variable"),
            py::arg("data_shape"))
        ;
}

}  // namespace Python.
} // Namespace Kratos

