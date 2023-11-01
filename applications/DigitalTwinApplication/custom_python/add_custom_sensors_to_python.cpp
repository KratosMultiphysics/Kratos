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
#include <pybind11/stl.h>

// Project includes

// Application includes
#include "custom_sensors/sensor_specification.h"
#include "custom_sensors/nodal_sensor_specification.h"
#include "custom_sensors/adjoint_sensor.h"
#include "custom_sensors/adjoint_displacement_sensor.h"

// Include base h
#include "custom_python/add_custom_sensors_to_python.h"

namespace Kratos::Python {

void  AddCustomSensorsToPython(pybind11::module& m)
{
    namespace py = pybind11;

    auto sensor_module = m.def_submodule("Sensors");

    // Add sensor specifications
    py::class_<SensorSpecification, SensorSpecification::Pointer, IndexedObject, DataValueContainer>(sensor_module, "SensorSpecification")
        .def("GetLocation", &SensorSpecification::GetLocation)
        .def("GetName", &SensorSpecification::GetName)
        .def("GetSensorValue", &SensorSpecification::GetSensorValue)
        .def("AddNodalExpression", &SensorSpecification::AddNodalExpression, py::arg("nodal_expression_name"), py::arg("nodal_expression"))
        .def("GetNodalExpression", &SensorSpecification::GetNodalExpression, py::arg("nodal_expression_name"))
        .def("GetNodalExpressionsMap", &SensorSpecification::GetNodalExpressionsMap)
        .def("AddConditionExpression", &SensorSpecification::AddConditionExpression, py::arg("condition_expression_name"), py::arg("condition_expression"))
        .def("GetConditionExpression", &SensorSpecification::GetConditionExpression, py::arg("condition_expression_name"))
        .def("GetConditionExpressionsMap", &SensorSpecification::GetConditionExpressionsMap)
        .def("AddElementExpression", &SensorSpecification::AddElementExpression, py::arg("element_expression_name"), py::arg("element_expression"))
        .def("GetElementExpression", &SensorSpecification::GetElementExpression, py::arg("element_expression_name"))
        .def("GetElementExpressionsMap", &SensorSpecification::GetElementExpressionsMap)
        .def("__str__", PrintObject<SensorSpecification>);
        ;

    py::class_<NodalSensorSpecification, NodalSensorSpecification::Pointer, SensorSpecification>(sensor_module, "NodalSensorSpecification")
        .def(py::init<const std::string&, const IndexType, const double, const double, const ModelPart::NodeType::Pointer>(), py::arg("sensor_name"), py::arg("sensor_id"), py::arg("sensor_value"), py::arg("sensor_weight"), py::arg("sensor_node"))
        .def("GetNode", &NodalSensorSpecification::GetNode)
        ;

    // Add sensor adjint responses
    py::class_<AdjointSensor, AdjointSensor::Pointer, AdjointResponseFunction>(sensor_module, "AdjointSensor")
        .def("SetSensorSpecification", &AdjointSensor::SetSensorSpecification, py::arg("sensor_specification"))
        ;

    py::class_<AdjointDisplacementSensor, AdjointDisplacementSensor::Pointer, AdjointSensor>(sensor_module, "AdjointDisplacementSensor")
        .def(py::init<Model&, Parameters>(), py::arg("model"), py::arg("sensor_parameters"))
        ;
}

} // namespace Kratos::Python
