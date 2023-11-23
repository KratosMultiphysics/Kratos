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
#include "response_functions/adjoint_response_function.h"

// Application includes
#include "custom_sensors/sensor.h"
#include "custom_sensors/sensor_view.h"
#include "custom_sensors/measurement_residual_response_function.h"
#include "custom_sensors/displacement_sensor.h"

// Include base h
#include "custom_python/add_custom_sensors_to_python.h"

namespace Kratos::Python {

void  AddCustomSensorsToPython(pybind11::module& m)
{
    namespace py = pybind11;

    auto sensor_module = m.def_submodule("Sensors");

    // Add sensor specifications
    py::class_<Sensor, Sensor::Pointer, AdjointResponseFunction, DataValueContainer>(sensor_module, "Sensor")
        .def("GetName", &Sensor::GetName)
        .def("GetLocation", &Sensor::GetLocation)
        .def("GetWeight", &Sensor::GetWeight)
        .def("GetSensorValue", &Sensor::GetSensorValue)
        .def("SetSensorValue", &Sensor::SetSensorValue)
        .def("GetSensorParameters", &Sensor::GetSensorParameters)
        .def("AddContainerExpression", &Sensor::AddContainerExpression<ModelPart::NodesContainerType>, py::arg("expression_name"), py::arg("nodal_expression"))
        .def("AddContainerExpression", &Sensor::AddContainerExpression<ModelPart::ConditionsContainerType>, py::arg("expression_name"), py::arg("condition_expression"))
        .def("AddContainerExpression", &Sensor::AddContainerExpression<ModelPart::ElementsContainerType>, py::arg("expression_name"), py::arg("element_expression"))
        .def("AddNodalExpression", &Sensor::AddNodalExpression, py::arg("nodal_expression_name"), py::arg("nodal_expression"))
        .def("GetNodalExpression", &Sensor::GetNodalExpression, py::arg("nodal_expression_name"))
        .def("GetNodalExpressionsMap", &Sensor::GetNodalExpressionsMap)
        .def("AddConditionExpression", &Sensor::AddConditionExpression, py::arg("condition_expression_name"), py::arg("condition_expression"))
        .def("GetConditionExpression", &Sensor::GetConditionExpression, py::arg("condition_expression_name"))
        .def("GetConditionExpressionsMap", &Sensor::GetConditionExpressionsMap)
        .def("AddElementExpression", &Sensor::AddElementExpression, py::arg("element_expression_name"), py::arg("element_expression"))
        .def("GetElementExpression", &Sensor::GetElementExpression, py::arg("element_expression_name"))
        .def("GetElementExpressionsMap", &Sensor::GetElementExpressionsMap)
        .def("GetDataVariableNames", &Sensor::GetDataVariableNames)
        .def("ClearNodalExpressions", &Sensor::ClearNodalExpressions)
        .def("ClearConditionExpressions", &Sensor::ClearConditionExpressions)
        .def("ClearElementExpressions", &Sensor::ClearElementExpressions)
        .def("__str__", PrintObject<Sensor>);
        ;

    using nodal_sensor_view = SensorView<ModelPart::NodesContainerType>;
    py::class_<nodal_sensor_view, nodal_sensor_view::Pointer>(sensor_module, "NodalSensorView")
        .def(py::init<Sensor::Pointer, const std::string&>(), py::arg("sensor"), py::arg("expression_name"))
        .def("GetSensor", &nodal_sensor_view::GetSensor)
        .def("GetContainerExpression", &nodal_sensor_view::GetContainerExpression)
        .def("__str__", PrintObject<nodal_sensor_view>);
        ;

    using condition_sensor_view = SensorView<ModelPart::ConditionsContainerType>;
    py::class_<condition_sensor_view, condition_sensor_view::Pointer>(sensor_module, "ConditionSensorView")
        .def(py::init<Sensor::Pointer, const std::string&>(), py::arg("sensor"), py::arg("expression_name"))
        .def("GetSensor", &condition_sensor_view::GetSensor)
        .def("GetContainerExpression", &condition_sensor_view::GetContainerExpression)
        .def("__str__", PrintObject<condition_sensor_view>);
        ;

    using element_sensor_view = SensorView<ModelPart::ElementsContainerType>;
    py::class_<element_sensor_view, element_sensor_view::Pointer>(sensor_module, "ElementSensorView")
        .def(py::init<Sensor::Pointer, const std::string&>(), py::arg("sensor"), py::arg("expression_name"))
        .def("GetSensor", &element_sensor_view::GetSensor)
        .def("GetContainerExpression", &element_sensor_view::GetContainerExpression)
        .def("__str__", PrintObject<element_sensor_view>);
        ;

    py::class_<MeasurementResidualResponseFunction, MeasurementResidualResponseFunction::Pointer, AdjointResponseFunction>(sensor_module, "MeasurementResidualResponseFunction")
        .def(py::init<>())
        .def("AddSensor", &MeasurementResidualResponseFunction::AddSensor, py::arg("sensor"))
        .def("Clear", &MeasurementResidualResponseFunction::Clear)
        .def("GetSensorsList", &MeasurementResidualResponseFunction::GetSensorsList)
        .def("__str__", PrintObject<MeasurementResidualResponseFunction>)
        ;

    py::class_<DisplacementSensor, DisplacementSensor::Pointer, Sensor>(sensor_module, "DisplacementSensor")
        .def(py::init<const std::string&,const Point&,const array_1d<double, 3>&,const Element&,const double>(),
            py::arg("name"),
            py::arg("location"),
            py::arg("direction"),
            py::arg("element"),
            py::arg("weight"))
        ;
}

} // namespace Kratos::Python
