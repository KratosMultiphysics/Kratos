//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         SystemIdentificationApplication/license.txt
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
#include "custom_sensors/displacement_sensor.h"
#include "custom_sensors/strain_sensor.h"
#include "custom_sensors/velocity_sensor.h"
#include "custom_sensors/pressure_sensor.h"

// Include base h
#include "custom_python/add_custom_sensors_to_python.h"

namespace Kratos::Python {

void  AddCustomSensorsToPython(pybind11::module& m)
{
    namespace py = pybind11;

    auto sensor_module = m.def_submodule("Sensors");

    // Add sensor specifications
    py::class_<Sensor, Sensor::Pointer, AdjointResponseFunction>(sensor_module, "Sensor")
        .def("GetName", &Sensor::GetName)
        .def("GetNode", &Sensor::GetNode)
        .def("GetWeight", &Sensor::GetWeight)
        .def("GetSensorValue", &Sensor::GetSensorValue)
        .def("SetSensorValue", &Sensor::SetSensorValue)
        .def("GetSensorParameters", &Sensor::GetSensorParameters)
        .def("AddContainerExpression", &Sensor::AddContainerExpression, py::arg("expression_name"), py::arg("expression"))
        .def("GetContainerExpression", &Sensor::GetContainerExpression, py::arg("expression_name"))
        .def("GetContainerExpressionsMap", &Sensor::GetContainerExpressionsMap)
        .def("ClearContainerExpressions", &Sensor::ClearContainerExpressions)
        .def("__str__", PrintObject<Sensor>);
        ;

    using nodal_sensor_view = SensorView<ModelPart::NodesContainerType>;
    py::class_<nodal_sensor_view, nodal_sensor_view::Pointer>(sensor_module, "NodalSensorView")
        .def(py::init<Sensor::Pointer, const std::string&>(), py::arg("sensor"), py::arg("expression_name"))
        .def("GetSensor", &nodal_sensor_view::GetSensor)
        .def("GetContainerExpression", &nodal_sensor_view::GetContainerExpression)
        .def("GetExpressionName", &nodal_sensor_view::GetExpressionName)
        .def("AddAuxiliaryExpression", &nodal_sensor_view::AddAuxiliaryExpression, py::arg("suffix"), py::arg("nodal_expression"))
        .def("GetAuxiliaryExpression", &nodal_sensor_view::GetAuxiliaryExpression, py::arg("suffix"))
        .def("GetAuxiliarySuffixes", &nodal_sensor_view::GetAuxiliarySuffixes)
        .def("__str__", PrintObject<nodal_sensor_view>);
        ;

    using condition_sensor_view = SensorView<ModelPart::ConditionsContainerType>;
    py::class_<condition_sensor_view, condition_sensor_view::Pointer>(sensor_module, "ConditionSensorView")
        .def(py::init<Sensor::Pointer, const std::string&>(), py::arg("sensor"), py::arg("expression_name"))
        .def("GetSensor", &condition_sensor_view::GetSensor)
        .def("GetContainerExpression", &condition_sensor_view::GetContainerExpression)
        .def("GetExpressionName", &condition_sensor_view::GetExpressionName)
        .def("AddAuxiliaryExpression", &condition_sensor_view::AddAuxiliaryExpression, py::arg("suffix"), py::arg("condition_expression"))
        .def("GetAuxiliaryExpression", &condition_sensor_view::GetAuxiliaryExpression, py::arg("suffix"))
        .def("GetAuxiliarySuffixes", &condition_sensor_view::GetAuxiliarySuffixes)
        .def("__str__", PrintObject<condition_sensor_view>);
        ;

    using element_sensor_view = SensorView<ModelPart::ElementsContainerType>;
    py::class_<element_sensor_view, element_sensor_view::Pointer>(sensor_module, "ElementSensorView")
        .def(py::init<Sensor::Pointer, const std::string&>(), py::arg("sensor"), py::arg("expression_name"))
        .def("GetSensor", &element_sensor_view::GetSensor)
        .def("GetContainerExpression", &element_sensor_view::GetContainerExpression)
        .def("GetExpressionName", &element_sensor_view::GetExpressionName)
        .def("AddAuxiliaryExpression", &element_sensor_view::AddAuxiliaryExpression, py::arg("suffix"), py::arg("element_expression"))
        .def("GetAuxiliaryExpression", &element_sensor_view::GetAuxiliaryExpression, py::arg("suffix"))
        .def("GetAuxiliarySuffixes", &element_sensor_view::GetAuxiliarySuffixes)
        .def("__str__", PrintObject<element_sensor_view>);
        ;

    py::class_<DisplacementSensor, DisplacementSensor::Pointer, Sensor>(sensor_module, "DisplacementSensor")
        .def(py::init<const std::string&,Node::Pointer,const array_1d<double, 3>&,const Element&,const double>(),
            py::arg("name"),
            py::arg("node"),
            py::arg("direction"),
            py::arg("element"),
            py::arg("weight"))
        .def_static("GetDefaultParameters", &DisplacementSensor::GetDefaultParameters)
        .def_static("Create", &DisplacementSensor::Create, py::arg("domain_model_part"), py::arg("sensor_model_part"), py::arg("sensor_id"), py::arg("sensor_parameters"))
        ;

    auto strain_sensor = py::class_<StrainSensor, StrainSensor::Pointer, Sensor>(sensor_module, "StrainSensor");
    py::enum_<StrainSensor::StrainType>(strain_sensor, "StrainType")
        .value("STRAIN_XX", StrainSensor::StrainType::STRAIN_XX)
        .value("STRAIN_YY", StrainSensor::StrainType::STRAIN_YY)
        .value("STRAIN_ZZ", StrainSensor::StrainType::STRAIN_ZZ)
        .value("STRAIN_XY", StrainSensor::StrainType::STRAIN_XY)
        .value("STRAIN_XZ", StrainSensor::StrainType::STRAIN_XZ)
        .value("STRAIN_YZ", StrainSensor::StrainType::STRAIN_YZ)
        .export_values();
    strain_sensor
        .def(py::init<const std::string&,Node::Pointer, const Variable<Matrix>&, const StrainSensor::StrainType&, const Element&,const double>(),
            py::arg("name"),
            py::arg("node"),
            py::arg("strain_variable"),
            py::arg("strain_type"),
            py::arg("element"),
            py::arg("weight"))
        .def_static("GetDefaultParameters", &StrainSensor::GetDefaultParameters)
        .def_static("Create", &StrainSensor::Create, py::arg("domain_model_part"), py::arg("sensor_model_part"), py::arg("sensor_id"), py::arg("sensor_parameters"))
        ;

    py::class_<VelocitySensor, VelocitySensor::Pointer, Sensor>(sensor_module, "VelocitySensor")
        .def(py::init<const std::string&,Node::Pointer,const array_1d<double, 3>&,const Element&,const double>(),
            py::arg("name"),
            py::arg("node"),
            py::arg("direction"),
            py::arg("element"),
            py::arg("weight"))
        .def_static("GetDefaultParameters", &VelocitySensor::GetDefaultParameters)
        .def_static("Create", &VelocitySensor::Create, py::arg("domain_model_part"), py::arg("sensor_model_part"), py::arg("sensor_id"), py::arg("sensor_parameters"))
        ;

    py::class_<PressureSensor, PressureSensor::Pointer, Sensor>(sensor_module, "PressureSensor")
        .def(py::init<const std::string&,Node::Pointer,const Element&,const double>(),
            py::arg("name"),
            py::arg("node"),
            py::arg("element"),
            py::arg("weight"))
        .def_static("GetDefaultParameters", &PressureSensor::GetDefaultParameters)
        .def_static("Create", &PressureSensor::Create, py::arg("domain_model_part"), py::arg("sensor_model_part"), py::arg("sensor_id"), py::arg("sensor_parameters"))
        ;

}

} // namespace Kratos::Python
