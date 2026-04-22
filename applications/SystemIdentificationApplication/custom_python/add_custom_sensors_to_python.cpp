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
        .def("AddTensorAdaptor", &Sensor::AddTensorAdaptor, py::arg("tensor_adaptor_name"), py::arg("tensor_adaptor"))
        .def("GetTensorAdaptor", &Sensor::GetTensorAdaptor, py::arg("tensor_adaptor_name"))
        .def("GetTensorAdaptorsMap", &Sensor::GetTensorAdaptorsMap)
        .def("ClearTensorAdaptors", &Sensor::ClearTensorAdaptors)
        .def("__str__", PrintObject<Sensor>);
        ;

    py::class_<SensorView, SensorView::Pointer>(sensor_module, "SensorView")
        .def(py::init<Sensor::Pointer, const std::string&>(), py::arg("sensor"), py::arg("tensor_adaptor_name"))
        .def("GetSensor", &SensorView::GetSensor)
        .def("GetTensorAdaptor", &SensorView::GetTensorAdaptor)
        .def("GetTensorAdaptorName", &SensorView::GetTensorAdaptorName)
        .def("AddAuxiliaryTensorAdaptor", &SensorView::AddAuxiliaryTensorAdaptor, py::arg("suffix"), py::arg("tensor_adaptor"))
        .def("GetAuxiliaryTensorAdaptor", &SensorView::GetAuxiliaryTensorAdaptor, py::arg("suffix"))
        .def("GetAuxiliarySuffixes", &SensorView::GetAuxiliarySuffixes)
        .def("__str__", PrintObject<SensorView>);
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
}

} // namespace Kratos::Python
