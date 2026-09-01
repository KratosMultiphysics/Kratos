//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         SystemIdentificationApplication/license.txt
//
//  Main authors:    Máté Kelemen
//

// --- External Includes ---
#include <pybind11/stl.h>

// --- SysId Includes ---
#include "custom_python/add_custom_adjoints_to_python.hpp"
#include "custom_adjoint/displacement_sensor_response.hpp"
#include "custom_adjoint/sensor_aggregate_response.hpp"


namespace Kratos::Python {


void AddCustomAdjointsToPython(pybind11::module& rModule) {
    pybind11::class_<
        DisplacementSensorResponse,
        DisplacementSensorResponse::Pointer,
        SensorResponse
    >(rModule, "DisplacementSensorResponse")
        .def(pybind11::init<>())
        .def(
            "Create",
            [] (const DisplacementSensorResponse& rThis,
                const ModelPart& rDomainModelPart,
                const ModelPart& rSensorModelPart,
                std::size_t Id,
                Parameters SensorParameters) -> SensorResponse::Pointer {
                    return rThis.Create(
                        rDomainModelPart,
                        rSensorModelPart,
                        Id,
                        SensorParameters);
                },
            pybind11::arg("rDomainModelPart"),
            pybind11::arg("rSensorModelPart"),
            pybind11::arg("Id"),
            pybind11::arg("SensorParameters"))
        ;

    pybind11::class_<
        SensorAggregateResponse,
        SensorAggregateResponse::Pointer,
        ResponseFunction
    >(rModule, "SensorAggregateResponse")
        .def(pybind11::init<>())
        .def(pybind11::init<std::size_t>())
        .def(
            "AddSensors",
            [] (SensorAggregateResponse& rThis, const std::vector<SensorResponse::Pointer>& rSensors) -> void {
                return rThis.AddSensors({rSensors.data(), rSensors.size()});
            })
        ;
}


} // namespace Kratos::Python
