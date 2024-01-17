//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ignasi de Pouplana
//

#if defined(KRATOS_PYTHON)

// System includes
#include <pybind11/pybind11.h>

// Project includes
#include "includes/define.h"
#include "includes/define_python.h"

// Application includes
#include "custom_python/add_custom_processes_to_python.h"
#include "laserdrilling_application.h"

namespace Kratos
{

namespace Python
{

namespace py = pybind11;

PYBIND11_MODULE(KratosLaserDrillingApplication,m)
{

    py::class_<KratosLaserDrillingApplication,
           KratosLaserDrillingApplication::Pointer,
           KratosApplication >(m,"KratosLaserDrillingApplication")
           .def(py::init<>())
           ;
    AddCustomProcessesToPython(m);

    // Registering variables in python
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, NO2)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, TRANSPORT_VELOCITY)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, VECTOR_VELOCITAT_VENT)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, IS_POTENTIAL_FLOW_STEP)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, VELOCITAT_VENT)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, DIRECCIO_VENT)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, METEO_DIRECCIO_VENT)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, STARTING_DATE)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SIMULATION_DURATION_IN_DAYS)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, WIND_AUTOMATIC_PROCESS)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, POLLUTANT_AUTOMATIC_PROCESS)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, CITY)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, CASE_ID)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, IN_PRODUCTION)

}

}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
