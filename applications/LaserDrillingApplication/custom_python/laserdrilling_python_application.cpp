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
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, THERMAL_ENERGY)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, THERMAL_ENERGY_PER_VOLUME)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, THERMAL_COUNTER)

}

}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
