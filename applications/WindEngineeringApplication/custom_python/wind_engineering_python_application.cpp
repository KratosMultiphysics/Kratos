//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Máté Kelemen
//

#ifdef KRATOS_PYTHON

// External includes
#include <pybind11/pybind11.h>

// Project includes
#include "includes/define.h"
#include "wind_engineering_application.h"
#include "add_custom_utilities_to_python.h"
#include "add_custom_processes_to_python.h"


namespace Kratos
{
namespace Python
{


PYBIND11_MODULE(KratosWindEngineeringApplication, module)
{
    pybind11::class_<KratosWindEngineeringApplication,
                     KratosWindEngineeringApplication::Pointer,
                     KratosApplication>(module, "KratosWindEngineeringApplication")
        .def(pybind11::init<>())
        ;

    AddCustomUtilitiesToPython(module);
    AddCustomProcessesToPython(module);

    // Register custom variables

} // PYBIND11_MODULE


} // namespace Python
} // namespace Kratos

#endif // KRATOS_PYTHON