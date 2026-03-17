//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    @{KRATOS_APP_AUTHOR}
//


// System includes

// External includes

// Project includes
#include "includes/define_python.h"
#include "sph_application.h"
#include "sph_application_variables.h"
#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_python/add_custom_processes_to_python.h"
#include "custom_python/add_custom_constitutive_laws_to_python.h"

namespace Kratos::Python {

PYBIND11_MODULE(KratosSPHApplication,m)
{
    namespace py = pybind11;

    py::class_<KratosSPHApplication,
        KratosSPHApplication::Pointer,
        KratosApplication>(m, "KratosSPHApplication")
        .def(py::init<>())
        ;

    AddCustomStrategiesToPython(m);
    AddCustomUtilitiesToPython(m);
    AddCustomProcessesToPython(m);
    AddCustomConstitutiveLawsToPython(m);

}

} // namespace Kratos::Python
