//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//


// System includes

// External includes

// Project includes
#include "includes/define_python.h"
#include "med_application.h"
#include "custom_python/add_custom_io_to_python.h"


namespace Kratos::Python {

PYBIND11_MODULE(KratosMedApplication,m)
{
    namespace py = pybind11;

    py::class_<KratosMedApplication,
        KratosMedApplication::Pointer,
        KratosApplication>(m, "KratosMedApplication")
        .def(py::init<>())
        ;

    AddCustomIOToPython(m);
}

} // namespace Kratos::Python
