// KRATOS  __  __          _    _                _ _           _   _
//        |  \/  | ___  __| |  / \   _ __  _ __ | (_) ___ __ _| |_(_) ___  _ ___
//        | |\/| |/ _ \/ _` | / _ \ | '_ \| '_ \| | |/ __/ _` | __| |/ _ \| '_  |
//        | |  | |  __/ (_| |/ ___ \| |_) | |_) | | | (_| (_| | |_| | (_) | | | |
//        |_|  |_|\___|\__,_/_/   \_\ .__/| .__/|_|_|\___\__,_|\__|_|\___/|_| |_|
//                                  |_|   |_|
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
#include "custom_python/add_custom_utilities_to_python.h"


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
    AddCustomUtilitiesToPython(m);
}

} // namespace Kratos::Python
