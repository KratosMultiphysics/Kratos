//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//

// System includes

#if defined(KRATOS_PYTHON)
// External includes
#include "pybind11/pybind11.h"


// Project includes
#include "includes/define_python.h"
#include "hdf5_application.h"
#include "hdf5_application_variables.h"
#include "custom_python/add_custom_io_to_python.h"
#include "custom_python/add_custom_processes_to_python.h"

namespace Kratos {
namespace Python {

PYBIND11_MODULE(KratosHDF5Application,m)
{
    namespace py = pybind11;

    py::class_<KratosHDF5Application,
        KratosHDF5Application::Pointer,
        KratosApplication >(m,"KratosHDF5Application")
        .def(py::init<>())
        ;

    AddCustomIOToPython(m);
    AddCustomProcessesToPython(m);
}


}  // namespace Python.
}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
