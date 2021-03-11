// KRATOS ___ ___  _  ___   __   ___ ___ ___ ___
//       / __/ _ \| \| \ \ / /__|   \_ _| __| __|
//      | (_| (_) | .` |\ V /___| |) | || _|| _|
//       \___\___/|_|\_| \_/    |___/___|_| |_|  APPLICATION
//
//  License:         BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Julio Marti
//

// System includes

#if defined(KRATOS_PYTHON)
// External includes

// Project includes
#include "includes/define_python.h"
#include "pfem_melting_application.h"
//#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"
//#include "custom_python/add_custom_response_functions_to_python.h"

namespace Kratos
{

namespace Python
{

namespace py = pybind11;

PYBIND11_MODULE(KratosPfemMeltingApplication,m)
{

    py::class_<KratosPfemMeltingApplication,
           KratosPfemMeltingApplication::Pointer,
           KratosApplication >(m,"KratosPfemMeltingApplication")
           .def(py::init<>())
           ;
    //AddCustomStrategiesToPython(m);
    AddCustomUtilitiesToPython(m);
    //AddCustomResponseFunctionsToPython(m);

    // Registering variables in python

}

}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
