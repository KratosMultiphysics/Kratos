// KRATOS  / ___|___/ ___|(_)_ __ ___  _   _| | __ _| |_(_) ___  _ ___
//        | |   / _ \___ \| | '_ ` _ \| | | | |/ _` | __| |/ _ \| '_  |
//        | |__| (_) |__) | | | | | | | |_| | | (_| | |_| | (_) | | | |
//         \____\___/____/|_|_| |_| |_|\__,_|_|\__,_|\__|_|\___/|_| |_|
//
//  License:		 BSD License
//					 license: CoSimulationApplication/license.txt
//
//  Main authors:    Aditya Ghantasala
//                   Philipp Bucher
//

// System includes

#if defined(KRATOS_PYTHON)
// External includes
#include <pybind11/pybind11.h>

// Project includes
#include "includes/define.h"
#include "co_simulation_application_variables.h"
#include "co_simulation_application.h"
#include "custom_python/add_custom_io_to_python.h"

namespace Kratos {
namespace Python {

PYBIND11_MODULE(KratosCoSimulationApplication,m)
{
    namespace py = pybind11;

    py::class_<KratosCoSimulationApplication,
        KratosCoSimulationApplication::Pointer,
        KratosApplication>(m, "KratosCoSimulationApplication")
        .def(py::init<>())
        ;

    AddCustomIOToPython(m);

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SCALAR_DISPLACEMENT );
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SCALAR_ROOT_POINT_DISPLACEMENT );
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SCALAR_REACTION );
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SCALAR_FORCE );
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SCALAR_VOLUME_ACCELERATION );


}

} // namespace Python.
} // namespace Kratos.

#endif // KRATOS_PYTHON defined
