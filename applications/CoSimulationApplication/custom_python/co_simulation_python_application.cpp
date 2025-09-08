// KRATOS  / ___|___/ ___|(_)_ __ ___  _   _| | __ _| |_(_) ___  _ ___
//        | |   / _ \___ \| | '_ ` _ \| | | | |/ _` | __| |/ _ \| '_  |
//        | |__| (_) |__) | | | | | | | |_| | | (_| | |_| | (_) | | | |
//         \____\___/____/|_|_| |_| |_|\__,_|_|\__,_|\__|_|\___/|_| |_|
//
//  License:         BSD License
//                   license: CoSimulationApplication/license.txt
//
//  Main authors:    Aditya Ghantasala 
//                   Philipp Bucher (https://github.com/philbucher)
//

// System includes

#if defined(KRATOS_PYTHON)
// External includes
#include <pybind11/pybind11.h>

// Project includes
#include "includes/define.h"
#include "co_simulation_application_variables.h"
#include "co_simulation_application.h"
#include "custom_python/add_co_sim_io_to_python.h"
#include "custom_python/add_custom_io_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_python/add_custom_processes_to_python.h"

namespace Kratos::Python {

PYBIND11_MODULE(KratosCoSimulationApplication,m)
{
    namespace py = pybind11;

    py::class_<KratosCoSimulationApplication,
        KratosCoSimulationApplication::Pointer,
        KratosApplication>(m, "KratosCoSimulationApplication")
        .def(py::init<>())
        ;

    AddCoSimIOToPython(m);
    AddCustomIOToPython(m);
    AddCustomUtilitiesToPython(m);
    AddCustomProcessesToPython(m);

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SCALAR_DISPLACEMENT );
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SCALAR_ROOT_POINT_DISPLACEMENT );
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SCALAR_REACTION );
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SCALAR_FORCE );
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SCALAR_VOLUME_ACCELERATION );

    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, RESULTANT_FORCE);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, RESULTANT_MOMENT);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, GLOBAL_DISPLACEMENT);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, GLOBAL_ROTATION);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, IMPOSED_FORCE);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, IMPOSED_MOMENT);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, IMPOSED_DISPLACEMENT);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, IMPOSED_ROTATION);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, EFFECTIVE_FORCE);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, EFFECTIVE_MOMENT);

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, COUPLING_ITERATION_NUMBER);
}

} // namespace Kratos::Python.

#endif // KRATOS_PYTHON defined
