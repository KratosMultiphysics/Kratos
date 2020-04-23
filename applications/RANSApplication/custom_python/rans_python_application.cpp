//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya (https://github.com/sunethwarna)
//


// System includes

#if defined(KRATOS_PYTHON)
// External includes
#include <pybind11/pybind11.h>


// Project includes
#include "includes/define.h"
#include "rans_application.h"
#include "rans_application_variables.h"
#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_python/add_custom_processes_to_python.h"


namespace Kratos {
namespace Python {

PYBIND11_MODULE(KratosRANSApplication,m)
{
    namespace py = pybind11;

    py::class_<KratosRANSApplication,
        KratosRANSApplication::Pointer,
        KratosApplication>(m, "KratosRANSApplication")
        .def(py::init<>())
        ;

    AddCustomStrategiesToPython(m);
    AddCustomUtilitiesToPython(m);
    AddCustomProcessesToPython(m);

    //registering variables in python
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, TURBULENT_KINETIC_ENERGY )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, TURBULENT_ENERGY_DISSIPATION_RATE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, TURBULENT_KINETIC_ENERGY_RATE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, TURBULENT_ENERGY_DISSIPATION_RATE_2 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, IS_CO_SOLVING_PROCESS_ACTIVE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, RANS_Y_PLUS )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, RANS_AUXILIARY_VARIABLE_1 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, RANS_AUXILIARY_VARIABLE_2 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, WALL_SMOOTHNESS_BETA )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, WALL_VON_KARMAN )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, TURBULENCE_RANS_C_MU )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, TURBULENCE_RANS_C1 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, TURBULENCE_RANS_C2 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, TURBULENT_KINETIC_ENERGY_SIGMA )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, TURBULENT_ENERGY_DISSIPATION_RATE_SIGMA )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, NUMBER_OF_NEIGHBOUR_CONDITIONS )

    //registering variables in python
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, RANS_NUT_SCALAR_PARTIAL_DERIVATIVES)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m , RANS_SCALAR_1_ADJOINT_1 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m , RANS_SCALAR_1_ADJOINT_2 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m , RANS_SCALAR_1_ADJOINT_3 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m , RANS_AUX_ADJOINT_SCALAR_1 )

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m , RANS_SCALAR_2_ADJOINT_1 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m , RANS_SCALAR_2_ADJOINT_2 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m , RANS_SCALAR_2_ADJOINT_3 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m , RANS_AUX_ADJOINT_SCALAR_2 )
}

} // namespace Python.
} // namespace Kratos.

#endif // KRATOS_PYTHON defined
