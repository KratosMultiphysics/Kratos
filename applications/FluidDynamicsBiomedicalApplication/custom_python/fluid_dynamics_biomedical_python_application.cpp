//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Eduardo Soudah
//

// System includes

#if defined(KRATOS_PYTHON)
// External includes
#include <pybind11/pybind11.h>

// Project includes
#include "includes/define_python.h"

// Application includes
#include "fluid_dynamics_biomedical_application.h"
#include "fluid_dynamics_biomedical_application_variables.h"
#include "custom_python/add_custom_utilities_to_python.h"

// Project includes

namespace Kratos
{

namespace Python
{

PYBIND11_MODULE(KratosFluidDynamicsBiomedicalApplication,m)
{
    namespace py = pybind11;

    py::class_<KratosFluidDynamicsBiomedicalApplication,
           KratosFluidDynamicsBiomedicalApplication::Pointer,
           KratosApplication >(m,"KratosFluidDynamicsBiomedicalApplication")
           .def(py::init<>())
           ;

    AddCustomUtilitiesToPython(m);

    // WSS variables
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, TAWSS);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, TWSS);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, ECAP);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, RRT);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, OSI);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, WSS);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, WALL_DISTANCE);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, WSS_TANGENTIAL_STRESS);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, WSS_NORMAL_STRESS);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, TEMPORAL_OSI);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, WALL_NORMAL);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, OUTLET_NORMAL);

}

}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_PYTHON defined