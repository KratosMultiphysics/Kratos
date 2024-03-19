//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//


// System includes

// External includes

// Project includes
#include "includes/define_python.h"
#include "fluid_dynamics_hydraulics_application.h"
#include "fluid_dynamics_hydraulics_application_variables.h"
#include "custom_python/add_custom_utilities_to_python.h"


namespace Kratos::Python {

PYBIND11_MODULE(KratosFluidDynamicsHydraulicsApplication,m)
{
    namespace py = pybind11;

    py::class_<KratosFluidDynamicsHydraulicsApplication,
        KratosFluidDynamicsHydraulicsApplication::Pointer,
        KratosApplication>(m, "KratosFluidDynamicsHydraulicsApplication")
        .def(py::init<>())
        ;

    AddCustomUtilitiesToPython(m);

    //registering variables in python

    //	KRATOS_REGISTER_IN_PYTHON_VARIABLE(NODAL_AREA);

}

} // namespace Kratos::Python
