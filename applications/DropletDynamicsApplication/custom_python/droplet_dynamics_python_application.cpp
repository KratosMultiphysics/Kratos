//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Mohammad R. Hashemi
//


// System includes

#if defined(KRATOS_PYTHON)
// External includes
#include <pybind11/pybind11.h>


// Project includes
#include "includes/define_python.h"
#include "droplet_dynamics_application.h"
#include "droplet_dynamics_application_variables.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_python/add_custom_processes_to_python.h"


namespace Kratos {
namespace Python {

PYBIND11_MODULE(KratosDropletDynamicsApplication,m)
{
    namespace py = pybind11;

    py::class_<KratosDropletDynamicsApplication,
        KratosDropletDynamicsApplication::Pointer,
        KratosApplication>(m, "KratosDropletDynamicsApplication")
        .def(py::init<>())
        ;

    AddCustomUtilitiesToPython(m);
    AddCustomProcessesToPython(m);

    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m,EXT_INT_FORCE)

    //registering variables in python
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, EPOTENTIAL )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PERMITTIVITYPOS )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PERMITTIVITYNEG )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, CONDUCTIVITYPOS )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, CONDUCTIVITYNEG )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, EFIELD )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, SCHARGE )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, EFORCE )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, EFIELDNEG )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, EFIELDPOS )

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, INV_K_ENRICH )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, BIJ_ENRICH_ROW )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, BJI_ENRICH_ROW )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, POS_GRAD_ENRICH )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, NEG_GRAD_ENRICH )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, RHS_ENRICH )


    //	KRATOS_REGISTER_IN_PYTHON_VARIABLE(NODAL_AREA);

}

} // namespace Python.
} // namespace Kratos.

#endif // KRATOS_PYTHON defined
