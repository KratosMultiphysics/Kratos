//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Author1 Fullname
//                   Author2 Fullname
//


// System includes

#if defined(KRATOS_PYTHON)

// System includes
#include <pybind11/pybind11.h>

// Project includes
#include "includes/define.h"
#include "includes/define_python.h"

#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "fluid_transport_application.h"


namespace Kratos
{

namespace Python
{

namespace py = pybind11;

PYBIND11_MODULE(KratosFluidTransportApplication, m)
{
	py::class_<KratosFluidTransportApplication,
	KratosFluidTransportApplication::Pointer,
    KratosApplication>(m, "KratosFluidTransportApplication")
	.def(py::init<>());

    AddCustomStrategiesToPython(m);
    AddCustomUtilitiesToPython(m);

	//registering variables in python
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,PECLET);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,THETA);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,PHI_THETA);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m,PHI_GRADIENT);
}


}  // namespace Python.
}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
