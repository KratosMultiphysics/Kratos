//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

// System includes

#if defined(KRATOS_PYTHON)
// External includes
#include <pybind11/pybind11.h>


// Project includes
#include "includes/define.h"
#include "compressible_potential_flow_application.h"
#include "compressible_potential_flow_application_variables.h"
#include "custom_python/add_custom_processes_to_python.h"


namespace Kratos {
namespace Python {

PYBIND11_MODULE(KratosCompressiblePotentialFlowApplication,m)
{
	namespace py = pybind11;

	py::class_<KratosCompressiblePotentialFlowApplication,
		KratosCompressiblePotentialFlowApplication::Pointer,
		KratosApplication >(m,"KratosCompressiblePotentialFlowApplication")
		.def(py::init<>())
		;

	AddCustomProcessesToPython(m);

	//registering variables in python
	KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m,VELOCITY_INFINITY);


}

}  // namespace Python.
}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
