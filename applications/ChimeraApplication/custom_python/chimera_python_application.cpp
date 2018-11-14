//
// ==============================================================================
//  ChimeraApplication
//
//  License:         BSD License
//                   license: ChimeraApplication/license.txt
//
//  Main authors:    Aditya Ghantasala, https://github.com/adityaghantasala
//                   Navaneeth K Narayanan
//
// ==============================================================================
// System includes

#if defined(KRATOS_PYTHON)
// External includes

// Project includes
#include "includes/define_python.h"
#include "chimera_application.h"
#include "chimera_application_variables.h"
#include "custom_python/add_custom_processes_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_python/add_custom_strategies_to_python.h"


namespace Kratos
{

namespace Python
{

  using namespace pybind11;

PYBIND11_MODULE(KratosChimeraApplication, m)
{
    class_<KratosChimeraApplication,
		KratosChimeraApplication::Pointer,
		KratosApplication >(m,"KratosChimeraApplication")
		.def(init<>())
		;

	AddCustomProcessesToPython(m);
	AddCustomUtilitiesToPython(m);
	AddCustomStrategiesToPython(m);

	//registering variables in python
	KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, CHIM_NEUMANN_COND )

	//KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, IS_WEAK )

	KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, BOUNDARY_NODE);
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, FLUX);
	KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, TRACTION);
	KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, SHEAR_FORCE);
	KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, PRESSURE_FORCE);
}

}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
