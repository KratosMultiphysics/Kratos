//
//  Main authors:    Miguel Angel Celigueta   maceli@cimne.upc.edu
//
//

// System includes

#if defined(KRATOS_PYTHON)
// External includes
#include <pybind11/pybind11.h>

// Project includes
#include "includes/define.h"
#include "../drilling_beam_application.h"
#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_python/add_custom_processes_to_python.h"


namespace Kratos
{

namespace Python
{

using namespace pybind11;


PYBIND11_MODULE(KratosDrillingBeamApplication,m)
{
	class_<KratosDrillingBeamApplication,
			KratosDrillingBeamApplication::Pointer,
			KratosApplication >(m,"KratosDrillingBeamApplication")
			.def(init<>());
	;

	AddCustomStrategiesToPython(m);
	AddCustomUtilitiesToPython(m);
	AddCustomProcessesToPython(m);

	//registering variables in python
	//KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, NODAL_AREA)



  }


}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
