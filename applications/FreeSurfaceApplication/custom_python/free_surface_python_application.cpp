//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Antonia Larese
//


// System includes

#if defined(KRATOS_PYTHON)
// External includes
#include <pybind11/pybind11.h>


// Project includes
#include "includes/define_python.h"
#include "free_surface_application.h"
#include "free_surface_application_variables.h"
#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_python/add_custom_edgebased_levelset_solver_to_python.h"
#include "custom_python/add_custom_io_to_python.h"


namespace Kratos
{

namespace Python
{

  namespace py = pybind11;



  PYBIND11_MODULE(KratosFreeSurfaceApplication, pymodule)
  {

	  py::class_<KratosFreeSurfaceApplication,
			  KratosFreeSurfaceApplication::Pointer,
			  KratosApplication >(pymodule,"KratosFreeSurfaceApplication")
                          .def(py::init<>())
			;

	AddCustomStrategiesToPython(pymodule);
	AddCustomUtilitiesToPython(pymodule);
	AddCustomIOToPython(pymodule);
	AddCustomEdgeBasedLevelSetToPython(pymodule);

	//registering variables in python


//	KRATOS_REGISTER_IN_PYTHON_VARIABLE(NODAL_AREA);


  }


}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
