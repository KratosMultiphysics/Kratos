//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Author Julio Marti
//


// System includes

#if defined(KRATOS_PYTHON)
// External includes


// Project includes
#include "includes/define.h"
#include "radiation_application.h"
#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"


namespace Kratos
{

namespace Python
{

  namespace py = pybind11;



  PYBIND11_MODULE(KratosRadiationApplication,m)
  {

	  py::class_<KratosRadiationApplication,
			  KratosRadiationApplication::Pointer,
			  KratosApplication>(m,"KratosRadiationApplication").def(init<>())
			;

	AddCustomStrategiesToPython(m);
	AddCustomUtilitiesToPython(m);

	//registering variables in python
//	KRATOS_REGISTER_IN_PYTHON_VARIABLE(NODAL_AREA);


  }


}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
