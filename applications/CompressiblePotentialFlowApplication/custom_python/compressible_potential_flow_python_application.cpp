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
//#include <boost/python.hpp>
#include <pybind11/pybind11.h>


// Project includes
#include "includes/define_python.h"
#include "compressible_potential_flow_application.h"
#include "compressible_potential_flow_application_variables.h"
#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_python/add_custom_processes_to_python.h"


namespace Kratos
{

namespace Python
{

  using namespace pybind11;



  PYBIND11_MODULE(KratosCompressiblePotentialFlowApplication,m)
  {

	  class_<KratosCompressiblePotentialFlowApplication,
			  KratosCompressiblePotentialFlowApplication::Pointer,
			  KratosApplication >(m,"KratosCompressiblePotentialFlowApplication")
			  .def(init<>())
			;

	AddCustomStrategiesToPython(m);
	AddCustomUtilitiesToPython(m);
	AddCustomProcessesToPython(m);
	//registering variables in python

//	KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,NODAL_AREA);
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, UPPER_SURFACE )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, LOWER_SURFACE )
  KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, WAKE_NORMAL )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PROJECTION_MATRIX )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, UPPER_PROJECTION )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, LOWER_PROJECTION )

  }


}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
