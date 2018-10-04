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

    AddCustomProcessesToPython(m);
	//registering variables in python

	KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,POSITIVE_POTENTIAL);
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,NEGATIVE_POTENTIAL);
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,WAKE_DISTANCE);
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,WAKE_ELEMENTAL_DISTANCES);
	KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m,VELOCITY_INFINITY);


  }


}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
