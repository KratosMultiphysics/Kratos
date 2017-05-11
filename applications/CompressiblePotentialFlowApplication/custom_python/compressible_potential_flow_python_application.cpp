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
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "compressible_potential_flow_application.h"
#include "compressible_potential_flow_application_variables.h"
#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_python/add_custom_processes_to_python.h"


namespace Kratos
{

namespace Python
{

  using namespace boost::python;



  BOOST_PYTHON_MODULE(KratosCompressiblePotentialFlowApplication)
  {

	  class_<KratosCompressiblePotentialFlowApplication,
			  KratosCompressiblePotentialFlowApplication::Pointer,
			  bases<KratosApplication>, boost::noncopyable >("KratosCompressiblePotentialFlowApplication")
			;

	AddCustomStrategiesToPython();
	AddCustomUtilitiesToPython();
        AddCustomProcessesToPython();
	//registering variables in python

//	KRATOS_REGISTER_IN_PYTHON_VARIABLE(NODAL_AREA);


  }


}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
