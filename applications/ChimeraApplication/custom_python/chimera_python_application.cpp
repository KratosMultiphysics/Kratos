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
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "chimera_application.h"
#include "chimera_application_variables.h"
#include "custom_python/add_custom_processes_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_python/add_custom_strategies_to_python.h"




namespace Kratos
{

namespace Python
{

  using namespace boost::python;



  BOOST_PYTHON_MODULE(KratosChimeraApplication)
  {

	  class_<KratosChimeraApplication,
			  KratosChimeraApplication::Pointer,
			  bases<KratosApplication>, boost::noncopyable >("KratosChimeraApplication")
			;
		AddCustomProcessesToPython();			
		AddCustomUtilitiesToPython();
		AddCustomStrategiesToPython();

			//registering variables in python
  		KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( CHIM_NEUMANN_COND )

		//KRATOS_REGISTER_IN_PYTHON_VARIABLE( IS_WEAK )
	
		KRATOS_REGISTER_IN_PYTHON_VARIABLE(BOUNDARY_NODE);
	    KRATOS_REGISTER_IN_PYTHON_VARIABLE(FLUX);
	    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(TRACTION);
}


}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
