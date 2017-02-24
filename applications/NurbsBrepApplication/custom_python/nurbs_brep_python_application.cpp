//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/NurbsBrepApplication/license.txt
//
//  Main authors:    Tobias Tescheamacher
//                   Riccardo Rossi
//

// System includes

#if defined(KRATOS_PYTHON)
// External includes
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "nurbs_brep_application.h"
#include "nurbs_brep_application_variables.h"
#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"


namespace Kratos
{

namespace Python
{

  using namespace boost::python;



  BOOST_PYTHON_MODULE(KratosNurbsBrepApplication)
  {

	  class_<KratosNurbsBrepApplication,
			  KratosNurbsBrepApplication::Pointer,
			  bases<KratosApplication>, boost::noncopyable >("KratosNurbsBrepApplication")
			;

	AddCustomStrategiesToPython();
	AddCustomUtilitiesToPython();

	//registering variables in python

	//To enhance weighting to nodes:
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(CONTROL_POINT_WEIGHT)


	//KRATOS_REGISTER_IN_PYTHON_VARIABLE(SHAPE_FUNCTION_VALUES)
	//KRATOS_REGISTER_IN_PYTHON_VARIABLE(SHAPE_FUNCTION_LOCAL_DERIVATIVES)
	//KRATOS_REGISTER_IN_PYTHON_VARIABLE(SHAPE_FUNCTION_LOCAL_DERIVATIVES_MASTER)
	//KRATOS_REGISTER_IN_PYTHON_VARIABLE(SHAPE_FUNCTION_LOCAL_DERIVATIVES_SLAVE)

	//KRATOS_REGISTER_IN_PYTHON_VARIABLE(TANGENTS)

	//KRATOS_REGISTER_IN_PYTHON_VARIABLE(PENALTY_FACTOR)
	//
	//KRATOS_REGISTER_IN_PYTHON_VARIABLE(DISPLACEMENT_ROTATION_FIX)

	//KRATOS_REGISTER_IN_PYTHON_VARIABLE(LOAD_TYPE)
	//KRATOS_REGISTER_IN_PYTHON_VARIABLE(DISTRIBUTED_LOAD_FACTOR)

  }


}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
