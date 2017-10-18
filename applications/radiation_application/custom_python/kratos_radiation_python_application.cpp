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
#include <boost/python.hpp>


// Project includes 
#include "includes/define.h"
#include "radiation.h"
#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"

 
namespace Kratos
{

namespace Python
{

  using namespace boost::python;


  
  BOOST_PYTHON_MODULE(KratosRadiationApplication)
  {

	  class_<KratosRadiationApplication, 
			  KratosRadiationApplication::Pointer, 
			  bases<KratosApplication>, boost::noncopyable >("KratosRadiationApplication")
			;
		AddCustomStrategiesToPython();
		AddCustomUtilitiesToPython();

		//registering variables in python

		KRATOS_REGISTER_IN_PYTHON_VARIABLE( TEMP_CONV_PROJ)
		KRATOS_REGISTER_IN_PYTHON_VARIABLE( EMISSIVITY)
		KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(CONVECTION_VELOCITY)
			 
  }
  
  
}  // namespace Python.
  
}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
