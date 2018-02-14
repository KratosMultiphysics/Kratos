//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Author1 Fullname
//                   Author2 Fullname 
//


// System includes

#if defined(KRATOS_PYTHON)
// External includes
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "my_laplacian_application.h"
#include "my_laplacian_application_variables.h"
#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"


namespace Kratos
{

namespace Python
{

  using namespace boost::python;



  BOOST_PYTHON_MODULE(KratosMyLaplacianApplication)
  {

	  class_<KratosMyLaplacianApplication,
			  KratosMyLaplacianApplication::Pointer,
			  bases<KratosApplication>, boost::noncopyable >("KratosMyLaplacianApplication")
			;

	AddCustomStrategiesToPython();
	AddCustomUtilitiesToPython();

	//registering variables in python
  KRATOS_REGISTER_IN_PYTHON_VARIABLE( MY_SCALAR )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE( MY_BOOL )
  KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( MY_VECTOR )

//	KRATOS_REGISTER_IN_PYTHON_VARIABLE(NODAL_AREA);


  }


}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
