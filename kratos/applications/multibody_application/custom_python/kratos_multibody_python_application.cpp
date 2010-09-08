//   
//   Project Name:        Kratos       
//   Last modified by:    $Author: anonymous $
//   Date:                $Date: 2008-11-05 17:07:44 $
//   Revision:            $Revision: 1.4 $
//
//

// System includes 

#if defined(KRATOS_PYTHON)
// External includes 
#include <boost/python.hpp>


// Project includes 
#include "includes/define.h"
#include "kratos_multibody_application.h"
#include "custom_python/add_custom_processes_to_python.h"


 
namespace Kratos
{

namespace Python
{

  using namespace boost::python;


  
  BOOST_PYTHON_MODULE(KratosMultibodyApplication)
  {

	  class_<KratosMultibodyApplication, 
			  KratosMultibodyApplication::Pointer, 
			  bases<KratosApplication>, boost::noncopyable >("KratosMultibodyApplication")
			;
		AddProcessesToPython();
		
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(IS_BODY)
		
  }
  
  
}  // namespace Python.
  
}  // namespace Kratos.

#endif // KRATOS_PYTHON defined



