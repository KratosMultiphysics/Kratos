//   
//   Project Name:        Kratos       
//   Last modified by:    $Author: janosch $
//   Date:                $Date: 2008-07-23 14:56:07 $
//   Revision:            $Revision: 1.1 $
//
//

// System includes 

#if defined(KRATOS_PYTHON)
// External includes 
#include <boost/python.hpp>


// Project includes 
#include "includes/define.h"
#include "mkl_solvers_application.h"
#include "custom_python/add_linear_solvers_to_python.h"

 
namespace Kratos
{

namespace Python
{

  using namespace boost::python;


  
  BOOST_PYTHON_MODULE(KratosMKLSolversApplication)
  {

	  class_<KratosMKLSolversApplication, 
			  KratosMKLSolversApplication::Pointer, 
			  bases<KratosApplication>, boost::noncopyable >("KratosMKLSolversApplication")
			;

	AddLinearSolversToPython();


  }
  
  
}  // namespace Python.
  
}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
