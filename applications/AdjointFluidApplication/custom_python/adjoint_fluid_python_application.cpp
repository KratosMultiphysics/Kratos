//
//   Project Name:        KratosAdjointFluidApplication $
//   Last modified by:    $Author: michael.andre@tum.de $
//   Date:                $Date:          February 2015 $
//   Revision:            $Revision:                0.0 $
//
//

// System includes
#if defined(KRATOS_PYTHON)
// External includes
#include <boost/python.hpp>

// Application includes
#include "adjoint_fluid_application.h"
#include "custom_python/add_custom_schemes_to_python.h"
#include "custom_python/add_custom_response_functions_to_python.h"

namespace Kratos
{
  
namespace Python
{

using namespace boost::python;
  
BOOST_PYTHON_MODULE(KratosAdjointFluidApplication)
{
  class_<KratosAdjointFluidApplication,
	 KratosAdjointFluidApplication::Pointer,
	 bases<KratosApplication>, boost::noncopyable >("KratosAdjointFluidApplication");

  AddCustomSchemesToPython();
  AddCustomResponseFunctionsToPython();

}

} // namespace Python

} // namespace Kratos

#endif // KRATOS_PYTHON defined
