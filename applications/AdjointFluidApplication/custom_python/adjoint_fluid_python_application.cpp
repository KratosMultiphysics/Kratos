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
#include "custom_python/add_custom_processes_to_python.h"

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
  AddCustomProcessesToPython();

  // Moved to Kratos Core for trilinos_application
  //KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( ADJOINT_VELOCITY );
  //KRATOS_REGISTER_IN_PYTHON_VARIABLE( ADJOINT_PRESSURE );
  //KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( PRIMAL_VELOCITY );
  //KRATOS_REGISTER_IN_PYTHON_VARIABLE( PRIMAL_PRESSURE );
  //KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( SHAPE_SENSITIVITY );
  //KRATOS_REGISTER_IN_PYTHON_VARIABLE( NORMAL_SENSITIVITY );
}

} // namespace Python

} // namespace Kratos

#endif // KRATOS_PYTHON defined
