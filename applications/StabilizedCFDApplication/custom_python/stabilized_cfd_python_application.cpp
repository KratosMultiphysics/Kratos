//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//

// System includes

#if defined(KRATOS_PYTHON)
// External includes
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "stabilized_cfd_application.h"
#include "stabilized_cfd_application_variables.h"
#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"


namespace Kratos
{

namespace Python
{

  using namespace boost::python;



  BOOST_PYTHON_MODULE(KratosStabilizedCFDApplication)
  {

	  class_<KratosStabilizedCFDApplication,
			  KratosStabilizedCFDApplication::Pointer,
			  bases<KratosApplication>, boost::noncopyable >("KratosStabilizedCFDApplication")
			;

	AddCustomStrategiesToPython();
	AddCustomUtilitiesToPython();

	//registering variables in python
  KRATOS_REGISTER_IN_PYTHON_VARIABLE( FIC_BETA )
  KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( DIRECTIONAL_BETA )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE( RECORDED_STEPS )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE( MEAN_KINETIC_ENERGY )
  KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( MEAN_VELOCITY )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE( MEAN_PRESSURE )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE( VELOCITY_COVARIANCES )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE( TURBULENCE_STATISTICS )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE( TRACE_XI )
  KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( DIV_XI )
  KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( MOMENTUM_PROJECTION )
  KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( MOMENTUM_PROJECTION_RHS )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE( MASS_PROJECTION )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE( MASS_PROJECTION_RHS )

//	KRATOS_REGISTER_IN_PYTHON_VARIABLE(NODAL_AREA);


  }


}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
