//
//   Project Name:        KratosContactMechanicsApplication $
//   Created by:          $Author:              JMCarbonell $
//   Last modified by:    $Co-Author:                       $
//   Date:                $Date:                  July 2016 $
//   Revision:            $Revision:                    0.0 $
//
//


// System includes

#if defined(KRATOS_PYTHON)
// External includes
#include <boost/python.hpp>


// Project includes
#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_python/add_custom_processes_to_python.h"
#include "custom_python/add_custom_bounding_to_python.h"
#include "custom_python/add_custom_modelers_to_python.h"
#include "custom_python/add_custom_friction_laws_to_python.h"

#include "contact_mechanics_application.h"

namespace Kratos
{

namespace Python
{

  using namespace boost::python;



  BOOST_PYTHON_MODULE(KratosContactMechanicsApplication)
  {

	  class_<KratosContactMechanicsApplication,
		 KratosContactMechanicsApplication::Pointer,
		 bases<KratosApplication>, boost::noncopyable >("KratosContactMechanicsApplication")
	    ;

	AddCustomStrategiesToPython();
	AddCustomUtilitiesToPython();
	AddCustomProcessesToPython();
        AddCustomBoundingToPython();
	AddCustomModelersToPython();
	AddCustomFrictionLawsToPython();

	//registering variables in python
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( FRICTION_ACTIVE )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( PENALTY_PARAMETER )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( LAGRANGE_MULTIPLIER_NORMAL )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( LAGRANGE_MULTIPLIER_NORMAL_REACTION )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( LAGRANGE_MULTIPLIER_TANGENTIAL )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( LAGRANGE_MULTIPLIER_TANGENTIAL_REACTION )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( TAU_STAB )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( MU_STATIC )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( MU_DYNAMIC )
        KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( CONTACT_STRESS )
        KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( EFFECTIVE_CONTACT_STRESS )
        KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( EFFECTIVE_CONTACT_FORCE )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( CONTACT_ADHESION )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( CONTACT_FRICTION_ANGLE )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( TANGENTIAL_PENALTY_RATIO )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( CONTACT_PLASTIC_SLIP )

	//KRATOS_REGISTER_IN_PYTHON_VARIABLE(NODAL_AREA)
	
	KRATOS_REGISTER_IN_PYTHON_VARIABLE( FRICTION_LAW )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( FRICTION_LAW_NAME )

  }


}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
