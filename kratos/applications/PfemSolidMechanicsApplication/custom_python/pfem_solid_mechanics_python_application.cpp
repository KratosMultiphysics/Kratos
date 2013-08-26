//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Last modified by:    $Author:                JMCarbonell $
//   Date:                $Date:                    July 2013 $
//   Revision:            $Revision:                      0.0 $
//
//

// System includes 
#if defined(KRATOS_PYTHON)
// External includes 
#include <boost/python.hpp>

// Project includes 
#include "includes/define.h"
#include "pfem_solid_mechanics_application.h"
#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_python/add_custom_processes_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_python/add_custom_conditions_to_python.h"
#include "custom_python/add_custom_constitutive_laws_to_python.h"
#include "custom_python/add_custom_modelers_to_python.h"
 
namespace Kratos
{

  namespace Python
  {

    using namespace boost::python;


  
    BOOST_PYTHON_MODULE(KratosPfemSolidMechanicsApplication)
    {

      class_<KratosPfemSolidMechanicsApplication, 
	     KratosPfemSolidMechanicsApplication::Pointer, 
	     bases<KratosApplication>, boost::noncopyable >("KratosPfemSolidMechanicsApplication")
	  ;

      AddCustomConditionsToPython(); 
      AddCustomProcessesToPython();
      AddCustomUtilitiesToPython();
      AddCustomStrategiesToPython();
      AddCustomConstitutiveLawsToPython();
      AddCustomModelersToPython();
      
      //registering variables in python ( if must to be seen from python )
      KRATOS_REGISTER_IN_PYTHON_VARIABLE( NUMBER_OF_ACTIVE_CONTACTS );
      KRATOS_REGISTER_IN_PYTHON_VARIABLE( NUMBER_OF_STICK_CONTACTS );
      KRATOS_REGISTER_IN_PYTHON_VARIABLE( NUMBER_OF_SLIP_CONTACTS );
 
      KRATOS_REGISTER_IN_PYTHON_VARIABLE( MEAN_ERROR );
 
      KRATOS_REGISTER_IN_PYTHON_VARIABLE( RIGID_WALL );
      KRATOS_REGISTER_IN_PYTHON_VARIABLE( WALL_TIP_RADIUS );
      KRATOS_REGISTER_IN_PYTHON_VARIABLE( WALL_REFERENCE_POINT );
      KRATOS_REGISTER_IN_PYTHON_VARIABLE( WALL_VELOCITY );

      KRATOS_REGISTER_IN_PYTHON_VARIABLE( OFFSET );

      KRATOS_REGISTER_IN_PYTHON_VARIABLE( FORCE_CONTACT_NORMAL );
      KRATOS_REGISTER_IN_PYTHON_VARIABLE( FORCE_CONTACT_TANGENT );

      KRATOS_REGISTER_IN_PYTHON_VARIABLE( SHRINK_FACTOR );

    }
  
  
  }  // namespace Python.
  
}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
