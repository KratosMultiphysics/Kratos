//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                    July 2013 $
//   Revision:            $Revision:                      0.0 $
//
//

// System includes 
#if defined(KRATOS_PYTHON)
// External includes 
#include <boost/python.hpp>

// Project includes 
#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_python/add_custom_processes_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_python/add_custom_constitutive_laws_to_python.h"

#include "pfem_solid_mechanics_application.h"
 
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
	     //bases<KratosSolidMechanicsApplication>, boost::noncopyable >("KratosPfemSolidMechanicsApplication")
	  ;

      AddCustomProcessesToPython();
      AddCustomUtilitiesToPython();
      AddCustomStrategiesToPython();
      AddCustomConstitutiveLawsToPython();

      KRATOS_REGISTER_IN_PYTHON_VARIABLE( JACOBIAN )
      KRATOS_REGISTER_IN_PYTHON_VARIABLE( REACTION_JACOBIAN )
     
      //registering variables in python ( if must to be seen from python )
      KRATOS_REGISTER_IN_PYTHON_VARIABLE( TOTAL_CAUCHY_STRESS )
      KRATOS_REGISTER_IN_PYTHON_VARIABLE( WATER_PRESSURE )
      KRATOS_REGISTER_IN_PYTHON_VARIABLE( REACTION_WATER_PRESSURE )
      KRATOS_REGISTER_IN_PYTHON_VARIABLE( IMPOSED_WATER_PRESSURE )
      KRATOS_REGISTER_IN_PYTHON_VARIABLE( DARCY_FLOW )
  
      KRATOS_REGISTER_IN_PYTHON_VARIABLE( WALL_TIP_RADIUS )
      KRATOS_REGISTER_IN_PYTHON_VARIABLE( WALL_REFERENCE_POINT )
      KRATOS_REGISTER_IN_PYTHON_VARIABLE( WALL_VELOCITY )

      KRATOS_REGISTER_IN_PYTHON_VARIABLE( PRECONSOLIDATION )
      KRATOS_REGISTER_IN_PYTHON_VARIABLE( STRESS_INV_P )
      KRATOS_REGISTER_IN_PYTHON_VARIABLE( STRESS_INV_J2 )
      KRATOS_REGISTER_IN_PYTHON_VARIABLE( STRESS_INV_THETA )
      KRATOS_REGISTER_IN_PYTHON_VARIABLE( VOLUMETRIC_PLASTIC )
      KRATOS_REGISTER_IN_PYTHON_VARIABLE( INCR_SHEAR_PLASTIC )
      KRATOS_REGISTER_IN_PYTHON_VARIABLE( M_MODULUS )
      KRATOS_REGISTER_IN_PYTHON_VARIABLE( POROSITY )
      KRATOS_REGISTER_IN_PYTHON_VARIABLE( INITIAL_POROSITY )
      KRATOS_REGISTER_IN_PYTHON_VARIABLE( COHESION )
      KRATOS_REGISTER_IN_PYTHON_VARIABLE( INTERNAL_FRICTION_ANGLE )
      KRATOS_REGISTER_IN_PYTHON_VARIABLE( INTERNAL_DILATANCY_ANGLE )
      KRATOS_REGISTER_IN_PYTHON_VARIABLE( DENSITY_WATER )
      KRATOS_REGISTER_IN_PYTHON_VARIABLE( K0 )
      KRATOS_REGISTER_IN_PYTHON_VARIABLE( PERMEABILITY )

    }
  
  
  }  // namespace Python.
  
}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
