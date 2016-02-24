//
//   Project Name:        KratosPfemFluidDynamicsApplication $
//   Created by:          $Author:               JMCarbonell $
//   Last modified by:    $Co-Author:                        $
//   Date:                $Date:               February 2016 $
//   Revision:            $Revision:                     0.0 $
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
#include "custom_python/add_custom_modelers_to_python.h"

#include "pfem_fluid_dynamics_application.h"
 
namespace Kratos
{

  namespace Python
  {

    using namespace boost::python;


  
    BOOST_PYTHON_MODULE(KratosPfemFluidDynamicsApplication)
    {

      class_<KratosPfemFluidDynamicsApplication, 
	     KratosPfemFluidDynamicsApplication::Pointer, 
	     bases<KratosApplication>, boost::noncopyable >("KratosPfemFluidDynamicsApplication")
	     //bases<KratosFluidDynamicsApplication>, boost::noncopyable >("KratosPfemFluidDynamicsApplication")
	  ;

      AddCustomProcessesToPython();
      AddCustomUtilitiesToPython();
      AddCustomStrategiesToPython();
      AddCustomConstitutiveLawsToPython();
      AddCustomModelersToPython();
      
      //registering variables in python ( if must to be seen from python )

    }
  
  
  }  // namespace Python.
  
}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
