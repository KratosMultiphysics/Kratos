// ==============================================================================
//  KratosTopologyOptimizationApplication
//
//  License:         BSD License
//                   license: TopologyOptimizationApplication/license.txt
//
//  Main authors:    Baumgärtner Daniel, https://github.com/dbaumgaertner
//                   Octaviano Malfavón Farías
//                   Eric Gonzales
//
// ==============================================================================

// System includes 

#if defined(KRATOS_PYTHON)
// External includes 
#include <boost/python.hpp>


// Project includes 
#include "includes/define.h"
#include "topology_optimization_application.h"
#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"

 
namespace Kratos
{

namespace Python
{

  using namespace boost::python;


  
  BOOST_PYTHON_MODULE(KratosTopologyOptimizationApplication)
  {

	  class_<KratosTopologyOptimizationApplication, 
                          KratosTopologyOptimizationApplication::Pointer, 
                          bases<KratosApplication>, boost::noncopyable >("KratosTopologyOptimizationApplication")
                         ;

    AddCustomStrategiesToPython();
    AddCustomUtilitiesToPython();

    //Registering variables in python
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( E_MIN )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( E_0 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( PENAL )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( X_PHYS )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( X_PHYS_OLD )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( DCDX )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( DVDX )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( SOLID_VOID )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( LOCAL_STRAIN_ENERGY )


  }
  
  
}  // namespace Python.
  
}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
