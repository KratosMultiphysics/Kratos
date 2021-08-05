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



// Project includes 
#include "includes/define_python.h"
#include "topology_optimization_application.h"
#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_python/add_custom_constitutive_laws_to_python.h"
#include "custom_python/add_custom_response_functions_to_python.h"


namespace Kratos
{

namespace Python
{

  


  
  PYBIND11_MODULE(KratosTopologyOptimizationApplication,m)
  {
    namespace py = pybind11;

	  py::class_<KratosTopologyOptimizationApplication, 
                          KratosTopologyOptimizationApplication::Pointer, 
                          KratosApplication >(m, "KratosTopologyOptimizationApplication")
                          .def(py::init<>())
                         ;

    AddCustomStrategiesToPython(m);
    AddCustomUtilitiesToPython(m);

    //Registering variables in python
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, E_MIN )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, E_0 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PENAL )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, X_PHYS )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, X_PHYS_OLD )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, X_PHYS_OLD_1 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, X_PHYS_OLD_2 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, DCDX )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, DCDX_OLD )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, DCDX_OLD_2 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, DVDX )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SOLID_VOID )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, LOCAL_STRAIN_ENERGY )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, LOW )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, UPP )


  }
  
  
}  // namespace Python.
  
}  // namespace Kratos.

#endif // KRATOS_PYTHON defined