//
//   Project Name:        KratosPfemApplication     $
//   Created by:          $Author:      JMCarbonell $
//   Last modified by:    $Co-Author:               $
//   Date:                $Date:      February 2016 $
//   Revision:            $Revision:            0.0 $
//
//

// System includes
#if defined(KRATOS_PYTHON)

// External includes 

// Project includes 
#include "custom_python/add_custom_processes_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_python/add_custom_modelers_to_python.h"
#include "custom_python/add_custom_bounding_to_python.h"

#include "pfem_application.h"
 
namespace Kratos
{

namespace Python
{

using namespace pybind11;
  
PYBIND11_MODULE(KratosPfemApplication,m)
{

  class_<KratosPfemApplication, 
         KratosPfemApplication::Pointer, 
         KratosApplication>(m,"KratosPfemApplication")
      .def(init<>())
      ;

  AddCustomProcessesToPython(m);
  AddCustomUtilitiesToPython(m);
  AddCustomModelersToPython(m);
  AddCustomBoundingToPython(m);
      
  //registering variables in python ( if must to be seen from python )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, INITIALIZED_DOMAINS )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MESHING_STEP_PERFORMED )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, RIGID_WALL )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MEAN_ERROR )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, OFFSET )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SHRINK_FACTOR )
  }
  
  
}  // namespace Python.
  
}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
