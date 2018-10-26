//
//   Project Name:        KratosDelaunayMeshingApplication $
//   Created by:          $Author:             JMCarbonell $
//   Last modified by:    $Co-Author:                      $
//   Date:                $Date:                April 2018 $
//   Revision:            $Revision:                   0.0 $
//
//

// System includes
#if defined(KRATOS_PYTHON)

// External includes

// Project includes
#include "custom_python/add_custom_processes_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_python/add_custom_meshers_to_python.h"
#include "custom_python/add_custom_bounding_to_python.h"

#include "delaunay_meshing_application.h"

namespace Kratos
{

namespace Python
{

using namespace pybind11;

PYBIND11_MODULE(KratosDelaunayMeshingApplication,m)
{

  class_<KratosDelaunayMeshingApplication,
         KratosDelaunayMeshingApplication::Pointer,
         KratosApplication>(m,"KratosDelaunayMeshingApplication")
      .def(init<>())
      ;

  AddCustomProcessesToPython(m);
  AddCustomUtilitiesToPython(m);
  AddCustomMeshersToPython(m);
  AddCustomBoundingToPython(m);

  //registering variables in python ( if must to be seen from python )

  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, INITIALIZED_DOMAINS )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MESHING_STEP_TIME )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, RIGID_WALL )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MEAN_ERROR )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, OFFSET )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SHRINK_FACTOR )

}

}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
