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


#include "pfem_application.h"

namespace Kratos
{

namespace Python
{

namespace py = pybind11;;

PYBIND11_MODULE(KratosPfemApplication,m)
{

  py::class_<KratosPfemApplication,
         KratosPfemApplication::Pointer,
         KratosApplication>(m,"KratosPfemApplication")
      .def(py::init<>())
      ;

  AddCustomProcessesToPython(m);

  //registering variables in python ( if must to be seen from python )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PROPERTIES_VECTOR )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MATERIAL_PERCENTAGE )
}


}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
