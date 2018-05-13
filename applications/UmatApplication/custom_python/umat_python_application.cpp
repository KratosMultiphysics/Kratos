//
//   Project Name:        KratosUmatApplication      $
//   Created by:          $Author:      JMCarbonell  $
//   Last modified by:    $Co-Author:                $
//   Date:                $Date:      September 2017 $
//   Revision:            $Revision:             0.0 $
//
//

// System includes
#if defined(KRATOS_PYTHON)

// External includes

// Project includes
#include "custom_python/add_custom_constitutive_laws_to_python.h"

#include "umat_application.h"

namespace Kratos
{

namespace Python
{

using namespace pybind11;



PYBIND11_MODULE(KratosUmatApplication,m)
{

  class_<KratosUmatApplication,
         KratosUmatApplication::Pointer,
         KratosApplication>(m,"KratosUmatApplication")
      .def(init<>())
      ;

  AddCustomConstitutiveLawsToPython(m);

  //registering variables in python

}


}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
