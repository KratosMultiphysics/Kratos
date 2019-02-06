//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

// System includes
#if defined(KRATOS_PYTHON)

// External includes

// Project includes
#include "custom_python/add_custom_processes_to_python.h"
#include "custom_python/add_custom_constitutive_laws_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"

#include "constitutive_models_application.h"


namespace Kratos
{

namespace Python
{

namespace py = pybind11;


PYBIND11_MODULE(KratosConstitutiveModelsApplication,m)
{

  py::class_<KratosConstitutiveModelsApplication,
         KratosConstitutiveModelsApplication::Pointer,
         KratosApplication>(m,"KratosConstitutiveModelsApplication")
      .def(py::init<>())
      ;

  AddCustomConstitutiveLawsToPython(m);
  AddCustomUtilitiesToPython(m);
  AddCustomProcessesToPython(m);

      //registering variables in python


  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PROPERTIES_LAYOUT )
      
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, RHOS )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, RHOT )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, KSIM )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, CHIS )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, CHIT )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, REFERENCE_PRESSURE )

  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PS )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PT )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PM )

  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PLASTIC_VOL_DEF )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, NONLOCAL_PLASTIC_VOL_DEF )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PLASTIC_VOL_DEF_ABS )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, NONLOCAL_PLASTIC_VOL_DEF_ABS )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PLASTIC_DEV_DEF )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, NONLOCAL_PLASTIC_DEV_DEF )

}


}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
