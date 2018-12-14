//
//   Project Name:        KratosContactMechanicsApplication $
//   Created by:          $Author:              JMCarbonell $
//   Last modified by:    $Co-Author:                       $
//   Date:                $Date:                  July 2016 $
//   Revision:            $Revision:                    0.0 $
//
//

// System includes
#if defined(KRATOS_PYTHON)

// External includes

// Project includes
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_python/add_custom_processes_to_python.h"
#include "custom_python/add_custom_bounding_to_python.h"
#include "custom_python/add_custom_meshers_to_python.h"
#include "custom_python/add_custom_friction_laws_to_python.h"

#include "contact_mechanics_application.h"

namespace Kratos
{

namespace Python
{

namespace py = pybind11;


PYBIND11_MODULE(KratosContactMechanicsApplication,m)
{

  py::class_<KratosContactMechanicsApplication,
         KratosContactMechanicsApplication::Pointer,
         KratosApplication>(m,"KratosContactMechanicsApplication")
      .def(py::init<>())
      ;

  AddCustomUtilitiesToPython(m);
  AddCustomProcessesToPython(m);
  AddCustomBoundingToPython(m);
  AddCustomMeshersToPython(m);
  AddCustomFrictionLawsToPython(m);

  //registering variables in python
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, FRICTION_ACTIVE )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PENALTY_PARAMETER )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, LAGRANGE_MULTIPLIER_NORMAL )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, LAGRANGE_MULTIPLIER_NORMAL_REACTION )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, LAGRANGE_MULTIPLIER_TANGENTIAL )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, LAGRANGE_MULTIPLIER_TANGENTIAL_REACTION )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, TAU_STAB )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MU_STATIC )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MU_DYNAMIC )
  KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, CONTACT_STRESS )
  KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, EFFECTIVE_CONTACT_STRESS )
  KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, EFFECTIVE_CONTACT_FORCE )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, CONTACT_ADHESION )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, CONTACT_FRICTION_ANGLE )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, TANGENTIAL_PENALTY_RATIO )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, CONTACT_PLASTIC_SLIP )

  //KRATOS_REGISTER_IN_PYTHON_VARIABLE(NODAL_AREA)

  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, FRICTION_LAW )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, FRICTION_LAW_NAME )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, CONTACT_STEP_TIME )
  }


}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
