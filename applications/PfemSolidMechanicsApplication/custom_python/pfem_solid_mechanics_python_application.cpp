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

// Project includes
#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_python/add_custom_processes_to_python.h"
#include "custom_python/add_custom_constitutive_laws_to_python.h"

#include "pfem_solid_mechanics_application.h"

namespace Kratos
{

namespace Python
{

using namespace pybind11;



PYBIND11_MODULE(KratosPfemSolidMechanicsApplication,m)
{

  class_<KratosPfemSolidMechanicsApplication,
         KratosPfemSolidMechanicsApplication::Pointer,
         KratosApplication>(m,"KratosPfemSolidMechanicsApplication")
      .def(init<>())
      ;

  AddCustomProcessesToPython(m);
  AddCustomStrategiesToPython(m);
  AddCustomConstitutiveLawsToPython(m);

  KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, WATER_DISPLACEMENT )
  KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, WATER_VELOCITY )
  KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, WATER_ACCELERATION )
  KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, WATER_DISPLACEMENT_REACTION )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, WATER_PRESSURE_VELOCITY )

  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, JACOBIAN )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, REACTION_JACOBIAN )

  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, TOTAL_CAUCHY_STRESS )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, WATER_PRESSURE )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, REACTION_WATER_PRESSURE )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, DARCY_FLOW )

  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, WALL_TIP_RADIUS )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, WALL_REFERENCE_POINT )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, WALL_VELOCITY )

  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PRECONSOLIDATION )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, STRESS_INV_P )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, STRESS_INV_J2 )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, STRESS_INV_THETA )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, VOLUMETRIC_PLASTIC )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, INCR_SHEAR_PLASTIC )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, M_MODULUS )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, POROSITY )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, INITIAL_POROSITY )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, INTERNAL_FRICTION_ANGLE )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, DENSITY_WATER )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, K0 )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PERMEABILITY )

  }


}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
