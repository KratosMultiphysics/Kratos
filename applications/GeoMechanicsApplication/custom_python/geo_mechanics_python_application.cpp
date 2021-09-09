// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Ignasi de Pouplana,
//                   Vahid Galavi
//

// System includes

#if defined(KRATOS_PYTHON)

// External includes
#include <pybind11/pybind11.h>


// Project includes
#include "includes/define_python.h"


// Application includes
#include "geo_mechanics_application.h"
#include "geo_mechanics_application_variables.h"
#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_python/add_custom_constitutive_laws_to_python.h"
#include "custom_python/add_custom_processes_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"


namespace Kratos {
namespace Python {

  using namespace pybind11;

  PYBIND11_MODULE(KratosGeoMechanicsApplication,m)
  {

  class_<KratosGeoMechanicsApplication,
         KratosGeoMechanicsApplication::Pointer,
         KratosApplication>(m, "KratosGeoMechanicsApplication")
        .def(init<>());



    AddCustomStrategiesToPython(m);
    AddCustomUtilitiesToPython(m);
    AddCustomConstitutiveLawsToPython(m);
    AddCustomProcessesToPython(m);

    //Registering variables in python
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, VELOCITY_COEFFICIENT )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, DT_PRESSURE_COEFFICIENT )

    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, DT_WATER_PRESSURE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, NORMAL_FLUID_FLUX )

    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, HYDRAULIC_HEAD )

    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, HYDRAULIC_DISCHARGE )

    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( m, FLUID_FLUX_VECTOR )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( m, LOCAL_FLUID_FLUX_VECTOR )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( m, LOCAL_STRESS_VECTOR )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( m, LOCAL_RELATIVE_DISPLACEMENT_VECTOR )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, PERMEABILITY_MATRIX )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, LOCAL_PERMEABILITY_MATRIX )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, TOTAL_STRESS_TENSOR )

    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, DEGREE_OF_SATURATION )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, EFFECTIVE_SATURATION )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, BISHOP_COEFICIENT )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, DERIVATIVE_OF_SATURATION )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, RELATIVE_PERMEABILITY )

    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( m, TOTAL_DISPLACEMENT)

    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, IS_CONVERGED )

    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, ARC_LENGTH_LAMBDA )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, ARC_LENGTH_RADIUS_FACTOR )

    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, TIME_UNIT_CONVERTER )

    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, JOINT_WIDTH )

    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, NODAL_SMOOTHING )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, NODAL_CAUCHY_STRESS_TENSOR )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, NODAL_DAMAGE_VARIABLE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, NODAL_JOINT_AREA )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, NODAL_JOINT_WIDTH )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, NODAL_JOINT_DAMAGE )

    /* Reset displacement "flag" needed for GeoMechanicalApplication*/
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,RESET_DISPLACEMENTS)

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,CONSIDER_GAP_CLOSURE)

    // KRATOS_REGISTER_IN_PYTHON_VARIABLE(NODAL_AREA);


  }


} // namespace Python.
} // namespace Kratos.

#endif // KRATOS_PYTHON defined
