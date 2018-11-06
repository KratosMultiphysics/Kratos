//
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:            January 2016 $
//   Revision:            $Revision:                 1.0 $
//

#if defined(KRATOS_PYTHON)

// System includes
#include <pybind11/pybind11.h>

// Project includes
#include "includes/define.h"
#include "includes/define_python.h"

// Application includes
#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_python/add_custom_constitutive_laws_to_python.h"
#include "custom_python/add_custom_processes_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "poromechanics_application.h"

namespace Kratos
{

namespace Python
{

namespace py = pybind11;

PYBIND11_MODULE(KratosPoromechanicsApplication, m)
{
    py::class_<KratosPoromechanicsApplication,
    KratosPoromechanicsApplication::Pointer,
    KratosApplication>(m, "KratosPoromechanicsApplication")
    .def( py::init<>());

    AddCustomStrategiesToPython(m);
    AddCustomConstitutiveLawsToPython(m);
    AddCustomProcessesToPython(m);
    AddCustomUtilitiesToPython(m);

    //Registering variables in python
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, DT_WATER_PRESSURE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, NORMAL_FLUID_FLUX )

    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( m, FLUID_FLUX_VECTOR )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( m, LOCAL_FLUID_FLUX_VECTOR )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( m, LOCAL_STRESS_VECTOR )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( m, LOCAL_RELATIVE_DISPLACEMENT_VECTOR )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, PERMEABILITY_MATRIX )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, LOCAL_PERMEABILITY_MATRIX )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, TOTAL_STRESS_TENSOR )

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
}

}  // namespace Python.
}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
