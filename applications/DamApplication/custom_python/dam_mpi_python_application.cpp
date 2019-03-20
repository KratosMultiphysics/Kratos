//
//   Project Name:
//   Last modified by:    $Author:  $
//   Date:                $Date: $
//   Revision:            $Revision: $
//

#if defined(KRATOS_PYTHON)

// System includes
#include <pybind11/pybind11.h>

// Project includes
#include "includes/define.h"
#include "includes/define_python.h"

#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_python/add_custom_mpi_strategies_to_python.h"
#include "custom_python/add_custom_constitutive_laws_to_python.h"
#include "custom_python/add_custom_processes_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "dam_application.h"



namespace Kratos
{

namespace Python
{

namespace py = pybind11;

PYBIND11_MODULE(KratosDamApplication, m)
{
    py::class_<KratosDamApplication,
    KratosDamApplication::Pointer,
    KratosApplication>(m, "KratosDamApplication")
    .def(py::init<>());

    AddCustomStrategiesToPython(m);
    AddCustomMPIStrategiesToPython(m);
    AddCustomConstitutiveLawsToPython(m);
    AddCustomProcessesToPython(m);
    AddCustomUtilitiesToPython(m);

    //Registering variables in python
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, THERMAL_EXPANSION )

    // Thermal Variables
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, THERMAL_STRESS_TENSOR )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, MECHANICAL_STRESS_TENSOR )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, THERMAL_STRAIN_TENSOR )

    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, THERMAL_STRESS_VECTOR )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, MECHANICAL_STRESS_VECTOR )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, THERMAL_STRAIN_VECTOR )

    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, ALPHA_HEAT_SOURCE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, TIME_ACTIVATION )

    // Output Variables
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, Vi_POSITIVE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, Viii_POSITIVE )

    // Wave Eqaution
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, Dt_PRESSURE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, Dt2_PRESSURE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, VELOCITY_PRESSURE_COEFFICIENT )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, ACCELERATION_PRESSURE_COEFFICIENT )

    // Others
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, NODAL_YOUNG_MODULUS )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, ADDED_MASS )


}

}  // namespace Python.
}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
