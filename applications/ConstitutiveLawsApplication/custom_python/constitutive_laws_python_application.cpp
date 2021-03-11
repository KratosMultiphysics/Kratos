//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Alejandro Cornejo Velazquez
//                   Riccardo Rossi
//


// System includes

#if defined(KRATOS_PYTHON)
// External includes
#include <pybind11/pybind11.h>


// Project includes
#include "includes/define_python.h"
#include "constitutive_laws_application.h"
#include "constitutive_laws_application_variables.h"

#include "custom_python/add_custom_constitutive_laws_to_python.h"


namespace Kratos {
namespace Python {

PYBIND11_MODULE(KratosConstitutiveLawsApplication,m)
{
    namespace py = pybind11;

    py::class_<KratosConstitutiveLawsApplication,
        KratosConstitutiveLawsApplication::Pointer,
        KratosApplication>(m, "KratosConstitutiveLawsApplication")
        .def(py::init<>())
        ;

    AddCustomConstitutiveLawsToPython(m);

    // Constitutive laws variables
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, HIGH_CYCLE_FATIGUE_COEFFICIENTS)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, FATIGUE_REDUCTION_FACTOR)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, NUMBER_OF_CYCLES)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, LOCAL_NUMBER_OF_CYCLES)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, WOHLER_STRESS)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, REVERSION_FACTOR_RELATIVE_ERROR)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MAX_STRESS_RELATIVE_ERROR)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MAX_STRESS)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, THRESHOLD_STRESS)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, CYCLE_INDICATOR)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, CYCLES_TO_FAILURE)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, TIME_INCREMENT)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, DAMAGE_ACTIVATION)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PREVIOUS_CYCLE);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, CYCLE_PERIOD)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, ADVANCE_STRATEGY_APPLIED);
}

} // namespace Python.
} // namespace Kratos.

#endif // KRATOS_PYTHON defined
