//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Inigo Lopez and Riccardo Rossi
//

// System includes

#if defined(KRATOS_PYTHON)
// External includes
#include <pybind11/pybind11.h>

// Project includes
#include "includes/define.h"
#include "compressible_potential_flow_application.h"
#include "compressible_potential_flow_application_variables.h"
#include "custom_python/add_custom_processes_to_python.h"
#include "custom_python/add_custom_response_functions_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"

namespace Kratos
{
namespace Python
{

PYBIND11_MODULE(KratosCompressiblePotentialFlowApplication, m)
{
    namespace py = pybind11;

    py::class_<KratosCompressiblePotentialFlowApplication,
               KratosCompressiblePotentialFlowApplication::Pointer,
               KratosApplication>(m, "KratosCompressiblePotentialFlowApplication")
        .def(py::init<>());

    AddCustomResponseFunctionUtilitiesToPython(m);
    AddCustomProcessesToPython(m);
    AddCustomUtilitiesToPython(m);

    //registering variables in python
    // Degrees of freedom
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, VELOCITY_POTENTIAL);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, AUXILIARY_VELOCITY_POTENTIAL);

    //Embedded variables
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, GEOMETRY_DISTANCE)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, ROTATION_ANGLE);

    //Wake variables
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, WAKE_DISTANCE)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, WAKE_ELEMENTAL_DISTANCES)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, WAKE_ORIGIN)

    // Adjoint potential
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, ADJOINT_VELOCITY_POTENTIAL);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, ADJOINT_AUXILIARY_VELOCITY_POTENTIAL);

    // Flow field magnitudes
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, PERTURBATION_VELOCITY);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, VELOCITY_LOWER);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PRESSURE_LOWER);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, POTENTIAL_JUMP);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, ENERGY_NORM_REFERENCE);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, POTENTIAL_ENERGY_REFERENCE);

    // Free stream magnitudes
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, FREE_STREAM_VELOCITY);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, FREE_STREAM_DENSITY);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, FREE_STREAM_MACH);

    // Geometrical variables
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, REFERENCE_CHORD);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, WAKE_NORMAL);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, WING_SPAN_DIRECTION);

    // Solver parameters
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MACH_LIMIT);

    // Integral magnitudes
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, LIFT_COEFFICIENT);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, DRAG_COEFFICIENT);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MOMENT_COEFFICIENT);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, LIFT_COEFFICIENT_JUMP);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, LIFT_COEFFICIENT_FAR_FIELD);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, DRAG_COEFFICIENT_FAR_FIELD);

    // Markers
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, WAKE);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, KUTTA);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, WING_TIP);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, TRAILING_EDGE);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, UPPER_SURFACE);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, LOWER_SURFACE);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, UPPER_WAKE);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, LOWER_WAKE);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, AIRFOIL);

    // To be removed
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, TRAILING_EDGE_ELEMENT);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, DECOUPLED_TRAILING_EDGE_ELEMENT);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, DEACTIVATED_WAKE);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, ALL_TRAILING_EDGE);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, ZERO_VELOCITY_CONDITION);
}

} // namespace Python.
} // namespace Kratos.

#endif // KRATOS_PYTHON defined
