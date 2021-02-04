//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//

// System includes

#if defined(KRATOS_PYTHON)
// External includes
#include <pybind11/pybind11.h>


// Project includes
#include "includes/define_python.h"
#include "shallow_water_application.h"
#include "shallow_water_application_variables.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_python/add_custom_processes_to_python.h"
#include "custom_python/add_custom_strategies_to_python.h"


namespace Kratos
{

namespace Python
{

  namespace py = pybind11;

  PYBIND11_MODULE(KratosShallowWaterApplication, m)
  {
    py::class_<KratosShallowWaterApplication,
        KratosShallowWaterApplication::Pointer,
        KratosApplication>(m, "KratosShallowWaterApplication")
        .def(py::init<>())
        ;

    AddCustomUtilitiesToPython(m);
    AddCustomProcessesToPython(m);
    AddCustomStrategiesToPython(m);

    // Adding enums
    py::enum_<Framework>(m, "Framework")
    .value("EULERIAN_FRAMEWORK", EULERIAN_FRAMEWORK)
    .value("PFEM2_FRAMEWORK", PFEM2_FRAMEWORK)
    ;

    py::enum_<Formulation>(m, "Formulation")
    .value("REDUCED_VARIABLES", REDUCED_VARIABLES)
    .value("CONSERVED_VARIABLES", CONSERVED_VARIABLES)
    ;

    py::enum_<Variables>(m, "Variables")
    .value("FREE_SURFACE_VARIABLE", FREE_SURFACE_VARIABLE)
    .value("VELOCITY_VARIABLE", VELOCITY_VARIABLE)
    .value("FREE_SURFACE_AND_VELOCITY", FREE_SURFACE_AND_VELOCITY)
    ;

    // Registering variables in python
    // Primary variables
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, HEIGHT);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, FREE_SURFACE_ELEVATION);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, VERTICAL_VELOCITY);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, FLOW_RATE);

    // Physical variables
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, BATHYMETRY);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, TOPOGRAPHY);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, RAIN);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MANNING);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PERMEABILITY);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, ATMOSPHERIC_PRESSURE);

    // Auxiliary variables
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SHOCK_STABILIZATION_FACTOR);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, DRY_HEIGHT);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, RELATIVE_DRY_HEIGHT);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, DRY_DISCHARGE_PENALTY);

    // Post-process variables
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, TOPOGRAPHY_GRADIENT);

    // Specific variables for PFEM2
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,DELTA_SCALAR1)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,PROJECTED_SCALAR1)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m,DELTA_VECTOR1)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m,PROJECTED_VECTOR1)

    // Variables for Algebraic Flux Corrected Transport algorithm
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, POSITIVE_FLUX)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, NEGATIVE_FLUX)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, POSITIVE_RATIO)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, NEGATIVE_RATIO)

    // Benchmark variables
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, EXACT_HEIGHT)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, HEIGHT_ERROR)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, EXACT_FREE_SURFACE)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, FREE_SURFACE_ERROR)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, EXACT_VELOCITY)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, VELOCITY_ERROR)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, EXACT_MOMENTUM)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, MOMENTUM_ERROR)
  }

}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
