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
#include "custom_python/add_custom_utilities_to_python.h"


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

    // Registering variables in python
    // Shallow water variables
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,HEIGHT);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,BATHYMETRY);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,RAIN);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,FREE_SURFACE_ELEVATION);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,MANNING);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,TOPOGRAPHY_SLOPE);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,TOPOGRAPHY_GRADIENT);

    // Specific variables for PFEM2
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,MEAN_SIZE);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,DELTA_SCALAR1)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,PROJECTED_SCALAR1)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m,DELTA_VECTOR1)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m,PROJECTED_VECTOR1)

    // Units conversion
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,TIME_UNIT_CONVERTER)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,WATER_HEIGHT_UNIT_CONVERTER)
  }

}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
