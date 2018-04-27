// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
//
// ==============================================================================

#if defined(KRATOS_PYTHON)

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <pybind11/pybind11.h>

// ------------------------------------------------------------------------------
// External includes
// ------------------------------------------------------------------------------

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "includes/define_python.h"
#include "shape_optimization_application.h"
#include "custom_python/add_custom_utilities_to_python.h"

// ==============================================================================

namespace Kratos
{

namespace Python
{

using namespace pybind11;



  PYBIND11_MODULE(KratosShapeOptimizationApplication, m)
  {

	  class_<KratosShapeOptimizationApplication,
			  KratosShapeOptimizationApplication::Pointer,
			  KratosApplication >(m, "KratosShapeOptimizationApplication")
        .def(init<>())
			;

	AddCustomUtilitiesToPython(m);


	//registering variables in python

	// Geometry variables
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, NORMALIZED_SURFACE_NORMAL);

    // Optimization variables
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, OBJECTIVE_SENSITIVITY);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, OBJECTIVE_SURFACE_SENSITIVITY);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, MAPPED_OBJECTIVE_SENSITIVITY);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, CONSTRAINT_SENSITIVITY);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, CONSTRAINT_SURFACE_SENSITIVITY);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, MAPPED_CONSTRAINT_SENSITIVITY);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, SEARCH_DIRECTION);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, CONTROL_POINT_UPDATE);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, CONTROL_POINT_CHANGE);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, SHAPE_UPDATE);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, SHAPE_CHANGE);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, MESH_CHANGE);

    // For edge damping
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, DAMPING_FACTOR);

    // For mapping
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MAPPING_ID);
  }


}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
