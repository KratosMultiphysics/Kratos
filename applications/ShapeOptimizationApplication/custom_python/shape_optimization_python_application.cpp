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

// ------------------------------------------------------------------------------
// External includes
// ------------------------------------------------------------------------------
#include <boost/python.hpp>

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "includes/define.h"
#include "shape_optimization_application.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "shape_optimization_application.h"
 
// ==============================================================================

namespace Kratos
{

namespace Python
{

  using namespace boost::python;


  
  BOOST_PYTHON_MODULE(KratosShapeOptimizationApplication)
  {

	  class_<KratosShapeOptimizationApplication, 
			  KratosShapeOptimizationApplication::Pointer, 
			  bases<KratosApplication>, boost::noncopyable >("KratosShapeOptimizationApplication")
			;

	AddCustomUtilitiesToPython();

	//registering variables in python

	// Geometry variables
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(NORMALIZED_SURFACE_NORMAL);

    // Optimization variables
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(OBJECTIVE_SENSITIVITY);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(OBJECTIVE_SURFACE_SENSITIVITY);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(MAPPED_OBJECTIVE_SENSITIVITY);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(CONSTRAINT_SENSITIVITY);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(CONSTRAINT_SURFACE_SENSITIVITY);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(MAPPED_CONSTRAINT_SENSITIVITY);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(SEARCH_DIRECTION);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(CONTROL_POINT_UPDATE);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(CONTROL_POINT_CHANGE);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(SHAPE_UPDATE);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(SHAPE_CHANGE);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(MESH_CHANGE);    

    // For edge damping
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(DAMPING_FACTOR);

    // For Structure Sensitivity Analysis
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(STRAIN_ENERGY_SHAPE_GRADIENT);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(MASS_SHAPE_GRADIENT);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(ACTIVE_NODE_INDEX);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(DKDXU);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(DKDXU_X);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(DKDXU_Y);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(DKDXU_Z);

    // For mapping
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(MAPPING_ID);
  }
  
  
}  // namespace Python.
  
}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
