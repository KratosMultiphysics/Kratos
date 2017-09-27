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
#include "custom_python/add_custom_processes_to_python.h"
#include "custom_python/add_custom_schemes_to_python.h"

 
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
  AddCustomSchemesToPython();
  AddCustomProcessesToPython();

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
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(DESIGN_UPDATE);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(DESIGN_CHANGE_ABSOLUTE);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(SHAPE_UPDATE);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(SHAPE_CHANGE_ABSOLUTE);

    // For edge damping
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(DAMPING_FACTOR);

    // For Structure Sensitivity Analysis
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(STRAIN_ENERGY_SHAPE_GRADIENT);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(MASS_SHAPE_GRADIENT);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(EIGENFREQUENCY_SHAPE_GRADIENT);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(LOCAL_STRESS_GRADIENT);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(ACTIVE_NODE_INDEX);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(DKDXU);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(DKDXU_X);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(DKDXU_Y);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(DKDXU_Z);
    // Variables for output of sensitivities
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(CROSS_AREA_SENSITIVITY); 
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( POINT_LOAD_SENSITIVITY )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( IY_SENSITIVITY )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( IZ_SENSITIVITY )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( THICKNESS_SENSITIVITY )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( YOUNG_MODULUS_SENSITIVITY)
    // Variables to for computing parts of sensitivity analysis
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(TRACED_STRESS_TYPE); 
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(STRESS_DISP_DERIV_ON_NODE);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(STRESS_DISP_DERIV_ON_GP); 
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(STRESS_DV_DERIV_ON_NODE); 
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(STRESS_DV_DERIV_ON_GP); 
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(STRESS_ON_GP); 
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(STRESS_ON_NODE); 
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(DESIGN_VARIABLE_NAME); 

    // For mapping
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(MAPPING_ID);
  }
  
  
}  // namespace Python.
  
}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
