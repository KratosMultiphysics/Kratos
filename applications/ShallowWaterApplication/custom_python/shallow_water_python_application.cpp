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
#include <boost/python.hpp>

// Project includes 
#include "includes/define.h"
#include "shallow_water_application.h"
#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"

 
namespace Kratos
{

namespace Python
{

  using namespace boost::python;

  
  BOOST_PYTHON_MODULE(KratosShallowWaterApplication)
  {
    class_<KratosShallowWaterApplication, 
        KratosShallowWaterApplication::Pointer, 
        bases<KratosApplication>, boost::noncopyable >("KratosShallowWaterApplication")
        ;

    AddCustomStrategiesToPython();
    AddCustomUtilitiesToPython();

    // Registering variables in python
    // Shallow water variables
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(HEIGHT);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(BATHYMETRY);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(RAIN);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(FREE_SURFACE_ELEVATION);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(MANNING);

    // Specific variables for PFEM2
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(MEAN_SIZE);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(DELTA_SCALAR1)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(PROJECTED_SCALAR1)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(DELTA_VECTOR1)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(PROJECTED_VECTOR1)

    // Units conversion
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(TIME_UNIT_CONVERTER)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(WATER_HEIGHT_UNIT_CONVERTER)
  }
  
}  // namespace Python.
  
}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
