//
//   Project Name:
//   Last modified by:    $Author:  $
//   Date:                $Date: $
//   Revision:            $Revision: $
//

#if defined(KRATOS_PYTHON)

// External includes
#include <boost/python.hpp>

// Project includes
#include "includes/define.h"

#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_python/add_custom_constitutive_laws_to_python.h"
#include "custom_python/add_custom_processes_to_python.h"

#include "dam_application.h"

namespace Kratos
{

namespace Python
{

using namespace boost::python;

BOOST_PYTHON_MODULE(KratosDamApplication)
{
    class_<KratosDamApplication, KratosDamApplication::Pointer, bases<KratosApplication>, boost::noncopyable >("KratosDamApplication");

    AddCustomStrategiesToPython();
    AddCustomConstitutiveLawsToPython();
    AddCustomProcessesToPython();

    //Registering variables in python
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( NEWMARK_BETA )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( NEWMARK_GAMMA )

    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( IMPOSED_POINT_LOAD )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( IMPOSED_LINE_LOAD )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( IMPOSED_SURFACE_LOAD )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( IMPOSED_NORMAL_STRESS )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( IMPOSED_TANGENTIAL_STRESS )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( IMPOSED_TEMPERATURE )
    
    //Bofang and Hidrostatic variables for evolution changes
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( GRAVITY_DIRECTION )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( COORDINATE_BASE_DAM )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( SURFACE_TEMP )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( BOTTOM_TEMP )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( HEIGHT_DAM )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( AMPLITUDE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( FREQUENCY )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( DAY_MAXIMUM )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( SPECIFIC_WEIGHT )   
    
    // Thermal Variables
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( REFERENCE_TEMPERATURE );
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( THERMAL_STRESS_TENSOR )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( MECHANICAL_STRESS_TENSOR )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( THERMAL_STRAIN_TENSOR )
    
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( THERMAL_STRESS_VECTOR )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( MECHANICAL_STRESS_VECTOR )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( THERMAL_STRAIN_VECTOR )
}

}  // namespace Python.
}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
