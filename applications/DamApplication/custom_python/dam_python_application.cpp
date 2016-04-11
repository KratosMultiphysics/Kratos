//
//   Project Name:
//   Last modified by:    $Author:  $
//   Date:                $Date: $
//   Revision:            $Revision: $
//

#if defined(KRATOS_PYTHON)

// External includes
#include <boost/python.hpp>
#include "boost/python/detail/wrap_python.hpp"

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
    
    //Bofang and Hidrostatic variables for evolution changes
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( GRAVITY_DIRECTION )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( COORDINATE_BASE_DAM )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( SURFACE_TEMP )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( BOTTOM_TEMP )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( HEIGHT_DAM )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( AMPLITUDE )
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
