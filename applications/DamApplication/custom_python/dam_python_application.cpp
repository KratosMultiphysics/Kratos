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
#include "custom_python/add_custom_utilities_to_python.h"
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
    AddCustomUtilitiesToPython();

    //Registering variables in python
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( THERMAL_EXPANSION )    
 
    // Thermal Variables
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( THERMAL_STRESS_TENSOR )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( MECHANICAL_STRESS_TENSOR )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( THERMAL_STRAIN_TENSOR )
    
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( THERMAL_STRESS_VECTOR )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( MECHANICAL_STRESS_VECTOR )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( THERMAL_STRAIN_VECTOR )

    KRATOS_REGISTER_IN_PYTHON_VARIABLE( ALPHA_HEAT_SOURCE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( TIME_ACTIVATION )
    
    // Output Variables
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( Vi_POSITIVE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( Viii_POSITIVE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( NODAL_JOINT_WIDTH )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( NODAL_JOINT_AREA )
    
    // Wave Eqaution
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( Dt_PRESSURE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( Dt2_PRESSURE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( VELOCITY_PRESSURE_COEFFICIENT )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( ACCELERATION_PRESSURE_COEFFICIENT )    

    // Others
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( NODAL_YOUNG_MODULUS )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( ADDED_MASS )    
    
}

}  // namespace Python.
}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
