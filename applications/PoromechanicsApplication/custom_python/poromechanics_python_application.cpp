//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:            January 2016 $
//   Revision:            $Revision:                 1.0 $
//

#if defined(KRATOS_PYTHON)

// External includes 
#include <boost/python.hpp>

// Project includes 
#include "includes/define.h"

#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_python/add_custom_constitutive_laws_to_python.h"

#include "poromechanics_application.h"
 
namespace Kratos
{

namespace Python
{

using namespace boost::python;
  
BOOST_PYTHON_MODULE(KratosPoromechanicsApplication)
{
    class_<KratosPoromechanicsApplication, KratosPoromechanicsApplication::Pointer, bases<KratosApplication>, boost::noncopyable >("KratosPoromechanicsApplication");

    AddCustomStrategiesToPython();
    AddCustomConstitutiveLawsToPython();

    //Registering variables in python
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( BETA_NEWMARK )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( GAMMA_NEWMARK )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( THETA_NEWMARK )
    
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( CONSTITUTIVE_LAW_POINTER );
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( CONSTITUTIVE_LAW_NAME );

    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( IMPOSED_DISPLACEMENT )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( POINT_LOAD );
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( IMPOSED_POINT_LOAD )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( LINE_LOAD );
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( IMPOSED_LINE_LOAD )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( SURFACE_LOAD );
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( IMPOSED_SURFACE_LOAD )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( IMPOSED_NORMAL_STRESS )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( IMPOSED_TANGENTIAL_STRESS )

    KRATOS_REGISTER_IN_PYTHON_VARIABLE( DERIVATIVE_WATER_PRESSURE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( IMPOSED_FLUID_PRESSURE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( NORMAL_FLUID_FLUX )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( IMPOSED_NORMAL_FLUID_FLUX )

    KRATOS_REGISTER_IN_PYTHON_VARIABLE( VON_MISES_STRESS )
    
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( FLUID_FLUX )
}

}  // namespace Python.
}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
