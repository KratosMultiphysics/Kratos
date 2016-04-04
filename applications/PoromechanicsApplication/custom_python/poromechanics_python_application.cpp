//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:            January 2016 $
//   Revision:            $Revision:                 1.0 $
//

#if defined(KRATOS_PYTHON)

// External includes 
#include <boost/python.hpp>
#include "boost/python/detail/wrap_python.hpp"

// Project includes 
#include "includes/define.h"

// Application includes
#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_python/add_custom_constitutive_laws_to_python.h"
#include "custom_python/add_custom_processes_to_python.h"
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
    AddCustomProcessesToPython();

    //Registering variables in python
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( DT_WATER_PRESSURE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( NORMAL_FLUID_FLUX )

    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( FLUID_FLUX_VECTOR )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( LOCAL_FLUID_FLUX_VECTOR )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( LOCAL_STRESS_VECTOR )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( LOCAL_RELATIVE_DISPLACEMENT_VECTOR )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( PERMEABILITY_MATRIX )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( LOCAL_PERMEABILITY_MATRIX )
        
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( NO_CONVERGENCE )
    
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( TOTAL_STRESS_TENSOR )
}

}  // namespace Python.
}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
