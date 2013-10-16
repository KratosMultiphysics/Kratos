//
//   Project Name:        KratosSolidMechanicsApplication $
//   Last modified by:    $Author:            JMCarbonell $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes
#if defined(KRATOS_PYTHON)
// External includes
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_python/add_custom_constitutive_laws_to_python.h"

#include "solid_mechanics_application.h"

namespace Kratos
{

namespace Python
{

using namespace boost::python;



BOOST_PYTHON_MODULE(KratosSolidMechanicsApplication)
{

    class_<KratosSolidMechanicsApplication,
           KratosSolidMechanicsApplication::Pointer,
           bases<KratosApplication>, boost::noncopyable >("KratosSolidMechanicsApplication")
           ;

    AddCustomUtilitiesToPython();
    AddCustomStrategiesToPython();
    AddCustomConstitutiveLawsToPython();

    //registering variables in python ( if must to be seen from python )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( CONSTITUTIVE_LAW_NAME );
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( WRITE_ID );
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( PREVIOUS_DELTA_TIME );

    KRATOS_REGISTER_IN_PYTHON_VARIABLE( FORCE_INTERNAL );
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( FORCE_EXTERNAL);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( FORCE_DYNAMIC );
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( CROSS_AREA );

    KRATOS_REGISTER_IN_PYTHON_VARIABLE( IMPOSED_DISPLACEMENT );
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( REACTION_PRESSURE );

    KRATOS_REGISTER_IN_PYTHON_VARIABLE( VON_MISES_STRESS );
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( PLASTIC_STRAIN );
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( DELTA_PLASTIC_STRAIN );

}


}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
