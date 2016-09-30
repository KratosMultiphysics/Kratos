// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Vicente Mataix
//

// System includes

#if defined(KRATOS_PYTHON)
// External includes
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "contact_structural_mechanics_application_variables.h"
#include "contact_structural_mechanics_application.h"
#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_python/add_mappers_to_python.h"

namespace Kratos
{

namespace Python
{

using namespace boost::python;



BOOST_PYTHON_MODULE(KratosContactStructuralMechanicsApplication)
{

    class_<KratosContactStructuralMechanicsApplication,
           KratosContactStructuralMechanicsApplication::Pointer,
           bases<KratosApplication>, boost::noncopyable >("KratosContactStructuralMechanicsApplication")
           ;

    AddCustomStrategiesToPython();
    AddCustomUtilitiesToPython();
    AddCustomMappersToPython();

    //registering variables in python

    // CONDITIONS
    // CONTACT
    /* Mortar contact */
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( WEIGHTED_GAP )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( IS_ACTIVE_SET )
}


}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
