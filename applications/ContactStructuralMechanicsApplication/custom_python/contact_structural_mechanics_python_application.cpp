// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: StructuralMechanicsApplication/license.txt
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
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( INTEGRATION_ORDER_CONTACT )    // The integration order considered
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( ACTIVE_CHECK_FACTOR )          // The factor employed to search an active/inactive node
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( WEIGHTED_GAP )                 // The integrated gap employed in mortar formulation
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( WEIGHTED_SLIP )                // The integrated slip employed in mortar formulation
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( WEIGHTED_FRICTION )            // The integrated friction employed in mortar formulation
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( AUXILIAR_ACTIVE )              // Auxiliar boolean to check if the node is active or not
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( AUXILIAR_SLIP )                // Auxiliar boolean to check if the node is stick or not
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( ACTIVE_CHECK_FACTOR )          // The factor employed to search an active/inactive node
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( NORMAL_AUGMENTATION_FACTOR )   // The constant that is considered for the check of active or inactive (when 0 it doesn't accept traction)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( TANGENT_AUGMENTATION_FACTOR )  // The constant that is considered for the check if the node is slip/stick
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( PENALTY_FACTOR )               // The penalty factor needed for the ALM
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( SCALE_FACTOR )                 // The scale factor needed for the ALM
}


}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
