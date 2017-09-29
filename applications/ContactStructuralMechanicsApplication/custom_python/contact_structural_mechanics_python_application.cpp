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
#include "custom_python/add_custom_processes_to_python.h"
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
    AddCustomProcessesToPython();
    AddCustomMappersToPython();

    //registering variables in python

    // CONDITIONS
    // CONTACT
    /* Mortar contact */
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( INTEGRATION_ORDER_CONTACT )                   // The integration order considered
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( ACTIVE_CHECK_FACTOR )                         // The factor employed to search an active/inactive node
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( WEIGHTED_GAP )                                // The integrated gap employed in mortar formulation
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( WEIGHTED_SLIP )                               // The integrated slip employed in mortar formulation
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( WEIGHTED_SCALAR_RESIDUAL )                    // The integrated scalar residual  
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( WEIGHTED_VECTOR_RESIDUAL ) // The integrated vector residual    
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( ACTIVE_CHECK_FACTOR )                         // The factor employed to search an active/inactive node
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( AUGMENTED_NORMAL_CONTACT_PRESSURE )           // The resultant augmented pressure in the normal direction
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( AUGMENTED_TANGENT_CONTACT_PRESSURE )          // The resultant augmented pressure in the tangent direction
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( TANGENT_FACTOR )                              // The proportion between the tangent and normal penalty
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( CONSIDER_NORMAL_VARIATION )                   // A value used to check if consider normal variation or not
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( CONSIDER_PAIR_VARIATION )                     // A value used to check if consider variation or not in the active inactive pairs
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( ADAPT_PENALTY )                               // To set if the penalty is recalculated or not
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( MAX_GAP_FACTOR )                              // The factor between the nodal H and the max gap considered to recalculate the penalty
}


}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
