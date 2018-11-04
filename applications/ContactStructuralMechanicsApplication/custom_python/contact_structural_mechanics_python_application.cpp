// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

#if defined(KRATOS_PYTHON)
// External includes


// Project includes
#include "includes/define.h"
#include "includes/define_python.h"
#include "contact_structural_mechanics_application_variables.h"
#include "contact_structural_mechanics_application.h"
#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_python/add_custom_processes_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_python/add_mappers_to_python.h"
#include "custom_python/add_custom_linear_solvers_to_python.h"

namespace Kratos
{
namespace Python
{
namespace py = pybind11;

PYBIND11_MODULE(KratosContactStructuralMechanicsApplication, m)
{
    py::class_<KratosContactStructuralMechanicsApplication,
        KratosContactStructuralMechanicsApplication::Pointer,
        KratosApplication >(m, "KratosContactStructuralMechanicsApplication")
        .def(py::init<>())
        ;

    AddCustomStrategiesToPython(m);
    AddCustomUtilitiesToPython(m);
    AddCustomProcessesToPython(m);
    AddCustomMappersToPython(m);
    AddCustomLinearSolversToPython(m);

    // Adding enums
    py::enum_<NormalDerivativesComputation>(m, "NormalDerivativesComputation")
    .value("NO_DERIVATIVES_COMPUTATION",NO_DERIVATIVES_COMPUTATION)
    .value("ELEMENTAL_DERIVATIVES",ELEMENTAL_DERIVATIVES)
    .value("NODAL_ELEMENTAL_DERIVATIVES",NODAL_ELEMENTAL_DERIVATIVES)
    .value("NO_DERIVATIVES_COMPUTATION_WITH_NORMAL_UPDATE",NO_DERIVATIVES_COMPUTATION_WITH_NORMAL_UPDATE)
    ;

    //Registering variables in python

    // CONDITIONS
    // CONTACT
    /* Mortar contact */
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, AUXILIAR_COORDINATES )               // Auxiliar coordinates used to map
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, DELTA_COORDINATES )                  // Delta coordinates used to map
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, INTEGRATION_ORDER_CONTACT )                             // The integration order considered
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, DISTANCE_THRESHOLD )                                    // The distance threshold considered
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, ACTIVE_CHECK_FACTOR )                                   // The factor employed to search an active/inactive node
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, NORMAL_GAP )                                            // The normal gap employed in contact formulation
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, WEIGHTED_GAP )                                          // The integrated gap employed in mortar formulation
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, WEIGHTED_SLIP )                      // The integrated slip employed in mortar formulation
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, WEIGHTED_SCALAR_RESIDUAL )                              // The integrated scalar residual
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, WEIGHTED_VECTOR_RESIDUAL )           // The integrated vector residual
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, ACTIVE_CHECK_FACTOR )                                   // The factor employed to search an active/inactive node
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, ACTIVE_SET_CONVERGED )                                  // To know if the active set has converged
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, DYNAMIC_FACTOR )                                        // The factor considered for dynamic problems (in order to take intro account the gap evolution)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, LAGRANGE_MULTIPLIER_CONTACT_PRESSURE )                  // The lagrange multiplier for normal contact pressure
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, AUGMENTED_NORMAL_CONTACT_PRESSURE )                     // The resultant augmented pressure in the normal direction
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, AUGMENTED_TANGENT_CONTACT_PRESSURE ) // The resultant augmented pressure in the tangent direction
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, TANGENT_FACTOR )                                        // The proportion between the tangent and normal penalty
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, CONSIDER_NORMAL_VARIATION )                             // A value used to check if consider normal variation or not
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, ADAPT_PENALTY )                                         // To set if the penalty is recalculated or not
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MAX_GAP_FACTOR )                                        // The factor between the nodal H and the max gap considered to recalculate the penalty
}


}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
