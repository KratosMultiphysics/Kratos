//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                 May 2018 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_solvers/solution_local_flags.hpp"

namespace Kratos
{

///@name Type Definitions
///@{


/**
 * Flags for the solution control
 */
KRATOS_CREATE_LOCAL_FLAG( SolverLocalFlags, INITIALIZED,               0 );
KRATOS_CREATE_LOCAL_FLAG( SolverLocalFlags, CONVERGED,                 1 );
KRATOS_CREATE_LOCAL_FLAG( SolverLocalFlags, DOFS_INITIALIZED,          2 );
KRATOS_CREATE_LOCAL_FLAG( SolverLocalFlags, ELEMENTS_INITIALIZED,      3 );
KRATOS_CREATE_LOCAL_FLAG( SolverLocalFlags, CONDITIONS_INITIALIZED,    4 );
KRATOS_CREATE_LOCAL_FLAG( SolverLocalFlags, ADAPTIVE_SOLUTION,         5 );

/**
 * Flags for the solution options
 */
KRATOS_CREATE_LOCAL_FLAG( SolverLocalFlags, MOVE_MESH,                 0 );
KRATOS_CREATE_LOCAL_FLAG( SolverLocalFlags, UPDATE_VARIABLES,          1 );
KRATOS_CREATE_LOCAL_FLAG( SolverLocalFlags, REFORM_DOFS,               2 );
KRATOS_CREATE_LOCAL_FLAG( SolverLocalFlags, INCREMENTAL_SOLUTION,      3 );
KRATOS_CREATE_LOCAL_FLAG( SolverLocalFlags, COMPUTE_REACTIONS,         4 );
KRATOS_CREATE_LOCAL_FLAG( SolverLocalFlags, CONSTANT_SYSTEM_MATRIX,    5 );
KRATOS_CREATE_LOCAL_FLAG( SolverLocalFlags, RAYLEIGH_DAMPING,          6 );
KRATOS_CREATE_LOCAL_FLAG( SolverLocalFlags, IMPLEX,                    7 );


/**
 * Flags for the convergence criterion control
 */
KRATOS_CREATE_LOCAL_FLAG( CriterionLocalFlags, INITIALIZED,               0 );
KRATOS_CREATE_LOCAL_FLAG( CriterionLocalFlags, INCREMENTAL,               1 );
KRATOS_CREATE_LOCAL_FLAG( CriterionLocalFlags, CONVERGED,                 2 );
KRATOS_CREATE_LOCAL_FLAG( CriterionLocalFlags, AND,                       3 );
KRATOS_CREATE_LOCAL_FLAG( CriterionLocalFlags, OR,                        4 );
KRATOS_CREATE_LOCAL_FLAG( CriterionLocalFlags, UPDATE_RHS,                5 );
KRATOS_CREATE_LOCAL_FLAG( CriterionLocalFlags, SUPPLIED_DOF,              6 );

/**
 * Flags for the time integration options
 */
KRATOS_CREATE_LOCAL_FLAG( TimeIntegrationLocalFlags, PREDICT_PRIMARY_VARIABLE,   0 );

///@}
///@name Input and output
///@{

///@}

///@} addtogroup block

}  // namespace Kratos.

