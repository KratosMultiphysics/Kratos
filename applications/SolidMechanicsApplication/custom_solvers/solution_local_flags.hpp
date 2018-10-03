//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:            November 2017 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_SOLUTION_LOCAL_FLAGS_H_INCLUDED)
#define  KRATOS_SOLUTION_LOCAL_FLAGS_H_INCLUDED

// System includes

// External includes

// Project includes
#include "containers/flags.h"

namespace Kratos
{
///@addtogroup SolidMechanicsApplication
///@{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/** @brief Solver local flags class definition
 *  @details This is the base class for solver local flags
 */
class KRATOS_API(SOLID_MECHANICS_APPLICATION) SolverLocalFlags
{
 public:
  /// Flags for the solution control:
  KRATOS_DEFINE_LOCAL_FLAG( INITIALIZED );
  KRATOS_DEFINE_LOCAL_FLAG( CONVERGED );
  KRATOS_DEFINE_LOCAL_FLAG( DOFS_INITIALIZED );
  KRATOS_DEFINE_LOCAL_FLAG( ELEMENTS_INITIALIZED );
  KRATOS_DEFINE_LOCAL_FLAG( CONDITIONS_INITIALIZED );
  KRATOS_DEFINE_LOCAL_FLAG( ADAPTIVE_SOLUTION );

  /// Flags for the solution options:
  KRATOS_DEFINE_LOCAL_FLAG( MOVE_MESH );
  KRATOS_DEFINE_LOCAL_FLAG( UPDATE_VARIABLES );
  KRATOS_DEFINE_LOCAL_FLAG( REFORM_DOFS );
  KRATOS_DEFINE_LOCAL_FLAG( INCREMENTAL_SOLUTION );
  KRATOS_DEFINE_LOCAL_FLAG( COMPUTE_REACTIONS );
  KRATOS_DEFINE_LOCAL_FLAG( CONSTANT_SYSTEM_MATRIX );
  KRATOS_DEFINE_LOCAL_FLAG( RAYLEIGH_DAMPING );
  KRATOS_DEFINE_LOCAL_FLAG( IMPLEX );
};


/** @brief Solver local flags class definition
 *  @details This is the base class for criterion local flags
 */
class KRATOS_API(SOLID_MECHANICS_APPLICATION) CriterionLocalFlags
{
 public:
  /// Flags for the solution control:
  KRATOS_DEFINE_LOCAL_FLAG( INITIALIZED );
  KRATOS_DEFINE_LOCAL_FLAG( INCREMENTAL );
  KRATOS_DEFINE_LOCAL_FLAG( CONVERGED );
  KRATOS_DEFINE_LOCAL_FLAG( AND );
  KRATOS_DEFINE_LOCAL_FLAG( OR );
  KRATOS_DEFINE_LOCAL_FLAG( UPDATE_RHS );
  KRATOS_DEFINE_LOCAL_FLAG( SUPPLIED_DOF );
};


/** @brief Solver local flags class definition
 *  @details This is the base class for time integration local flags
 */

class KRATOS_API(SOLID_MECHANICS_APPLICATION) TimeIntegrationLocalFlags
{
 public:

  /// Flags for the solution options:
  KRATOS_DEFINE_LOCAL_FLAG( PREDICT_PRIMARY_VARIABLE );

  /// Flags for the solution control:
}; // Class TimeIntegrationLocalFlags


///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_SOLUTION_LOCAL_FLAGS_H_INCLUDED defined
