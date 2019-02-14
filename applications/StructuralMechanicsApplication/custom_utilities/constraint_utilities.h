// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_CONSTRAINT_UTILITIES)
#define KRATOS_CONSTRAINT_UTILITIES

// System includes

// External includes

// Project includes
#include "includes/model_part.h"

namespace Kratos
{
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
/**
 * @namespace ConstraintUtilities
 * @ingroup StructuralMechanicsApplication
 * @brief This namespace includes several utilities necessaries for the computation of the MPC
 * @author Vicente Mataix Ferrandiz
 */
namespace ConstraintUtilities
{
    /**
     * @brief This method resets the values of the slave dofs
     * @param rModelPart The model of the problem to solve
     */
    void KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) ResetSlaveDofs(ModelPart& rModelPart);

    /**
     * @brief This method resets the values of the slave dofs
     * @param rModelPart The model of the problem to solve
     */
    void KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) ApplyConstraints(ModelPart& rModelPart);

}; // namespace ConstraintUtilities
}  // namespace Kratos
#endif /* KRATOS_CONSTRAINT_UTILITIES defined */
