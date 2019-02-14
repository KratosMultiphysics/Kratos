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

// System includes

// External includes

// Project includes
#include "custom_utilities/constraint_utilities.h"

namespace Kratos
{
namespace ConstraintUtilities
{
void ResetSlaveDofs(ModelPart& rModelPart)
{
    KRATOS_TRY

    // The number of constraints
    const int number_of_constraints = static_cast<int>(rModelPart.MasterSlaveConstraints().size());

    // Setting to zero the slave dofs
    #pragma omp parallel
    {
        #pragma omp for schedule(guided, 512)
        for (int i_const = 0; i_const < number_of_constraints; ++i_const) {
            auto it_const = rModelPart.MasterSlaveConstraints().begin() + i_const;

            // Detect if the constraint is active or not. If the user did not make any choice the constraint
            // It is active by default
            bool constraint_is_active = true;
            if (it_const->IsDefined(ACTIVE))
                constraint_is_active = it_const->Is(ACTIVE);

            if (constraint_is_active) {
                const typename MasterSlaveConstraint::DofPointerVectorType& r_slave_equation_dofs = it_const->GetSlaveDofsVector();

                for (std::size_t i = 0; i < r_slave_equation_dofs.size(); ++i) {
                    #pragma omp atomic
                    r_slave_equation_dofs[i]->GetSolutionStepValue() *= 0.0;
                }
            }
        }
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void ApplyConstraints(ModelPart& rModelPart)
{
    KRATOS_TRY

    // The number of constraints
    const int number_of_constraints = static_cast<int>(rModelPart.MasterSlaveConstraints().size());

    // The current process info
    const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

    // Adding MPC contribution
    #pragma omp parallel
    {
        #pragma omp for schedule(guided, 512)
        for (int i_const = 0; i_const < number_of_constraints; ++i_const) {
            auto it_const = rModelPart.MasterSlaveConstraints().begin() + i_const;

            // Detect if the constraint is active or not. If the user did not make any choice the constraint
            // It is active by default
            bool constraint_is_active = true;
            if (it_const->IsDefined(ACTIVE))
                constraint_is_active = it_const->Is(ACTIVE);

            if (constraint_is_active) {
                it_const->Apply(r_current_process_info);
            }
        }
    }

    KRATOS_CATCH("")
}

} // namespace ConstraintUtilities
} // namespace Kratos
