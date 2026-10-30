//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    DamApplication contributors
//

#if !defined(KRATOS_DAM_INITIAL_STRESS_UTILITIES_H_INCLUDED)
#define KRATOS_DAM_INITIAL_STRESS_UTILITIES_H_INCLUDED

// System includes
#include <algorithm>

// Project includes
#include "includes/checks.h"
#include "includes/model_part.h"
#include "utilities/parallel_utilities.h"

// Application includes
#include "dam_application_variables.h"

namespace Kratos
{

/// @brief Utility functions for Dam initial-stress nodal storage.
class DamInitialStressUtilities
{
public:
    /**
     * @brief Ensures that every node stores INITIAL_STRESS_TENSOR as a 3x3 matrix.
     * @details Existing 3x3 tensors are preserved. Empty (0x0) tensors are
     * initialized to zero, and 2x2 tensors are copied into the top-left block
     * of a zero 3x3 tensor. Any other size is rejected.
     * The operation is idempotent.
     * @param rModelPart Model part whose nodes are initialized.
     */
    static void EnsureInitialStressTensorStorage(ModelPart& rModelPart)
    {
        block_for_each(rModelPart.Nodes(), [](Node& rNode) {
            Matrix& r_initial_stress = rNode.FastGetSolutionStepValue(INITIAL_STRESS_TENSOR);

            const std::size_t rows = r_initial_stress.size1();
            const std::size_t cols = r_initial_stress.size2();

            const bool is_empty = (rows == 0 && cols == 0);
            const bool is_2d = (rows == 2 && cols == 2);
            const bool is_3d = (rows == 3 && cols == 3);

            KRATOS_ERROR_IF_NOT(is_empty || is_2d || is_3d)
                << "DamInitialStressUtilities: unsupported INITIAL_STRESS_TENSOR size "
                << rows << "x" << cols << " at node " << rNode.Id()
                << ". Supported sizes are 0x0 (empty), 2x2 and 3x3." << std::endl;

            if (is_3d) {
                return;
            }

            Matrix resized_initial_stress = ZeroMatrix(3, 3);
            for (std::size_t i = 0; i < rows; ++i) {
                for (std::size_t j = 0; j < cols; ++j) {
                    resized_initial_stress(i, j) = r_initial_stress(i, j);
                }
            }
            r_initial_stress = resized_initial_stress;
        });
    }
};

} // namespace Kratos

#endif // KRATOS_DAM_INITIAL_STRESS_UTILITIES_H_INCLUDED
