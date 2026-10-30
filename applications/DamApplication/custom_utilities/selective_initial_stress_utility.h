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

#if !defined(KRATOS_SELECTIVE_INITIAL_STRESS_UTILITY_H_INCLUDED)
#define KRATOS_SELECTIVE_INITIAL_STRESS_UTILITY_H_INCLUDED

// Project includes
#include "includes/model_part.h"

namespace Kratos
{

/// @brief Assigns selected Gauss-point stresses from a source model part to the
///        constitutive initial state of an independent destination model part.
/// @details This is a same-mesh, stress-only operation. It does not transfer
/// displacements, reactions, time, step counters, loads or any other source
/// variable, and it does not restart the destination analysis. All elements are
/// validated before any destination constitutive law is modified.
class KRATOS_API(DAM_APPLICATION) SelectiveInitialStressUtility
{
public:
    /**
     * @brief Validates the source/destination pair and assigns one InitialState
     *        per integration point using the source Cauchy stress.
     * @param rSourceModelPart Source model part loaded from a restart checkpoint.
     * @param rDestinationModelPart Fresh model part owned by the new analysis.
     * @throws If the configuration is outside the supported validated scope.
     */
    static void AssignInitialStresses(
        ModelPart& rSourceModelPart,
        ModelPart& rDestinationModelPart);
};

} // namespace Kratos

#endif // KRATOS_SELECTIVE_INITIAL_STRESS_UTILITY_H_INCLUDED
