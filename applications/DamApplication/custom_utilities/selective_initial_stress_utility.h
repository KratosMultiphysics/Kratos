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
#include <map>
#include <vector>

#include "includes/model_part.h"

namespace Kratos
{

/// @brief Assigns selected Gauss-point stresses from a source model part to the
///        constitutive initial state of an independent destination model part.
/// @details This is a same-mesh, stress-only operation. It does not transfer
/// displacements, reactions, time, step counters, loads or any other source
/// variable, and it does not restart the destination analysis. All elements are
/// validated before any destination constitutive law is modified.
/// Supported geometries are the 2D T3/T6/Q4/Q8/Q9 and 3D T4/T10/P6/P15/H8/H20/H27
/// SmallDisplacement solid elements, with linear elastic or linear thermoelastic
/// constitutive laws. The actual element integration rule is used and must match
/// between source and destination.
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

    /**
     * @brief Returns the current Cauchy stress vector at every integration point.
     * @param rModelPart Model part whose elements are queried.
     * @return Stresses keyed by element ID.
     */
    static std::map<IndexType, std::vector<Vector>> GetCurrentCauchyStresses(
        ModelPart& rModelPart);

    /**
     * @brief Returns the stress stored in each constitutive-law InitialState.
     * @details Missing or malformed stored states are rejected; this method never
     *          falls back to the current Cauchy stress.
     * @param rModelPart Model part whose elements are queried.
     * @return Stored initial stresses keyed by element ID.
     */
    static std::map<IndexType, std::vector<Vector>> GetStoredInitialStresses(
        ModelPart& rModelPart);

    /**
     * @brief Returns the actual integration-method name of an element.
     */
    static std::string GetIntegrationMethodName(const Element& rElement);

    /**
     * @brief Returns the integration-point local coordinates of an element.
     */
    static std::vector<Vector> GetIntegrationPointCoordinates(const Element& rElement);

    /**
     * @brief Validates stored initial stresses against the destination laws
     *        without modifying any constitutive state.
     * @param rDestinationModelPart Fresh destination model part.
     * @param rStresses Imported stresses keyed by destination element ID.
     */
    static void ValidateStoredInitialStresses(
        ModelPart& rDestinationModelPart,
        const std::map<IndexType, std::vector<Vector>>& rStresses);

    /**
     * @brief Assigns imported per-integration-point stresses to the destination
     *        constitutive InitialState objects.
     * @param rDestinationModelPart Fresh destination model part.
     * @param rStresses Imported stresses keyed by destination element ID.
     */
    static void AssignStoredInitialStresses(
        ModelPart& rDestinationModelPart,
        const std::map<IndexType, std::vector<Vector>>& rStresses);
};

} // namespace Kratos

#endif // KRATOS_SELECTIVE_INITIAL_STRESS_UTILITY_H_INCLUDED
