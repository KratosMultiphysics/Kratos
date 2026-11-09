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

/// @brief Extracts final total Cauchy stresses and assigns imported Gauss-point
///        initial stresses for lightweight Dam state files.
/// @details All elements are validated before any destination constitutive law
/// is modified. Displacements, reactions, time, step counters, loads and any
/// other source variable are never transferred.
/// Supported geometries are the 2D T3/T6/Q4/Q8/Q9 and 3D T4/T10/P6/P15/H8/H20/H27
/// SmallDisplacement solid elements, with linear elastic or linear thermoelastic
/// constitutive laws. The actual element integration rule is used.
class KRATOS_API(DAM_APPLICATION) DamInitialStateTransferUtility
{
public:
    /**
     * @brief Returns the final total Cauchy stress vector at every integration point.
     * @param rModelPart Model part whose elements are queried.
     * @return Stresses keyed by element ID.
     */
    static std::map<IndexType, std::vector<Vector>> GetFinalCauchyStresses(
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
     * @brief Validates imported initial stresses against the destination laws
     *        without modifying any constitutive state.
     * @param rDestinationModelPart Fresh destination model part.
     * @param rStresses Imported stresses keyed by destination element ID.
     * @param CompensateInitialThermalStrain When true, an initial strain equal to
     *        minus the current thermoelastic strain is prepared for each Gauss
     *        point, so the imported stress is reproduced at zero displacement even
     *        when T != T_ref.
     */
    static void ValidateInitialStresses(
        ModelPart& rDestinationModelPart,
        const std::map<IndexType, std::vector<Vector>>& rStresses,
        const bool CompensateInitialThermalStrain = false);

    /**
     * @brief Assigns imported per-integration-point stresses to the destination
     *        constitutive InitialState objects.
     * @param rDestinationModelPart Fresh destination model part.
     * @param rStresses Imported stresses keyed by destination element ID.
     * @param CompensateInitialThermalStrain When true, the constitutive initial
     *        strain is set to minus the current thermoelastic strain. This makes
     *        the imported stress the initial total stress at zero displacement
     *        while preserving the imported temperature and reference fields.
     */
    static void AssignInitialStresses(
        ModelPart& rDestinationModelPart,
        const std::map<IndexType, std::vector<Vector>>& rStresses,
        const bool CompensateInitialThermalStrain = false);
};

} // namespace Kratos

#endif // KRATOS_SELECTIVE_INITIAL_STRESS_UTILITY_H_INCLUDED
