//
//   Project Name:        KratosDamApplication
//
//

#if !defined(KRATOS_DAM_NONLOCAL_DAMAGE_UTILITIES_H_INCLUDED)
#define KRATOS_DAM_NONLOCAL_DAMAGE_UTILITIES_H_INCLUDED

// Project includes
#include "includes/model_part.h"
#include "utilities/parallel_utilities.h"

// Application includes
#include "dam_application_variables.h"
#include "custom_constitutive/thermal_nonlocal_damage_3D_law.hpp"

namespace Kratos
{

/**
 * @class DamNonlocalDamageUtilities
 * @ingroup DamApplication
 * @brief Dam-only orchestration of the LOCAL_EQUIVALENT_STRAIN production for
 * the thermal nonlocal-damage workflow.
 *
 * The Poromechanics nonlocal strategy owns the full sequence
 *   LOCAL -> spatial averaging -> NONLOCAL.
 * Dam owns only
 *   current element kinematics -> current LOCAL.
 *
 * This utility therefore only invokes
 *   Element::CalculateOnIntegrationPoints(LOCAL_EQUIVALENT_STRAIN, ...)
 * on the applicable active elements; the element provides the kinematics and
 * the constitutive law computes and stores the LOCAL driving quantity (which
 * the Poromechanics averaging utility later reads via GetValue). It does NOT
 * compute B matrices, displacement gradients, thermal strain, nor nonlocal
 * averaging, and it never touches NONLOCAL_EQUIVALENT_STRAIN.
 */
class DamNonlocalDamageUtilities
{
public:
    /// True when the scheme-based LOCAL ownership is active (set by the Dam
    /// solvers from the 'nonlocal_damage' setting).
    static bool IsProcessBasedLocalOwnership(const ModelPart& rModelPart)
    {
        return rModelPart.GetProcessInfo().Has(USE_PROCESS_BASED_LOCAL_EQUIVALENT_STRAIN)
            && rModelPart.GetProcessInfo()[USE_PROCESS_BASED_LOCAL_EQUIVALENT_STRAIN];
    }

    /// Filters the elements that participate in the thermal nonlocal-damage
    /// LOCAL production. Only elements whose constitutive-law type belongs to
    /// the thermal nonlocal-damage family (ThermalNonlocalDamage3DLaw and its
    /// Simo-Ju / Modified-Mises / plane-strain / plane-stress derived classes)
    /// are applicable. The check is done on the constitutive law stored in the
    /// element Properties; no element-name string inspection is used.
    static bool IsApplicable(Element& rElement)
    {
        if (!rElement.GetProperties().Has(CONSTITUTIVE_LAW)) {
            return false;
        }
        const ConstitutiveLaw::Pointer& r_prototype_law = rElement.GetProperties().GetValue(CONSTITUTIVE_LAW);
        if (r_prototype_law == nullptr) {
            return false;
        }
        return dynamic_cast<const ThermalNonlocalDamage3DLaw*>(r_prototype_law.get()) != nullptr;
    }

    /**
     * @brief Recomputes the current LOCAL_EQUIVALENT_STRAIN on every applicable
     * active element through its generic integration-point interface.
     * @param rModelPart the model part whose elements are updated.
     * @note The operation is element-local (per integration-point constitutive
     * law) and performs no shared-node accumulation, so it can be executed in a
     * parallel element loop without nodal locks.
     */
    static void CalculateLocalEquivalentStrain(ModelPart& rModelPart)
    {
        const ProcessInfo& r_process_info = rModelPart.GetProcessInfo();
        block_for_each(rModelPart.Elements(), [&r_process_info](Element& r_element) {
            // Respect the ACTIVE convention; skip unrelated/inactive elements.
            if (!r_element.IsActive()) return;
            if (!IsApplicable(r_element)) return;
            // The element provides the kinematics; the law computes and stores
            // the LOCAL quantity. The returned vector is discarded (the law
            // stores the same value internally for the Poromechanics utility).
            std::vector<double> local_values;
            r_element.CalculateOnIntegrationPoints(
                LOCAL_EQUIVALENT_STRAIN, local_values, r_process_info);
        });
    }

private:
    DamNonlocalDamageUtilities() = delete;
};

} // namespace Kratos

#endif // KRATOS_DAM_NONLOCAL_DAMAGE_UTILITIES_H_INCLUDED defined
