// System includes
#include <string>
#include <iostream>
#include <cmath>

// Project includes
#include "DEM_KDEM_with_damage_parallel_bond_bilinear_CL.h"
#include "custom_elements/spheric_continuum_particle.h"

namespace Kratos {

    DEMContinuumConstitutiveLaw::Pointer DEM_KDEM_with_damage_parallel_bond_bilinear::Clone() const {
        DEMContinuumConstitutiveLaw::Pointer p_clone(new DEM_KDEM_with_damage_parallel_bond_bilinear(*this));
        return p_clone;
    }

    void DEM_KDEM_with_damage_parallel_bond_bilinear::SetConstitutiveLawInProperties(Properties::Pointer pProp, bool verbose) {
        KRATOS_INFO("DEM") << "Assigning DEM_KDEM_with_damage_parallel_bond to Properties " << pProp->Id() << std::endl;
        pProp->SetValue(DEM_CONTINUUM_CONSTITUTIVE_LAW_POINTER, this->Clone());
        this->Check(pProp);
        SetDebugPrintingOptionValue(pProp);
    }

    void DEM_KDEM_with_damage_parallel_bond_bilinear::Check(Properties::Pointer pProp) const {

        DEM_KDEM_with_damage_parallel_bond::Check(pProp);

        if (!pProp->Has(SIGMA_SLOPE_CHANGE_THRESHOLD)) {
            KRATOS_WARNING("DEM") << std::endl;
            KRATOS_WARNING("DEM") << "WARNING: Variable SIGMA_SLOPE_CHANGE_THRESHOLD was not found in the Properties when using DEM_KDEM_with_damage_parallel_bond_bilinear. A default value of 0.0 was assigned."<<std::endl;
            KRATOS_WARNING("DEM") << std::endl;
            pProp->GetValue(SIGMA_SLOPE_CHANGE_THRESHOLD) = 0.0;
        }

        if (!pProp->Has(INTERNAL_FRICTION_AFTER_THRESHOLD)) {
            KRATOS_WARNING("DEM") << std::endl;
            KRATOS_WARNING("DEM") << "WARNING: Variable INTERNAL_FRICTION_AFTER_THRESHOLD was not found in the Properties when using DEM_KDEM_with_damage_parallel_bond_bilinear. A default value of 0.0 was assigned."<<std::endl;
            KRATOS_WARNING("DEM") << std::endl;
            pProp->GetValue(INTERNAL_FRICTION_AFTER_THRESHOLD) = 0.0;
        }
    }

    void DEM_KDEM_with_damage_parallel_bond_bilinear::AdjustTauStrengthAndUpdatedMaxTauStrength(double& tau_strength, double& updated_max_tau_strength, const double internal_friction,
                                                                                                double contact_sigma, SphericContinuumParticle* element1, SphericContinuumParticle* element2) {
        KRATOS_TRY
        const double sigma_slope_change_threshold = 0.5 * (element1->GetProperties()[SIGMA_SLOPE_CHANGE_THRESHOLD] + element2->GetProperties()[SIGMA_SLOPE_CHANGE_THRESHOLD]);
        const double internal_friction_after_threshold = 0.5 * (element1->GetProperties()[INTERNAL_FRICTION_AFTER_THRESHOLD] + element2->GetProperties()[INTERNAL_FRICTION_AFTER_THRESHOLD]);

        if (contact_sigma <= sigma_slope_change_threshold) {
            tau_strength += (1.0 - mDamageTangential) * internal_friction * contact_sigma;
            updated_max_tau_strength += internal_friction * contact_sigma;
        } else {
            const double aux = (internal_friction * sigma_slope_change_threshold + internal_friction_after_threshold * (contact_sigma - sigma_slope_change_threshold));
            tau_strength += (1.0 - mDamageTangential) * aux;
            updated_max_tau_strength += aux;
        }
        KRATOS_CATCH("")
    }
} // namespace Kratos
