// System includes
#include <string>
#include <iostream>
#include <cmath>

// Project includes
#include "DEM_KDEM_soft_torque_with_noise_CL.h"
#include "custom_elements/spheric_continuum_particle.h"

namespace Kratos {

    void DEM_KDEM_soft_torque_with_noise::Initialize() {

        KRATOS_TRY
        KRATOS_CATCH("")
    }

    DEMContinuumConstitutiveLaw::Pointer DEM_KDEM_soft_torque_with_noise::Clone() const {
        DEMContinuumConstitutiveLaw::Pointer p_clone(new DEM_KDEM_soft_torque_with_noise(*this));
        return p_clone;
    }

    void DEM_KDEM_soft_torque_with_noise::SetConstitutiveLawInProperties(Properties::Pointer pProp, bool verbose) const {

        KRATOS_INFO("DEM") << "Assigning DEM_KDEM_soft_torque_with_noise to Properties " << pProp->Id() << std::endl;
        pProp->SetValue(DEM_CONTINUUM_CONSTITUTIVE_LAW_POINTER, this->Clone());
        this->Check(pProp);
    }

    void DEM_KDEM_soft_torque_with_noise::Check(Properties::Pointer pProp) const {

        DEM_KDEM_soft_torque::Check(pProp);

        if (!pProp->Has(KDEM_STANDARD_DEVIATION_TAU_ZERO)) {
            KRATOS_WARNING("DEM")<<std::endl;
            KRATOS_WARNING("DEM")<<"WARNING: Variable KDEM_STANDARD_DEVIATION_TAU_ZERO should be present in the properties when using DEM_KDEM_soft_torque_with_noise. 0.0 value assigned by default."<<std::endl;
            KRATOS_WARNING("DEM")<<std::endl;
            pProp->GetValue(KDEM_STANDARD_DEVIATION_TAU_ZERO) = 0.0;
        }
        if (!pProp->Has(KDEM_STANDARD_DEVIATION_FRICTION)) {
            KRATOS_WARNING("DEM")<<std::endl;
            KRATOS_WARNING("DEM")<<"WARNING: Variable KDEM_STANDARD_DEVIATION_FRICTION should be present in the properties when using DEM_KDEM_soft_torque_with_noise. 0.0 value assigned by default."<<std::endl;
            KRATOS_WARNING("DEM")<<std::endl;
            pProp->GetValue(KDEM_STANDARD_DEVIATION_FRICTION) = 0.0;
        }
    }

    double DEM_KDEM_soft_torque_with_noise::GetTauZero(SphericContinuumParticle* element1) {

        double id = element1->GetId();
        srand(id);
        return rand_normal(DEM_KDEM::GetTauZero(element1), element1->GetProperties()[KDEM_STANDARD_DEVIATION_TAU_ZERO]);
    }

    double DEM_KDEM_soft_torque_with_noise::GetInternalFricc(SphericContinuumParticle* element1) {

        double id = element1->GetId();
        srand(id);
        double internal_fricc = rand_normal(DEM_KDEM::GetInternalFricc(element1), element1->GetProperties()[KDEM_STANDARD_DEVIATION_FRICTION]);
        KRATOS_WATCH(id)
        KRATOS_WATCH(internal_fricc)
        return internal_fricc;
    }

    double DEM_KDEM_soft_torque_with_noise::rand_normal(const double mean, const double stddev) {

        KRATOS_TRY

        if (!stddev) return mean;
        double x, y, r;

        do {
            x = 2.0 * rand() / RAND_MAX - 1;
            y = 2.0 * rand() / RAND_MAX - 1;
            r = x*x + y*y;
        } while (r == 0.0 || r > 1.0);

        double d = sqrt(- 2.0 * log(r) / r);
        return x * d * stddev + mean;

        KRATOS_CATCH("")
    }

} // namespace Kratos
