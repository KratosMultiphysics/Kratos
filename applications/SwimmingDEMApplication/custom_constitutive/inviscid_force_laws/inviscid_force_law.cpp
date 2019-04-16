// Author: Guillermo Casas (gcasas@cimne.upc.edu)
// Date: February 2019
#include "swimming_DEM_application.h"
#include "inviscid_force_law.h"

namespace Kratos {

    void InviscidForceLaw::Initialize(const ProcessInfo& r_process_info) {}

    std::string InviscidForceLaw::GetTypeOfLaw() {
        std::string type_of_law = "Generic inviscid force law";
        return type_of_law;
    }

    void InviscidForceLaw::SetInviscidForceLawInProperties(Properties::Pointer pProp) const {
        pProp->SetValue(SDEM_INVISCID_FORCE_LAW_POINTER, this->Clone());
    }

    InviscidForceLaw::Pointer InviscidForceLaw::Clone() const {
        InviscidForceLaw::Pointer p_clone(new InviscidForceLaw(*this));
        return p_clone;
    }

    double InviscidForceLaw::ComputeParticleAccelerationNumber(const double particle_radius,
                                                               const array_1d<double, 3>& minus_slip_velocity,
                                                               const array_1d<double, 3>& minus_slip_acceleration)
    {
        const double norm_of_slip_vel = SWIMMING_MODULUS_3(minus_slip_velocity);
        return SWIMMING_POW_3(norm_of_slip_vel) / std::abs(2 * particle_radius * SWIMMING_INNER_PRODUCT_3(minus_slip_velocity, minus_slip_acceleration));
    }

} // namespace Kratos
