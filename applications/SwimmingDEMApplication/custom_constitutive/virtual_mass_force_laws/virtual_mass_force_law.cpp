// Author: Guillermo Casas (gcasas@cimne.upc.edu)
// Date: February 2019
#include "swimming_DEM_application.h"
#include "virtual_mass_force_law.h"

namespace Kratos {

    void VirtualMassForceLaw::Initialize(const ProcessInfo& r_process_info) {}

    std::string VirtualMassForceLaw::GetTypeOfLaw() {
        std::string type_of_law = "Generic virtual mass force law";
        return type_of_law;
    }

    void VirtualMassForceLaw::SetVirtualMassForceLawInProperties(Properties::Pointer pProp) const {
        pProp->SetValue(SDEM_VIRTUAL_MASS_FORCE_LAW_POINTER, this->Clone());
    }

    VirtualMassForceLaw::Pointer VirtualMassForceLaw::Clone() const {
        VirtualMassForceLaw::Pointer p_clone(new VirtualMassForceLaw(*this));
        return p_clone;
    }

    double VirtualMassForceLaw::ComputeParticleAccelerationNumber(const double particle_radius,
                                                               const array_1d<double, 3>& minus_slip_velocity,
                                                               const array_1d<double, 3>& minus_slip_acceleration)
    {
        const double norm_of_slip_vel = SWIMMING_MODULUS_3(minus_slip_velocity);
        return SWIMMING_POW_3(norm_of_slip_vel) / std::abs(2 * particle_radius * SWIMMING_INNER_PRODUCT_3(minus_slip_velocity, minus_slip_acceleration));
    }

} // namespace Kratos
