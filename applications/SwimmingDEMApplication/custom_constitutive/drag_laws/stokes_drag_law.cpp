// Author: Guillermo Casas (gcasas@cimne.upc.edu)
// Date: February 2019

#include "stokes_drag_law.h"

namespace Kratos {

    DragLaw::Pointer StokesDragLaw::Clone() const {
        StokesDragLaw::Pointer p_clone(new StokesDragLaw(*this));
        return p_clone;
    }

    void StokesDragLaw::Initialize(const ProcessInfo& r_process_info) {}

    std::string StokesDragLaw::GetTypeOfLaw() {
        std::string type_of_law = "Stokes drag law";
        return type_of_law;
    }

    void StokesDragLaw::ComputeForce(SphericParticle* p_particle,
                                     const double reynolds_number,
                                     double particle_radius,
                                     double fluid_density,
                                     double fluid_kinematic_viscosity,
                                     array_1d<double, 3>& minus_slip_velocity,
                                     array_1d<double, 3>& drag_force,
                                     const ProcessInfo& r_current_process_info)
    {
        double drag_coeff = 6.0 * Globals::Pi * fluid_kinematic_viscosity * fluid_density * particle_radius;
        noalias(drag_force) = drag_coeff * minus_slip_velocity;
    }
} // namespace Kratos
