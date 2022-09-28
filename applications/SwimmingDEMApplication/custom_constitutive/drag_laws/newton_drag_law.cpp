// Author: Guillermo Casas (gcasas@cimne.upc.edu)
// Date: February 2019
#include "swimming_DEM_application.h"
#include "newton_drag_law.h"

namespace Kratos {

    DragLaw::Pointer NewtonDragLaw::Clone() const {
        NewtonDragLaw::Pointer p_clone(new NewtonDragLaw(*this));
        return p_clone;
    }

    void NewtonDragLaw::Initialize(const ProcessInfo& r_process_info) {}

    std::string NewtonDragLaw::GetTypeOfLaw() {
        std::string type_of_law = "NewtonDragLaw";
        return type_of_law;
    }

    void NewtonDragLaw::ComputeForce(SphericParticle* p_particle,
                                       const double reynolds_number,
                                       double particle_radius,
                                       double fluid_density,
                                       double fluid_kinematic_viscosity,
                                       array_1d<double, 3>& minus_slip_velocity,
                                       array_1d<double, 3>& drag_force,
                                       const ProcessInfo& r_current_process_info)
    {
        double drag_coeff  = 0.5 * Globals::Pi * SWIMMING_POW_2(particle_radius) * fluid_density * SWIMMING_MODULUS_3(minus_slip_velocity);

        drag_coeff *= 0.44;
        noalias(drag_force) = drag_coeff * minus_slip_velocity;
    }
} // namespace Kratos
