#include "swimming_DEM_application.h"
#include "chien_drag_law.h"

namespace Kratos {

    DragLaw::Pointer ChienDragLaw::Clone() const {
        ChienDragLaw::Pointer p_clone(new ChienDragLaw(*this));
        return p_clone;
    }

    void ChienDragLaw::Initialize(const ProcessInfo& r_process_info) {}

    std::string ChienDragLaw::GetTypeOfLaw() {
        std::string type_of_law = "ChienDragLaw";
        return type_of_law;
    }

    void ChienDragLaw::ComputeForce(SphericParticle* p_particle,
                                       const double reynolds_number,
                                       double particle_radius,
                                       double fluid_density,
                                       double fluid_kinematic_viscosity,
                                       array_1d<double, 3>& minus_slip_velocity,
                                       array_1d<double, 3>& drag_force,
                                       const ProcessInfo& r_current_process_info)
    {
        Geometry<Node<3> > geometry = p_particle->GetGeometry();
        const double sphericity = geometry[0].FastGetSolutionStepValue(PARTICLE_SPHERICITY);
        double drag_coeff       = 0.5 * Globals::Pi * SWIMMING_POW_2(particle_radius) * fluid_density * SWIMMING_MODULUS_3(minus_slip_velocity);
        drag_coeff *= 30 / reynolds_number + 67.289 * exp(- 5.03 * sphericity);

        noalias(drag_force) = drag_coeff * minus_slip_velocity;
    }
} // namespace Kratos
