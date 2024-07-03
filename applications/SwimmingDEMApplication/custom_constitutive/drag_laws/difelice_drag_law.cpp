// Author: Joaquin Gonzalez-Usua (jgonzalez@cimne.upc.edu)
// Date: May 2023
//This drag law is explained and implemented in detail in the following paper:

#include "swimming_DEM_application.h"
#include "difelice_drag_law.h"

namespace Kratos {

    DragLaw::Pointer DiFeliceDragLaw::Clone() const {
        DiFeliceDragLaw::Pointer p_clone(new DiFeliceDragLaw(*this));
        return p_clone;
    }

    void DiFeliceDragLaw::Initialize(const ProcessInfo& r_process_info) {}

    std::string DiFeliceDragLaw::GetTypeOfLaw() {
        std::string type_of_law = "DiFelice drag law";
        return type_of_law;
    }

    void DiFeliceDragLaw::ComputeForce(SphericParticle* p_particle,
                                       const double reynolds_number,
                                       double particle_radius,
                                       double fluid_density,
                                       double fluid_kinematic_viscosity,
                                       array_1d<double, 3>& minus_slip_velocity,
                                       array_1d<double, 3>& drag_force,
                                       const ProcessInfo& r_current_process_info)
    {

        Geometry<Node >& r_geometry = p_particle->GetGeometry();
        const double eps = r_geometry[0].FastGetSolutionStepValue(FLUID_FRACTION_PROJECTED);
        const array_1d<double, 3> fluid_vel = r_geometry[0].FastGetSolutionStepValue(FLUID_VEL_PROJECTED);

        const array_1d<double, 3> truly_slip_velocity = r_geometry[0].FastGetSolutionStepValue(SLIP_VELOCITY);
        const double norm_minus_slip_velocity = SWIMMING_MODULUS_3(truly_slip_velocity);
        const double reynolds_particle = eps * norm_minus_slip_velocity * 2.0*particle_radius / fluid_kinematic_viscosity;
        const double chi = 3.7 - 0.65 * std::exp(-std::pow(1.5 - std::log10(reynolds_particle), 2) / 2);

        double drag_coeff = std::pow((2.654 / (1.0 + 3.213) + 4.8 / (std::sqrt(reynolds_particle))),2);
        double drag_force_magnitude = 1.0 / 8.0 * drag_coeff * fluid_density  * Globals::Pi * std::pow(2.0*particle_radius, 2) * norm_minus_slip_velocity *  std::pow(eps, 2.0 - chi);
        double& drag_magnitude = r_geometry[0].FastGetSolutionStepValue(DRAG_COEFFICIENT);
        drag_magnitude = drag_force_magnitude;
        noalias(drag_force) = drag_force_magnitude * minus_slip_velocity;
    }

    array_1d<double, 3> DiFeliceDragLaw::CalculateParticleVelocityMean(SphericParticle* p_particle){

        //array_1d<double, 3> particle_velocity = p_particle->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
        array_1d<double, 3> total_velocities = p_particle->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
        double total_particles = 1;

        for (unsigned int i = 0; i < p_particle->mNeighbourElements.size(); ++i){
            for (unsigned int d = 0; d < 3; ++d)
                total_velocities[d] += p_particle->mNeighbourElements[i]->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY)[d];

            total_particles += 1;
        }

        return total_velocities / total_particles;
    }

} // namespace Kratos